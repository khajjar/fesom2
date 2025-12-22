!==========================================================
MODULE mod_transit
! Parameters, variables and functions for transient tracer simulations.
! By mbutzin, 2019-2021
! By khajjar, 2025

  implicit none
  save


! Atmospheric pressure, local (dummy) variable, global-mean SLP, and local wind speed at 10 m
  real(kind=8) :: press_a, mean_slp = 1.01325e5, wind_2

! Atmospheric trace gas (dummy) values used in air-sea flux calculations
! Isotopic ratios are normalized and fractionation-corrected,
! volume mixing ratios are mole fractions in dry air.
  real(kind=8) :: xf11_a  = 0.0, &       ! CFC-11 (latitude dependent)
                  xf12_a  = 0.0, &       ! CFC-12 (latitude dependent)
                  xsf6_a  = 0.0          ! SF6 (latitude dependent)

! Transient values of atmospheric trace gases (1d-arrays of variable length to be specified in namelist.config -> length_transit)
  real(kind=8), allocatable, dimension(:) ::   xf11_nh, xf11_sh, &          ! CFC-11 conc, latitude-dependent
                                                xf12_nh, xf12_sh, &          ! CFC-12 conc, latitude-dependent
                                                xsf6_nh, xsf6_sh             ! SF6 conc, latitude-dependent
                                                                                                                                        
  real(kind=8) :: f11t_nh=0.0, f11t_sh=0.0, &      ! CFC-11 trend, latitude-dependent
                  f12t_nh=0.0, f12t_sh=0.0, &      ! CFC-12 trend, latitude-dependent
                  sf6t_nh=0.0, sf6t_sh=0.0         ! SF6 trend, latitude-dependent

  integer, allocatable, dimension(:)      :: year_ce                      ! current year in anthropenic runs (control output)
! Further internal parameters
! Latitude of atmospheric boundary conditions and latitudinal interpolation weight
  real(kind=8) :: y_abc, yy_nh
! Tracer indices of transient tracers
  integer ::      id_f11 , id_f12, id_sf6, id_pf11 , id_pf12, id_psf6
! Time index (=year) in transient simulations
  integer ::      ti_transit
  
! Namelist to modify default parameter settings
   NAMELIST /transit_param/ xf11_nh, xf11_sh, xf12_nh, xf12_sh, xsf6_nh, xsf6_sh, f11t_nh, f11t_sh, f12t_nh, f12t_sh, sf6t_nh, sf6t_sh



  contains


    function iso_flux(strVarGas, tempC, sal, wind_2, iceFraction, p_atm, x_gas, r_air, r_sea, c_surf)
!     Calculate isotopic air-sea exchange fluxes in 1 / (m**2 * s) assuming local solubility equilibrium
!     for the abundant isotopologue. Positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: iso_flux
!     Input parameters
      character(len=3), intent(in) :: strVarGas    ! trace gas name

      real(kind=8), intent(in) :: tempC, sal, &    ! SST (deg C) and SSS ("PSU" or permil)
                                  wind_2, &        ! wind speed at 10 m heigth squared
                                  iceFraction, &   ! sea-ice fractional coverage
                                  p_atm, &         ! total atmospheric pressure (Pa)
                                  x_gas, &         ! atmospheric mole fraction of the abundant isotope
                                  r_air, r_sea, &  ! isotopic ratios in atmosphere and ocean
                                  c_surf           ! surface water concentration of the abundant isotope (mol / m**3)

      iso_flux = gasTransferVelocity(strVarGas, tempC, wind_2, iceFraction) * &
                 solub(strVarGas, tempC, sal) * (p_atm / 1.01325e5 * x_gas) * &
                 (r_air - r_sea) / c_surf
      return
    end function iso_flux


    function gas_flux(strVarGas, tempC, sal, wind_2, iceFraction, p_atm, x_gas, c_surf)
!     Computes air-sea exchange gas fluxes in mol / (m**2 * s) , positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: gas_flux
!     Input parameters
      character(len=3), intent(in) :: strVarGas        ! trace gas name
      real(kind=8), intent(in) :: tempC, sal, &        ! SST (deg C) and SSS ("PSU" or permil)
                                  wind_2,     &        ! wind speed at 10 m heigth squared
                                  iceFraction, &       ! sea-ice fractional coverage
                                  p_atm, &             ! total atmospheric pressure (Pa)
                                  x_gas, &             ! atmospheric mole fraction 
                                  c_surf               ! marine surface water concentration (mol / m**3)
!     Internal variables
      real(kind=8) :: c_sat                      ! marine saturation concentration (mol / m**3)
      c_sat = solub(strVarGas, tempC, sal) * (p_atm / 1.01325e5) * x_gas
      gas_flux = gasTransferVelocity(strVarGas, tempC, wind_2, iceFraction) * (c_sat - c_surf)

      return
    end function gas_flux
    

    function solub(strVarGas, tempC, sal)
!     Computes the solubility of trace gases in seawater.
      implicit none
      real(kind=8) :: solub                      ! solubility ((p)mol / (m**3 * atm))
!     Input parameters
      character(len=3), intent(in) :: strVarGas  ! tracer name
      real(kind=8), intent(in) :: tempC, &       ! temperature (deg C)
                                  sal            ! salinity ("PSU" or permil)
      real(kind=8) :: a1, a2, a3, a4, &          ! polynomial coefficients of the
                      b1, b2, b3, b4, c1, &      ! solubility function
                      tempK, &                   ! water temperature in K
                      con2con                    ! concentration units conversion factor
      integer ::      pow                        ! power in solubility function

      tempK = tempC + 273.15

      select case (strVarGas)
      case ("f11")
!       CFC-11 in mol / (L * atm) (Warner & Weiss 1985, doi:10.1016/0198-0149(85)90099-8, Table 5)
        a1 = -229.9261;  a2 = 319.6552;   a3 = 119.4471;   a4 = -1.39165; pow = 2
        b1 = -0.142382;  b2 = 0.091459;   b3 = -0.0157274
        con2con = 1000. ! convert to mol / (m**3 * atm)
      case ("f12") 
!       CFC-12 in mol / (L * atm) (Warner & Weiss 1985, doi:10.1016/0198-0149(85)90099-8, Table 5)
        a1 = -218.0971;  a2 = 298.9702;   a3 = 113.8049;   a4 = -1.39165; pow = 2
        b1 = -0.143566;  b2 = 0.091015;   b3 = -0.0153924
        con2con = 1000. ! convert to mol / (m**3 * atm)
      case ("sf6") 
!       SF6 in mol / (L * atm) (Bullister et al. 2002, doi:10.1016/S0967-0637(01)00051-6, Table 3)
        a1 = -80.0343;   a2 = 117.232;    a3 = 29.5817;    a4 = 0.;       pow = 2
        b1 =  0.0335183; b2 = -0.0373942; b3 = 0.00774862
        con2con = 1000. ! convert to mol / (m**3 * atm)
      end select

      solub = exp(a1 + a2 * (100.0/tempK) + a3 * log(tempK/100.0) + a4 * (tempK/100.0)**pow + & 
                  sal * (b1 + b2 * (tempK/100.0) + b3 * (tempK/100.0)**pow))
      solub = solub * con2con

      return
    end function solub


    function schmidtNumber(strVarGas, tempC)
!     Schmidt numbers of trace gases in sea water 
!     Coefficients are from Wanninkhof 2014, tab. 1.
      implicit none

!     Result
      real(kind=8) :: schmidtNumber                ! Schmidt number
!     Input parameters
      character(len=3), intent(in) :: strVarGas    ! tracer name
      real(kind=8), intent(in) :: tempC            ! temperature (deg C)
!     Internal parameters and/or variables
      real(kind=8) :: as, bs, cs, ds, es           ! polynomial coefficients

      select case (strVarGas)
      case ("f11") ! CFC-11
        as = 3579.2; bs = -222.63; cs = 7.5749; ds = -0.14595; es = 0.0011874
      case ("f12") ! CFC-12
        as = 3828.1; bs = -249.86; cs = 8.7603; ds = -0.171600; es = 0.0014080
      case ("sf6") ! SF6
        as = 3177.5; bs = -200.57; cs = 6.8865; ds = -0.133350; es = 0.0010877
      end select
      
      schmidtNumber = (as + bs *tempC + cs * tempC**2 + ds * tempC**3 + es * tempC**4)
      
      return
    end function schmidtNumber


    function gasTransferVelocity(strVarGas, tempC, wind_2, iceFraction)
!     Compute gas transfer velocities of / for tracers
      implicit none
!     Result
      real(kind=8) :: gasTransferVelocity                  ! transfer velocity (m / s)
!     Input parameters
      character(len=3), intent(in) :: strVarGas            ! tracer name
      real(kind=8), intent(in) :: tempC, &                 ! temperature (deg C)
                                  wind_2, &                ! wind speed squared at 10 m height (m / s)
                                  iceFraction              ! ice fraction

      gasTransferVelocity = 6.97e-7 * (schmidtNumber(strVarGas, tempC)/660)**(-0.5) * wind_2 * (1. - iceFraction)

      return
    end function gasTransferVelocity


    function speed_2(windstr_x, windstr_y)
!     Computes the square of wind speed at 10 m height from wind stress fields
!     in coupled simulations as long as it is not provided by the AGCM / OASIS.
!     We follow Peixoto & Oort (1992, Eq. (10.28), (10,29)) and Charnock (1955); 
!     also see MPI report 349 (2003), Eq. (5.7).
      implicit none
      real(kind=8) :: speed_2

!     Input
      real(kind=8), intent(in) :: windstr_x, windstr_y

!     Internal variables and parameters
!     Zonal and meridional velocities at 10 m height
      real(kind=8) :: u_10, v_10
!     Zonal and meridional friction velocities
      real(kind=8) :: u_fric, v_fric
!     Zonal and meridional roughness lengths
      real(kind=8) :: l_rough_x, l_rough_y
!     Inverse von-Karman constant (0.4), Charnock constant (0.018) divided by g, inverse density of air (1.3), log(10)
      real(kind=8), parameter :: inv_karm = 2.5, charn_g = 0.00173, inv_dens_air = 0.76923, log_10 = 2.30258
     
!     Calculate friction velocities (Peixoto & Oort, 1992, Eq. (10.28))
      u_fric = sqrt(abs(windstr_x) * inv_dens_air)
      v_fric = sqrt(abs(windstr_y) * inv_dens_air)

!     Calculate roughness lengths (MPI report 349, 2003, Eq. (5.7), quoting Charnock, 1955)
      l_rough_x = max((charn_g * u_fric**2), 1.5e-5)
      l_rough_y = max((charn_g * v_fric**2), 1.5e-5)

!     Calculate wind speed at 10 m (Peixoto & Oort, 1992, Eq. (10.29))
      u_10 = inv_karm * u_fric * (log_10 - log(l_rough_x))
      v_10 = inv_karm * v_fric * (log_10 - log(l_rough_y))
     
      speed_2 = u_10**2 + v_10**2
      
      return
    end function speed_2
    

    subroutine read_transit_input
!     Read atmospheric input of isoCO2 and / or other tracers
      use g_config, only: anthro_transit, paleo_transit, length_transit
      implicit none

!     Internal variables
      integer :: jj
!      real(kind=8), allocatable, dimension(:) :: d14c_nh, d14c_tz, d14c_sh, d14c_ti, d13c_dummy


      if (anthro_transit .or. paleo_transit) then
!       Anthropogenic input for 1850 - 2015 CE
	allocate(year_ce(length_transit))
	allocate(xf11_nh(length_transit))
	allocate(xf11_sh(length_transit))
	allocate(xf12_nh(length_transit))
	allocate(xf12_sh(length_transit))
	allocate(xsf6_nh(length_transit))
	allocate(xsf6_sh(length_transit))

!       Skip header lines
        do jj = 1,3
          read (20, fmt=*)
        end do
!       Read input values
        do jj = 1, length_transit
          read (20, fmt=*) year_ce(jj), &
                           xf11_nh(jj), xf11_sh(jj), & 
                           xf12_nh(jj), xf12_sh(jj), & 
                           xsf6_nh(jj), xsf6_sh(jj)
        end do

                
      else
!       Read constant parameter values from namelist.oce. 
	allocate(xf11_nh(1))
	allocate(xf11_sh(1))
	allocate(xf12_nh(1))
	allocate(xf12_sh(1))
	allocate(xsf6_nh(1))
	allocate(xsf6_sh(1))
	read (20,NML=transit_param)
      end if

!     Convert volume mixing ratios
      xf11_nh = xf11_nh * 1.e-12
      xf11_sh = xf11_sh * 1.e-12
      xf12_nh = xf12_nh * 1.e-12
      xf12_sh = xf12_sh * 1.e-12
      xsf6_nh = xsf6_nh * 1.e-12
      xsf6_sh = xsf6_sh * 1.e-12
      f11t_nh = f11t_nh * 1.e-12
      f11t_sh = f11t_sh * 1.e-12
      f12t_nh = f12t_nh * 1.e-12
      f12t_sh = f12t_sh * 1.e-12
      sf6t_nh = sf6t_nh * 1.e-12
      sf6t_sh = sf6t_sh * 1.e-12

      return
    end subroutine read_transit_input
    
END MODULE mod_transit
!==========================================================


