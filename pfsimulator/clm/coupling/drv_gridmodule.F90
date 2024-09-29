!#include <misc.h>

module drv_gridmodule 

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Module for grid space variable specification.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!=========================================================================

  use MOD_Precision
  !use clm_varpar, only : max_nlevsoi
  use MOD_Vars_Global, only : nl_soil
  implicit none
  public griddec

  type griddec

!=== GRID SPACE User-Defined Parameters ==================================

     integer :: patchclass
     real(r8):: patchlatr
     real(r8):: patchlonr
     real(r8):: vf_quartz(nl_soil)  !soil, nc readin ##not used at all
     real(r8):: vf_gravels(nl_soil) !soil, nc readin ##used in soil_hcap_cond, maybe setup in 2nd file
     real(r8):: vf_om(nl_soil)      !soil, nc readin ##used in soil_hcap_cond, maybe setup in 2nd file 
     real(r8):: vf_sand(nl_soil)    !soil, nc readin             
     real(r8):: wf_gravels(nl_soil) !soil, nc readin
     real(r8):: wf_sand(nl_soil)    !soil, nc readin   
     !real(r8):: porsl(1:,i)      !soil, nc readin, from pf 
     real(r8):: psi0(nl_soil)       !soil, nc readin              
     real(r8):: bsw(nl_soil)        !soil, nc readin 
     !real(r8):: theta_r(nl_soil)    !soil, nc readin                                                
     real(r8):: alpha_vgm(nl_soil)  !soil, nc readin 
     real(r8):: n_vgm(nl_soil)      !soil, nc readin 
     real(r8):: L_vgm(nl_soil)      !soil, nc readin  
     real(r8):: hksati(nl_soil)     !soil, nc readin 
     real(r8):: csol(nl_soil)       !soil, nc readin     
     real(r8):: k_solids(nl_soil)   !soil, nc readin   
     real(r8):: dksatu(nl_soil)     !soil, nc readin            
     real(r8):: dksatf(nl_soil)     !soil, nc readin     
     real(r8):: dkdry(nl_soil)      !soil, nc readin                                           
     real(r8):: BA_alpha(nl_soil)   !soil, nc readin  
     real(r8):: BA_beta(nl_soil)    !soil, nc readin                                        
     !real(r8):: htoplc              !Forest_Height.nc 
     real(r8):: OM_density(nl_soil)   !soil, nc readin 
     real(r8):: BD_all(nl_soil)       !soil, nc readin   


! Leaf constants

     !real(r8) :: dewmx     ! Maximum allowed dew [mm]

! Roughness lengths

     !real(r8) :: zlnd      ! Roughness length for soil [m]
     !real(r8) :: zsno      ! Roughness length for snow [m]
     !real(r8) :: csoilc    ! Drag coefficient for soil under canopy [-]

! Soil parameters

     !real(r8) :: wtfact    ! Fraction of model area with high water table
     !real(r8) :: trsmx0    ! Max transpiration for moist soil+100% veg. [mm/s]
     !real(r8) :: scalez    ! Soil layer thickness discretization (m)
     !real(r8) :: hkdepth   ! length scale for Ksat decrease(m)

! Land surface parameters

     !real(r8) :: latdeg    ! Latitude in Degrees 
     !real(r8) :: londeg    ! Longitude in Degrees

     !real(r8) :: sand(max_nlevsoi) ! Percent sand in soil 
     !real(r8) :: clay(max_nlevsoi) ! Percent clay in soil 

     !real(r8), pointer :: fgrd(:) !Fraction of vegetation class in grid     
     !integer , pointer :: pveg(:) !Predominance of vegetation class in grid

     !integer :: mask      ! Land=1, Not Land (i.e. not modeled)=0

!=== CLM Forcing parameters

     !real(r8) :: forc_hgt_u ! Observational height of wind [m]
     !real(r8) :: forc_hgt_t ! Observational height of temperature [m]
     !real(r8) :: forc_hgt_q ! Observational height of humidity [m] 

!=== Land Surface Fluxes

     !real(r8) :: qflx_evap_tot   ! evapotranspiration from canopy height to atmosphere [mm/s]
     !real(r8) :: eflx_sh_tot     ! sensible heat from canopy height to atmosphere [W/m2]
     !real(r8) :: eflx_lh_tot     ! latent heat flux from canopy height to atmosphere [W/2]
     !real(r8) :: eflx_lwrad_out  ! outgoing long-wave radiation from ground+canopy
     !real(r8) :: t_ref2m         ! 2 m height air temperature [K]
     !real(r8) :: t_rad           ! radiative temperature [K]

!=== CLM Vegetation parameters

     !real(r8) :: rootfr            ! Root Fraction (depth average)

!=== CLM Soil parameters

     !real(r8) :: smpmax       ! Wilting point potential in mm
     !integer  :: isoicol      ! Soil color index

!=== Numerical finite-difference

     !real(r8) :: capr         ! Tuning factor to turn first layer T into surface T
     !real(r8) :: cnfac        ! Crank Nicholson factor between 0 and 1
     !real(r8) :: smpmin       ! Restriction for min of soil poten. (mm)
     !real(r8) :: ssi          ! Irreducible water saturation of snow
     !real(r8) :: wimp         ! Water impremeable if porosity < wimp
     !real(r8) :: pondmx       ! Ponding depth (mm)

     !integer  :: tilei        ! Tile index at x,y (or c,r); used to convert tile data into grid data;
                              ! works only if there is a single land cover type per grid cell!!

!=== End Variable List ===================================================

  end type griddec

end module drv_gridmodule




