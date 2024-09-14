#include <define.h>

SUBROUTINE CLM_LSM(pressure,saturation,evap_trans,topo,porosity,pf_dz_mult,istep_pf,dt,time,           &
   start_time,pdx,pdy,pdz,ix,iy,nx,ny,nz,nx_f,ny_f,nz_f,nz_rz,ip,npp,npq,npr,gnx,gny,rank,sw_pf,lw_pf,    &
   prcp_pf,tas_pf,u_pf,v_pf,patm_pf,qatm_pf,lai_pf,sai_pf,z0m_pf,displa_pf,                               &
   slope_x_pf,slope_y_pf,                                                                                 &
   eflx_lh_pf,eflx_lwrad_pf,eflx_sh_pf,eflx_grnd_pf,                                                     &
   qflx_tot_pf,qflx_grnd_pf,qflx_soi_pf,qflx_eveg_pf,qflx_tveg_pf,qflx_in_pf,swe_pf,t_g_pf,               &
   t_soi_pf,clm_dump_interval,clm_1d_out,clm_forc_veg,clm_output_dir,clm_output_dir_length,clm_bin_output_dir,         &
   write_CLM_binary,slope_accounting_CLM,beta_typepf,veg_water_stress_typepf,wilting_pointpf,field_capacitypf,                 &
   res_satpf,irr_typepf, irr_cyclepf, irr_ratepf, irr_startpf, irr_stoppf, irr_thresholdpf,               &
   qirr_pf,qirr_inst_pf,irr_flag_pf,irr_thresholdtypepf,soi_z,clm_next,clm_write_logs,                    &
   clm_last_rst,clm_daily_rst, pf_nlevsoi, pf_nlevlak)
!-----------------------------------------------------------------------------
! Description:
!   This is the main program for the Common Land Model (CoLM)
!
!   @Copyright Yongjiu Dai Land Modeling Grop at the School of Atmospheric Sciences
!   of the Sun Yat-sen University, Guangdong, CHINA.
!   All rights reserved.
!
! Initial : Yongjiu Dai, 1998-2014
! Revised : Hua Yuan, Shupeng Zhang, Nan Wei, Xingjie Lu, Zhongwang Wei, Yongjiu Dai
!           2014-2024
!-----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_Const_Physical
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   !USE MOD_Vars_2DForcing
   USE MOD_Vars_1DFluxes
   !USE MOD_Vars_1DAccFluxes
   !USE MOD_Forcing
   !USE MOD_Hist
   USE MOD_TimeManager
   !USE MOD_RangeCheck

   !USE MOD_Block
   !USE MOD_Pixel
   !USE MOD_Mesh
   !USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   !USE MOD_LandPatch
#ifdef URBAN_MODEL
   USE MOD_LandUrban
   USE MOD_Urban_LAIReadin
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_ElmVector
#endif
#ifdef CATCHMENT
   USE MOD_HRUVector
#endif
#if(defined CaMa_Flood)
   USE MOD_CaMa_colmCaMa ! whether cama-flood is used
#endif
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

#if (defined CatchLateralFlow)
   USE MOD_Catch_LateralFlow
#endif

   !USE MOD_Ozone, only: init_ozone_data, update_ozone_data

   !USE MOD_SrfdataRestart
   !USE MOD_LAIReadin

#ifdef BGC
   USE MOD_NitrifData
   USE MOD_NdepData
   USE MOD_FireData
   USE MOD_LightningData
#endif

#ifdef LULCC
   USE MOD_Lulcc_Driver
#endif

#ifdef CoLMDEBUG
   USE MOD_Hydro_SoilWater
#endif

   ! SNICAR
   !USE MOD_SnowSnicar, only: SnowAge_init, SnowOptics_init
   !USE MOD_Aerosol, only: AerosolDepInit, AerosolDepReadin

#ifdef DataAssimilation
   USE MOD_DataAssimilation
#endif

#ifdef USEMPI
   USE MOD_HistWriteBack
#endif

   IMPLICIT NONE

   character(len=256) :: nlfile
   character(len=256) :: casename
   character(len=256) :: dir_landdata
   character(len=256) :: dir_forcing
   character(len=256) :: dir_hist
   character(len=256) :: dir_restart
   character(len=256) :: fsrfdata

   real(r8) :: deltim       ! time step (senconds)
   integer  :: sdate(3)     ! calendar (year, julian day, seconds)
   integer  :: idate(3)     ! calendar (year, julian day, seconds)
   integer  :: edate(3)     ! calendar (year, julian day, seconds)
   integer  :: pdate(3)     ! calendar (year, julian day, seconds)
   integer  :: jdate(3)     ! calendar (year, julian day, seconds), year beginning style
   logical  :: greenwich    ! greenwich time

   logical :: doalb         ! true => start up the surface albedo calculation
   logical :: dolai         ! true => start up the time-varying vegetation paramter
   logical :: dosst         ! true => update sst/ice/snow

   integer :: Julian_1day_p, Julian_1day
   integer :: Julian_8day_p, Julian_8day
   integer :: s_year, s_month, s_day, s_seconds, s_julian
   integer :: e_year, e_month, e_day, e_seconds, e_julian
   integer :: p_year, p_month, p_day, p_seconds, p_julian
   integer :: lc_year, lai_year
   integer :: month, mday, year_p, month_p, mday_p
   integer :: spinup_repeat, istep

   type(timestamp) :: ststamp, itstamp, etstamp, ptstamp

   integer*8 :: start_time, end_time, c_per_sec, time_used

   integer  :: numpatch

   !@CY: below from parflow=================================
   integer  :: pf_nlevsoi                         ! number of soil levels, passed from PF
   integer  :: pf_nlevlak                         ! number of lake levels, passed from PF
  
   !=== Local Variables =====================================================
 
   ! basic indices, counters
   integer  :: t                                   ! tile space counter
   integer  :: l,ll                                   ! layer counter 
   integer  :: r,c                                 ! row,column indices
   !integer  :: ierr                                ! error output 
 
   ! values passed from parflow
   integer  :: nx,ny,nz,nx_f,ny_f,nz_f,nz_rz
   integer  :: soi_z                               ! NBE: Specify layer shold be used for reference temperature
   real(r8) :: pressure((nx+2)*(ny+2)*(nz+2))     ! pressure head, from parflow on grid w/ ghost nodes for current proc
   real(r8) :: saturation((nx+2)*(ny+2)*(nz+2))   ! saturation from parflow, on grid w/ ghost nodes for current proc
   real(r8) :: evap_trans((nx+2)*(ny+2)*(nz+2))   ! ET flux from CLM to ParFlow on grid w/ ghost nodes for current proc
   real(r8) :: topo((nx+2)*(ny+2)*(nz+2))         ! mask from ParFlow 0 for inactive, 1 for active, on grid w/ ghost nodes for current proc
   real(r8) :: porosity((nx+2)*(ny+2)*(nz+2))     ! porosity from ParFlow, on grid w/ ghost nodes for current proc
   real(r8) :: pf_dz_mult((nx+2)*(ny+2)*(nz+2))   ! dz multiplier from ParFlow on PF grid w/ ghost nodes for current proc
   real(r8) :: dt                                 ! parflow dt in parflow time units not CLM time units
   real(r8) :: time                               ! parflow time in parflow units
   real(r8) :: start_time                         ! starting time in parflow units
   real(r8) :: pdx,pdy,pdz                        ! parflow DX, DY and DZ in parflow units
   integer  :: istep_pf                           ! istep, now passed from PF
   integer  :: ix                                 ! parflow ix, starting point for local grid on global grid
   integer  :: iy                                 ! parflow iy, starting point for local grid on global grid
   integer  :: ip                               
   integer  :: npp,npq,npr                        ! number of processors in x,y,z
   integer  :: gnx, gny                           ! global grid, nx and ny
   integer  :: rank                               ! processor rank, from ParFlow
 
   integer :: clm_next                           ! NBE: Passing flag to sync outputs
   integer :: d_stp                              ! NBE: Dummy for CLM restart
   integer :: clm_write_logs                     ! NBE: Enable/disable writing of the log files
   integer :: clm_last_rst                       ! NBE: Write all the CLM restart files or just the last one
   integer :: clm_daily_rst                      ! NBE: Write daily restart files or hourly
 
   ! surface fluxes & forcings
   real(r8) :: eflx_lh_pf((nx+2)*(ny+2)*3)        ! e_flux   (lh)    output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: eflx_lwrad_pf((nx+2)*(ny+2)*3)     ! e_flux   (lw)    output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: eflx_sh_pf((nx+2)*(ny+2)*3)        ! e_flux   (sens)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: eflx_grnd_pf((nx+2)*(ny+2)*3)      ! e_flux   (grnd)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_tot_pf((nx+2)*(ny+2)*3)       ! h2o_flux (total) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_grnd_pf((nx+2)*(ny+2)*3)      ! h2o_flux (grnd)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_soi_pf((nx+2)*(ny+2)*3)       ! h2o_flux (soil)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_eveg_pf((nx+2)*(ny+2)*3)      ! h2o_flux (veg-e) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_tveg_pf((nx+2)*(ny+2)*3)      ! h2o_flux (veg-t) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: qflx_in_pf((nx+2)*(ny+2)*3)        ! h2o_flux (infil) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: swe_pf((nx+2)*(ny+2)*3)            ! swe              output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: t_g_pf((nx+2)*(ny+2)*3)            ! t_grnd           output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
   real(r8) :: t_soi_pf((nx+2)*(ny+2)*(pf_nlevsoi+2))!tsoil             output var to send to ParFlow, on grid w/ ghost nodes for current proc, but nz=10 (3D)
   real(r8) :: sw_pf((nx+2)*(ny+2)*3)             ! SW rad, passed from PF
   real(r8) :: lw_pf((nx+2)*(ny+2)*3)             ! LW rad, passed from PF
   real(r8) :: prcp_pf((nx+2)*(ny+2)*3)           ! Precip, passed from PF
   real(r8) :: tas_pf((nx+2)*(ny+2)*3)            ! Air temp, passed from PF
   real(r8) :: u_pf((nx+2)*(ny+2)*3)              ! u-wind, passed from PF
   real(r8) :: v_pf((nx+2)*(ny+2)*3)              ! v-wind, passed from PF
   real(r8) :: patm_pf((nx+2)*(ny+2)*3)           ! air pressure, passed from PF
   real(r8) :: qatm_pf((nx+2)*(ny+2)*3)           ! air specific humidity, passed from PF
   real(r8) :: lai_pf((nx+2)*(ny+2)*3)            ! BH: lai, passed from PF
   real(r8) :: sai_pf((nx+2)*(ny+2)*3)            ! BH: sai, passed from PF
   real(r8) :: z0m_pf((nx+2)*(ny+2)*3)            ! BH: z0m, passed from PF
   real(r8) :: displa_pf((nx+2)*(ny+2)*3)         ! BH: displacement height, passed from PF
   real(r8) :: irr_flag_pf((nx+2)*(ny+2)*3)       ! irrigation flag for deficit-based scheduling -- 1 = irrigate, 0 = no-irrigate
   real(r8) :: qirr_pf((nx+2)*(ny+2)*3)           ! irrigation applied above ground -- spray or drip (2D)
   real(r8) :: qirr_inst_pf((nx+2)*(ny+2)*(pf_nlevsoi+2))! irrigation applied below ground -- 'instant' (3D)
 
   real(r8) :: slope_x_pf((nx+2)*(ny+2)*3)        ! Slope in x-direction from PF
   real(r8) :: slope_y_pf((nx+2)*(ny+2)*3)        ! Slope in y-direction from PF
 
   ! output keys
   integer :: clm_dump_interval                  ! dump inteval for CLM output, passed from PF, always in interval of CLM timestep, not time
   integer  :: clm_1d_out                         ! whether to dump 1d output 0=no, 1=yes
   integer  :: clm_forc_veg                       ! BH: whether vegetation (LAI, SAI, z0m, displa) is being forced 0=no, 1=yes
   integer  :: clm_output_dir_length              ! for output directory
   integer  :: clm_bin_output_dir                 ! output directory
   integer  :: write_CLM_binary                   ! whether to write CLM output as binary 
   integer  :: slope_accounting_CLM               ! account for slope is solar zenith angle calculations
   character (LEN=clm_output_dir_length) :: clm_output_dir ! output dir location
 
   ! ET keys
   integer  :: beta_typepf                        ! beta formulation for bare soil Evap 0=none, 1=linear, 2=cos
   integer  :: veg_water_stress_typepf            ! veg transpiration water stress formulation 0=none, 1=press, 2=sm
   real(r8) :: wilting_pointpf                    ! wilting point in m if press-type, in saturation if soil moisture type
   real(r8) :: field_capacitypf                   ! field capacity for water stress same as units above
   real(r8) :: res_satpf                          ! residual saturation from ParFlow
 
   ! irrigation keys
   integer  :: irr_typepf                         ! irrigation type flag (0=none,1=spray,2=drip,3=instant)
   integer  :: irr_cyclepf                        ! irrigation cycle flag (0=constant,1=deficit)
   real(r8) :: irr_ratepf                         ! irrigation application rate for spray and drip [mm/s]
   real(r8) :: irr_startpf                        ! irrigation daily start time for constant cycle
   real(r8) :: irr_stoppf                         ! irrigation daily stop tie for constant cycle
   real(r8) :: irr_thresholdpf                    ! irrigation threshold criteria for deficit cycle (units of soil moisture content)
   integer  :: irr_thresholdtypepf                ! irrigation threshold criteria type -- top layer, bottom layer, column avg
 
   ! local indices & counters
   integer  :: i,j,k,k1,j1,l1                     ! indices for local looping
   !integer  :: bj,bl                              ! indices for local looping !BH
 
   integer  :: j_incr,k_incr                      ! increment for j and k to convert 1D vector to 3D i,j,k array
   integer, allocatable :: counter(:,:) 
   real(r8) :: total
   character*100 :: RI
   real(r8) :: u         ! Tempoary UNDEF Variable  

   save
   !@CY: we don't want to initialize everything every timestep

      !@CY:
      numpatch = nx*ny

if (time == start_time) then !initialization
#ifdef USEMPI
      CALL spmd_init ()
#endif

      CALL getarg (1, nlfile)

      CALL read_namelist (nlfile)

#ifdef USEMPI
      IF (DEF_HIST_WriteBack) THEN
         CALL spmd_assign_writeback ()
      ENDIF

      !IF (p_is_writeback) THEN
      !   CALL hist_writeback_daemon ()
      !ELSE
#endif

      !IF (p_is_master) THEN
      !   CALL system_clock (start_time)
      !ENDIF

      casename     = DEF_CASE_NAME
      dir_landdata = DEF_dir_landdata
      dir_forcing  = DEF_dir_forcing
      dir_hist     = DEF_dir_history
      dir_restart  = DEF_dir_restart

#ifdef SinglePoint
      fsrfdata = trim(dir_landdata) // '/srfdata.nc'
#ifndef URBAN_MODEL
      CALL read_surface_data_single (fsrfdata, mksrfdata=.false.)
#else
      CALL read_urban_surface_data_single (fsrfdata, mksrfdata=.false., mkrun=.true.)
#endif
#endif

      deltim    = DEF_simulation_time%timestep
      greenwich = DEF_simulation_time%greenwich
      s_year    = DEF_simulation_time%start_year
      s_month   = DEF_simulation_time%start_month
      s_day     = DEF_simulation_time%start_day
      s_seconds = DEF_simulation_time%start_sec
      e_year    = DEF_simulation_time%end_year
      e_month   = DEF_simulation_time%end_month
      e_day     = DEF_simulation_time%end_day
      e_seconds = DEF_simulation_time%end_sec
      p_year    = DEF_simulation_time%spinup_year
      p_month   = DEF_simulation_time%spinup_month
      p_day     = DEF_simulation_time%spinup_day
      p_seconds = DEF_simulation_time%spinup_sec

      spinup_repeat = DEF_simulation_time%spinup_repeat

      CALL initimetype(greenwich)
      CALL monthday2julian(s_year,s_month,s_day,s_julian)
      CALL monthday2julian(e_year,e_month,e_day,e_julian)
      CALL monthday2julian(p_year,p_month,p_day,p_julian)

      sdate(1) = s_year; sdate(2) = s_julian; sdate(3) = s_seconds
      edate(1) = e_year; edate(2) = e_julian; edate(3) = e_seconds
      pdate(1) = p_year; pdate(2) = p_julian; pdate(3) = p_seconds

      CALL Init_GlobalVars
      CAll Init_LC_Const
      CAll Init_PFT_Const

      !CALL pixel%load_from_file    (dir_landdata)
      !CALL gblock%load_from_file   (dir_landdata)

#ifdef LULCC
      lc_year = s_year
#else
      lc_year = DEF_LC_YEAR
#endif

      !CALL mesh_load_from_file (dir_landdata, lc_year)

      !CALL pixelset_load_from_file (dir_landdata, 'landelm'  , landelm  , numelm  , lc_year)

#ifdef CATCHMENT
      CALL pixelset_load_from_file (dir_landdata, 'landhru'  , landhru  , numhru  , lc_year)
#endif

      !CALL pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch, lc_year)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL pixelset_load_from_file (dir_landdata, 'landpft'  , landpft  , numpft  , lc_year)
      CALL map_patch_to_pft
#endif

#ifdef URBAN_MODEL
      CALL pixelset_load_from_file (dir_landdata, 'landurban', landurban, numurban, lc_year)
      CALL map_patch_to_urban
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
      CALL elm_vector_init ()
#ifdef CATCHMENT
      CALL hru_vector_init ()
#endif
#endif

      CALL adj2end(sdate)
      CALL adj2end(edate)
      CALL adj2end(pdate)

      ststamp = sdate
      etstamp = edate
      ptstamp = pdate

      ! date in beginning style
      jdate = sdate
      CALL adj2begin(jdate)

      !IF (ptstamp <= ststamp) THEN
      !   spinup_repeat = 0
      !ELSE
      !   spinup_repeat = max(0, spinup_repeat)
      !ENDIF

      ! ----------------------------------------------------------------------
      ! Read in the model time invariant constant data
      CALL allocate_TimeInvariants (numpatch)
      !CALL READ_TimeInvariants (lc_year, casename, dir_restart)

      ! Read in the model time varying data (model state variables)
      CALL allocate_TimeVariables  (numpatch)
      !CALL READ_TimeVariables (jdate, lc_year, casename, dir_restart)

      ! Read in SNICAR optical and aging parameters
      !CALL SnowOptics_init( DEF_file_snowoptics ) ! SNICAR optical parameters
      !CALL SnowAge_init( DEF_file_snowaging )     ! SNICAR aging   parameters

      ! ----------------------------------------------------------------------
      doalb = .true.
      dolai = .true.
      dosst = .false.

      ! Initialize meteorological forcing data module
      CALL allocate_1D_Forcing (numpatch)
      !CALL forcing_init (dir_forcing, deltim, ststamp, lc_year, etstamp)
      !CALL allocate_2D_Forcing (gforc)

      ! Initialize history data module
      ! CALL hist_init (dir_hist)
      CALL allocate_1D_Fluxes (numpatch)


#if(defined CaMa_Flood)
      CALL colm_CaMa_init !initialize CaMa-Flood
#endif

      IF(DEF_USE_OZONEDATA)THEN
         !CALL init_Ozone_data (sdate)
      ENDIF

      ! Initialize aerosol deposition forcing data
      IF (DEF_Aerosol_Readin) THEN
         !CALL AerosolDepInit ()
      ENDIF

#ifdef BGC
      IF (DEF_USE_NITRIF) THEN
         CALL init_nitrif_data (sdate)
      ENDIF

      IF (DEF_NDEP_FREQUENCY==1)THEN ! Initial annual ndep data readin
         CALL init_ndep_data_annually (sdate(1))
      ELSEIF(DEF_NDEP_FREQUENCY==2)THEN ! Initial monthly ndep data readin
         CALL init_ndep_data_monthly (sdate(1),s_month) ! sf_add
      ELSE
         write(6,*) 'ERROR: DEF_NDEP_FREQUENCY should be only 1-2, Current is:', &
                     DEF_NDEP_FREQUENCY
         CALL CoLM_stop ()
      ENDIF

      IF (DEF_USE_FIRE) THEN
         CALL init_fire_data (sdate(1))
         CALL init_lightning_data (sdate)
      ENDIF
#endif

#if (defined CatchLateralFlow)
      CALL lateral_flow_init (lc_year)
#endif

#ifdef DataAssimilation
      CALL init_DataAssimilation ()
#endif

      idate   = sdate
      itstamp = ststamp

endif

      ! ======================================================================
      ! begin time stepping loop
      ! ======================================================================

      istep   = istep_pf


      !TIMELOOP : DO WHILE (itstamp < etstamp)

         CALL julian2monthday (jdate(1), jdate(2), month_p, mday_p)

         year_p = jdate(1)

         IF (p_is_master) THEN
            !IF (itstamp < ptstamp) THEN
            !   write(*, 99) istep, jdate(1), month_p, mday_p, jdate(3), spinup_repeat
            !ELSE
               write(*,100) istep, jdate(1), month_p, mday_p, jdate(3)
            !ENDIF
         ENDIF


         Julian_1day_p = int(calendarday(jdate)-1)/1*1 + 1
         Julian_8day_p = int(calendarday(jdate)-1)/8*8 + 1

         ! Read in the meteorological forcing
         ! ----------------------------------------------------------------------
         ! CALL read_forcing (jdate, dir_forcing)

         IF(DEF_USE_OZONEDATA)THEN
            !CALL update_Ozone_data(itstamp, deltim)
         ENDIF
#ifdef BGC
         IF(DEF_USE_FIRE)THEN
            CALL update_lightning_data (itstamp, deltim)
         ENDIF
#endif

         ! Read in aerosol deposition forcing data
         IF (DEF_Aerosol_Readin) THEN
            !CALL AerosolDepReadin (jdate)
         ENDIF

         ! Calendar for NEXT time step
         ! ----------------------------------------------------------------------
         CALL TICKTIME (deltim,idate)
         itstamp = itstamp + int(deltim)
         jdate = idate
         CALL adj2begin(jdate)

         CALL julian2monthday (jdate(1), jdate(2), month, mday)

#ifdef BGC
         IF(DEF_USE_NITRIF) THEN
            IF (month /= month_p) THEN
               CALL update_nitrif_data (month)
            ENDIF
         ENDIF

         IF (DEF_NDEP_FREQUENCY==1)THEN ! Read Annual Ndep data
            IF (jdate(1) /= year_p) THEN
               CALL update_ndep_data_annually (idate(1), iswrite = .true.)
            ENDIF
         ELSEIF(DEF_NDEP_FREQUENCY==2)THEN! Read Monthly Ndep data
            IF (jdate(1) /= year_p .or. month /= month_p) THEN  !sf_add
               CALL update_ndep_data_monthly (jdate(1), month, iswrite = .true.) !sf_add
            ENDIF
         ELSE
            write(6,*) 'ERROR: DEF_NDEP_FREQUENCY should be only 1-2, Current is:',&
                        DEF_NDEP_FREQUENCY
            CALL CoLM_stop ()
         ENDIF

         IF(DEF_USE_FIRE)THEN
            IF (jdate(1) /= year_p) THEN
               CALL update_hdm_data (idate(1))
            ENDIF
         ENDIF
#endif


         ! Call colm driver
         ! ----------------------------------------------------------------------
         IF (p_is_worker) THEN
            CALL CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oroflag,numpatch)
         ENDIF


#if (defined CatchLateralFlow)
         CALL lateral_flow (deltim)
#endif

#if(defined CaMa_Flood)
         CALL colm_CaMa_drv(idate(3)) ! run CaMa-Flood
#endif

#ifdef DataAssimilation
         CALL do_DataAssimilation (idate, deltim)
#endif

         ! Write out the model variables for restart run and the histroy file
         ! ----------------------------------------------------------------------
         !CALL hist_out (idate, deltim, itstamp, etstamp, ptstamp, dir_hist, casename)

         ! DO land USE and land cover change simulation
         ! ----------------------------------------------------------------------
#ifdef LULCC
         IF ( isendofyear(idate, deltim) ) THEN
            CALL deallocate_1D_Forcing (numpatch)
            CALL deallocate_1D_Fluxes (numpatch)

            CALL LulccDriver (casename,dir_landdata,dir_restart,&
                              idate,greenwich)

            CALL allocate_1D_Forcing (numpatch)
            CALL forcing_init (dir_forcing, deltim, itstamp, jdate(1))
            !CALL deallocate_acc_fluxes
            CALL hist_init (dir_hist)
            CALL allocate_1D_Fluxes (numpatch)
         ENDIF
#endif

         ! Get leaf area index
         ! ----------------------------------------------------------------------
#if(defined DYN_PHENOLOGY)
         ! Update once a day
         dolai = .false.
         Julian_1day = int(calendarday(jdate)-1)/1*1 + 1
         IF(Julian_1day /= Julian_1day_p)THEN
            dolai = .true.
         ENDIF
#else
         ! READ in Leaf area index and stem area index
         ! ----------------------------------------------------------------------
         ! Hua Yuan, 08/03/2019: read global monthly LAI/SAI data
         ! zhongwang wei, 20210927: add option to read non-climatological mean LAI
         ! Update every 8 days (time interval of the MODIS LAI data)
         ! Hua Yuan, 06/2023: change namelist DEF_LAI_CLIM to DEF_LAI_MONTHLY
         ! and add DEF_LAI_CHANGE_YEARLY for monthly LAI data
         !
         ! NOTES: Should be caution for setting DEF_LAI_CHANGE_YEARLY to ture in non-LULCC
         ! case, that means the LAI changes without condisderation of land cover change.

         IF (DEF_LAI_CHANGE_YEARLY) THEN
            lai_year = jdate(1)
         ELSE
            lai_year = DEF_LC_YEAR
         ENDIF

         IF (DEF_LAI_MONTHLY) THEN
            IF ((itstamp < etstamp) .and. (month /= month_p)) THEN
               !CALL LAI_readin (lai_year, month, dir_landdata)
#ifdef URBAN_MODEL
               CALL UrbanLAI_readin(lai_year, month, dir_landdata)
#endif
            ENDIF
         ELSE
            ! Update every 8 days (time interval of the MODIS LAI data)
            Julian_8day = int(calendarday(jdate)-1)/8*8 + 1
            IF ((itstamp < etstamp) .and. (Julian_8day /= Julian_8day_p)) THEN
               !CALL LAI_readin (jdate(1), Julian_8day, dir_landdata)
               ! 06/2023, yuan: or depend on DEF_LAI_CHANGE_YEARLY nanemlist
               !CALL LAI_readin (lai_year, Julian_8day, dir_landdata)
            ENDIF
         ENDIF
#endif

         IF (save_to_restart (idate, deltim, itstamp, ptstamp)) THEN
#ifdef LULCC
            CALL WRITE_TimeVariables (jdate, jdate(1), casename, dir_restart)
#else
            !CALL WRITE_TimeVariables (jdate, lc_year,  casename, dir_restart)
#endif
#if(defined CaMa_Flood)
            IF (p_is_master) THEN
               CALL colm_cama_write_restart (jdate, lc_year,  casename, dir_restart)
            ENDIF
#endif
         ENDIF
#ifdef RangeCheck
         CALL check_TimeVariables ()
#endif
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CoLMDEBUG
         CALL print_VSF_iteration_stat_info ()
#endif


         !IF (p_is_master) THEN
         !   CALL system_clock (end_time, count_rate = c_per_sec)
         !   time_used = (end_time - start_time) / c_per_sec
         !   IF (time_used >= 3600) THEN
         !      write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
         !   ELSEIF (time_used >= 60) THEN
         !      write(*,102) time_used/60, mod(time_used,60)
         !   ELSE
         !      write(*,103) time_used
         !   ENDIF
         !ENDIF

         !IF ((spinup_repeat > 1) .and. (ptstamp <= itstamp)) THEN
         !   spinup_repeat = spinup_repeat - 1
         !   idate   = sdate
         !   jdate   = sdate
         !   itstamp = ststamp
         !   CALL adj2begin(jdate)
         !   !CALL forcing_reset ()
         !ENDIF

         !istep = istep + 1

      !ENDDO TIMELOOP

      !CALL deallocate_TimeInvariants (numpatch)
      !CALL deallocate_TimeVariables  (numpatch)
      !CALL deallocate_1D_Forcing     (numpatch)
      !CALL deallocate_1D_Fluxes      (numpatch)

#if (defined CatchLateralFlow)
      CALL lateral_flow_final ()
#endif

      !CALL forcing_final ()
      !CALL hist_final    ()

#ifdef SinglePoint
      CALL single_srfdata_final ()
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#if(defined CaMa_Flood)
      CALL colm_cama_exit ! finalize CaMa-Flood
#endif

#ifdef DataAssimilation
      CALL final_DataAssimilation ()
#endif

      IF (p_is_master) THEN
         write(*,'(/,A25)') 'CoLM Execution Completed.'
      ENDIF

      99  format(/, 'TIMESTEP = ', I0, ' | DATE = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5, ' Spinup (', I0, ' repeat left)')
      100 format(/, 'TIMESTEP = ', I0, ' | DATE = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5)
      101 format(/, 'Time elapsed : ', I4, ' hours', I3, ' minutes', I3, ' seconds.')
      102 format(/, 'Time elapsed : ', I3, ' minutes', I3, ' seconds.')
      103 format(/, 'Time elapsed : ', I3, ' seconds.')

#ifdef USEMPI
      !ENDIF

      IF (DEF_HIST_WriteBack) THEN
         CALL hist_writeback_exit ()
      ENDIF

      CALL spmd_exit
#endif

END SUBROUTINE CLM_LSM
! ---------- EOP ------------
