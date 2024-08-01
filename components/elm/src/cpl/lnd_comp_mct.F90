module lnd_comp_mct
  
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Interface of the active land model component of CESM the ELM (E3SM Land Model)
  !  with the main E3SM driver. This is a thin interface taking E3SM driver information
  !  in MCT (Model Coupling Toolkit) format and converting it to use by ELM.
  !
  ! !uses:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_sys_mod      , only : shr_sys_flush
  use mct_mod          , only : mct_avect, mct_gsmap
  use decompmod        , only : bounds_type, ldecomp
  use lnd_import_export
  !
  ! !public member functions:
  implicit none
  save
  private                     ! by default make data private
  !
  ! !public member functions:
  public :: lnd_init_mct      ! elm initialization
  public :: lnd_run_mct       ! elm run phase
  public :: lnd_final_mct     ! elm finalization/cleanup
  !
  ! !private member functions:
  private :: lnd_setgsmap_mct ! set the land model mct gs map
  private :: lnd_domain_mct   ! set the land model domain information
  !---------------------------------------------------------------------------

contains

  !====================================================================================

  subroutine lnd_init_mct( EClock, cdata_l, x2l_l, l2x_l, NLFilename )
#if defined(CLDERA_PROFILING)
   use iso_c_binding, only: c_loc
   use cldera_interface_mod, only: cldera_init, cldera_set_log_unit, &
                                   cldera_set_masterproc, max_str_len, &
                                   cldera_add_partitioned_field, &
                                   cldera_set_field_part_extent, &
                                   cldera_set_field_part_data, &
                                   cldera_commit_all_fields,   &
                                   cldera_commit_field, cldera_switch_context
    use elm_varpar      , only : nlevgrnd
    use elm_instMod     , only : solarabs_vars, surfrad_vars, &
                                 veg_es, veg_wf, &
                                 canopystate_vars, photosyns_vars
    use GridcellType   , only : grc_pp
    use ColumnType     , only : col_pp
    use ColumnDataType , only : col_ws
    use VegetationType , only : veg_pp
#endif
    !
    ! !DESCRIPTION:
    ! Initialize land surface model and obtain relevant atmospheric model arrays
    ! back from (i.e. albedos, surface temperature and snow cover over land).
    !
    ! !USES:
    use abortutils       , only : endrun
    use shr_kind_mod     , only : SHR_KIND_CL
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, set_nextsw_cday
    use elm_initializeMod, only : initialize1, initialize2, initialize3
    use elm_instMod      , only : lnd2atm_vars, lnd2glc_vars
    use elm_instance     , only : elm_instance_init
    use elm_varctl       , only : finidat,single_column, elm_varctl_set, iulog, noland
    use elm_varctl       , only : inst_index, inst_suffix, inst_name, precip_downscaling_method
    use elm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds
    use domainMod        , only : ldomain
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod     , only : shr_file_getLogUnit, shr_file_getLogLevel
    use shr_file_mod     , only : shr_file_getUnit, shr_file_setIO
    use shr_taskmap_mod  , only : shr_taskmap_write
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use seq_comm_mct     , only : info_taskmap_comp
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use seq_flds_mod     , only : seq_flds_x2l_fields, seq_flds_l2x_fields
    use spmdMod          , only : masterproc, npes, spmd_init
    use elm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use elm_cpl_indices  , only : elm_cpl_indices_set
    use perf_mod         , only : t_startf, t_stopf
    use mct_mod
    use ESMF
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock           ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_l          ! Input land-model driver data
    type(mct_aVect),            intent(inout) :: x2l_l, l2x_l     ! land model import and export states
    character(len=*), optional, intent(in)    :: NLFilename       ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    integer                          :: LNDID	     ! Land identifyer
    integer                          :: mpicom_lnd   ! MPI communicator
    type(mct_gsMap),         pointer :: GSMap_lnd    ! Land model MCT GS map
    type(mct_gGrid),         pointer :: dom_l        ! Land model domain
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer  :: lsz                                  ! size of attribute vector
    integer  :: g,i,j                                ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_elm                            ! elm time-step
    logical  :: exists                               ! true if file exists
    logical  :: verbose_taskmap_output               ! true then use verbose task-to-node mapping format
    logical  :: atm_aero                             ! Flag if aerosol data sent from atm model
    logical  :: atm_present                          ! Flag if atmosphere model present
    real(r8) :: scmlat                               ! single-column latitude
    real(r8) :: scmlon                               ! single-column longitude
    real(r8) :: nextsw_cday                          ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=8)           :: c_inst_index       ! instance number           
    character(len=8)           :: c_npes             ! number of pes
    integer :: nsrest                                ! elm restart type
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: nstep
#if defined(CLDERA_PROFILING)
    character(len=max_str_len) :: fname
    integer :: p, c, nfields, idx, rank, icmp, nparts, part_dim, ipart, fsize, ncols, icall, tag_loop
    integer :: nlcols,irank,part_alloc_size
    integer :: dims(3)
    character(len=max_str_len) :: dimnames(3)
    ! arrays for copied data
    integer, allocatable :: cols_gids(:)
    real(r8), allocatable :: cols_area(:)
    real(r8), allocatable :: subgrid_lat(:)
    real(r8), allocatable :: subgrid_lon(:)
    real(r8), allocatable :: subgrid_wts_sum(:)
    ! pointers for data arrays
    real(r8), pointer :: field1d(:), field2d(:,:), field3d(:,:,:)
    integer, pointer :: intfield1d(:)
    character(len=5) :: int_str
    ! bounds for different data types
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
#endif
    type(bounds_type) :: bounds                      ! bounds
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    !-----------------------------------------------------------------------

    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    ! Set and save LNDID for easy access by other modules

    call elm_instance_init( LNDID )

    ! Determine attriute vector indices

    call elm_cpl_indices_set()

    ! Initialize elm MPI communicator 

    call spmd_init( mpicom_lnd, LNDID )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_mct:start::',lbnum)
    endif
#endif                      

    inst_name   = seq_comm_name(LNDID)
    inst_index  = seq_comm_inst(LNDID)
    inst_suffix = seq_comm_suffix(LNDID)

    ! Initialize io log unit

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "ELM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Identify SMP nodes and process/SMP mapping for this instance
    ! (Assume that processor names are SMP node names on SMP clusters.)
    write(c_inst_index,'(i8)') inst_index

    if (info_taskmap_comp > 0) then

       if (info_taskmap_comp == 1) then
          verbose_taskmap_output = .false.
       else
          verbose_taskmap_output = .true.
       endif

       write(c_npes,'(i8)') npes

       if (masterproc) then
          write(iulog,'(/,3A)') &
             trim(adjustl(c_npes)), &
             ' pes participating in computation of ELM instance #', &
             trim(adjustl(c_inst_index))
          call shr_sys_flush(iulog)
       endif

       call t_startf("shr_taskmap_write")
       call shr_taskmap_write(iulog, mpicom_lnd,                    &
                              'LND #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output        )
       call t_stopf("shr_taskmap_write")

    endif

    ! Use infodata to set orbital values

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Consistency check on namelist filename	

    call control_setNL("lnd_in"//trim(inst_suffix))

    ! Initialize elm
    ! initialize1 reads namelist, grid and surface data (need this to initialize gsmap) 
    ! initialize2 performs rest of initialization	

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, case_name=caseid,    &
                              case_desc=ctitle, single_column=single_column,    &
                              scmlat=scmlat, scmlon=scmlon,                     &
                              brnch_retain_casename=brnch_retain_casename,      &
                              start_type=starttype, model_version=version,      &
                              hostname=hostname, username=username )
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod)
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

#if defined(CLDERA_PROFILING)
    ! Initialize CLDERA profiling before elm_init, but after time init
    call t_startf('cldera_init')
    call cldera_init("elm",mpicom_lnd,start_ymd,start_tod,ref_ymd,ref_tod,stop_ymd,stop_tod)
    call cldera_set_log_unit (iulog)
    call cldera_set_masterproc (masterproc)
    call t_stopf('cldera_init')
#endif

    call elm_varctl_set(caseid_in=caseid, ctitle_in=ctitle,                     &
                        brnch_retain_casename_in=brnch_retain_casename,         &
                        single_column_in=single_column, scmlat_in=scmlat,       &
                        scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                        hostname_in=hostname, username_in=username)

    ! Read namelist, grid and surface data

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland ) then
       call seq_infodata_PutData( infodata, lnd_present   =.false.)
       call seq_infodata_PutData( infodata, lnd_prognostic=.false.)
       return
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    call seq_infodata_GetData(infodata, atm_present=atm_present)
    call seq_infodata_GetData(infodata, atm_aero=atm_aero )
    !DMR 6/12/15 - remove this requirement (CPL_BPYASS mode uses SATM)
    if ( .not. atm_aero .and. atm_present )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to ELM' )
    end if

    ! Initialize elm gsMap, elm domain and elm attribute vectors

    call get_proc_bounds( bounds )

    call lnd_SetgsMap_mct( bounds, mpicom_lnd, LNDID, gsMap_lnd ) 	
    lsz = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)

    call lnd_domain_mct( bounds, lsz, gsMap_lnd, dom_l )

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsz)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsz)
    call mct_aVect_zero(l2x_l)

    ! Finish initializing elm

    call initialize2()
    call initialize3()

#if defined(CLDERA_PROFILING)
    ! this region needs to be after initialize3, which sets elm_varpar%nlevgrnd, needed for H2OSOI
    call t_startf('cldera_add_fields')
    call cldera_switch_context("elm")
    begp = bounds%begp; endp = bounds%endp ! these are the bounds for data defined on subgrid cells
    begc = bounds%begc; endc = bounds%endc ! likely unused
    begl = bounds%begl; endl = bounds%endl ! likely unused
    begg = bounds%begg; endg = bounds%endg ! local bounds for columns
    
    ! All fields are partitioned over cols index, which is the first
    part_dim = 1
    nparts = 1

    !
    ! Start with data on the coarse grid (begg:endg) for lat/lon/gid/area
    !
    ncols = endg - begg + 1
    dimnames(1) = "ncol"
    dims(1) = ncols
    part_alloc_size = -1 ! let cldera_tools automatically size arrays

    if (masterproc) then
       write(int_str,'(I5)'), ncols
       write(iulog,*)'GH LND_INIT_MCT with '//int_str//' columns'
       write(iulog,*)'[cldera profiling] Start adding ELM fields'
       call shr_sys_flush(iulog)
    endif
   
    ! 1d latitude
    field1d => grc_pp%latdeg
    call cldera_add_partitioned_field("lat",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("lat",1,ncols)
    call cldera_set_field_part_data("lat",1,field1d)
    call cldera_commit_field("lat")

    ! 1d longitude
    field1d => grc_pp%londeg
    call cldera_add_partitioned_field("lon",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("lon",1,ncols)
    call cldera_set_field_part_data("lon",1,field1d)
    call cldera_commit_field("lon")

    ! Add GIDs and areas
    allocate(cols_gids(ncols))
    allocate(cols_area(ncols))
    allocate(subgrid_wts_sum(ncols))

    ! NOTE: .false. is to declare the field as a copy of input data, rather than a view
    call cldera_add_partitioned_field("col_gids",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.,"int")
    call cldera_add_partitioned_field("col_area",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
    call cldera_set_field_part_extent("col_gids",1,ncols)
    call cldera_set_field_part_extent("col_area",1,ncols)
    call cldera_commit_field("col_gids")
    call cldera_commit_field("col_area")

    do g = 1,ncols
      ! multiply area by land fraction; otherwise data will be biased toward coastal land
      cols_area(g) = ldomain%area(g) * ldomain%frac(g)
      ! initialize column weights
      cols_gids(g) = 1
      subgrid_wts_sum(g) = 0
    enddo
    
    call cldera_set_field_part_data("col_gids",1,cols_gids)
    call cldera_set_field_part_data("col_area",1,cols_area)

    !
    ! Now work on data defined on the subgrid
    !
    ncols = endp - begp + 1
    dims(1) = ncols

    ! Add lat/lon data for subgrid cells (this is the lat/lon of the corresponding coarse cell)
    ! lat
    allocate(subgrid_lat(ncols))
    field1d => grc_pp%latdeg
    do p = 1,ncols
      g = veg_pp%gridcell(p)
      subgrid_lat(p) = field1d(g)
    enddo
    call cldera_add_partitioned_field("subgrid_lat",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
    call cldera_set_field_part_extent("subgrid_lat",1,ncols)
    call cldera_commit_field("subgrid_lat")
    call cldera_set_field_part_data("subgrid_lat",1,subgrid_lat)
    ! lon
    allocate(subgrid_lon(ncols))
    field1d => grc_pp%londeg
    do p = 1,ncols
      g = veg_pp%gridcell(p)
      subgrid_lon(p) = field1d(g)
    enddo
    call cldera_add_partitioned_field("subgrid_lon",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
    call cldera_set_field_part_extent("subgrid_lon",1,ncols)
    call cldera_commit_field("subgrid_lon")
    call cldera_set_field_part_data("subgrid_lon",1,subgrid_lon)

    ! Add subgrid to column mapping data (performed on local ids, "lids")
    intfield1d => veg_pp%gridcell
    call cldera_add_partitioned_field("col_lids",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("col_lids",1,ncols)
    call cldera_set_field_part_data("col_lids",1,intfield1d)
    call cldera_commit_field("col_lids")

    ! Add subgrid weights relative to the individual column
    field1d => veg_pp%wtcol
    call cldera_add_partitioned_field("subgrid_wts",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("subgrid_wts",1,ncols)
    call cldera_set_field_part_data("subgrid_wts",1,field1d)
    call cldera_commit_field("subgrid_wts")
    
    ! Validate subgrid weights
    do p = 1,ncols
      g = veg_pp%gridcell(p)
      subgrid_wts_sum(g) = subgrid_wts_sum(g) + field1d(p)
    enddo

    !
    ! Add miscellaneous column metadata fields that may be important
    !

    ! logical: whether a column is active, meaning ELM performs calculations on the column
    ! field1d => col_pp%active
    ! call cldera_add_partitioned_field("active",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    ! call cldera_set_field_part_extent("active",1,ncols)
    ! call cldera_set_field_part_data("active",1,field1d)
    ! call cldera_commit_field("active")

    !
    ! Add data for QOIs
    !

    ! Add SolarRadiation variables
    field1d => surfrad_vars%fsds_vis_d_patch
    call cldera_add_partitioned_field("FSDSVD",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("FSDSVD",1,ncols)
    call cldera_set_field_part_data("FSDSVD",1,field1d)
    call cldera_commit_field("FSDSVD")

    field1d => surfrad_vars%fsds_vis_i_patch
    call cldera_add_partitioned_field("FSDSVI",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("FSDSVI",1,ncols)
    call cldera_set_field_part_data("FSDSVI",1,field1d)
    call cldera_commit_field("FSDSVI")

    ! Add SolarAbsorbedType variables
    field1d => solarabs_vars%fsds_nir_d_patch
    call cldera_add_partitioned_field("FSDSND",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("FSDSND",1,ncols)
    call cldera_set_field_part_data("FSDSND",1,field1d)
    call cldera_commit_field("FSDSND")

    field1d => solarabs_vars%fsds_nir_i_patch
    call cldera_add_partitioned_field("FSDSNI",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("FSDSNI",1,ncols)
    call cldera_set_field_part_data("FSDSNI",1,field1d)
    call cldera_commit_field("FSDSNI")

    ! Add VegetationDataType variables
    field1d => veg_es%t_veg
    call cldera_add_partitioned_field("TV",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("TV",1,ncols)
    call cldera_set_field_part_data("TV",1,field1d)
    call cldera_commit_field("TV")

    field1d => veg_wf%qflx_tran_veg
    call cldera_add_partitioned_field("QVEGT",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("QVEGT",1,ncols)
    call cldera_set_field_part_data("QVEGT",1,field1d)
    call cldera_commit_field("QVEGT")

    ! Add CanopyStateType variables
    field1d => canopystate_vars%tlai_hist_patch
    call cldera_add_partitioned_field("TLAI",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("TLAI",1,ncols)
    call cldera_set_field_part_data("TLAI",1,field1d)
    call cldera_commit_field("TLAI")

    ! Add PhotosynthesisType variables
    field1d => photosyns_vars%fpsn_patch
    call cldera_add_partitioned_field("FPSN",1,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("FPSN",1,ncols)
    call cldera_set_field_part_data("FPSN",1,field1d)
    call cldera_commit_field("FPSN")

    !
    ! 2d fields
    !
    dims(2) = nlevgrnd
    dimnames(2) = 'dz'
    field2d => col_pp%dz
    call cldera_add_partitioned_field("dz",2,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("dz",1,ncols)
    call cldera_set_field_part_data("dz",1,field2d)
    call cldera_commit_field("dz")

    field2d => col_ws%h2osoi_vol
    call cldera_add_partitioned_field("H2OSOI",2,dims,dimnames,nparts,part_dim,part_alloc_size)
    call cldera_set_field_part_extent("H2OSOI",1,ncols)
    call cldera_set_field_part_data("H2OSOI",1,field2d)
    call cldera_commit_field("H2OSOI")

    if (masterproc) then
      write(iulog,*)'[cldera profiling] Finished adding ELM fields'
    endif

    if (masterproc) then
      write(iulog, *) 'begp', begp
      write(iulog, *) 'endp', endp
      write(iulog, *) 'begc', begc
      write(iulog, *) 'endc', endc
      write(iulog, *) 'begl', begl
      write(iulog, *) 'endl', endl
      write(iulog, *) 'begg', begg
      write(iulog, *) 'endg', endg

      field1d => grc_pp%latdeg
      write(iulog, *) 'lat', field1d(:)
      field1d => grc_pp%londeg
      write(iulog, *) 'lon', field1d(:)
      write(iulog, *) 'cols_area', cols_area(:)
      write(iulog, *) 'cols_gids', cols_gids(:)
      write(iulog, *) 'subgrid_lat', subgrid_lat(:)
      write(iulog, *) 'subgrid_lon', subgrid_lon(:)
      write(iulog, *) 'subgrid_wts_sum', subgrid_wts_sum(:)
      call shr_sys_flush(iulog)
    end if

    call cldera_commit_all_fields()
    call t_stopf('cldera_add_fields')
#endif

    ! Check that elm internal dtime aligns with elm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_elm = get_step_size()
    if (masterproc) then
       write(iulog,*)'dtime_sync= ',dtime_sync,&
            ' dtime_elm= ',dtime_elm,' mod = ',mod(dtime_sync,dtime_elm)
    end if
    if (mod(dtime_sync,dtime_elm) /= 0) then
       write(iulog,*)'elm dtime ',dtime_elm,' and Eclock dtime ',&
            dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state 

    if (atm_present) then 
      call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, l2x_l%rattr)
    endif

    ! Fill in infodata settings

    call seq_infodata_PutData(infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData(infodata, lnd_nx=ldomain%ni, lnd_ny=ldomain%nj, precip_downscaling_method = precip_downscaling_method)

    ! Get infodata info

    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )
    call set_nextsw_cday(nextsw_cday)

    if (.not. atm_present) then 
      !Calculate next radiation calendar day (since atm model did not run to set
      !this)
      !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
      nstep = get_nstep()
      nextsw_cday = mod((nstep/(86400._r8/dtime_elm))*1.0_r8,365._r8)+1._r8
      call set_nextsw_cday( nextsw_cday )
    end if

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine lnd_init_mct

  !====================================================================================

  subroutine lnd_run_mct(EClock, cdata_l, x2l_l, l2x_l)
#if defined(CLDERA_PROFILING)
    use cldera_interface_mod , only: cldera_switch_context
#endif
    !
    ! !DESCRIPTION:
    ! Run elm model
    !
    ! !USES:
    use shr_kind_mod    ,  only : r8 => shr_kind_r8
    use elm_instMod     , only : lnd2atm_vars, atm2lnd_vars, lnd2glc_vars, glc2lnd_vars
    use elm_driver      ,  only : elm_drv
    use clm_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use clm_time_manager,  only : advance_timestep, set_nextsw_cday,update_rad_dtime
    use decompMod       ,  only : get_proc_bounds
    use abortutils      ,  only : endrun
    use elm_varctl      ,  only : iulog
    use elm_varorb      ,  only : eccen, obliqr, lambm0, mvelpp
    use shr_file_mod    ,  only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod    ,  only : shr_file_getLogUnit, shr_file_getLogLevel
    use seq_cdata_mod   ,  only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,  only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
    use seq_timemgr_mod ,  only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use seq_infodata_mod,  only : seq_infodata_type, seq_infodata_GetData
    use spmdMod         ,  only : masterproc, mpicom
    use perf_mod        ,  only : t_startf, t_stopf, t_barrierf
    use shr_orb_mod     ,  only : shr_orb_decl
    use mct_mod
    use ESMF
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    !
    ! !LOCAL VARIABLES:
    integer      :: ymd_sync             ! Sync date (YYYYMMDD)
    integer      :: yr_sync              ! Sync current year
    integer      :: mon_sync             ! Sync current month
    integer      :: day_sync             ! Sync current day
    integer      :: tod_sync             ! Sync current time of day (sec)
    integer      :: ymd                  ! ELM current date (YYYYMMDD)
    integer      :: yr                   ! ELM current year
    integer      :: mon                  ! ELM current month
    integer      :: day                  ! ELM current day
    integer      :: tod                  ! ELM current time of day (sec)
    integer      :: dtime                ! time step increment (sec)
    integer      :: nstep                ! time step index
    logical      :: rstwr_sync           ! .true. ==> write restart file before returning
    logical      :: rstwr                ! .true. ==> write restart file before returning
    logical      :: nlend_sync           ! Flag signaling last time-step
    logical      :: nlend                ! .true. ==> last time-step
    logical      :: dosend               ! true => send data back to driver
    logical      :: doalb                ! .true. ==> do albedo calculation on this time step
    real(r8)     :: nextsw_cday          ! calday from clock of next radiation computation
    real(r8)     :: caldayp1             ! elm calday plus dtime offset
    integer      :: shrlogunit,shrloglev ! old values for share log unit and log level
    integer      :: lbnum                ! input to memory diagnostic
    integer      :: g,i,lsz              ! counters
    real(r8)     :: calday               ! calendar day for nstep
    real(r8)     :: declin               ! solar declination angle in radians for nstep
    real(r8)     :: declinp1             ! solar declination angle in radians for nstep+1
    real(r8)     :: eccf                 ! earth orbit eccentricity factor
    real(r8)     :: recip                ! reciprical
    logical,save :: first_call = .true.  ! first call work
    logical      :: atm_present
    type(seq_infodata_type),pointer :: infodata             ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_l                ! Land model domain data
    type(bounds_type)               :: bounds               ! bounds
    character(len=32)               :: rdate                ! date char string for restart file names
    character(len=32), parameter    :: sub = "lnd_run_mct"
    !---------------------------------------------------------------------------

    ! Determine processor bounds

    call get_proc_bounds(bounds)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    call seq_cdata_setptrs(cdata_l, infodata=infodata, dom=dom_l)
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )
    dtime = get_step_size()

    call seq_infodata_GetData(infodata, atm_present=atm_present)
    if (.not. atm_present) then 
      !Calcualte next radiation calendar day (since atm model did not run to set this)
      !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
      nstep = get_nstep()
      nextsw_cday = mod((nstep/(86400._r8/dtime))*1.0_r8,365._r8)+1._r8 
      call set_nextsw_cday( nextsw_cday )
    end if
 
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    ! Map MCT to land data type
    ! Perform downscaling if appropriate

    
    ! Map to elm (only when state and/or fluxes need to be updated)

    call t_startf ('lc_lnd_import')
    call lnd_import( bounds, x2l_l%rattr, atm2lnd_vars, glc2lnd_vars, lnd2atm_vars)
    call t_stopf ('lc_lnd_import')

    ! Use infodata to set orbital values if updated mid-run

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

#if defined(CLDERA_PROFILING)
    call cldera_switch_context("elm")
#endif

    ! Loop over time steps in coupling interval

    dosend = .false.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       nstep = get_nstep()
       caldayp1 = get_curr_calday(offset=dtime)
       if (nstep == 0) then
	  doalb = .false. 	
       else if (nstep == 1) then 
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8) 
       else
          doalb = (nextsw_cday >= -0.5_r8) 
       end if
       call update_rad_dtime(doalb)

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run elm 

       call t_barrierf('sync_elm_run1', mpicom)
       call t_startf ('elm_run')
       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')
       call elm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('elm_run')

       ! Create l2x_l export state - add river runoff input to l2x_l if appropriate

#ifndef CPL_BYPASS       
       call t_startf ('lc_lnd_export')
       call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, l2x_l%rattr)
       call t_stopf ('lc_lnd_export')
#endif

       ! Advance elm time step
       
       call t_startf ('lc_elm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_elm2_adv_timestep')

    end do

    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' elm ymd=',ymd     ,'  elm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: ELM clock not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call  = .false.

  end subroutine lnd_run_mct

  !====================================================================================

  subroutine lnd_final_mct( EClock, cdata_l, x2l_l, l2x_l)
    !
    ! !DESCRIPTION:
    ! Finalize land surface model

    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
    use seq_timemgr_mod ,only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use mct_mod
    use esmf
    use elm_finalizeMod, only : final
#if defined(CLDERA_PROFILING)
    use cldera_interface_mod, only: cldera_clean_up, cldera_switch_context
    use perf_mod            , only: t_startf, t_stopf
#endif
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    !---------------------------------------------------------------------------

    ! fill this in
    call final()

#if defined(CLDERA_PROFILING)
    call t_startf('cldera_clean_up')
    call cldera_switch_context("elm")
    call cldera_clean_up ()
    call t_stopf('cldera_clean_up')
#endif
  end subroutine lnd_final_mct

  !====================================================================================

  subroutine lnd_setgsmap_mct( bounds, mpicom_lnd, LNDID, gsMap_lnd )
    !
    ! !DESCRIPTION:
    ! Set the MCT GS map for the land model
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use domainMod    , only : ldomain
    use mct_mod      , only : mct_gsMap, mct_gsMap_init
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds     ! bounds
    integer           , intent(in)  :: mpicom_lnd ! MPI communicator for the elm land model
    integer           , intent(in)  :: LNDID      ! Land model identifyer number
    type(mct_gsMap)   , intent(out) :: gsMap_lnd  ! Resulting MCT GS map for the land model
    !
    ! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: i, j, n, gi            ! Indices
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code
    !---------------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    allocate(gindex(bounds%begg:bounds%endg),stat=ier)

    ! number the local grid

    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1
    gsize = ldomain%ni * ldomain%nj

    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

  !====================================================================================

  subroutine lnd_domain_mct( bounds, lsz, gsMap_l, dom_l )
    !
    ! !DESCRIPTION:
    ! Send the land model domain information to the coupler
    !
    ! !USES:
    use elm_varcon  , only: re
    use domainMod   , only: ldomain
    use spmdMod     , only: iam
    use mct_mod     , only: mct_gsMap, mct_gGrid, mct_gGrid_importIAttr
    use mct_mod     , only: mct_gGrid_importRAttr, mct_gGrid_init, mct_gsMap_orderedPoints
    use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other
    !
    ! !ARGUMENTS: 
    type(bounds_type), intent(in)  :: bounds  ! bounds
    integer        , intent(in)    :: lsz     ! land model domain data size
    type(mct_gsMap), intent(inout) :: gsMap_l ! Output land model MCT GS map
    type(mct_ggrid), intent(out)   :: dom_l   ! Output domain information for land model
    !
    ! Local Variables
    integer :: g,i,j              ! index
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !---------------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_l, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsz )
    !
    ! Allocate memory
    !
    allocate(data(lsz))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsz)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsz)
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsz)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"frac",data,lsz)

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct

end module lnd_comp_mct
