! Module for X-ray cluster fgas data
! (c) 2013, 2020 Adam Mantz, amantz@stanford.edu

module fgas
  use constants, only : G, eV, m_p, Mpc, mass_ratio_He_H, c, const_pi
  use Precision
  use settings
  use likelihood
  use Likelihood_Cosmology
  use CosmologyTypes
  use Calculator_Cosmology
  use CosmoTheory
  implicit none
  private

  type, extends(TCosmoCalcLikelihood) :: FgasLikelihood
     integer :: numdatasets
     integer, allocatable, dimension(:) :: num_clusters
     real(dl), allocatable, dimension(:) :: z, dpar
     ! these could be made allocatable
     real(dl), dimension(19+12) :: mpar ! magic number (see wrapper.hpp)
     real(dl), dimension(5) :: const ! magic number (see wrapper.hpp)
   contains
     procedure :: LogLike => Fgas_LnLike
  end type FgasLikelihood

  Type(FgasLikelihood), pointer :: like

  !! options
  ! generate a fake dataset at the redshifts of the input data file
  logical :: fgas_simulate = .false.
  ! mass pivot for fgas(M) scaling relation
  real(dl) :: massPivot = 3.0e14
  ! use lensing data to calibrate X-ray masses
  logical :: calibrate = .true.
  ! use the fgas data
  logical :: do_fgas = .true.

  ! enumerate cluster model parameters, for clarity
  integer, parameter :: U_baryon_0    =  1
  integer, parameter :: U_baryon_1    =  2
  integer, parameter :: fgas_scatter  =  3
  integer, parameter :: fgas_rslope   =  4
  integer, parameter :: Mtot_rslope   =  5
  integer, parameter :: Mscaling      =  6
  integer, parameter :: cl_cal        =  7
  integer, parameter :: cl_calev      =  8
  integer, parameter :: cl_calscat    =  9
  integer, parameter :: cl_lenssys    = 10
  integer, parameter :: cl_lenssysev  = 11
  integer, parameter :: cl_masses     = 12

  real(dl), parameter :: Msun = 1.9891e30 ! solar mass in kg

  public FgasLikelihood, FgasLikelihood_Add, like
contains

  function rhocrit(this, z)
    Class(FgasLikelihood) :: this
    real(dl) :: z
    real(dl) :: rhocrit
    real(dl) :: H
    H = (c * this%Calculator%Hofz(z)) ! m/s/Mpc
    rhocrit = 3 * H**2  * Mpc / (8 * const_pi * G * Msun) ! Msun/Mpc^3
  end function rhocrit
    
  subroutine FgasLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    Type(TSettingIni) :: ini

    integer :: fclwrapperinit, fclwrapperloaddata, fclwrapperdatasetsize, fclwrappernumdatapars, fclwrappernumclusters
    real(dl) :: fclwrappergetredshift
    external fclwrapperinit, fclwrapperloaddata, fclwrapperdatasetsize, fclwrappernumdatapars, fclwrappergetredshift, fclwrappersetdataoptions, fclwrappersetconst, fclwrappernumclusters, fclwrappersetmasspivot

    integer :: numdatasets, i, j, k

    if (Ini%Read_Logical('use_fgas', .false.)) then

       ! read data
       numdatasets = Ini%Read_Int('fgas_numdatasets', 0)
       if (numdatasets<1) call MpiStop('use_fgas but fgas_numdatasets = 0')
       if (Ini%Haskey('fgas_simulate')) then
          Fgas_simulate = Ini%Read_Logical('fgas_simulate', fgas_simulate)
       end if
       if (Ini%Haskey('fgas_massPivot')) then
          massPivot = Ini%Read_Double('fgas_massPivot', massPivot)
       end if
       if (Ini%Haskey('fgas_do_lensing_cal')) then
          calibrate = Ini%Read_Logical('fgas_do_lensing_cal', calibrate)
       end if
       if (Ini%Haskey('fgas_do_fgas')) then
          do_fgas = Ini%Read_Logical('fgas_do_fgas', do_fgas)
       end if
       allocate(like)
       like%LikelihoodType = 'Fgas'
       like%name = 'Fgas'
       like%needs_background_functions = .true.
       if (.not. fclwrapperinit('')) call MpiStop('Error initializing fgas module')
       call fclwrappersetmasspivot(massPivot)
       i = 0
       if (calibrate) i = 1
       j = 0
       if (do_fgas) j = 1
       call fclwrappersetdataoptions(i, j)
       allocate(like%num_clusters( numdatasets ))
       do i= 1, numdatasets
          if (.not. fclwrapperloaddata( trim(Ini%Read_String(numcat('fgas_dataset', i))) )) call MpiStop('Error reading in a fgas dataset: ' // trim(Ini%Read_String(numcat('fgas_dataset', i))))
          like%num_clusters(i) = fclwrapperdatasetsize(i-1)
       end do
       allocate(like%z( fclwrappernumclusters() ))
       allocate(like%dpar( fclwrappernumdatapars() ))

       ! get redshifts
       k = 1
       do i = 1, numdatasets
          do j = 1, like%num_clusters(i)
             like%z(k) = fclwrappergetredshift(i-1, j-1)
             k = k + 1
          end do
       end do

       ! set physical contants for consistency
       like%const(1) = G / (1.0e3*eV) / Mpc * Msun**2 ! G in keV Mpc/Msun^2
       like%const(2) = m_p / Msun ! proton mass in Msun
       like%const(3) = Mpc * 100.0 ! 1 Mpc in cm
       like%const(4) = mass_ratio_He_H ! ratio of alpha to proton mass
       like%const(5) = c**2/G * Mpc / Msun ! c^2/G in Msun/Mpc
       call fclwrappersetconst(like%const)

       ! load cluster model parameters
       call Like%loadParamNames('fgas/cosmomc/fgas.paramnames')

       like%numdatasets = numdatasets

       like%needs_background_functions = .true.
       call LikeList%Add(like)
       if (Feedback > 0) write(*,*) 'read fgas datasets'
    end if
    
  end subroutine FgasLikelihood_Add


  function Fgas_LnLike(this, CMB, Theory, DataParams)
    Class(CMBParams) CMB
    Class(FgasLikelihood) :: this
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) Fgas_LnLike

    integer :: fclwrapperloadparameters, fclwrappernumclusters
    real(dl) :: fclwrapperlnp
    external fclwrapperloadparameters, fclwrapperlnp, fclwrappersimulate, fclwrappernumclusters

    integer :: i, j, k, l, k1
    real(dl) :: junk

    if (Feedback > 2) write(*,*) 'cluster parameters  = ', DataParams

    ! load parameters
    this%mpar(1)  = DataParams(U_baryon_0)
    this%mpar(2)  = DataParams(U_baryon_1)
    this%mpar(3)  = 0.0 !exp(DataParams(U_coldgas_0)) ! NB sampling the log
    this%mpar(4)  = 0.0 ! DataParams(U_coldgas_1)
    this%mpar(5)  = DataParams(fgas_scatter)
    this%mpar(6)  = DataParams(fgas_rslope)
    this%mpar(7)  = DataParams(Mtot_rslope)
    this%mpar(8)  = DataParams(Mscaling)
    this%mpar(9)  = CMB%omb / (CMB%omb + CMB%omdm)
    this%mpar(10) = CMB%YHe
    this%mpar(11) = -1. ! using YHe rather than ne_nH,
    this%mpar(12) = -1. ! ntot_ne,
    this%mpar(13) = -1. ! mu
    this%mpar(14) = 1.0 !DataParams(cl_nonth)
    this%mpar(15) = DataParams(cl_cal)
    this%mpar(16) = DataParams(cl_calev)
    this%mpar(17) = DataParams(cl_calscat)
    this%mpar(18) = DataParams(cl_lenssys)
    this%mpar(19) = DataParams(cl_lenssysev)
    do i=20, size(this%mpar)
       this%mpar(i) = DataParams(cl_masses+i-20)
    end do
    if (.not. fclwrapperloadparameters( this%mpar )) call MpiStop('Error setting model parameters in fgas module')

    ! calculate and load data parameters
    k1 = 1
    l = 1
    do i = 1, this%numdatasets
       k = k1
       do j = 1, this%num_clusters(i)
          this%dpar(l) = this%Calculator%AngularDiameterDistance(this%z(k))
          k = k + 1
          l = l + 1
       end do
       k = k1
       do j = 1, this%num_clusters(i)
          this%dpar(l) = this%Calculator%LuminosityDistance(this%z(k))
          k = k + 1
          l = l + 1
       end do
       k = k1
       do j = 1, this%num_clusters(i)
          this%dpar(l) = rhocrit(this, this%z(k))
          k = k + 1
          l = l + 1
       end do
       k1 = k
    end do

    ! optionally simulate data
    if (fgas_simulate) then
       junk = 0.0
       call fclwrappersimulate(this%dpar, junk, junk)
       fgas_simulate = .false. ! only simulate once
    end if

    ! calculate likelihood
    fgas_LnLike = -fclwrapperlnp(this%dpar)

    if (Feedback > 2) write(*,*) 'fgas -ln(like) = ', fgas_LnLike

  end function fgas_LnLike
    
end module fgas


function fangulardistance2(z1, z2)
  use Precision
  use fgas, only : like
  implicit none
  real(dl) fangulardistance2
  real(dl), intent(in) :: z1, z2
  fangulardistance2 = like%Calculator%AngularDiameterDistance2(z1, z2)
end function fangulardistance2
