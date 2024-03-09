module list_particle_module_omp

!use kdtree2_precision_module
!use kdtree2_module

use quadtree_module
implicit none


! ==============================================================================
!                              GLOBAL VARIABLES
! ==============================================================================

integer,    parameter :: ndim=2, nspec=2
integer,    parameter :: spkind = kind(1.0), dpkind = kind(1.0d0)
! choose single or double precision, below, or kind 'wp' from the quadtree code
integer,    parameter :: pDki = wp       ! Inherited from prior code - will keep this crappy notation
real(pDki), parameter :: pi = 4.0_pDki * atan(1.0_pDki)
real(pDki), parameter :: eps = 10**(-15)

!!!  Create any derived types here! 

! derived type for particles
type ParticleType
    real(pDki) :: loc(ndim) ! real-valued spatial location vector (x,y, ...)
    real(pDki) :: mass(nspec) ! real-valued mass vector (A,B,...)
    integer    :: box ! which box the particle is in 
    logical    :: active ! whether a given particle is currently active
!    logical    :: ghost(8) ! indicates whether a particle is a ghost particle in a corresponding direction.
!    logical    :: jumped(8) ! indicates whether the particle jumped subdomains and which way it went.

                            ! The logicals for both ghost(8) and jumped(8) corresopnd to a direction based on the 
                            ! value of the 1 in the vector. The positions of the 1 correspond to the following send
                            ! directions: L, R, D, U, DL, DR, UL, UR
end type ParticleType

type SparseMatType
    integer    :: row
    integer    :: col
    real       :: val
end type SparseMatType

! holds the results of the kD tree fixed radius search
type kDRS_ResultType
    integer                 :: num    ! number of nearest neighbors (size of idx and rad)
    integer,    allocatable :: idx(:) ! indices of nearest neighbors
    real(pDki), allocatable :: rad(:) ! distances to nearest neighbors
end type kDRS_ResultType

!type(kdtree2_result), allocatable  :: results(:)

!!$OMP THREADPRIVATE (results)

! all the subroutines are below
contains


subroutine mt_all_dts(A,constnx,dt,dtot,kappa,xmin,xmax,np,ttot,concfactor,beta)

character(80) :: blah
integer :: ntsteps,Stot,constnx,np,ncount, bucket
integer :: nxmt(ndim),ijks(3)
integer :: i,j,k,l,ll,m,n
integer :: numlocal,numghost,numnearest,rownow,partitionsize, Nclose

type(ParticleType), intent(inout)       :: A(:) ! particle array
! type(ParticleType), allocatable       :: X(:) ! subset of particle array
type(kDRS_ResultType), allocatable      :: neighbors(:) ! holds the results of the kD tree fixed radius search
type(SparseMatType), allocatable        :: Emat(:) ! sparse distance matrix
!type(kdtree2), pointer                  :: tree ! the KD tree
type(qtree)                             :: tree ! the quad-tree


real(pDki), allocatable                 :: locations(:,:) 
integer,allocatable                     :: start(:), finish(:)
real(wp)                                :: xmin(ndim), xmax(ndim), dxmt(ndim)
real(wp)                                :: sep(ndim),pad(ndim)
real(wp)                                :: ttot,dt,dtot,kappa,dmt,drw,sigma,concfactor
real(wp)                                :: beta,cutdist2,cutdist,denom,ds,total_volume
integer                                 :: TID, OMP_GET_THREAD_NUM, NTHREADS, OMP_GET_NUM_THREADS

!  Making allocatable arrays decreases demands on the stack, I think? Maybe allocate within parallel do loop!

real,allocatable                        :: rvec(:)
real(wp),allocatable                    :: tmpmass(:,:),tmpmass2(:,:),ith_part_mass(:)
real,allocatable                        :: xlo(:,:),xhi(:,:),ones(:)
integer,allocatable                     :: idx(:),idxghost(:),idxlocal(:),idxnearest(:),idxActive(:),startpart(:),endpart(:)

allocate(locations(np,ndim))
allocate(tmpmass(np,nspec),tmpmass2(np,nspec))

!allocate(mask(np),masklocal(np),maskghost(np))
!logical,allocatable :: mask(:),masklocal(:),maskghost(:),masknearest(:)


ntsteps=ceiling(ttot/dt)
dt=ttot/real(ntsteps)
drw=(1.0-kappa)*dtot
dmt=kappa*dtot          ! Diffusion stuff
idx=(/(i,i=1,np)/)      ! for packing into indexing 
 
!!! Do the mass-transfer
!!! First chop up space: ideally each chunk is bigger than 3 sigma
  sigma = sqrt(4.0*dmt*dt)    ! Char. diffusion dist.
  Stot=1
!  do i=1,ndim
!    pad(i)=3.0*sigma 
 
!    nxmt(i)=constnx                     ! same # of subdomains in each direction
!    dxmt(i)=(xmax(i)-xmin(i))/nxmt(i)   ! for constant dx
!    Stot=Stot*nxmt(i)
!  enddo
!  print*,'Ratio of pad to dx (assuming constants) = ',pad(1)/dxmt(1)

  total_volume=1.0 
  do i=1,ndim
     total_volume=total_volume*(xmax(i)-xmin(i))
  enddo
  ds=total_volume/real(np)


!  allocate(xlo(ndim,maxval(nxmt)),xhi(ndim,maxval(nxmt)))  ! This is public
  
!  do i=1,ndim
!    do m=1,nxmt(i)    ! This can allow variable dxs later (i.e. kdtree)
!        xlo(i,m)=xmin(i)+(m-1)*dxmt(i);
!        xhi(i,m)=xmin(i)+m*dxmt(i);
!    enddo
!  enddo
 
 bucket = max(10,int(real(np)/10000.0))
 bucket = int(0.3*(real(np)**0.5))

!$OMP PARALLEL
   NTHREADS=OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
 
 print*, 'Qtree bucket size = ', bucket, ';    Number of threads: ', NTHREADS
 print*, ' '


 do n=1,ntsteps    !  This is the time loop !!!!!!!!!!!!!!!!

  print*,'number of subdomains is ',Stot,'; tstep ',n,' of ',ntsteps,'; mass = ',concfactor*sum(A%mass(1))
  
  !  Build a kd-tree of all particle locations
  do l=1,ndim
 !   locs1 = real(p%loc(1), kdkind)
 !   locs2 = real(p%loc(2), kdkind)
     locations(1:np,l) = (A%loc(l))
 !   locs(2,:) = locs2
  enddo

   do i=1,nspec
        tmpmass(:,i)=A(:)%mass(i)     ! Pull masses from A; push to tmpmass!  Public to all
    enddo
    tmpmass2=tmpmass

!!!!!!!!  Divy up the particles among the Stot partitions
!    allocate(startpart(Stot),endpart(Stot))
!    partitionsize=floor(real(np)/real(Stot))
!    startpart(1)=1
!    endpart(1)=partitionsize
!    do m=2,Stot
!        startpart(m)=1+endpart(m-1)
!        endpart(m)=startpart(m)+partitionsize-1
!    enddo
!    endpart(Stot)=min(np,endpart(Stot))         ! Final partition might not be full
!!!!!!!!!!!!!!!!

!     Make a quad-tree of the particle locations !!!!!!
!  Arguments: Center_x, center_y, width, height of bounds, bucket size

  call tree%init((xmax(1)+xmin(1))/2.0,(xmax(2)+xmin(2))/2.0,xmax(1)-xmin(1),xmax(2)-xmin(2),bucket)    

  call tree%populate(locations)

  ! print*, 'made it 2'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    cutdist2=(3.0*sqrt(4.0*dmt*dt))**2.0   !  Seems pretty big, like too big.
    cutdist=3.0*sqrt(4.0*dmt*dt)   !  Seems pretty big, like too big.
    denom=-4.0*dmt*dt/beta

!$OMP PARALLEL PRIVATE(i,numlocal,Nclose,neighbors,Emat,TID,ith_part_mass)
 
!$OMP DO SCHEDULE(DYNAMIC,1) 
    
    do i=1,np     ! Each particle gets a calculation

        numlocal = 1                    ! Number of particles assigned to this thread
        TID = OMP_GET_THREAD_NUM()      !    Obtain thread number for debugging
    
        allocate (neighbors(numlocal),ith_part_mass(nspec)) 
  
        call FRsearch(tree, cutdist, neighbors, TID, i, i, 1, np, locations)

        Nclose = sum(neighbors%num) 

        if (Nclose>0) then                    !!! At least 1 real particle to change?  

            ! build the sparse mass transfer matrix (all in one subroutine!!)
 
            call build_MT_matrix_sparse(i, i, Emat, Nclose, TID, neighbors, numlocal, denom, beta, ds)

     ! print *, i, TID, Emat

            ! conduct the mass transfers with sparse matrix-vector multiplication
            ! only change the temporary masses

            ! print *,'Here 1',sum(tmpmass2(:,1))

            call SP_matVecMult(Emat, tmpmass, Nclose,ith_part_mass)

             ! print *,'Here 2',sum(tmpmass2(:,1))

           tmpmass2(i,:)=ith_part_mass(:)          
           deallocate(Emat,ith_part_mass)
 
        endif

        deallocate(neighbors)
 

  enddo !!!!! cell (subdivision) loop ******************************************

 !$OMP END DO NOWAIT

 !$OMP END PARALLEL    

  ! call kdtree2_destroy(tree)


 
    do i=1,nspec
        A(:)%mass(i)=tmpmass2(:,i)    ! Push mass back up to A matrix at end of cell loop 
 !    print *, tmpmass2(:,i)

    enddo

  
!!!!  Do random walks here ...
  if(kappa<1) then
  allocate(rvec(np))
  	  do i=1,ndim
        call random_number(rvec);
        A(:)%loc(i)=A(:)%loc(i)+sqrt(2*DRW*dt)*rvec;
        where (A(:)%loc(i)<xmin(i))
           A(:)%loc(i)=2*xmin(i)+A(:)%loc(i)
        end where
        where (A(:)%loc(i)>xmax(i)) 
           A(:)%loc(i)=2*xmax(i)-A(:)%loc(i)
        end where
      enddo
    deallocate(rvec)
  endif

 enddo  ! Timesteps loop
 end subroutine mt_all_dts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine builds the (squared) distance matrix, then the mass-transfer matrix all in one

subroutine build_MT_matrix_sparse(startpart, endpart, Emat, Nclose, TID, neighbors, numlocal, denom, beta, ds)

    integer,                          intent(in   ) :: startpart, endpart, TID, numlocal ! indices of local particles
 !   real(pDki),                       intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
 !   type(ParticleType),               intent(in   ) :: p(:) ! particle array
    type(SparseMatType), allocatable, intent(  out) :: Emat(:) ! sparse distance matrix
    type(kDRS_ResultType),             intent(in) :: neighbors(numlocal) ! holds the results of the kD tree fixed radius search
 !   integer, allocatable,             start(:), finish(:) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    integer,                          intent(  out) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)
    real(pDki),          intent(in   ) :: beta, denom, ds  ! SPH bandwidth, denominator of the exponential in the co-location PDF
 
    ! ========================== LOCAL VARIABLES ===============================

    integer, allocatable               :: start(:), finish(:) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat  
    integer                            :: i,j                   ! loop iterator
    integer                            :: tmpstart ! temporary variable to handle the n+1^th calculation of start in the loop
    integer                            :: rownow   ! holds the actual row of the ith particle (same as the full tree neighbors!)
    real(pDki)                         :: prefactor 
    real(pDki), allocatable             :: deficit(:)
 
    ! allocate Distmat to have length = total number of neighbors found by the kD search
    Nclose = sum(neighbors%num)


 !   nlocal=endpart-startpart+1

 !   allocate(start(nlocal),finish(nlocal))
    allocate(start(numlocal),finish(numlocal))
    allocate(Emat(Nclose))
    allocate(deficit(numlocal))

    ! fill in Distmat
    prefactor=ds/((-denom * pi)**(real(ndim)/2.0))
    tmpstart = 1
    do i = 1, numlocal
        rownow=startpart+i-1
        start(i) = tmpstart
        finish(i) = start(i) - 1 + neighbors(i)%num
        Emat(start(i) : finish(i))%col = neighbors(i)%idx
        Emat(start(i) : finish(i))%row = rownow
        Emat(start(i) : finish(i))%val = prefactor*exp( (neighbors(i)%rad)/denom )
        tmpstart = finish(i) + 1
        deficit(i)=1.0-sum(Emat(start(i) : finish(i))%val) 
    enddo

    do i = 1,numlocal
        do j = start(i), finish(i)
            if (Emat(j)%row == Emat(j)%col) then
              Emat(j)%val = Emat(j)%val + deficit(i)  ! Add deficit to main diag.
              !  rowsum is now 1.0, but some of these main diagonals may be less than ZERO!!!   
            endif
        enddo
    enddo

    Emat%val = beta*Emat%val

    do i = 1,numlocal
        do j = start(i), finish(i)
            if (Emat(j)%row == Emat(j)%col) then
                Emat(j)%val = Emat(j)%val + 1.0_pDki - beta
            endif
        enddo
    enddo
   !    print *,startpart, Emat(start(1):finish(1))%val, sum(Emat(start(1) : finish(1))%val)

    deallocate(start,finish)


end subroutine build_MT_matrix_sparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solves y = A * x, when A is in sparse coordinate format
! y is [n x nspec], A is [n x n], x is [n x nspec]
! source: https://www.it.uu.se/education/phd_studies/phd_courses/pasc/lecture-1


subroutine SP_matVecMult(A, x, Nclose, y)
    type(SparseMatType), intent(in   ) :: A(:)
    real(pDki),          intent(in   ) :: x(:,:)
    integer,             intent(in   ) :: Nclose ! number of entries in A (= Nclose)
 !   real(pDki),          intent(out  ) :: y(:,:)
    real(pDki),          intent(out  ) :: y(:)
    integer                            :: i, startpart, endpart

    y(1:nspec) = 0.0_pDki
!   y(A(1)%row,1:nspec) = 0.0_pDki

    do i = 1,Nclose
        y(1:nspec) = y(1:nspec) + A(i)%val * x(A(i)%col,1:nspec)
!        y(A(i)%row,1:nspec) = y(A(i)%row,1:nspec) + A(i)%val * x(A(i)%col,1:nspec)
!        print*,i,A(i)%row, A(i)%col, A(i)%val, y(A(i)%row,1:nspec),x(A(i)%col,1:nspec)
    enddo
!print*, 'Finished matrix multiply'

end subroutine SP_matVecMult
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! this searches an already built quad-tree  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FRsearch(tree, cutdist, neighbors,TID,startpart,endpart,nlocal,np,points)
    integer,                intent(in   ) :: TID,startpart,endpart,nlocal,np ! number of particles
    type(qtree),            intent(in   ) :: tree
 !   type(kdtree2), pointer, intent(in   ) :: tree ! the KD tree
    real(wp),               intent(in   ) :: cutdist !  search radius
    type(kDRS_ResultType),  intent(out  ) :: neighbors(nlocal) ! holds the results of the kD tree fixed radius search
    real(wp),               intent(in   ) :: points(np,ndim)
    ! ========================== LOCAL VARIABLES ===============================
    integer                               :: i, particle_now, numnear, nalloc
    real(wp)                              :: locnow(2)
    integer, allocatable                  :: idxs(:)
 
    ! allocate results as big as it could possibly be
    ! theres probably a more memory-efficient way to do this

    ! loop over all particles in this partition
    do i = 1,nlocal
    
        particle_now=i+startpart-1
        ! the type of search used here finds all the particles within
        ! squared distance r2 from the i^th particle in the list
        ! the hard-coded 0 is the 'correlation time' of the search
        ! as far as i can tell, this means a particle, itself, is included the
        ! idx list, while 1 would leave the particle, itself, out of the list

 
        locnow(1)=points(particle_now,1)
        locnow(2)=points(particle_now,2)

        idxs = tree%query_radius(locnow,cutdist)
        numnear = size(idxs)
     
 
    ! print *,'thread ',TID,' I am here 1, particle # ',particle_now,' numnear: ', numnear


        ! allocate these based on how many nearby particles were found
        allocate (neighbors(i)%idx(numnear), neighbors(i)%rad(numnear))
     !print*, 'I am here 1',i,particle_now,neighbors(i)%idx
  
        ! put the results in the neighbors array
        neighbors(i)%num = numnear
        neighbors(i)%idx = idxs(:)
        neighbors(i)%rad = (points(idxs,1)-locnow(1))**2.0_wp + (points(idxs,2)-locnow(2))**2.0_wp 
 
     !print*, 'I am here 2',i,particle_now,neighbors(i)%idx
     

    enddo
    
end subroutine FRsearch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!subroutine find_particles(numlocal,numghost,idxlocal,idxghost,A,m,xlo,xhi,np,nxmt,pad,TID)
!    integer :: nxmt(ndim),ijks(3), TID
!    integer numlocal, numghost
!    real,allocatable :: xlo(:,:),xhi(:,:)
!    integer,allocatable :: idx(:),idxghost(:),idxlocal(:)
!    logical,allocatable :: masklocal(:),maskghost(:)
!    integer :: i,np,m
!    real :: pad(ndim)
!    type(ParticleType), intent(in)   :: A(:) ! particle array

!    allocate(idx(np),idxghost(np),idxlocal(np))
!    allocate(masklocal(np),maskghost(np))
!    allocate(xlo(ndim,maxval(nxmt)),xhi(ndim,maxval(nxmt)))  

!    idx=(/(i,i=1,np)/)      ! for packing into indexing 
!    maskghost=.false.
!    masklocal=.false.
 
! !  figure out the i,j,k indices    
!      if(ndim==1) then
!        ijks(1)=m
!      elseif(ndim==2) then
!        ijks(1)=int(0.5+amod(-0.00001+m,float(nxmt(1)) ) )
!        ijks(2)=int(0.99999+float(m)/float(nxmt(1)))
!      else 
!        ijks(1)=int(0.5+amod(-0.1+(-0.1+ float(mod(m,nxmt(1)*nxmt(2))) ),float(nxmt(1))))
!        ijks(2)=int(0.99999 +  (amod(-0.000001+m,float(nxmt(1)*nxmt(2))) / nxmt(1) ))
!        ijks(3)=int(0.99999 + float(m)/(float(nxmt(1)*nxmt(2))))
!      endif


!  ! figure out which particles are ghost and which are local and true
!      if(ndim==1) then
!        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)))  maskghost=.true.
!        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) )  masklocal=.true.       
!      elseif(ndim==2) then
!        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)) .and.  &
!               A(1:np)%loc(2)>(xlo(2,ijks(2))-pad(2)) .and. A(1:np)%loc(2)<(xhi(2,ijks(2))+pad(2))  ) maskghost=.true.
!        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) .and.  &
!               A(1:np)%loc(2)>xlo(2,ijks(2)) .and. A(1:np)%loc(2)<xhi(2,ijks(2)) )  masklocal=.true.       
!      elseif(ndim==3) then
!        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)) .and.  &
!               A(1:np)%loc(2)>(xlo(2,ijks(2))-pad(2)) .and. A(1:np)%loc(2)<(xhi(2,ijks(2))+pad(2)) .and.  &
!               A(1:np)%loc(3)>(xlo(3,ijks(3))-pad(3)) .and. A(1:np)%loc(3)<(xhi(3,ijks(3))+pad(3)) ) maskghost=.true.
!        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) .and.  &
!               A(1:np)%loc(2)>xlo(2,ijks(2)) .and. A(1:np)%loc(2)<xhi(2,ijks(2)) .and.  &
!               A(1:np)%loc(3)>xlo(3,ijks(3)) .and. A(1:np)%loc(3)<xhi(3,ijks(3)) )  masklocal=.true.       
!      endif 
!
!      numlocal=count(masklocal)
!      numghost=count(maskghost)
!      
!      idxghost=pack(idx,maskghost)
!      idxlocal=pack(idx,masklocal)     


! end subroutine find_particles      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


end module list_particle_module_omp

