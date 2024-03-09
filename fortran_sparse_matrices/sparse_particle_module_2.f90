module sparse_particle_module_2

!use kdtree2_precision_module
!use kdtree2_module
!implicit none
use quadtree_module

! ==============================================================================
!                              GLOBAL VARIABLES
! ==============================================================================

integer,    parameter :: ndim=2, nspec=2
integer,    parameter :: spkind = kind(1.0), dpkind = kind(1.0d0)
! choose single or double precision, below
integer,    parameter :: pDki = wp
real(pDki), parameter :: pi = 4.0_pDki * atan(1.0_pDki)
! integer               :: ndim = 2  ! Number of spatial dimensions
real(pDki), parameter :: eps = 10**(-15)

!!!  Create any derived types here! 

! derived type for particles
type ParticleType
    real(pDki) :: loc(ndim) ! real-valued spatial location vector (x,y, ...)
    real(pDki) :: mass(nspec) ! real-valued mass vector (A,B,...)
    integer    :: box ! which box the particle is in 
    logical    :: active ! whether a given particle is currently active
end type ParticleType

type SparseMatType
    integer    :: row
    integer    :: col
    real       :: val
end type SparseMatType

! holds the results of the tree fixed radius search
type kDRS_ResultType
    integer                 :: num    ! number of nearest neighbors (size of idx and rad)
    integer,    allocatable :: idx(:) ! indices of nearest neighbors
    real(pDki), allocatable :: rad(:) ! distances to nearest neighbors
end type kDRS_ResultType

! all the subroutines are below
contains


subroutine mt_all_dts(A,constnx,dt,dtot,kappa,xmin,xmax,np,ttot,concfactor,beta)

character(80) :: blah
integer :: ntsteps,Stot,constnx,np,ncount
integer :: nxmt(ndim),ijks(3)
integer :: i,j,k,l,ll,m,n
integer :: numlocal,numghost,numnearest,rownow

type(ParticleType), intent(inout)   :: A(:) ! particle array
type(ParticleType), allocatable     :: X(:) ! subset of particle array

real(wp)    :: xmin(ndim), xmax(ndim), dxmt(ndim)
real(pDki)  :: sep(ndim),pad(ndim)
real(pDki)  :: ttot,dt,dtot,kappa,dmt,drw,sigma,concfactor,beta,cutdist,denom,ds,volume
integer     :: TID, OMP_GET_THREAD_NUM,NTHREADS,OMP_GET_NUM_THREADS

!  Making allocatable arrays decreases demands on the stack, I think? Maybe allocate within parallel do loop!

real(pDki),allocatable      :: rvec(:)
real(pDki),allocatable      :: massnew(:,:),masstemp(:,:)
real(pDki),allocatable      :: xlo(:,:),xhi(:,:),ones(:)
integer,allocatable :: idx(:),idxghost(:),idxlocal(:),idxnearest(:),idxActive(:);
!logical,allocatable :: mask(:),masklocal(:),maskghost(:),masknearest(:)

allocate(massnew(np,nspec),masstemp(np,nspec))

ntsteps=ceiling(ttot/dt)
dt=ttot/real(ntsteps)
drw=(1.0-kappa)*dtot
dmt=kappa*dtot          ! Diffusion stuff
idx=(/(i,i=1,np)/)      ! for packing into indexing 
 
!!! Do the mass-transfer
!!! First chop up space: ideally each chunk is bigger than 3 sigma
  sigma = sqrt(4.0*dmt*dt)    ! Char. diffusion dist.
  Stot=1
  volume=1.0
  do i=1,ndim
    pad(i)=3.0*sigma 
 
    nxmt(i)=constnx                     ! same # of subdomains in each direction
    dxmt(i)=(xmax(i)-xmin(i))/nxmt(i)   ! for constant dx
    Stot=Stot*nxmt(i)
    volume=volume*(xmax(i)-xmin(i))
  enddo
  ds=volume/real(np)
  print*,'Ratio of pad to dx (assuming constants) = ',pad(1)/dxmt(1)

  allocate(xlo(ndim,maxval(nxmt)),xhi(ndim,maxval(nxmt)))  ! This is public
  
 
  do i=1,ndim
    do m=1,nxmt(i)    ! This can allow variable dxs in x and y later 
        xlo(i,m)=xmin(i)+(m-1)*dxmt(i);
        xhi(i,m)=xmin(i)+m*dxmt(i);
    enddo
  enddo
 
 do n=1,ntsteps    !  This is the time loop !!!!!!!!!!!!!!!!

  print*,'number of subdomains is ',Stot,'; tstep ',n,' of ',ntsteps,'; mass = ',concfactor*sum(A%mass(1))
  
  
  ! Build sparse matrix of separations including ghosts from adjoining cells
  
    do i=1,nspec
        massnew(:,i)=A(:)%mass(i)     ! Pull masses from A; push to massnew!  Public to all
    enddo
    masstemp=massnew      ! This holds temporary masses (incl. ghosts): Public to all

!$OMP PARALLEL PRIVATE(m,X,numlocal,numghost,idxghost,idxlocal,TID,treelocal,diag,ijks)
   NTHREADS=OMP_GET_NUM_THREADS()
  !   print*,'    Number of threads: ', NTHREADS
 
!$OMP DO SCHEDULE(DYNAMIC,1) 
     do m=1,Stot            ! go to every subdomain (cell loop)  ******* Need to OMP parallel this loop! *******

     !     Obtain thread number for debugging
       TID = OMP_GET_THREAD_NUM()
     !  PRINT *, 'Hello from thread = ', TID, '   subdiv. = ',m  ! For debugging

    call find_particles(numlocal,numghost,idxlocal,idxghost,A,m,xlo,xhi,np,nxmt,pad,TID,ijks)

    denom=-4.0*dmt*dt/beta   ! This is for 2-dimensional Gaussian 

    if(numlocal>0) then       !!! At least 1 real particle to change?  

!  Make a subset X of the full particle array for this subdomain
    
       allocate(X(numghost))
       X(1:numghost)=A(idxghost(1:numghost)) 
      
       call massTrans_octree_Mat(numghost, numlocal, xlo, xhi, m, pad, denom, X, beta, TID, ds, ijks)
     
      
       do i=1,nspec
          masstemp(idxghost(1:numghost),i)=X(1:numghost)%mass(i)                ! update new mass on all real and ghosts
          massnew(idxlocal(1:numlocal),i)=masstemp(idxlocal(1:numlocal),i)      ! push only new real particle masses up
       enddo
       deallocate(X)

    endif  !!!!!!! more than 1 particle if statement
 
   deallocate(idxlocal,idxghost)

     enddo !!!!! cell (subdivision) loop ******************************************

 !$OMP END DO NOWAIT

 !$OMP END PARALLEL    

    do i=1,nspec
        A(:)%mass(i)=massnew(:,i)    ! Push mass back up to A matrix at end of cell loop 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine find_particles(numlocal,numghost,idxlocal,idxghost,A,m,xlo,xhi,np,nxmt,pad,TID,ijks)
    integer,                intent(in   ) :: TID, nxmt(ndim)
    real(pDki),allocatable, intent(in   ) :: xlo(:,:),xhi(:,:)
    integer,                intent(  out) :: ijks(3)
    integer,allocatable :: idx(:),idxghost(:),idxlocal(:)
    logical,allocatable :: masklocal(:),maskghost(:)
    integer :: i,np,m
    integer :: numlocal, numghost
 
    real(pDki) :: pad(ndim)
    type(ParticleType), intent(in)   :: A(:) ! particle array

    allocate(idx(np),idxghost(np),idxlocal(np))
    allocate(masklocal(np),maskghost(np))

    idx=(/(i,i=1,np)/)      ! for packing into indexing 
    maskghost=.false.
    masklocal=.false.
 
 !  figure out the i,j,k indices    
      if(ndim==1) then
        ijks(1)=m
      elseif(ndim==2) then
        ijks(1)=int(0.5+amod(-0.00001+m,float(nxmt(1)) ) )
        ijks(2)=int(0.99999+float(m)/float(nxmt(1)))
      else 
        ijks(1)=int(0.5+amod(-0.1+(-0.1+ float(mod(m,nxmt(1)*nxmt(2))) ),float(nxmt(1))))
        ijks(2)=int(0.99999 +  (amod(-0.000001+m,float(nxmt(1)*nxmt(2))) / nxmt(1) ))
        ijks(3)=int(0.99999 + float(m)/(float(nxmt(1)*nxmt(2))))
      endif


  ! figure out which particles are ghost and which are local and true
      if(ndim==1) then
        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)))  maskghost=.true.
        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) )  masklocal=.true.       
      elseif(ndim==2) then
        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)) .and.  &
               A(1:np)%loc(2)>(xlo(2,ijks(2))-pad(2)) .and. A(1:np)%loc(2)<(xhi(2,ijks(2))+pad(2))  ) maskghost=.true.
        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) .and.  &
               A(1:np)%loc(2)>xlo(2,ijks(2)) .and. A(1:np)%loc(2)<xhi(2,ijks(2)) )  masklocal=.true.       
      elseif(ndim==3) then
        where( A(1:np)%loc(1)>(xlo(1,ijks(1))-pad(1)) .and. A(1:np)%loc(1)<(xhi(1,ijks(1))+pad(1)) .and.  &
               A(1:np)%loc(2)>(xlo(2,ijks(2))-pad(2)) .and. A(1:np)%loc(2)<(xhi(2,ijks(2))+pad(2)) .and.  &
               A(1:np)%loc(3)>(xlo(3,ijks(3))-pad(3)) .and. A(1:np)%loc(3)<(xhi(3,ijks(3))+pad(3)) ) maskghost=.true.
        where( A(1:np)%loc(1)>xlo(1,ijks(1)) .and. A(1:np)%loc(1)<xhi(1,ijks(1)) .and.  &
               A(1:np)%loc(2)>xlo(2,ijks(2)) .and. A(1:np)%loc(2)<xhi(2,ijks(2)) .and.  &
               A(1:np)%loc(3)>xlo(3,ijks(3)) .and. A(1:np)%loc(3)<xhi(3,ijks(3)) )  masklocal=.true.       
      endif 

      numlocal=count(masklocal)
      numghost=count(maskghost)
      
      idxghost=pack(idx,maskghost)
      idxlocal=pack(idx,masklocal)     


 end subroutine find_particles      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


! this is the sparse matrix-based mass transfer algorithm (SPaM)
 
subroutine massTrans_octree_Mat(nghost, nlocal, xlo, xhi, m, pad, denom, X, beta, TID, ds, ijks)
    integer,            intent(in   ) :: nghost  ! number of real and ghost particles in *local partition* particle array X
    integer,            intent(in   ) :: nlocal  ! number of real particles in partition
    real(pDki),         intent(in   ) :: xlo(ndim,m), xhi(ndim,m)  ! coordinates of the partitions
    integer,            intent(in   ) :: m  ! number of this partition
    real(pDki),         intent(in   ) :: pad(ndim) ! cutoff distance for the oct-tree fixed-radius search
    real(pDki),         intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    type(ParticleType), intent(inout) :: X(nghost) ! particle array for this partition
    real(pDki),         intent(in   ) :: beta, ds ! bandwidth of the mass-transfer kernel, particle support vol.
    integer,            intent(in   ) :: ijks(3)
    ! ========================== LOCAL VARIABLES ===============================
    type(SparseMatType), allocatable :: Emat(:) ! mass transfer matrix
    ! integer, allocatable             :: start(:), finish(:)
    integer                          :: Nclose ! total number of entries in Emat matrix
    real(pDki), allocatable          :: tmpmass(:,:), tmpmass2(:,:) ! temporary array for holding particle masses
    integer                          :: i, numactive, TID, bucket
 
    real(wp),allocatable             :: locations(:,:)
    real(wp)                         :: xmin, xmax, ymin, ymax  ! For the quad-tree
    type(qtree)                      :: treelocal                             ! the quad-tree
    type(kDRS_ResultType),allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search


     allocate(tmpmass(nghost,nspec), tmpmass2(nghost,nspec)) 
     allocate(locations(nghost,ndim),neighbors(nghost))

    !  figure out the dimensions of this particular partition ...
    xmin = xlo(1,ijks(1))-pad(1)
    xmax = xhi(1,ijks(1))+pad(1)
    ymin = xlo(2,ijks(2))-pad(2)
    ymax = xhi(2,ijks(2))+pad(2)

    ! This is an empirical bucket size based on testing this oct-tree
    bucket = max(10,int(0.3*(real(nghost)**0.5)))

    do l=1,ndim
      locations(1:nghost,l) = X(1:nghost)%loc(l)
    enddo
    
    call maketree(treelocal, nghost, locations, xmin, xmax, ymin, ymax, bucket)

    call FRsearch(treelocal, pad, neighbors, TID, nlocal, nghost, locations)

    call treelocal%kill

    ! Build the matrix of co-location probabilities then mass transfer coefficients
    call build_MT_matrix_sparse(nghost, Nclose, Emat, TID, neighbors, denom, beta, ds)
   
    do i=1,nspec
        tmpmass(1:nghost,i) = X(1:nghost)%mass(i)
    enddo
    tmpmass2=tmpmass

! conduct the mass transfers with sparse matrix-vector multiplication
    call SP_matVecMult(Emat, tmpmass2, Nclose, tmpmass)

    ! change the masses of all particles in this partition (main subroutine will choose real particles)

    do i=1,nspec
      X(1:nghost)%mass(i) = tmpmass(1:nghost,i)
    enddo
    
    deallocate(Emat,tmpmass,tmpmass2,locations)
    deallocate(neighbors)

end subroutine massTrans_octree_Mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine builds the (squared) distance matrix, then the mass-transfer matrix all in one
! Some things are commented out that do the newer list-based mass transfer scheme 

subroutine build_MT_matrix_sparse(nghost, Nclose, Emat, TID, neighbors, denom, beta, ds)

    integer,                          intent(in   ) :: TID, nghost ! indices of ALL local particles
    integer,                          intent(  out) :: Nclose  ! Number of total entries in Emat
    type(SparseMatType), allocatable, intent(  out) :: Emat(:) ! sparse distance matrix
    type(kDRS_ResultType),            intent(inout)    :: neighbors(nghost) ! holds the results of the kD tree fixed radius search
    real(pDki),                       intent(in   ) :: beta, denom, ds  ! SPH bandwidth, denominator of the exponential in the co-location PDF
 
    ! ========================== LOCAL VARIABLES ===============================

   ! integer, allocatable               :: start(:), finish(:) ! indices (in the Emat vectors) for the start and finish of each column of Emat  
    integer                            :: i,j                   ! loop iterator
    integer                            :: tmpstart ! temporary variable to handle the n+1^th calculation of start in the loop
    integer                            :: rownow   ! holds the actual row of the ith particle (same as the full tree neighbors!)
    integer                            :: startnow,finishnow  ! Counters for EMat actual size
    integer, allocatable               :: diag(:) ! linear indices of the diagonal elements of Pmat
    real(pDki)                         :: prefactor 
  !  real(pDki), allocatable            :: deficit(:)
    real(pDki), allocatable            :: rowsum(:), colsum(:) ! arrays holding the row/col sums of Pmat

 
    ! allocate Distmat to have length = total number of neighbors found by the kD search
    Nclose = sum(neighbors%num)
  
 !   nlocal=endpart-startpart+1

 !   allocate(start(nlocal),finish(nlocal))
 !   allocate(start(nlocal),finish(nlocal))
    allocate(Emat(Nclose))
    allocate(diag(nghost))
 !   allocate(deficit(nlocal))
    allocate(rowsum(Nclose),colsum(Nclose))



    ! fill in probability density matrix (start of Emat)
    prefactor=ds/((-denom * pi)**(real(ndim)/2.0))
    rownow = 1
    startnow=1
    do i = 1, nghost
        if (neighbors(i)%num > 0) then
           finishnow = startnow - 1 + neighbors(i)%num
           Emat(startnow : finishnow)%row = i
           Emat(startnow : finishnow)%col = neighbors(i)%idx
           Emat(startnow : finishnow)%val = prefactor*exp( (neighbors(i)%rad)/denom )
           startnow = finishnow + 1
    !       deficit(i)=1.0-sum(Emat(start(i) : finish(i))%val) 
        endif
    enddo

! compute the rowsums of Pmat
    rowsum = 0.0_pDki
    colsum = 0.0_pDki
    do i = 1, Nclose
        rowsum(Emat(i)%row) = rowsum(Emat(i)%row) + Emat(i)%val
        colsum(Emat(i)%col) = colsum(Emat(i)%col) + Emat(i)%val
         if (Emat(i)%row == Emat(i)%col) then
            diag(Emat(i)%row) = i
        endif
    enddo

    do i = 1, Nclose
        Emat(i)%val = Emat(i)%val / ((rowsum(Emat(i)%row) + colsum(Emat(i)%col)) / 2.0_pDki)
    enddo

    rowsum = 0.0_pDki
    do i = 1, Nclose
        rowsum(Emat(i)%row) = rowsum(Emat(i)%row) + Emat(i)%val
    enddo

    ! this is the I - beta * diag(P * 1) step
    rowsum = 1.0_pDki - beta * rowsum

    ! this is the beta * \vec Pmat step
    Emat%val = beta * Emat%val

    ! finally, add them together
    Emat(diag)%val =  Emat(diag)%val + rowsum


!  NEW ALGORITHM
!
!    do i = 1,numlocal
!        do j = start(i), finish(i)
!            if (Emat(j)%row == Emat(j)%col) then
!              Emat(j)%val = Emat(j)%val + deficit(i)  ! Add deficit to main diag.
!              !  rowsum is now 1.0, but some of these main diagonals may be less than ZERO!!!   
!            endif
!        enddo
!    enddo

!   Emat%val = beta*Emat%val

!    do i = 1,numlocal
!        do j = start(i), finish(i)
!            if (Emat(j)%row == Emat(j)%col) then
!                Emat(j)%val = Emat(j)%val + 1.0_pDki - beta
!            endif
!        enddo
!    enddo
   !    print *,startpart, Emat(start(1):finish(1))%val, sum(Emat(start(1) : finish(1))%val)

!    deallocate(start,finish)
   deallocate(diag,rowsum,colsum)
 !  deallocate(neighbors)

end subroutine build_MT_matrix_sparse


! solves y = A * x, when A is in sparse coordinate format
! y is [n x nspec], A is [n x n], x is [n x nspec]
! source: https://www.it.uu.se/education/phd_studies/phd_courses/pasc/lecture-1
  subroutine SP_matVecMult(A, x, n, y)
    type(SparseMatType), intent(in   ) :: A(:)
    real(pDki),          intent(in   ) :: x(:,:)
    integer,             intent(in   ) :: n ! number of entries in A
    real(pDki),          intent(out  ) :: y(:,:)
    integer                            :: i

    y = 0.0_pDki

    do i = 1, n
        y(A(i)%row,1:nspec) = y(A(i)%row,1:nspec) + A(i)%val * x(A(i)%col,1:nspec)
    enddo

end subroutine SP_matVecMult


! this searches an already built quad-tree  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FRsearch(treelocal, pad, neighbors,TID,nlocal,nghost,points)
    integer,                intent(in   ) :: TID,nlocal,nghost ! number of particles
    type(qtree),            intent(in   ) :: treelocal
 !   type(kdtree2), pointer, intent(in   ) :: tree ! the KD tree
    real(wp),               intent(in   ) :: pad(ndim) !  search radius
    type(kDRS_ResultType),  intent(inout) :: neighbors(nghost) ! holds the results of the kD tree fixed radius search
    real(wp),               intent(in   ) :: points(nghost,ndim)
    ! ========================== LOCAL VARIABLES ===============================
    integer                               :: i, particle_now, numnear, nalloc
    real(wp)                              :: locnow(2)
    integer, allocatable                  :: idxs(:)
    real(wp)                              :: cutdist    ! for now an isotropic search
  !  type(kDRS_ResultType),allocatable     :: neighbors(:)

    ! allocate results as big as it could possibly be
    ! theres probably a more memory-efficient way to do this

    allocate(idxs(nghost))
   ! allocate(neighbors(nghost))
    cutdist = pad(1)

   ! cutdist=10.0_wp
   ! locnow(1)=5.0_wp
   ! locnow(2)=5.0_wp
   ! idxs = treelocal%query_radius(locnow,cutdist)
   ! numnear = size(idxs,1)
   ! do i=1,numnear
   !    print*,points(i,1),points(i,2)
   ! enddo


    ! loop over all particles in this partition
    do i = 1,nghost
    
        particle_now=i
 
        locnow(1)=points(particle_now,1)
        locnow(2)=points(particle_now,2)

        idxs = treelocal%query_radius(locnow,cutdist)
        numnear = size(idxs,1)
     
 
  !    print *,'thread ',TID,' I am here 1, particle # ',particle_now,'loc:',locnow,' numnear: ', numnear

        ! allocate these based on how many nearby particles were found
        allocate (neighbors(i)%idx(numnear), neighbors(i)%rad(numnear))
     ! print*, 'I am here 1',TID,i,particle_now,neighbors(i)%idx
  
        ! put the results in the neighbors array
        neighbors(i)%num = numnear
        neighbors(i)%idx = idxs(:)
        neighbors(i)%rad = (points(idxs,1)-locnow(1))**2.0_wp + (points(idxs,2)-locnow(2))**2.0_wp 
 
    enddo

    deallocate(idxs)
 
    
end subroutine FRsearch


! this builds a quadtree
subroutine maketree(treelocal, n, locs, xmin, xmax, ymin, ymax, bucket)
    type(qtree),            intent(  out) :: treelocal         ! the quad-tree

    integer,                intent(in   ) :: n            !  number of particles
    real(pDki),             intent(in   ) :: locs(n,ndim) ! location array for particles, with dimension n x d (number of particles x number of spatial dimensions)
    real(wp), intent(in)                  :: xmin, xmax, ymin, ymax
    integer, intent(in)                   :: bucket
    integer j


!  Make a quad-tree of the particle locations !!!!!!
!  Arguments: Center_x, center_y, width, height of bounds, bucket size

  call treelocal%init((xmax+xmin)/2.0,(ymax+ymin)/2.0,xmax-xmin,ymax-ymin,bucket)    

  call treelocal%populate(locs)

end subroutine maketree



end module sparse_particle_module_2

