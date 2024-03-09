
program main_mt
! use OMP_LIB
! define numerical and physical parameters
use list_particle_module_omp 
use quadtree_module

character(80) :: fnin,fnout,speedfile,speedout,dummy

integer     :: np, nsx, nens            ! # particles, # of discretizations, # in ensemble 
integer     :: ntsteps                  ! # of timesteps
integer     :: i, j, k, n, nn, m        ! counters of various use
integer     :: nsource, constnx         ! # of sources, discretization value
integer, allocatable :: nxmt(:), sx(:)  ! 
!integer, allocatable :: idx(:),idxnear(:),
integer,dimension(8) :: start_time_values,stop_time_values,clock_elapsed

real(wp)      :: ttot, dt, dtot, dmt, drw, kappa, concfactor, beta
real      :: start_time, finish_time, clock_elapsed_sec

type(ParticleType), allocatable :: A(:)

real(wp), allocatable :: xmin(:), xmax(:), randvec(:), pad(:)
real(wp), allocatable :: dxmt(:), speed(:,:), sep(:), analytic(:)

print*,' Enter input filename:'
read(*,'(A)')fnin
print*,' Enter output speeds filename:'
read(*,'(A)')speedfile
print*,' Enter general output filename:'
read(*,'(A)')fnout

open(unit=30,file=trim(fnin))
read(30,*)np, i, i, nsx, dummy                      ! ndim and nspec deprecated here 
read(30,*)ttot, dt, dtot, kappa, beta, dummy
print*, '# particles= ',np,'# species= ',nspec,'# dimensions= ',ndim,'# of discretizations=',nsx
print*,'Total time= ',ttot,'dt= ', dt,'D_total= ', Dtot
print*,'kappa= ', kappa,'beta= ', beta
allocate(xmin(ndim), xmax(ndim), sx(nsx), speed(nsx,2))
do i=1,ndim
   read(30,*)xmin(i),xmax(i), dummy   ! Define domain boundaries
   print*,'for dimension ',i,' min= ',xmin(i),' max= ',xmax(i)
enddo
!do i=1,nsx
  read(30,*)(sx(i),i=1,nsx),dummy
!enddo
print*,'discretizations = ',(sx(i),i=1,nsx)
!!!!!  Output file names
write(speedout,'(A,A1,f4.1,A1,i0)')trim(speedfile),'_',dt,'_',np
open(unit=32,file=trim(speedout))
open(unit=31,file=trim(fnout))

allocate(A(np))  
allocate(randvec(np))

!idx=(/(i,i=1,np)/)      ! for packing into indexing 

drw=(1.0-kappa)*dtot; dmt=kappa*dtot;  ! Diffusion stuff

do nn=1,nsx       ! Go through a bunch of discretizations
   constnx=sx(nn)

do i = 1,np
    do j=1,ndim
       A(i)%loc(j)=0.0
      enddo
    do j=1,nspec
       A(i)%mass(j)=0.0
    enddo
    A(i)%active=.true.
enddo   


! Define ICs %%%%%%%%%%%
! Random positions:

do i=1,ndim
	call random_number(randvec)
   A(1:np)%loc(i)=(xmax(i)-xmin(i))*randvec+xmin(i)
enddo


nsource=1
! locations of A and B ICs;
concfactor=1.0
do i=1,ndim
   concfactor=concfactor*(xmax(i)-xmin(i))
enddo
concfactor=concfactor/real(np)
do i=1,ndim
    A(1:nsource)%loc(i)=0.501*(xmax(i)-xmin(i)) 
enddo

A(1:nsource)%mass(1)=1.0/real(nsource)/concfactor
A(1:nsource)%mass(2)=1.0/real(nsource)/concfactor


!!!!!!!!!  Call the mass transfer subroutine inside a timer !!!!!!!!!!!!!!
call date_and_time(VALUES=start_time_values)
call cpu_time(start_time)

call mt_all_dts(A,constnx,dt,dtot,kappa,xmin,xmax,np,ttot,concfactor,beta)

call date_and_time(VALUES=stop_time_values)
clock_elapsed=stop_time_values-start_time_values
clock_elapsed_sec = 86400.*clock_elapsed(3)+3600.*clock_elapsed(5)+60.*clock_elapsed(6)+ & 
                   1.*clock_elapsed(7)+0.001*clock_elapsed(8)
call cpu_time(finish_time)

speed(nn,1)=float(constnx**ndim)
speed(nn,2)=finish_time-start_time
write(32,*)speed(nn,1),speed(nn,2)
print*,' '
print*,'   CPU time = ',speed(nn,2), 's'
print*,' Clock time = ',clock_elapsed_sec, 's'
print*,' '
enddo   ! End the various discretizations

!  Compare to analytic soln at particle points (only works for debugging right now)
allocate(analytic(np),sep(ndim))
analytic=0.0
do n=1,nsource
do m=1,np
   sep(1:ndim)=A(n)%loc-A(m)%loc
   sepsq=dot_product(sep,sep);
   analytic(m)=analytic(m)+(4*pi*Dtot*ttot)**(-float(ndim)/2.0)*exp(sepsq/(-4*Dtot*ttot));
enddo
enddo

do m=1,np
   if ( abs( A(m)%loc(1)-0.501*(xmax(1)-xmin(1)) ) < 0.01 ) then
      write(31,*)A(m)%loc(2),A(m)%mass(1),analytic(m)
   endif
enddo

RMSE=sqrt( (1./real(np)) *  sum( ( A(1:np)%mass(1)-analytic(1:np) )**2.0))
print*,'RMSE = ', RMSE
end program main_mt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  END OF MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


