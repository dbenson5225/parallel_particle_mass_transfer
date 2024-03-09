module quadtree_module

! Adapted from the tutorial at https://scipython.com/blog/quadtrees-2-implementation-in-python/

  implicit none

  integer, parameter :: wp = kind(1.0d0)
 ! integer, parameter :: wp = kind(1.0)

  type :: qtnode
    real(wp) :: c(2)   ! Center of a rectangle (first node is total bounding box)
    real(wp) :: w, h   ! Width and height of a rectangle 
    real(wp) :: bnd(4)
    integer :: max_points = 0
    integer :: depth = 0
    integer :: n = 0
    integer, allocatable :: points(:)
    type(qtnode), pointer :: children(:) => null()
    logical :: divided = .false. ! could also be called has_children
  end type

  type :: qtree
    type(qtnode) :: root
    integer :: n = 0
    real(wp), pointer :: points(:,:) => null()
  contains
    procedure :: init => qtree_init
    procedure :: populate => qtree_populate
    procedure :: output_gnuplot => qtree_output_gnuplot
    procedure :: size => qtree_size
    procedure :: query => qtree_query
    procedure :: query_radius => qtree_query_radius
    procedure :: kill => qtree_destroy
  end type

  integer, parameter :: DEFAULT_MAX_POINTS = 5000
 
  interface size
    module procedure qtree_size
  end interface

  interface operator(.contains.)
    module procedure :: contains
  end interface

  interface operator(.intersects.)
    module procedure :: intersects
  end interface

  ! sequential colormap
  integer, parameter :: colors(0:6) = [int(z'eff3ff'),&
                                       int(z'c6dbef'),&
                                       int(z'9ecae1'),&
                                       int(z'6baed6'),&
                                       int(z'4292c6'),&
                                       int(z'2171b5'),&
                                       int(z'084594')]

  ! qualitative color map
  ! integer, parameter :: colors(0:6) = [int(z'7fc97f'),&
  !                                      int(z'beaed4'),&
  !                                      int(z'fdc086'),&
  !                                      int(z'ffff99'),&
  !                                      int(z'386cb0'),&
  !                                      int(z'f0027f'),&
  !                                      int(z'bf5b17')]


contains

  subroutine qtree_init(self,cx,cy,w,h,max_points)
    class(qtree), intent(out) :: self
    real(wp), intent(in) :: cx, cy, w, h
    integer, intent(in), optional :: max_points

    call qtnode_init(self%root,cx,cy,w,h,max_points)
  end subroutine

  subroutine qtree_populate(self,points)
    class(qtree), intent(inout) :: self
    real(wp), intent(in), target :: points(:,:)

    real(wp) :: p(2)
    integer :: n, i
    ! assert dimensions ...

    self%points => points
    n = 0
    do i = 1, size(self%points,dim=1)
      p = self%points(i,:)
      if (qtnode_insert_point(self%root,p,i)) then
        n = n + 1
      end if
    end do
    self%n = n
  end subroutine

  subroutine qtree_output_gnuplot(self,unit,color)
    use iso_fortran_env, only: output_unit
    class(qtree), intent(in) :: self
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: color
    integer :: unit_

    unit_ = output_unit
    if (present(unit)) unit_ = unit

    call qtnode_print_boundaries(self%root,unit_,color)

  end subroutine

subroutine qtree_output_matlabplot(self,unit,color)
    use iso_fortran_env, only: output_unit
    class(qtree), intent(in) :: self
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: color
    integer :: unit_

    unit_ = output_unit
    if (present(unit)) unit_ = unit

    call qtnode_print_boundaries(self%root,unit_,color)

  end subroutine



  recursive subroutine qtnode_print_points(self,points,unit,color)
    class(qtnode), intent(in) :: self
    real(wp), intent(in) :: points(:,:)
    integer, intent(in) :: unit
    logical, intent(in), optional :: color

    integer :: i

    do i = 1, self%n
      write(unit,'(G0.6,1X,G0.6,1X,I0)') points(self%points(i),1:2), colors(self%depth)
    end do

    if (self%divided) then
      do i = 1, 4
        call qtnode_print_points(self%children(i),points,unit,color)
      end do
    end if      
  end subroutine

  recursive subroutine qtnode_print_boundaries(self,unit,color)
    class(qtnode), intent(in) :: self
    integer, intent(in) :: unit
    logical, intent(in), optional :: color

    logical :: color_
    character(200) :: fmt
    character(9) :: hex_string
    integer :: i
    character, parameter :: quote = achar(39)

    color_ = .false.
    !if (present(color)) color_ = color

    if (color_) then
      fmt = "('set object rectangle center ',(G0.6),',',(G0.6),' size ',G0.6,',',G0.6,' fs empty border lc rgb ',A, ' lw 2')"
      if (self%depth < 7) then
        write(hex_string,'(A1,"#",Z6,A1)') quote,colors(self%depth),quote
      else
        write(hex_string,'(A1,"#",Z6,A1)') quote,colors(6),quote
      end if
      write(unit,trim(fmt)) self%c(1), self%c(2), self%w, self%h, hex_string
    else
  !     fmt = "('set object rectangle center ',(G0.6),',',(G0.6),' size ',G0.6,',',G0.6,' fs empty border lc 2')"
  !     write(unit,trim(fmt)) self%bnd(4), self%bnd(2), self%bnd(3), self%bnd(1)
      fmt = "((G0.6),',',(G0.6),','(G0.6),',',(G0.6))"
      write(unit,trim(fmt)) self%c(1)-self%w/2.0_wp, self%c(2)-self%h/2.0_wp, self%w, self%h
    end if

    if (self%divided) then
      do i = 1, 4
        call qtnode_print_boundaries(self%children(i),unit,color_)
      end do
    end if
  end subroutine

  pure function boundaries(cx,cy,w,h) result(bnd)
    real(wp), intent(in) :: cx, cy, w, h
    real(wp) :: bnd(4)
    ! north
    bnd(1) = cy + h/2.0_wp
    ! south
    bnd(2) = cy - h/2.0_wp
    ! east
    bnd(3) = cx + w/2.0_wp 
    ! west
    bnd(4) = cx - w/2.0_wp
  end function


  subroutine qtnode_init(self,cx,cy,w,h,max_points,depth)
    class(qtnode), intent(out) :: self
    real(wp), intent(in) :: cx, cy, w, h
    integer, intent(in), optional :: max_points
    integer, intent(in), optional :: depth

    self%c(1) = cx
    self%c(2) = cy
    self%w = w
    self%h = h
    self%bnd = boundaries(cx,cy,w,h)

    self%max_points = DEFAULT_MAX_POINTS
    if (present(max_points)) self%max_points = max_points

    ! print*, max_points

    self%depth = 0
    if (present(depth)) self%depth = depth

    allocate(self%points(self%max_points))
    self%points = -1
  end subroutine

  subroutine qtnode_divide(self)
    class(qtnode), intent(inout) :: self

    real(wp) :: cx, cy, w, h
    
    cx = self%c(1)
    cy = self%c(2)
    
    w = self%w/2.0_wp
    h = self%h/2.0_wp

    allocate(self%children(4))

    ! north west
    call qtnode_init(self%children(1),cx - w/2.0_wp,cy + h/2.0_wp,w,h, &
      self%max_points,self%depth + 1)
    ! north east
    call qtnode_init(self%children(2),cx + w/2.0_wp,cy + h/2.0_wp,w,h, &
      self%max_points,self%depth + 1)
    ! south east
    call qtnode_init(self%children(3),cx + w/2.0_wp,cy - h/2.0_wp,w,h, &
      self%max_points,self%depth + 1)
    ! south west
    call qtnode_init(self%children(4),cx - w/2.0_wp,cy - h/2.0_wp,w,h, &
      self%max_points,self%depth + 1)

    self%divided = .true.
  end subroutine

  recursive logical function qtnode_insert_point(self,p,idx) result(success)
    class(qtnode), intent(inout) :: self
    real(wp), intent(in) :: p(2)
    integer, intent(in) :: idx
    logical :: nw, ne, se, sw

    if (.not. (self%bnd .contains. p)) then
      success = .false.
      return
    end if

    if (self%n < self%max_points) then
      self%n = self%n + 1
      self%points(self%n) = idx
      success = .true.
      return
    end if

    if (.not. self%divided) then
      call qtnode_divide(self)
    end if

    nw = qtnode_insert_point(self%children(1),p,idx) 
    ne = qtnode_insert_point(self%children(2),p,idx) 
    se = qtnode_insert_point(self%children(3),p,idx) 
    sw = qtnode_insert_point(self%children(4),p,idx) 

    success = (nw .or. ne) .or. (se .or. sw)
  end function

  pure recursive integer function qtnode_size(self) result(npoints)
    class(qtnode), intent(in) :: self
    integer :: i
    npoints = self%n
    if (self%divided) then
      do i = 1, 4
        npoints = npoints + qtnode_size(self%children(i))
      end do
    end if
  end function

  pure integer function qtree_size(self) result(npoints)
    class(qtree), intent(in) :: self
    npoints = qtnode_size(self%root)
  end function

  pure function distance(p,q) result(d)
    real(wp), intent(in) :: p(2)
    real(wp), intent(in) :: q(2)
    real(wp) :: d
    d = hypot(p(1)-q(1),p(2)-q(2))
  end function

  pure logical function contains(bnd,p)
    real(wp), intent(in) :: bnd(4)
    real(wp), intent(in) :: p(2)

    real(wp) :: x, y, ne, se, ee, we
    
    x = p(1)
    y = p(2)
    ne = bnd(1)
    se = bnd(2)
    ee = bnd(3)
    we = bnd(4)

    contains = (x >= we .and. x < ee .and. &
                y >= se .and. y < ne)
  end function

  logical function intersects(bnd1,bnd2)
    real(wp), intent(in) :: bnd1(4)
    real(wp), intent(in) :: bnd2(4)
    associate(ne1 => bnd1(1), se1 => bnd1(2), ee1 => bnd1(3), we1 => bnd1(4), &
              ne2 => bnd2(1), se2 => bnd2(2), ee2 => bnd2(3), we2 => bnd2(4))
      intersects = .not. ((we2 > ee1 .or. ee2 < we1) .or. (se2 > ne1 .or. ne2 < se1))
    end associate
    ! print *, "hello", intersects
  end function


  recursive integer function qtnode_count_blocks(self,bnd) result(nblocks)
    class(qtnode), intent(in) :: self
    real(wp), intent(in) :: bnd(4)
    integer :: i
    
    nblocks = 0

    if (self%bnd .intersects. bnd) then
      nblocks = nblocks + 1
    end if

    if (self%divided) then
      do i = 1, 4
        nblocks = nblocks + qtnode_count_blocks(self%children(i),bnd)
      end do
    end if
  end function


  recursive subroutine qtnode_query_rect(self,points,bnd,nfound,idxs)
    class(qtnode), intent(in) :: self
    real(wp), intent(in) :: bnd(4)
    real(wp), intent(in) :: points(:,:)
    integer, intent(inout) :: nfound
    integer, intent(inout) :: idxs(:)

    real(wp) :: p(2)
    integer :: i

    if (.not. (self%bnd .intersects. bnd)) then
      return
    end if

    do i = 1, self%n
      p = points(self%points(i),1:2)
      if (bnd .contains. p) then
        nfound = nfound + 1
        idxs(nfound) = self%points(i)
      end if
    end do

    if (self%divided) then
      do i = 1, 4
        call qtnode_query_rect(self%children(i),points,bnd,nfound,idxs)
      end do
    end if
  end subroutine

  function qtree_query(self,boundary) result(idxs)
    class(qtree), intent(in) :: self
    real(wp), intent(in) :: boundary(4)
    integer, allocatable :: idxs(:)

    integer :: nblocks, nfound
    integer, allocatable :: temp_idxs(:)

    nblocks = qtnode_count_blocks(self%root,boundary)

    allocate(temp_idxs(nblocks*self%root%max_points)) ! worst case scenario if all blocks are full

    nfound = 0
    call qtnode_query_rect(self%root,self%points,boundary,nfound,temp_idxs)

    allocate(idxs,source=temp_idxs(1:nfound))

  end function


  recursive subroutine qtnode_query_circle(self,points,bnd,center,radius,nfound,idxs)
    class(qtnode), intent(in) :: self
    real(wp), intent(in) :: points(:,:)
    real(wp), intent(in) :: bnd(4)
    real(wp), intent(in) :: center(2)
    real(wp), intent(in) :: radius
    integer, intent(inout) :: nfound
    integer, intent(inout) :: idxs(:)

    integer :: i
    real(wp) :: p(2)

    if (.not. (self%bnd .intersects. bnd)) then
      return
    end if

    do i = 1, self%n
      p = points(self%points(i),1:2)
      if ((bnd .contains. p) .and. distance(p,center) <= radius) then
        nfound = nfound + 1
        idxs(nfound) = self%points(i)
      end if
    end do

    if (self%divided) then
      do i = 1, 4
        call qtnode_query_circle(self%children(i),points,bnd,center,radius,nfound,idxs)
      end do
    end if

  end subroutine

  function qtree_query_radius(self,center,radius) result(idxs)
    class(qtree), intent(in) :: self
    real(wp), intent(in) :: center(2)
    real(wp), intent(in) :: radius
    integer,  allocatable :: idxs(:)

    integer :: nblocks, nfound
    real(wp) :: boundary(4)
    integer, allocatable :: temp_idxs(:)

    boundary = boundaries(center(1),center(2),2*radius,2*radius)
    nblocks = qtnode_count_blocks(self%root,boundary)
    allocate(temp_idxs(nblocks*self%root%max_points)) ! worst case scenario if all blocks are full

    nfound = 0
    call qtnode_query_circle(self%root,self%points,boundary,center,radius,nfound,temp_idxs)

    allocate(idxs,source=temp_idxs(1:nfound))
  end function

  subroutine qtree_destroy(self)
    class(qtree), intent(inout) :: self
    if (self%n > 0) then
      ! if the qtree is not empty
      call qtree_destroy_node(self%root)
      self%root%max_points = 0
      self%root%depth = 0
      self%root%n = 0
    end if
    nullify(self%points)

  end subroutine

!  contains
    recursive subroutine qtree_destroy_node(self)
      class(qtnode), intent(inout) :: self
      integer :: i
      if (self%divided) then
        do i = 1, 4
          call qtree_destroy_node(self%children(i))
        end do
        deallocate(self%children)
        self%divided = .false.
      end if
      deallocate(self%points)
    end subroutine
!  end subroutine

end module

