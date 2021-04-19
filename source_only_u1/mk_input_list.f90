program main
  implicit none
  integer :: model   ! 1: square lattice, 2: triangular lattice
  integer :: lx, ly  ! length of lattices
  real(8) :: jxx, jz ! xxz interaction
  integer :: i, j, p, pd, nos, nod, nohz, noxxz
  character(100) :: filexxz
  namelist /list_input/ model, lx, ly, jxx, jz, filexxz
  read(*,list_input)
  write(*,list_input)
  open(10, file=trim(filexxz),status='unknown')
  if(model == 1) then
    call bond_x
    call bond_y
  elseif(model == 2)then
    call bond_x
    call bond_y
    call bond_diag
  end if
  close(10)
  stop
contains
  subroutine bond_x
    p = 1
    do i = 1, lx-1
      do j = 1, ly
        write(10,*) min(p, p+ly), max(p, p+ly), jxx, jz
        p = p + 1
      end do
    end do
    pd = 1
    do j = 1, ly
      write(10,*) min(p, pd), max(p, pd), jxx, jz
      p = p + 1
      pd= pd+ 1
    end do
  end subroutine bond_x
  !
  subroutine bond_y
    p = 1
    do i = 1, lx
      do j = 1, ly-1
        write(10,*) min(p, p+1), max(p, p+1), jxx, jz
        p = p + 1
      end do
      write(10,*) min(p, p-ly+1), max(p, p-ly+1), jxx, jz
      p = p + 1
    end do    
  end subroutine bond_y
  !
  subroutine bond_diag
    p = 1
    do i = 1, lx-1
      write(10,*) min(p, p+2*ly-1), max(p, p+2*ly-1), jxx, jz
      p = p + 1
      do j = 2, ly
        write(10,*) min(p, p+ly-1), max(p, p+ly-1), jxx, jz
        p = p + 1
      end do
    end do
    pd = -ly+1
    write(10,*) min(p, pd+2*ly-1), max(p, pd+2*ly-1), jxx, jz
    p  = p  + 1
    pd = pd + 1
    do j = 2, ly
      write(10,*) min(p, pd+ly-1), max(p, pd+ly-1), jxx, jz
      p = p + 1
      pd= pd+ 1
    end do
  end subroutine bond_diag
end program main
