Program main
  implicit none
  !for square-cluster (Lx=Ly)
  integer::i,j,k,l,L_N,NOring,count
  integer,parameter::Lx=10
  integer,parameter::Ly=10
  integer,parameter::Lz=10
  real(8),parameter::J1_val=1.00d0
  !real(8),parameter::J2_val=0.00d0!0.00d0!0.89d0
  !real(8),parameter::J_ring=0.00d0!0.5d0!0.50d0
  integer,parameter::num_J1=Lx*Ly*Lz*3!Lx*Lx*2     !2020/6/25-modified
  integer::site_J1(1:2,1:num_J1)        !2020/1/30-add
  real(8)::dum1,dum2                    !2020/1/30-add
  character(4)::size_name
  logical::ref,J2_int



  L_N=Lx*Ly*Lz

  if(L_N>=1000)then
    write(size_name,'(i4)')L_N
  else if(L_N>=100)then
    write(size_name,'(i3)')L_N
  else if(L_N<100)then
    write(size_name,'(i2)')L_N
  end if

  open(10,file='list_xxz_term_'//trim(adjustl(size_name))//'.dat',position='append')

  !J1-val (x-direction)
  if(Lx>1)then
    do i=1, L_N
      if(mod(i,Lx).eq.0)then
        write(10,'(2I5,2es23.15)') i-Lx+1, i, J1_val, J1_val
      else 
        write(10,'(2I5,2es23.15)') i, i+1, J1_val, J1_val
      end if
    end do
  end if

  !J1-val (y-direction)
  if(Ly>1)then
    count=0
    do k=1, Lz
      do j=1, Ly
        do i=1, Lx
          count=count+1
          if(j.ne.Ly)then
            write(10,'(2I5,2es23.15)') count, count+Lx, J1_val, J1_val
          else 
            write(10,'(2I5,2es23.15)') count-Lx*(Ly-1), count, J1_val, J1_val
          end if
        end do
      end do
    end do
  end if

  !J1-val (z-direction)
  if(Lz>1)then
    count=0
    do k=1, Lz
      do j=1, Ly
        do i=1, Lx
          count=count+1
          if(k.ne.Lz)then
            write(10,'(2I5,2es23.15)') count, count+Lx*Ly, J1_val, J1_val
          else 
            write(10,'(2I5,2es23.15)') count-Lx*Ly*(Lz-1), count, J1_val, J1_val
          end if
        end do
      end do
    end do
  end if

  !do i=1, L_N
  ! if(i>(L_N-Lx))then
  ! if(mod(i,Lx*(Ly-1))>=1 .and. mod(i,Lx*(Ly-1))<=Lx .and. i>Lx )then
  !  write(10,'(2I4,2es23.15)') i-Lx*(Ly-1), i, J1_val, J1_val
  ! else
  !  write(10,'(2I4,2es23.15)') i, i+Lx, J1_val, J1_val
  ! end if
  !end do

  !local mag
  open(10,file='list_local_mag.dat',position='append')
  do i=1, L_N
    write(10,'(I5)')i
  end do
  close(10)

  !two-point-corr
  open(10,file='list_cf_ss.dat',position='append')
  do i=1, L_N
    do j=1, L_N
      write(10,'(2I5)') i,j
    end do
  end do
  close(10)


  !read J1-bond information
  open(10,file='list_xxz_term_'//trim(size_name)//'.dat',status='old')
  do i=1, num_J1
    read(10,*) site_J1(1,i), site_J1(2,i), dum1, dum2
  end do
  close(10)


  open(10,file='list_site_position_'//trim(size_name)//'_type1.dat',position='append')
  do k=0, Lz-1
    do j=0, Ly-1
      do i=0, Lx-1
        write(10,'(3es23.15)') dble(i), dble(j), dble(k)  
      end do
    end do
  end do
  close(10)


  open(10,file='condition.dat',position='append')
  write(10,*) L_N*3, "# of J1 (=NOxxz)"
  write(10,*) L_N,   "# of localMag (=NOLM)"
  write(10,*) L_N*L_N, "# of 2point corr (=NOCF)"
  close(10)   

end program main


