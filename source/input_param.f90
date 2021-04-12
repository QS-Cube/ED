module input_param
  implicit none
  integer :: NOS       ! number of sites
  integer :: NOD       ! number of down spins
  integer :: NOV       ! Number of lowest eigenvectors used for computing static physical quantities.
  integer :: THS       ! total hilbelt space of S^z_{tot} = NOS/2 - NOD: THS   = combination(NOS,NOD)
  integer :: NOLM      ! numver of expectation values <Sz_i>
  integer :: NOCF      ! numver of expectation values <Sz_i sz_j> or <S+_i S-_j>
  integer :: NOxxz     ! number of XXZ interaction
  integer :: ALG        !1:Conventional Lanczos,  2:Thick-Restarted Lanczos, 3:Full diagonalization
  integer, allocatable :: p_xxz(:,:)  ! pair of xxz interaction: p_xxz(1:2,1:NOxxz) !2020/7/17-add
  real(8), allocatable :: Jint(:,:)   ! XXZ term: Jint(1:NOxxz, 2:3)
  real(8), allocatable :: sJint(:,:)  ! scaled XXZ term: sJint(1:NOxxz, 2) = Jint(1:NOxxz, 2)/2, 
                                      !                  sJint(1:NOxxz, 3) = Jint(1:NOxxz, 2)/4, 
  character(100) :: FILExxz, FILElm, FILECF, FILEwf, FILEp, FILEp2 !2020/1/5-modified, 2020/1/12-add
  character(100) :: OUTDIR

  ! ! path of file containing interaction list of xxz interactions
  integer :: random, check_p, wr_wf, re_wf, cal_lm, cal_cf, cal_dsf   ! flags

  ! !for DSSF/DQSF
  integer ::spsmsz                      ! chose DSF, 1:S+(q,w), 2:S-(q,w), 3:Sz(q,w), 4:Q--(q,w)
  integer ::itr_dsf                     ! iteration for DSF, itr_dsf=0 then we do not calculate DSF
  real(8) ::qx, qy, qz, kvec_calc(1:3)  ! wave vectors
  real(8) :: rfield                     ! Magnetic field value
  character(100) ::FILEpos              ! site position file

  ! !for wavevector decomposition 
  integer :: LX, LY, LZ       
  integer :: KX, KY, KZ  
  real(8) :: PI, rkx, rky, rkz
  complex(8),allocatable :: explist(:,:,:)
  integer, allocatable :: list_s(:)                         
  real(8), allocatable :: list_r(:)
  integer, allocatable :: shift_x_SQ(:), shift_y_SQ(:), shift_z_SQ(:) 

  !for DSF under wavevector decomposition 
  integer::NOD_new
  integer::THS_new
  real(8), allocatable :: list_r_new(:)
  integer, allocatable::list_s_new(:)
  real(8)::rkx_new,rky_new,rkz_new
  complex(8),allocatable::explist_new(:,:,:)

  !
  namelist /input_parameters/ NOS, NOD, LX,LY,LZ,KX,KY,KZ,NOxxz, &
    &          ALG, cal_lm, cal_cf, cal_dsf, wr_wf, re_wf, FILExxz, FILEwf, OUTDIR

  namelist /input_static/ NOV, NOLM, NOCF, FILElm, FILECF

  namelist /input_dynamic/ spsmsz, itr_dsf, qx, qy, qz, rfield, FILEpos
  !
  real(8), allocatable :: ene(:)
  complex(8),allocatable :: psi(:,:) 
  !
contains
  !
  subroutine read_ip
    integer :: i,j,k
    integer :: tmp
    !
    write(*,'(" ### Read input_parameters. ")')
    read(*,input_parameters)
    write(*,input_parameters)
    !
    write(*,'(" ### Read input_static. ")')
    read(*,input_static)
    write(*,input_static)
    !
    write(*,'(" ### Read input_dynamic. ")')
    read(*,input_dynamic)
    write(*,input_dynamic)
    !
    write(*,'(" ### Set wavevectors. ")')
    PI = acos(-1.0d0)     
    rkx = 2.0d0 * PI * dble(KX)/dble(LX) 
    rky = 2.0d0 * PI * dble(KY)/dble(LY) 
    rkz = 2.0d0 * PI * dble(KZ)/dble(LZ) 
    !
    write(*,'(" ### Set random_seed. ")')
    call random_seed
    !
    write(*,'(" ### Allocation arrays. ")')
    allocate(p_xxz(2,NOxxz),Jint(NOxxz,2:3),sJint(NOxxz,2:3),shift_x_SQ(NOS),shift_y_SQ(NOS),&
      shift_z_SQ(NOS),explist(0:LX,0:LY,0:LZ))
    !
    write(*,'(" ### Set phase factors. ")')
    forall(i=0:LX,j=0:LY,k=0:LZ) explist(i,j,k) = &
      exp((0.0d0,1.0d0)*(rkx*dble(i)+rky*dble(j)+rkz*dble(k)))
    !
    write(*,'(" ### Read FILExxz. ")')
    open(10, file=trim(adjustl(FILExxz)),status='old')
    do i = 1, NOxxz
      read(10,*) p_xxz(1,i), p_xxz(2,i), Jint(i,2), Jint(i,3)
      sJint(i,2) = Jint(i,2) * 0.5d0
      sJint(i,3) = Jint(i,3) * 0.25d0
      if( p_xxz(1,i) > p_xxz(2,i) )then
        tmp = p_xxz(1,i)
        p_xxz(1,i) = p_xxz(2,i)
        p_xxz(2,i) = tmp
      end if
    end do
    close(10)
    !
    write(*,'(" ### write ouput/FILExxz.dat. ")')
    open(10,file=trim(adjustl(OUTDIR))//'input_xxz.dat',position='append')
    write(10,'("************************************************************************************")')
    write(10,*) "         i         j          J^{2}_{ij}           J^{3}_{ij}"
    write(10,'("************************************************************************************")')
    do i = 1, NOxxz
      write(10,*) p_xxz(1,i), p_xxz(2,i), Jint(i,2), Jint(i,3)
    end do
    write(10,'("************************************************************************************")')
    close(10)    
    !
    return
  end subroutine read_ip

  subroutine write_wf(dim,i_vec_min,i_vec_max, NOE, ene) 
    integer,intent(in)::dim
    integer,intent(in)::i_vec_min,i_vec_max
    integer,intent(in)::NOE
    real(8),intent(in)::ene(1:NOE)
    character(200)::file_name
    character(500)::file_name1
    integer::i_vec,i
    if(wr_wf /= 1) return

    do i_vec=i_vec_min, i_vec_max
      write(file_name,'(A100,A2,I1,A4)')trim(adjustl(FILEwf)),'wf',i_vec,'.bin'
      write(*,*) " Write wave functions (",i_vec,"-th state) "
      print *, file_name
      open(10, file=trim(adjustl(file_name)), status='unknown', action='write', form='unformatted')
      write(10) psi(1:dim,i_vec)
      close(10)
    end do

    write(file_name1,'(A300,A10)')trim(adjustl(FILEwf)),'energy.dat'
    open(10, file=trim(adjustl(file_name1)), form='formatted')
    do i=1, NOE
      write(10,'(2ES24.15)') dble(i), ene(i) 
    end do
    close(10)

    write(*,'("************************************************************************************")')
    return
  end subroutine write_wf

  subroutine read_wf(dim,i_vec_min, i_vec_max,NOE)
    integer,intent(in)::dim
    integer,intent(in)::i_vec_min, i_vec_max
    integer,intent(in)::NOE
    character(200)::file_name
    character(500)::file_name1
    integer::i_vec, i
    real(8)::dummy
    if(re_wf /= 1) return

    allocate(psi(1:dim,i_vec_min:i_vec_max))
    do i_vec=i_vec_min, i_vec_max
      write(file_name,'(A100,A2,I1,A4)')trim(adjustl(FILEwf)),'wf',i_vec,'.bin'
      write(*,*) " read wave functions (",i_vec,"-th state) "
      open(10, file=trim(adjustl(file_name)), form='unformatted')
      read(10) psi(1:dim,i_vec)
      close(10)
      write(*,*) "write psi(1:10,i_vec)"
      write(*,*) psi(1:10,i_vec)
      write(*,*) " finish reading wave function!"
    end do

    write(file_name1,'(A300,A10)')trim(adjustl(FILEwf)),'energy.dat'
    allocate(ene(1:NOE))
    ene(1:NOE)=0.0d0
    open(10, file=trim(adjustl(file_name1)),status='old')
    do i=1, NOE
      read(10,*) dummy, ene(i) 
      write(*,*) dummy, ene(i), "reading!"
    end do
    close(10)
    write(*,*) "finish reading energy!"

    write(*,'("************************************************************************************")')
    return
  end subroutine read_wf

end module input_param
