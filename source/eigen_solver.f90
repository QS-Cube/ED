module eigen_solver
  use ham2vec
  implicit none
  integer :: NOE,NOK, NOM, maxitr, ite_lancz, i_vec_min, i_vec_max
  integer :: dim
  real(8) :: lnc_ene_conv
  namelist /input_TRLan/ NOE, NOK, NOM, maxitr, lnc_ene_conv, i_vec_min, i_vec_max
contains
  subroutine my_trlanczos_routines_complex(ene,psi,THS)
    use input_param, only : OUTDIR
    real(8), allocatable, intent(out) :: ene(:)
    complex(8),allocatable, intent(out) :: psi(:,:)
    integer, intent(in) :: THS 
    integer :: i
    read(*,input_TRLan)
    write(*,input_TRLan)

    if(NOE > NOK) stop "ne must be eaqual or smaller than NOK!"
    if(NOM < 2*NOK) NOM = 2*NOK
    dim = THS
    if(dim < NOE .or. dim < NOE .or. dim < NOK .or. dim < NOM) stop 999
    allocate( psi(dim,NOM+1), ene(NOM) )
    write(*,*) "Reduced THS=", THS 
    write(*,'("************************************************************************************")')
    write(*,*) "Eigen solver: my thick restart Lanczos"
    call lancz_eigen_val_complex(ene, psi) 
    write(*,'("***********************************************************************************")')
    write(*,*) "total lanczos step:", ite_lancz
    write(*,'("***********************************************************************************")')
    write(*,'("    i           E(i)")')
    do i = 1, NOE
      write(*,'(i5, 1es25.15)') i, ene(i)
    end do

    open(10,file=trim(adjustl(OUTDIR))//'energy.dat',position='append')
    do i = 1, NOE
      write(10,'(i5, 1es25.15)') i, ene(i)
    end do
    close(10)

    write(*,'("***********************************************************************************")')
    return
  end subroutine my_trlanczos_routines_complex

  subroutine read_input_TRLan
    read(*,input_TRLan)
    write(*,input_TRLan)
  end subroutine read_input_TRLan

  subroutine Full_diag_routines(ene,psi,THS)
    real(8),allocatable,intent(out)::ene(:)
    complex(8),allocatable,intent(out)::psi(:,:)
    complex(8),allocatable::Ham(:,:)
    integer,intent(in)::THS
    dim=THS
    allocate(ene(dim),psi(dim,dim),Ham(dim,dim))
    Ham(:,:) = (0.0d0,0.0d0)
    ene(:)=0.0d0
    psi(:,:)=(0.0d0,0.0d0)

    print *, "make ham"
    call make_full_hamiltonian(dim,Ham)
    print *, "diag starts"
    call fulldia_routines_with_momentum(dim,Ham,psi,ene)

  end subroutine Full_diag_routines

  subroutine lancz_eigen_val_complex(ene, psi) 
    use input_param, only: NOD, list_s, list_r, explist 
    real(8), intent(inout) :: ene(1:NOM)
    complex(8), intent(inout) :: psi(1:dim,1:NOM+1)
    real(8), allocatable :: alpha(:), beta(:), coef(:,:), conv(:)
    complex(8),allocatable:: T(:,:), Qk(:,:), overlap(:),calpha(:),cbeta(:),ccoef(:,:)
    integer, allocatable :: isuppz(:), iwork_lanc(:), ifail_lanc(:)
    complex(8), allocatable :: work_lanc(:)
    real(8),allocatable::rwork_lanc(:)
    integer :: count_ne, info, lwork, lrwork,liwork, tritr, flag
    integer :: i, j, ell
    real(8) :: eps, dlamch, dznrm2, abstol, max_conv,tmp(2)
    integer t1,t2,t3,t4,t_rate,t_max
    complex(8) :: zdotc
    complex(8) :: prdct 
    flag = 0
    allocate(conv(NOK), alpha(NOM), beta(0:NOM), coef(NOM,NOM), isuppz(2*NOM), &
      overlap(NOM), T(NOM, NOM), Qk(dim,NOK),ccoef(NOM,NOM), calpha(NOM), cbeta(0:NOM))

    !*********************************
    !***** initial lanczos step ******
    !*********************************
    !$OMP parallel do private(i,tmp)
    do i = 1, dim
      call random_number(tmp)
      psi(i,1) = cmplx(tmp(1),tmp(2),8)
      psi(i,1) = psi(i,1) - (0.5d0, 0.5d0)
    end do
    !
    beta(0) = dznrm2( dim, psi(1,1), 1 )
    do i = 1, NOM
      call zdscal( dim, 1.0d0/beta(i-1), psi(1,i), 1 )
      call system_clock(t1)
      call ham_to_vec_wave_vector(Qk(:,1),psi(:,i),dim,NOD,list_s,list_r,explist) 
      call system_clock(t2,t_rate, t_max)
      !write(*,*) "ham_to_vec", (t2-t1)/dble(t_rate)
      alpha(i) = dble(zdotc(dim,psi(1,i),1,Qk(:,1),1))
      if(i > 1)then
        !$OMP parallel do private(ell)
        do ell = 1, dim
          psi(ell,i+1) = Qk(ell,1) - alpha(i)*psi(ell,i)-beta(i-1)*psi(ell,i-1)
        end do
        !
        call full_reorthogonalization_complex 
      else
        !$OMP parallel do private(ell)
        do ell = 1, dim
          psi(ell,i+1) = Qk(ell,1) - alpha(i)*psi(ell,i)
        end do
        !
      end if
      beta(i) = dznrm2( dim, psi(1,i+1), 1 )
    end do
    !
    call zdscal( dim, 1.0d0/beta(NOM), psi(1,NOM+1), 1 )
    !
    abstol = 2.0d0*dlamch('S')
    allocate(iwork_lanc(5*NOM), ifail_lanc(NOM), rwork_lanc(5*NOM))
    call dstevx('V', 'A', NOM, alpha, beta(1), 0.0d0, 1.0d0, 1, NOM, eps, &
      count_ne, ene, coef, NOM, rwork_lanc, iwork_lanc, ifail_lanc, info)
    if(info /= 0)then
      write(*,*) "error in dstevx in eigen_solver.f90. info =", info; stop
    end if
    deallocate(iwork_lanc, ifail_lanc, rwork_lanc)
    allocate(ifail_lanc(2*NOM))
    !
    !print *, "First lanczos step over"
    !do i=1,NOM
    !  print *, i, ene(i)
    !end do
    !
    !*********************************
    !***** restart lanczos step ******
    !*********************************
    ccoef = cmplx(coef,0.0d0,8)
    calpha = cmplx(alpha,0.0d0,8)
    cbeta = cmplx(beta,0.0d0,8)
    deallocate(coef,alpha,beta)
    do tritr = 1, maxitr
      call system_clock(t3)
      !
      ccoef = conjg(ccoef)
      call zgemm('N', 'N', dim, NOK, NOM, (1.0d0,0.0d0), &
        psi, dim, ccoef, NOM, (0.0d0,0.0d0), Qk, dim)
      call zcopy(dim*NOK,Qk,1,psi,1)
      conv(1:NOK) = abs(calpha(1:NOK)/ene(1:NOK) - 1.0d0)
      max_conv = maxval(conv)
      if(max_conv < lnc_ene_conv)then
        flag = flag + 1
      else
        flag = 0
      end if
      write(*,'("tritr, ene(NOK), max_conv, flag: ", i5, es25.15, es11.1, i5)') &
        tritr, ene(NOK), max_conv, flag
      if(flag == 2) exit
      if(tritr == maxitr)then
        write(*,*) "tritr = maxitr in trlanczos!"; exit
      end if
      call zcopy(dim,psi(1,NOM+1),1,psi(1,NOK+1),1)
      calpha(1:NOK) = cmplx(ene(1:NOK),0.0d0,8)
      !
      forall(j=1:NOK) cbeta(j) = cbeta(NOM) * ccoef(NOM,j)
      i = NOK+1
      call ham_to_vec_wave_vector(Qk(:,1),psi(:,i),dim,&
        NOD,list_s,list_r,explist) 
      calpha(i) = zdotc(dim,psi(1,i),1,Qk(1,1),1)
      do j=1, NOK
        call zaxpy(dim, -cbeta(j),psi(1,j),1,Qk(1,1),1)
      end do
      !$OMP parallel do private(ell)
      do ell = 1, dim
        psi(ell,i+1) = Qk(ell,1) - calpha(i)*psi(ell,i)
      end do
      !
      call full_reorthogonalization_complex
      cbeta(i) = cmplx(dznrm2( dim, psi(1,i+1), 1 ),0.0d0,8)
      do i = NOK+2, NOM
        call zscal( dim, (1.0d0,0.0d0)/cbeta(i-1), psi(1,i), 1 )
        call ham_to_vec_wave_vector(Qk(:,1),psi(:,i),dim,&
          NOD,list_s,list_r,explist)
        calpha(i) = zdotc(dim,psi(1,i),1,Qk(:,1),1)
        !$OMP parallel do private(ell)
        do ell = 1, dim
          psi(ell,i+1) = Qk(ell,1) - calpha(i)*psi(ell,i)-cbeta(i-1)*psi(ell,i-1)
        end do
        !
        call full_reorthogonalization_complex 
        cbeta(i) = cmplx(dznrm2( dim, psi(1,i+1), 1 ),0.0d0,8)
      end do
      call zscal( dim, (1.0d0,0.0d0)/cbeta(NOM), psi(1,NOM+1), 1 )
      T = (0.0d0,0.0d0)
      do j = 1, NOM
        T(j,j) = calpha(j)
      end do
      call zcopy(NOK,cbeta(1),1,T(1,NOK+1),1)
      do j = NOK+1, NOM-1
        T(j,j+1) = cbeta(j)
      end do
      if(tritr == 1) call get_work_mem_for_zheevr 
      call zheevr('V', 'A', 'U', NOM, T, NOM, 0.0d0, 1.0d0, 1, 1, abstol, &
        count_ne, ene, ccoef, NOM, isuppz, work_lanc, lwork, &
        rwork_lanc, lrwork, iwork_lanc, liwork, info)
      if(info /= 0)then
        write(*,*) "error in dxyevr in eigen_solver.f90. info =", info; stop 4
      end if

      !print *, "Temporal energy", "  tritr=", tritr
      !do j=1, NOE
      !  print *, j, ene(j)
      !end do

      call system_clock(t4,t_rate, t_max)
      !write(*,*) "time for 1 iteration", (t4-t3)/dble(t_rate), "sec"

    end do
    ite_lancz = tritr*NOM

    !check accuracy: calcu <Phi|H|Phi> 2020/7/30-add
    print *, "Check accuracy <Phi|H|Phi>"
    do j=1, NOE
      call ham_to_vec_wave_vector(Qk(1:dim,j),psi(1:dim,j),dim,&
        NOD,list_s,list_r,explist)
      prdct=(0.0d0,0.0d0)
      do i=1, dim
        prdct=prdct+conjg(psi(i,j))*Qk(i,j)
      end do
      print *, j, prdct  
    end do

    ! end if !Thick restart Lanczos ends

    return
  contains
    !
    subroutine full_reorthogonalization_complex 
      do j = 1, i-1 
        overlap(j) = zdotc(dim,psi(1,j),1,psi(:,i+1),1) 
        call zaxpy(dim,-overlap(j),psi(1,j),1,psi(:,i+1),1) 
      end do
    end subroutine full_reorthogonalization_complex
    !
    subroutine get_work_mem_for_zheevr 
      abstol = dlamch('S')
      lwork = -1
      lrwork = -1
      liwork= -1
      allocate(iwork_lanc(1), work_lanc(1), rwork_lanc(1))
      call zheevr('V', 'A', 'U', NOM, T, NOM, 0.0d0, 1.0d0, 1, 1, abstol, &
        count_ne, ene, ccoef, NOM, ifail_lanc, work_lanc, lwork, &
        rwork_lanc, lrwork, iwork_lanc, liwork, info)
      if(info /= 0)then
        write(*,*) "error in dxyevr in eigen_solver.f90. info =", info
        stop 3
      end if
      lwork = nint(dble(work_lanc(1)))
      lrwork = nint(rwork_lanc(1))
      liwork= iwork_lanc(1)
      deallocate(iwork_lanc, work_lanc, rwork_lanc)
      allocate(work_lanc(lwork), iwork_lanc(liwork), rwork_lanc(lrwork))
    end subroutine get_work_mem_for_zheevr
    !
  end subroutine lancz_eigen_val_complex

  subroutine fulldia_routines_with_momentum(THS,ham_complex,psi0,ene0)
    use input_param, only : OUTDIR
!!! Use blas_interfaces, Only: zscal
!!! Use lapack_example_aux, Only: nagf_file_print_matrix_complex_gen
!!! Use lapack_interfaces, Only: ddisna, zheev
!!! Use nag_library, Only: x04daf !ddisna, nag_wp, x02ajf, x04daf, zheev, zscal
    integer, intent(in)::THS
    real(8),intent(inout)::ene0(THS)
    complex(8),intent(in)::ham_complex(THS,THS)
    complex(8),intent(inout)::psi0(THS,THS)
    !complex(8), allocatable ::v0(:), v1(:)
    integer :: k
    !LAPACK
    !integer :: nin, nout, ifail, j
    !real(8) :: value
    integer :: nb, nmax
    integer :: lda, lwork
    real(8) :: eerrbd, eps
    integer :: i, info, lwkopt, n
    complex(8), allocatable :: a(:, :), work(:)
    complex(8) :: dummy(1)
    real(8), allocatable :: rcondz(:),rwork(:), w(:), zerrbd(:)
    real(8) :: dlamch
    real(8), external :: ddot
    complex(8), external :: zdotc
    External dlamch
    External ddisna, dsyev, dgemv, zheev, x04daf
    !Intrinsic :: abs, cmplx, conjg, epsilon, max, maxloc, nint, real

    nb=THS
    n=THS
    nmax=THS
    lda=nmax
    !allocate(v0(1:THS), v1(1:THS))
    allocate(a(lda,nmax))
    allocate(rcondz(nmax))
    allocate(w(nmax))
    allocate(zerrbd(nmax))
    allocate(rwork(3*nmax-2))
    a=(0.0d0,0.0d0)
    rcondz=0.0d0
    w=0.0d0
    zerrbd=0.0d0
    rwork=0.0d0

    !     Use routine workspace query to get optimal workspace.
    lwork = -1
    Call zheev('V', 'U', n, a, lda, w, dummy, lwork, rwork, info)

    !     Make sure that there is enough workspace for block size nb.
    lwork = max((nb+1)*n, nint(real(dummy(1))))
    allocate (work(lwork))

    !     Solve the Hermitian eigenvalue problem
    a=ham_complex
    Call zheev('V', 'U', n, a, lda, w, work, lwork, rwork, info)

    lwkopt = nint(dble(work(1)))
    if(info>0)then
      print *, "info is not 0, for dsyev"
      stop
    end if

    if(info==0)then
      !print *, "           "
      print *, "Eigenvalues by Full Diagonalization"
      !do k=1, THS
      !  print *, w(k)
      !end do
      open(111,file=trim(adjustl(OUTDIR))//'full_eigenvalues.dat',position='append')
      do k=1, THS
        write(111,'(es23.15)') w(k)
      end do
      close(111)

      ene0(1:THS)=w(1:THS) !2020/7/28-add

      open(112, file=trim(adjustl(OUTDIR))//'full_eigenvectors.dat',position='append')
      do k=1, THS
        write(112,*) a(k,1)
      end do
      close(112)

      psi0(1:THS,1:THS)  = a(1:THS,1:THS) 

    else 
      print *, "Failure in ZHEEV. INFO =", info
    end if

    !check the acuracy by <v|H|v>-E
    !print *, "check the acuracy"
    !do k=1, THS
    !  v0 = a(1:THS,k)
    !  v1 = 0.0d0
    !  !calcu H|v0> -> v1
    !  call zgemv('n',THS,THS,1.0d0,ham_complex,THS,v0,1,0.0d0,v1,1) 
    !end do

  end subroutine fulldia_routines_with_momentum

  subroutine calcu_FM_energy(ene0)
    use input_param, only:NOxxz, Jint
    real(8),intent(inout)::ene0

    ene0=sum(Jint(1:NOxxz,3))/4.0d0
    return
  end subroutine calcu_FM_energy

end module eigen_solver
