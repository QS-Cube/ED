module lanczos
  use my_l1_blas
  use ham2vec
  implicit none
  integer :: ite_lancz, maxitr, dim, minitr, itrint
  real(8) :: lnc_ene_conv
  namelist /input_lancz/ lnc_ene_conv, maxitr, minitr, itrint
contains
  subroutine read_lanczos_para
    read(*,input_lancz)
    write(*,input_lancz)
  end subroutine read_lanczos_para
  !
  subroutine lanczos_routines(ene,THS,psi)
    use input_param, only: NOD 
    real(8), allocatable, intent(out) :: ene(:)
    integer, intent(in) :: THS
    complex(8), allocatable, intent(out) :: psi(:,:)
    complex(8), allocatable :: xpsi(:,:)
    complex(8) :: zdotc
    real(8) :: alpha(2)
    dim = THS
    allocate(ene(2))
    ene = 99.0d0
    read(*,input_lancz)
    write(*,input_lancz)
    write(*,'("************************************************************************************")')
    write(*,*) "Eigen solver: Lanczos"
    allocate(psi(dim,1))
    call lancz_eigen_val(ene(1),psi)
    write(*,'("************************************************************************************")')
    write(*,*) "total lanczos step:", ite_lancz
    write(*,'("************************************************************************************")')
    write(*,'("E_0 = ",es25.15)') ene(1)
    write(*,'("************************************************************************************")')
    !check_wave_vector
    allocate(xpsi(dim,2))
    call ham_to_vec_wave_vector(xpsi(:,1),psi(:,1),dim,NOD) 
    call ham_to_vec_wave_vector(xpsi(:,2),xpsi(:,1),dim,NOD) 
    alpha(1) = dble(zdotc(dim,psi(1,1),1,xpsi(1,1),1))
    alpha(2) = dble(zdotc(dim,psi(1,1),1,xpsi(1,2),1))
    write(*,'("************************************************************************************")')
    write(*,'("quality of eigenvector")')
    write(*,'("                   <GS| H |GS>  = ", es25.15)') alpha(1)
    write(*,'("|1 - (<GS|H|GS>)^2/<GS|H^2|GS>| = ", es25.15)') 1.0d0 - (alpha(1)**2)/alpha(2)
    write(*,'("************************************************************************************")')
    return
  end subroutine lanczos_routines
  !
  subroutine lancz_eigen_val(ene,psi_io)
    use input_param, only: NOD
    real(8), intent(out) :: ene
    complex(8), intent(out) :: psi_io(dim)
    complex(8), allocatable :: psi(:,:)
    real(8), allocatable :: alpha(:), beta(:), z(:,:)
    integer :: i, m, info, ell, il = 1, iu = 1, mitr
    real(8) :: tmp(2), dznrm2, abstol, dlamch, conv, vl = 0.0d0, vu = 0.0d0
    integer, allocatable :: iwork(:), ifail(:)
    real(8), allocatable :: work(:), w(:), d(:), e(:)
    complex(8) :: zdotc
    abstol = 2.0d0*dlamch('S')
    allocate( psi(dim,3), alpha(maxitr), beta(0:maxitr) )
    !$omp parallel do private(i,tmp)
    do i = 1, dim
      call random_number(tmp)
      psi(i,1) = cmplx(tmp(1),tmp(2),8)
      psi(i,1) = psi(i,1) - (0.5d0, 0.5d0)
    end do
    beta(0) = dznrm2( dim, psi(1,1), 1 )
    !call zdscal( dim, 1.0d0 / beta(0), psi(1,1), 1 )
    call my_zdscal( dim, 1.0d0 / beta(0), psi(1,1))
    call zcopy(dim,psi(1,1),1,psi_io,1)
    tmp = 0.0d0
    mitr= min(maxitr,dim)
    do i = 1, mitr
      call ham_to_vec_wave_vector(psi(:,2),psi(:,1),dim,NOD) 
      alpha(i) = dble(zdotc(dim,psi(1,1),1,psi(1,2),1))
      if(i > 1)then
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1) - beta(i-1)*psi(ell,3)
        end do
      else
        !call zaxpy( dim, (-1.0d0,0.0d0)*alpha(i), psi(1,1), 1, psi(1,2), 1 )
        call my_zaxpy( dim, (-1.0d0,0.0d0)*alpha(i), psi(1,1), psi(1,2) )
      end if
      beta(i) = dznrm2( dim, psi(1,2), 1 )
      call zcopy(dim,psi(1,1),1,psi(1,3),1)
      !$omp parallel do
      do ell = 1, dim
        psi(ell,1) = psi(ell,2) / beta(i)
      end do
      if( i==mitr .or. i >= minitr .and. mod(i,itrint) == 0)then
        allocate(w(i),d(i),e(i),iwork(5*i), ifail(i), work(5*i), z(i,1))
        if(i == 1)then
           write(*,*) "### dim == 1"
           ene = alpha(1)
           psi_io(1) = (1.0d0,0.0d0)
           conv = 0.0d0
           write(*,'("itr, ene, conv: ", i5, es25.15, es11.1)') i, ene, conv
           return
        end if
        d = alpha(1:i); e = beta(1:i)
        call dstevx('V', 'I', i, d, e, vl, vu, il, iu, abstol, &
          m, w, z, size(z,1), work, iwork, ifail, info)
        if(info /= 0)then
          write(*,*) "error in dstevx in eigen_solver.f90. info =", info; stop
        end if
        ene = w(1)
        conv = abs(tmp(1)/ene - 1.0d0)
        write(*,'("itr, ene, conv: ", i5, es25.15, es11.1)') i, ene, conv
        if( i==maxitr .or. conv < lnc_ene_conv )then
          exit
        end if
        tmp(1) = ene
        deallocate(w, d, e, iwork, ifail, work, z)
      end if
    end do
    ite_lancz = i
    !
    ! calcu. eigenvector
    !
    call zcopy(dim,psi_io,1,psi(1,1),1)
    !call zdscal(dim, z(1,1), psi_io, 1)
    call my_zdscal(dim, z(1,1), psi_io)
    tmp = 0.0d0
    do i = 1, ite_lancz-1
      call ham_to_vec_wave_vector(psi(:,2),psi(:,1),dim,NOD) 
      if(i > 1)then
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1) - beta(i-1)*psi(ell,3)
        end do
      else
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1)
        end do
      end if
      call zcopy(dim,psi(1,1),1,psi(1,3),1)
      !$omp parallel do
      do ell = 1, dim
        psi(ell,1) = psi(ell,2) / beta(i)
      end do
      !call zaxpy( dim, (1.0d0,0.0d0)*z(i+1,1), psi(1,1), 1, psi_io, 1 )
      call my_zaxpy( dim, (1.0d0,0.0d0)*z(i+1,1), psi(1,1), psi_io)
    end do
    beta(0) = dznrm2( dim, psi_io, 1 )
    !call zdscal( dim, 1.0d0 / beta(0), psi_io, 1 )    
    call my_zdscal( dim, 1.0d0 / beta(0), psi_io)    
    return
  end subroutine lancz_eigen_val
  !
end module lanczos
