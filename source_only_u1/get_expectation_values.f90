module get_expectation_values
  use my_l1_blas
  implicit none
  integer, allocatable :: list_expe_sz(:), list_expe_ss(:,:)
  real(8), allocatable :: expe_sz(:,:)
  real(8), allocatable :: expe_szsz(:,:)
  complex(8), allocatable :: expe_spsm(:,:)
  real(8), allocatable :: Spos(:,:)          
contains
  !
  subroutine allocate_expe_mem(NOS,NOV,NOLM,NOCF,FILElm,FILECF,FILEpos) 
    integer, intent(in) :: NOS, NOV, NOLM, NOCF           
    character(*), intent(in) :: FILElm, FILECF, FILEpos 
    integer :: i
    allocate(list_expe_sz(NOLM), list_expe_ss(2,NOCF), &
      expe_sz(NOLM,NOV), expe_szsz(NOCF,NOV), expe_spsm(NOCF,NOV))
    open(10, file=trim(adjustl(FILElm)),status='old')
    do i = 1, NOLM
      read(10,*) list_expe_sz(i)
    end do
    close(10)
    !
    write(*,*) "### Read FILECF [list_expe_ss(2:1:-1,i)]"
    open(10, file=trim(adjustl(FILECF)),status='old')
    do i = 1, NOCF
      read(10,*) list_expe_ss(2,i), list_expe_ss(1,i)
    end do
    close(10)

    allocate(Spos(1:3,1:NOS))
    open(10, file=trim(adjustl(FILEpos)),status='old')
    do i=1, NOS
      read(10,*) Spos(1,i), Spos(2,i), Spos(3,i)
    end do
    close(10)

    return
  end subroutine allocate_expe_mem
  !
  subroutine get_lm_2_wave_vector(psi,NOD,dim,nvec,NOS)
    use input_param, only: OUTDIR, NOLM
    use ham2vec, only: calcu_lm
    integer, intent(in) :: NOD, nvec, NOS
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1:nvec)
    integer :: j, k
    !
    expe_sz=0.0d0    
    do k = 1, nvec
      do j = 1, NOLM
        call calcu_lm(psi(1,k),dim,NOD,NOS,expe_sz(j,k),list_expe_sz(j))
      end do
    end do

    write(*,'(" ### write ouput/local_mag.dat. ")')
    open(10,file=trim(adjustl(OUTDIR))//'local_mag.dat',position='append')
    do j = 1, NOLM
      write(10,'(i8,10000es23.15)') list_expe_sz(j), ( expe_sz(j,k), k=1, nvec )
    end do
    close(10)

    return

  end subroutine get_lm_2_wave_vector

  subroutine get_cf_2_wave_vector(psi,NOD,dim,nvec,NOS)
    use input_param, only: OUTDIR, NOCF
    use ham2vec, only: calcu_cf
    integer, intent(in) :: NOD,nvec,NOS
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1:nvec)
    integer :: j, k

    expe_szsz=0.0d0
    expe_spsm=0.0d0

    do k = 1, nvec
      do j = 1, NOCF
        call calcu_cf(psi(1,k),dim,NOD,NOS,list_expe_ss(1,j),expe_szsz(j,k),expe_spsm(j,k))
      end do
    end do
    !
    write(*,'(" ### write ouput/SzSz-corr.dat. ")')
    open(10,file=trim(adjustl(OUTDIR))//'SzSz-corr.dat',position='append')
    do j = 1, NOCF
      write(10,'(2i8,10000es25.15)') list_expe_ss(2,j), list_expe_ss(1,j), ( expe_szsz(j,k), k=1, nvec )
    end do
    close(10)    

    write(*,'(" ### write ouput/SpSm-corr.dat. ")')
    open(10,file=trim(adjustl(OUTDIR))//'SpSm-corr.dat',position='append')
    do j = 1, NOCF
      write(10,'(2i8,10000es25.15)') list_expe_ss(2,j), list_expe_ss(1,j), ( expe_spsm(j,k), k=1, nvec )
    end do
    close(10)

    return
  end subroutine get_cf_2_wave_vector

  subroutine calcu_Sq_a(spsmsz,psi,fn,kvec_calc,dim,dim_new,Spos)
    use input_param, only: NOS,NOD,NOD_new
    use state_lists, only: list_fly, search_Spi_state_SQ, search_Smi_state_SQ
    implicit none
    integer,intent(in)::dim,dim_new
    integer, intent(in)::spsmsz
    real(8),intent(in)::kvec_calc(1:3),Spos(1:3,1:NOS)
    complex(8),intent(in)::psi(1:dim)  
    complex(8),intent(inout)::fn(1:dim_new)
    logical :: fi
    integer:: i,l
    integer:: lp
    real(8) Siz
    integer, allocatable :: st_list(:)
    complex(8), allocatable :: exp_list_KS(:)
    allocate(st_list(NOD))
    allocate(exp_list_KS(NOS))

    fn = (0.0d0,0.0d0)

    do i = 1, NOS
      exp_list_KS(i) = exp((0.0d0,-1.0d0)*(dot_product(kvec_calc(1:3),Spos(1:3,i))))
    end do

    if(spsmsz==1)then
      do i=1, NOS
        !$omp parallel do private(l,lp,fi,st_list)
        do l=1, dim
          lp=0
          st_list = list_fly(l,NOD,NOS)
          fi=any(st_list==i)
          if(fi)then !If i-site has down spin
            call search_Spi_state_SQ(i,NOD,NOD_new,st_list,lp)
            !$omp atomic
            fn(lp)=fn(lp)+ psi(l)*exp_list_KS(i)
          end if
        end do
        !$omp end parallel do
      end do
      !
    else if(spsmsz==2)then
      do i=1, NOS
        !$omp parallel do private(l,lp,fi,st_list)
        do l=1, dim
          lp=0  
          st_list = list_fly(l,NOD,NOS)
          fi=any(st_list==i)
          if(.not. fi)then !If i-site has up spin
            call search_Smi_state_SQ(i,NOD,NOD_new,st_list,lp)  
            !$omp atomic
            fn(lp)=fn(lp)+ psi(l)*exp_list_KS(i)
            !
          end if
        end do
        !$omp end parallel do 
      end do
      !
    else if(spsmsz==3)then
      do i=1, NOS
        !$omp parallel do private(l,lp,fi,Siz,st_list)
        do l=1, dim
          lp=0  
          st_list = list_fly(l,NOD,NOS)
          fi=any(st_list==i)
          if(fi)then !If i-site has up spin
            Siz=-0.50d0
          else
            Siz=0.50d0
          end if
          lp = l
          fn(lp)=fn(lp)+ Siz*psi(l)*exp_list_KS(i)
        end do
        !$omp end parallel do
      end do
    end if
    !
    !call zdscal(dim_new,1.0d0/sqrt(dble(NOS)),fn,1)
    call my_zdscal(dim_new,1.0d0/sqrt(dble(NOS)),fn)
    !
  end subroutine calcu_Sq_a

  subroutine calcu_DSF_wavevector(psi,ene,kvec_calc,rfield,NOS,NOD,itr_dsf,dim) 
    use input_param, only: NOD_new,THS_new,spsmsz,OUTDIR
    use state_lists, only: combination
    use ham2vec, only:  ham_to_vec_wave_vector
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1)
    real(8),intent(in):: ene(1:2),rfield
    integer, intent(in) :: NOS, NOD,itr_dsf

    !for DSF
    integer::dim_new
    integer::ex_i,i,j,itr=1,itr_num !k,
    complex(8),allocatable :: fn(:)    !|f_n>
    complex(8),allocatable :: fnp1(:)  !|f_n+1>
    complex(8),allocatable :: fnm1(:)  !|f_n-1>
    real(8),intent(in)  :: kvec_calc(1:3)
    real(8) :: norm1
    real(8) :: eps, sztot
    integer, parameter :: la=4000
    real(8) :: alpha_dsf(la), beta_dsf(la)
    !
    ! alpha_dsf(i) and beta_dsf(i+1) correspond to 
    ! \alpha_i and \beta_i for Eq.~(23) in Sec.3.9 of arXiv:****.****
    !
    character(2) :: ex_index
    character(10) :: sztot_index
    real(8)::dznrm2,sztot_new
    complex(8) ::zdotc
    !
    print *, "DSF calculation starts"
    !
    write(*,'(" ### Set epsilon. ")')
    eps=1.0d-8
    !
    if(NOD==NOS/2 .and. spsmsz==2)then !if tsz=0
      print *, "You don't need to calculate S(q,w)^-"
      print *, "because S(q,w)^-=S(q,w)^+ in sztot=0 case."
      print *, "S(q,w)^x is obtained by (1/2)S(q,w)^+"
      print *, "But note that you should take rfield=0.0d0 when you calculate S(q,w)^+"
      print *, "Or just calculate S(q,w)^z because S(q,w)^x,y,z, are all same result in sztot=0."
      stop
    end if
    !
    if(NOD==0 .or. NOD==1 .and. spsmsz==1)then !if FM state or sztot=NOS/2-1 state
      print *, "You can not calculate S+(q,w) on the FM state and sztot=NOS/2-1 state!"  
      stop
    end if
    !
    print *, dznrm2(dim,psi(1,1),1), "psi"
    !
    write(*,'(" ### Set NOD_new. ")')
    if(spsmsz==1)then !Sp(q,w)
      NOD_new = NOD-1
    else if(spsmsz==2)then !Sm(q,w)
      NOD_new = NOD+1
    else if(spsmsz==3)then !Sz(q,w)
      NOD_new = NOD
    end if
    !
    write(*,'(" ### Allocate and Set arrays for intermediate state_lists. ")')
    THS_new = combination(NOS,NOD_new)  !Dimension before wave vector decomposition
    write(*,*) "  THS_new   = ", THS_new
    dim_new = THS_new
    !
    write(*,'(" ### Allocate work arrays. ")')
    allocate(fn(dim_new),fnp1(dim_new),fnm1(dim_new))
    !
    write(*,'(" ### Initializations: alpha, beta, fns ")')
    !$omp parallel 
    !$omp do private(i) 
    do i = 1, la
      alpha_dsf(i) = 1.0d300 !
      beta_dsf(i) = 0.0d0
    end do
    !$omp end do
    !$omp do private(i) 
    do i = 1, dim_new
      fn(i) = (0.0d0,0.0d0)
      fnp1(i) = (0.0d0,0.0d0)
      fnm1(i) = (0.0d0,0.0d0)
    end do
    !$omp end do
    !$omp end parallel
    !
    if(abs(ene(1)-ene(2))<eps)then
      print *, "Caution!!"
      print *, "E(2) and E(1) are degenerate"
    end if
    !
    write(*,'(" ### Calculate Sq_a. ")')
    !calc norm1=<psi|[(Sq^a)^+](Sq^a)|psi> (a=+,-,z)
    if(dim==1 .and. spsmsz==2)then !S^-(q,w) for FM state (must consider later...)2020/8/18-memo
      norm1=1.0d0             !<FM|\sum_{i,j} exp[-Q.(Ri-Rj)] Si^+Sj^-|FM>/N=1
    else
      !calcu Sq^a|GS> and |f1>
      call calcu_Sq_a(spsmsz,psi(1:dim,1),fn(1:dim_new),kvec_calc,dim,dim_new,Spos(1:3,1:NOS))
      norm1=dble(zdotc(dim_new,fn(1),1,fn(1),1))
      !call zdscal(dim_new,1.0d0/sqrt(norm1),fn,1)
      call my_zdscal(dim_new,1.0d0/sqrt(norm1),fn)
    end if

    write(*,'(" ### check <f_1|f_1>=1. ")')
    if((dznrm2(dim_new,fn(1),1)-1.0d0)>eps)then
      print *, "<f1|f1> is not 1"
      stop
    end if

    !, alpha(1), beta(1), and |f_2>=H|f1>-alpha(1)|f_1>-beta(1)|f_0>, beta(2)
    write(*,'(" ### calculate H|f_1>. ")')
    call  ham_to_vec_wave_vector(fnp1(1:dim_new),fn(1:dim_new),dim_new,NOD_new)
    ! 
    write(*,'(" ### consider zeeman term. ")')
    sztot_new=dble(NOS)/2.0d0-dble(NOD_new)
    !
    write(*,'(" ### now treating NOD +- 1 basis for S+-(q, w) or NOD_new=NOD for Sz(q,w). ")')
    !call zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,1,fnp1,1)
    call my_zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,fnp1)
    !
    write(*,'(" ### <f1|H|f1>/<f1|f1>, now <f1|f1>=1. ")')
    alpha_dsf(1)=dble(zdotc(dim_new,fn(1),1,fnp1(1),1))
    beta_dsf(1)=0.0d0
    !
    !output
    ex_i=1
    write(ex_index,'(i1)')ex_i-1
    itr_num=1
    sztot=dble(NOS)/2.0d0-dble(NOD) !present total Sz (not new)
    !
    open(162,file=trim(adjustl(OUTDIR))//'alpha_beta_values_E'//trim(adjustl(ex_index))//'-tmp.d', position='append')
    write(162, '(3I6, 7e23.15, I15)') ex_i, itr_num, itr, alpha_dsf(itr_num), beta_dsf(itr_num), norm1,&
      ene(ex_i)-rfield*sztot, sztot
    close(162)
    !
    if(dim_new > 1)then
      !
      write(*,'(" ### calculate |f_2>. ")')
      !call zaxpy(dim_new,(1.0d0,0.0d0)*(-alpha_dsf(1)),fn,1,fnp1,1)
      call my_zaxpy(dim_new,(1.0d0,0.0d0)*(-alpha_dsf(1)),fn,fnp1)
      !
      write(*,'(" ### beta_dsf(2)=beta(2)**2=<f2|f2>/<f1|f1>, <f1|f1>=1. ")')
      beta_dsf(2)=dble(zdotc(dim_new,fnp1(1),1,fnp1(1),1))
      !                                                    
      write(*,'(" ### normalize |f_2>. ")')
      !call zdscal(dim_new,1.0d0/sqrt(beta_dsf(2)),fnp1,1)
      call my_zdscal(dim_new,1.0d0/sqrt(beta_dsf(2)),fnp1)
      !
      write(*,'(" ### Start iterative calculations ")')
      do itr_num=2, itr_dsf
        !
        !$omp parallel do private(i)
        do i = 1, dim_new
          fnm1(i) = fn(i)
          fn(i) = fnp1(i)
          fnp1(i) = (0.0d0,0.0d0)
        end do
        !
        !calcu <f_n|H|f_n> and calcu H|f_n> for getting |f_n+1>
        call ham_to_vec_wave_vector(fnp1(1:dim_new),fn(1:dim_new),dim_new,NOD_new)
        !
        !consider zeeman term 
        sztot_new=dble(NOS)/2.0d0-dble(NOD_new) 
        !
        !now treating NOD +- 1 basis for S+-(q, w) or NOD_new=NOD for Sz(q,w)
        !call zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,1,fnp1,1)
        call my_zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,fnp1)
        !
        alpha_dsf(itr_num)=dble(zdotc(dim_new,fn(1),1,fnp1(1),1)) 
        !<f_n|H|f_n>/<f_n|f_n>, now |f_n> is normalized.
        !
        !if(mod(itr_num,2)==0)then
        !  call printProgressBar(itr_num,itr_dsf)
        !end if
        !
        !output temporal data
        write(ex_index,'(i1)')ex_i-1
        open(162,file=trim(adjustl(OUTDIR))//'alpha_beta_values_E'//trim(adjustl(ex_index))//'-tmp.d', position='append')
        write(162, '(3I6, 8e23.15)') ex_i, itr_num, itr, alpha_dsf(itr_num), beta_dsf(itr_num), norm1, &
          ene(ex_i)-rfield*sztot, sztot 
        close(162)
        !
        if(itr_num==dim_new) exit
        !
        !calcu |f_n+1>=H|f_n>-alpha(n)*|f_n>-beta(n)*|f_n-1>. Note that beta_dsf is equal to beta**2
        !call zaxpy(dim_new,(1.0d0,0.0d0)*(-alpha_dsf(itr_num)),fn,1,fnp1,1)
        !call zaxpy(dim_new,(1.0d0,0.0d0)*(-sqrt(beta_dsf(itr_num))),fnm1,1,fnp1,1)
        !$omp parallel do
        do i = 1, dim_new
          fnp1(i) = fnp1(i) -alpha_dsf(itr_num)*fn(i) -sqrt(beta_dsf(itr_num))*fnm1(i)
        end do
        !
        !calcu beta(n+1)**2 Note that |f_n> is already normalized
        beta_dsf(itr_num+1)=dble(zdotc(dim_new,fnp1(1),1,fnp1(1),1))  !<f_n+1|f_n+1>/<f_n|f_n>
        !
        !normalize |f_n+1>
        !call zdscal(dim_new,1.0d0/sqrt(beta_dsf(itr_num+1)),fnp1,1)
        call my_zdscal(dim_new,1.0d0/sqrt(beta_dsf(itr_num+1)),fnp1)
        !
      end do
    end if
    !output
    write(sztot_index,'(a,i4.4)')'sztot',int(sztot*10)
    write(ex_index,'(i1)')ex_i-1
    open(162,file=trim(adjustl(OUTDIR))//'alpha_beta_values_'//trim(adjustl(sztot_index))//'_E'//&
      trim(adjustl(ex_index))//'.d', position='append')
    do j=1, itr_dsf
      write(162, '(3I6, 8e23.15)') ex_i, j, itr, alpha_dsf(j), beta_dsf(j), norm1, &
        ene(ex_i)-rfield*sztot, sztot
    end do
    close(162)
    !
  end subroutine calcu_DSF_wavevector
  !
end module get_expectation_values
