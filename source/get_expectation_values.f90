module get_expectation_values
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
  subroutine representative_lm(ell,n) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(inout)::n(1) 
    integer,intent(out) :: ell(3)
    integer :: i,j,k
    integer :: nd(1)
    ell = 0
    nd = n
    do k = 1, LZ
      nd = shift_z_SQ(nd)
      do j = 1, LY
        nd = shift_y_SQ(nd)
        do i = 1, LX
          nd = shift_x_SQ(nd)
          !
          if(nd(1)<n(1))then
            ell=(/i,j,k/)
            n = nd
          end if
          !
        end do
      end do
    end do
  end subroutine representative_lm
  !
  subroutine mk_list_n(list_n,non,NOS,n) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)  :: n, NOS 
    integer,intent(out) :: list_n(NOS), non
    integer :: i,j,k
    integer :: nd
    non = 1
    list_n(non) = n
    !
    nd = n
    do k = 1, LZ
      nd = shift_z_SQ(nd)
      do j = 1, LY
        nd = shift_y_SQ(nd)
        do i = 1, LX
          nd = shift_x_SQ(nd)
          !
          if(nd<n)then
            non = 0
            return
          else if(nd>n)then
            non = non + 1
            list_n(non) = nd
          else
            return
          end if
          !
        end do
      end do
    end do
  end subroutine mk_list_n
  !
  subroutine mk_list_ij(list_ij,noij,NOS,r) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)  :: r(2), NOS 
    integer,intent(out) :: list_ij(2,NOS), noij
    integer :: i,j,k,l
    integer :: rd(2)
    noij = 1
    list_ij(:,noij) = r(:)
    !
    rd = r
    do k = 1, LZ
      rd = shift_z_SQ(rd)
      do j = 1, LY
        rd = shift_y_SQ(rd)
        do i = 1, LX
          rd = shift_x_SQ(rd)
          !
          do l = 2, 1, -1
            if(rd(l)<r(l))then
              noij = 0
              return
            else if(rd(l)>r(l))then
              noij = noij + 1
              list_ij(:,noij) = rd
              exit
!            else 
!              return
            end if
          end do
          if(l==0) return !*** modified 2021/02/03
          !
        end do
      end do
    end do
  end subroutine mk_list_ij
  !
  subroutine representative_szsz(ell,n) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(inout)::n(2) 
    integer,intent(out) :: ell(3)
    integer :: i,j,k,l
    integer :: nd(2), itmp
    ell = 0
    nd = n
    do k = 1, LZ
      nd = shift_z_SQ(nd)
      do j = 1, LY
        nd = shift_y_SQ(nd)
        do i = 1, LX
          nd = shift_x_SQ(nd)
          if(nd(2) < nd(1))then
            itmp=nd(2); nd(2)=nd(1); nd(1)=itmp
          end if
          !
          ! list_expe_ss(1,j)/(2,j) is the upper/lower index
          do l = 1, 2
            if(nd(l)>n(l))then
              exit
            else if(nd(l)<n(l))then
              ell=(/i,j,k/)
              n(1:l) = nd(1:l)
              exit
            end if
          end do
          !
        end do
      end do
    end do
  end subroutine representative_szsz
  !
  subroutine representative_spsm(ell,n)
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(inout)::n(2) 
    integer,intent(out) :: ell(3)
    integer :: i,j,k,l
    integer :: nd(2)
    ell = 0
    nd = n
    do k = 1, LZ
      nd = shift_z_SQ(nd)
      do j = 1, LY
        nd = shift_y_SQ(nd)
        do i = 1, LX
          nd = shift_x_SQ(nd)
          !
          do l = 2, 1, -1
            if(nd(l)>n(l))then
              exit
            else if(nd(l)<n(l))then
              ell=(/i,j,k/)
              n(1:l) = nd(1:l)
              exit
            end if
          end do
          !
        end do
      end do
    end do
  end subroutine representative_spsm
  !
  subroutine get_lm_2_wave_vector(psi,NOD,dim,nvec,NOS)
    use input_param, only: OUTDIR
    use ham2vec, only: calcu_lm_trans_2
    integer, intent(in) :: NOD, nvec, NOS
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1:nvec)
    integer :: j, k
    integer :: non, list_n(NOS), no_tot
    !
    expe_sz=0.0d0    
    do k = 1, nvec
      no_tot = 0
      do j = 1, NOS
        if(no_tot == NOS)then
          exit
        end if
        call mk_list_n(list_n,non,NOS,j)
        if(non==0)then
          cycle
        else
          no_tot = no_tot + non
          call calcu_lm_trans_2(psi(1,k),dim,NOD,expe_sz(1,k),NOS,list_n,non)
        end if
      end do
    end do

    write(*,'(" ### write ouput/local_mag.dat. ")')
    open(10,file=trim(adjustl(OUTDIR))//'local_mag.dat',position='append')
    do j = 1, NOS
      write(10,'(i8,10000es23.15)') j, ( expe_sz(j,k), k=1, nvec )
    end do
    close(10)

    return

  end subroutine get_lm_2_wave_vector

  subroutine get_cf_2_wave_vector(psi,NOCF,NOD,dim,nvec,list_r,list_s,NOS)
    use input_param, only: explist, OUTDIR
    use ham2vec, only: calcu_cf_trans_2
    integer, intent(in) :: NOCF,NOD,nvec,NOS
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1:nvec)
    real(8), intent(in) :: list_r(1:dim)
    integer, intent(in)::list_s(1:dim)
    integer :: j, k
    integer :: jmax, list_ij(2,NOS), noij, no_tot

    jmax=NOS**2
    expe_szsz=0.0d0
    expe_spsm=0.0d0

    do k = 1, nvec
      no_tot=0
      do j = 1, jmax
        if(no_tot==jmax) exit
        call mk_list_ij(list_ij,noij,NOS,list_expe_ss(:,j))
        if(noij==0)then
          cycle
        else
          no_tot = no_tot + noij
          call calcu_cf_trans_2(psi(1,k),dim,NOD,list_s,list_r,explist,NOS,list_ij,noij,&
            expe_szsz(1,k),expe_spsm(1,k))
        end if
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

  subroutine calcu_Sq_a(spsmsz,psi,fn,kvec_calc,rkx_new,rky_new,rkz_new,dim,dim_new,Spos)
    use input_param, only: list_r,list_r_new,list_s_new,NOS,NOD,NOD_new,list_s,explist_new
    use state_lists, only: search_representative_Spi_state_SQ, search_representative_Smi_state_SQ, &
      search_representative_Szi_state_SQ,list_fly
    implicit none
    integer,intent(in)::dim,dim_new
    integer, intent(in)::spsmsz
    real(8),intent(in)::kvec_calc(1:3),Spos(1:3,1:NOS),rkx_new,rky_new,rkz_new
    complex(8),intent(in)::psi(1:dim)  
    complex(8),intent(inout)::fn(1:dim_new)
    logical :: fi
    integer:: i,l,ell(3)
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
        !$omp parallel do private(l,ell,lp,fi,st_list)
        do l=1, dim
          ell=0
          lp=0
          st_list = list_fly(list_s(l),NOD,NOS)
          fi=any(st_list==i)
          if(fi)then !If i-site has down spin
            call search_representative_Spi_state_SQ(i,NOD,NOD_new,dim_new,&
              list_s_new(1:dim_new),st_list,lp,ell)
            if(lp>0)then 
              !$omp atomic
              fn(lp)=fn(lp)+ psi(l)*list_r_new(lp)/list_r(l) &
                * exp_list_KS(i) * conjg( explist_new(ell(1),ell(2),ell(3)) )
              !
            end if
          end if
        end do
        !$omp end parallel do
      end do
      !
    else if(spsmsz==2)then
      do i=1, NOS
        !$omp parallel do private(l,ell,lp,fi,st_list)
        do l=1, dim
          ell=0
          lp=0  
          st_list = list_fly(list_s(l),NOD,NOS)
          fi=any(st_list==i)
          if(.not. fi)then !If i-site has up spin
            call search_representative_Smi_state_SQ(i,NOD,NOD_new,dim_new,&
              list_s_new(1:dim_new),st_list,lp,ell)  
            if(lp>0)then 
              !$omp atomic
              fn(lp)=fn(lp)+ psi(l)*list_r_new(lp)/list_r(l) &
                * exp_list_KS(i) * conjg( explist_new(ell(1),ell(2),ell(3)) )
              !
            end if
          end if
        end do
        !$omp end parallel do 
      end do
      !
    else if(spsmsz==3)then
      do i=1, NOS
        !$omp parallel do private(l,ell,lp,fi,Siz,st_list)
        do l=1, dim
          ell=0
          lp=0  
          st_list = list_fly(list_s(l),NOD,NOS)
          fi=any(st_list==i)
          if(fi)then !If i-site has up spin
            Siz=-0.50d0
          else
            Siz=0.50d0
          end if
          call search_representative_Szi_state_SQ(NOD,NOD_new,dim_new,&
            list_s_new(1:dim_new),st_list,lp,ell)  
          if(lp>0)then
            fn(lp)=fn(lp)+ Siz*psi(l)*list_r_new(lp)/list_r(l) &
              * exp_list_KS(i) * conjg( explist_new(ell(1),ell(2),ell(3)) )
          end if
        end do
        !$omp end parallel do
      end do
    end if
    !
    call zdscal(dim_new,1.0d0/sqrt(dble(NOS)),fn,1)
    !
  end subroutine calcu_Sq_a

  subroutine calcu_DSF_wavevector(psi,ene,kvec_calc,rfield,NOS,NOD,itr_dsf,dim) 
    use input_param, only:LX,LY,LZ,spsmsz,qx,qy,qz,rkx,rky,rkz,rkx_new,rky_new,rkz_new,NOD_new,&
      list_s_new,list_r_new,THS_new,explist_new,OUTDIR
    use state_lists, only: combination,allocate_lists_omp_SQ_new
    use ham2vec, only:  ham_to_vec_wave_vector
    integer, intent(in) :: dim
    complex(8), intent(in) :: psi(1:dim,1)
    real(8),intent(in):: ene(1:2),rfield
    integer, intent(in) :: NOS, NOD,itr_dsf

    !for DSF
    integer::dim_new
    integer::ex_i,i,j,k,itr=1,itr_num
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
    write(*,'(" ### Set rk_a_new with a \in {x,y,z}. ")')
    rkx_new = rkx+qx
    rky_new = rky+qy
    rkz_new = rkz+qz
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
    call allocate_lists_omp_SQ_new      !THS_new, list_r_new,list_s_new are updated
    dim_new = THS_new
    write(*,*) "  dim_new   = ", dim_new
    !
    write(*,'(" ### Set phase factors for intermediate states. ")')
    allocate(explist_new(0:LX,0:LY,0:LZ))
    forall(i=0:LX,j=0:LY,k=0:LZ) explist_new(i,j,k) = &
      exp((0.0d0,1.0d0)*(rkx_new*dble(i)+rky_new*dble(j)+rkz_new*dble(k)))
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
      call calcu_Sq_a(spsmsz,psi(1:dim,1),fn(1:dim_new),&
        kvec_calc,rkx_new,rky_new,rkz_new,dim,dim_new,Spos(1:3,1:NOS)) 
      norm1=dble(zdotc(dim_new,fn(1),1,fn(1),1))
      call zdscal(dim_new,1.0d0/sqrt(norm1),fn,1)
    end if

    write(*,'(" ### check <f_1|f_1>=1. ")')
    if((dznrm2(dim_new,fn(1),1)-1.0d0)>eps)then
      print *, "<f1|f1> is not 1"
      stop
    end if

    !, alpha(1), beta(1), and |f_2>=H|f1>-alpha(1)|f_1>-beta(1)|f_0>, beta(2)
    write(*,'(" ### calculate H|f_1>. ")')
    call  ham_to_vec_wave_vector(fnp1(1:dim_new),fn(1:dim_new),dim_new,NOD_new,&
      list_s_new,list_r_new, explist_new)
    ! 
    write(*,'(" ### consider zeeman term. ")')
    sztot_new=dble(NOS)/2.0d0-dble(NOD_new)
    !
    write(*,'(" ### now treating NOD +- 1 basis for S+-(q, w) or NOD_new=NOD for Sz(q,w). ")')
    call zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,1,fnp1,1)
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
      call zaxpy(dim_new,(1.0d0,0.0d0)*(-alpha_dsf(1)),fn,1,fnp1,1)
      !
      write(*,'(" ### beta_dsf(2)=beta(2)**2=<f2|f2>/<f1|f1>, <f1|f1>=1. ")')
      beta_dsf(2)=dble(zdotc(dim_new,fnp1(1),1,fnp1(1),1))
      !                                                    
      write(*,'(" ### normalize |f_2>. ")')
      call zdscal(dim_new,1.0d0/sqrt(beta_dsf(2)),fnp1,1)
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
        call ham_to_vec_wave_vector(fnp1(1:dim_new),fn(1:dim_new),dim_new,NOD_new,&
          list_s_new,list_r_new,explist_new)
        !
        !consider zeeman term 
        sztot_new=dble(NOS)/2.0d0-dble(NOD_new) 
        !
        !now treating NOD +- 1 basis for S+-(q, w) or NOD_new=NOD for Sz(q,w)
        call zaxpy(dim_new,(1.0d0,0.0d0)*(-rfield*sztot_new),fn,1,fnp1,1)
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
        call zaxpy(dim_new,(1.0d0,0.0d0)*(-alpha_dsf(itr_num)),fn,1,fnp1,1)
        call zaxpy(dim_new,(1.0d0,0.0d0)*(-sqrt(beta_dsf(itr_num))),fnm1,1,fnp1,1)
        !
        !calcu beta(n+1)**2 Note that |f_n> is already normalized
        beta_dsf(itr_num+1)=dble(zdotc(dim_new,fnp1(1),1,fnp1(1),1))  !<f_n+1|f_n+1>/<f_n|f_n>
        !
        !normalize |f_n+1>
        call zdscal(dim_new,1.0d0/sqrt(beta_dsf(itr_num+1)),fnp1,1)
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
