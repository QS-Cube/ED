module state_lists
  implicit none
  integer, allocatable :: combination(:,:)
contains

  subroutine get_combination(NOS,NOD)
    integer, intent(in) :: NOS, NOD
    integer :: i, j
    allocate(combination(NOS,NOD))
    do j = 1, NOD
      do i = 1, NOS
        combination(i,j) = f_combination(i,j)
      end do
    end do
  end subroutine get_combination

  integer function f_combination(n,k) 
    integer, intent(in) :: n, k
    integer :: i
    integer ::nCk

    if(n < k)then
      f_combination = 0
      return
    end if

    nCk=1
    do i=1,k
      nCk=nCk*(n-k+i)
      nCk=nCk/i
    end do

    f_combination = nCk
  end function f_combination

  !   !for Lx*Ly*Lz lattice (x=[1,LX], y=[1,Ly], z=[1,Lz]) <=> x + LX*(y-1) + LX*LY*(z-1)
  function mk_shift_x_SQ(LX,LY,LZ) result(shift_x_SQ) 
    integer, intent(in) :: LX, LY, LZ
    integer :: shift_x_SQ(LX*LY*LZ)
    integer :: x, y, z
    forall(x=1:LX,y=1:LY,z=1:LZ)
      shift_x_SQ( x + LX*((y-1) + LY*(z-1)) ) = &
        (mod(x, LX)+1) + LX*(y-1) + LX*LY*(z-1)
    end forall
  end function mk_shift_x_SQ

  function mk_shift_y_SQ(LX,LY,LZ) result(shift_y_SQ) 
    integer, intent(in) :: LX, LY, LZ
    integer :: shift_y_SQ(LX*LY*LZ)
    integer :: x, y, z
    forall(x=1:LX,y=1:LY,z=1:LZ)
      shift_y_SQ( x + LX*((y-1) + LY*(z-1)) ) = &
        x + LX*mod(y, LY) + LX*LY*(z-1)
    end forall
  end function mk_shift_y_SQ

  function mk_shift_z_SQ(LX,LY,LZ) result(shift_z_SQ) 
    integer, intent(in) :: LX, LY, LZ
    integer :: shift_z_SQ(LX*LY*LZ)
    integer :: x, y, z
    forall(x=1:LX,y=1:LY,z=1:LZ)
      shift_z_SQ( x + LX*((y-1) + LY*(z-1)) ) = &
        x + LX*(y-1) + LX*LY*mod(z, LZ)
    end forall
  end function mk_shift_z_SQ

  subroutine checkstate_SQ(s,ro,NOD,pkx,pky,pkz)
    use input_param, only: LX,LY,LZ,NOS,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer, intent(in) :: s
    integer, intent(in) :: NOD
    real(8), intent(in) :: pkx,pky,pkz
    real(8), intent(out) :: ro
    complex(8) :: coef
    integer :: i, j, k, l, count
    integer, allocatable :: n(:), np(:)
    real(8) :: p1, p2
    allocate(n(NOD),np(NOD))
    n = list_fly(s,NOD,NOS)
    np = n
    !
    ro = 0.0d0
    count = 0
    coef = (0.0d0,0.0d0)
    !
    do k = 1, LZ
      np = shift_z_SQ(np)
      p1 = pkz * dble(k)
      do j = 1, LY
        np = shift_y_SQ(np)
        p2 = pky * dble(j) + p1
        do i = 1, LX
          np = shift_x_SQ(np)
          call insertion_sort(np,NOD)
          !
          do l = NOD, 1, -1
            if(np(l) < n(l))then
              return
            else if(np(l) > n(l))then
              exit
            end if
          end do
          if(l==0)then
            coef = coef + exp( (0.0d0,-1.0d0) * ( pkx * dble(i) + p2) )
            count = count + 1
          end if
          !
        end do
      end do
    end do
    !
    ro = (abs(coef)**2) *dble(NOS)/dble(count)
    !NOS/count='Ra', where, T^{Ra}|s>=|s>, and  ro='Na', see the eq. (118) in the Sandvik ref.
  end subroutine checkstate_SQ


  pure subroutine insertion_sort(a,NOD) 
    integer, intent(in) :: NOD
    integer, intent(inout) :: a(NOD)
    integer :: temp, i, j
    do i = 2, NOD
      j = i - 1
      temp = a(i)
      do while (a(j)>temp) 
        a(j+1) = a(j)
        j = j - 1
        if(j==0) exit
      end do
      a(j+1) = temp
    end do
  end subroutine insertion_sort

  integer function inv_list(ni,NOD) 
    integer, intent(in) :: NOD
    integer, intent(in) :: ni(NOD)     
    integer :: i
    inv_list = ni(1)
    do i = 2, NOD
      inv_list = inv_list + combination(ni(i)-1,i)
    end do
  end function inv_list

  recursive subroutine qsort_w_order(a, o, first, last)
    integer :: a(*), o(*), first, last
    integer :: x, t8, i, j
    integer :: t
    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
        i=i+1
      end do
      do while (x < a(j))
        j=j-1
      end do
      if (i >= j) exit
      t8 = a(i);  a(i) = a(j);  a(j) = t8
      t  = o(i);  o(i) = o(j);  o(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i - 1) call qsort_w_order(a, o, first, i - 1)
    if (j + 1 < last)  call qsort_w_order(a, o, j + 1, last)
  end subroutine qsort_w_order

  subroutine calcu_lm_trans(NOD,psi_l,list_r_l,site,n0,pkx,pky,pkz,lm_tmp) !see Notability
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)::NOD,site,n0(1:NOD)
    complex(8),intent(in)::psi_l
    real(8),intent(in)::list_r_l
    real(8),intent(in)::pkx,pky,pkz
    integer, allocatable :: n(:),np(:)
    integer :: r_x, r_y, r_z, rp_x, rp_y, rp_z
    integer :: lm_int
    complex(8),intent(out) :: lm_tmp
    real(8) :: p1, p2, c

    allocate(n(1:NOD),np(1:NOD))

    lm_int = 0
    lm_tmp = (0.0d0,0.0d0)
    c = dble(conjg(psi_l)*psi_l*0.50d0/(list_r_l**2))

    n = n0 !n0 means |a>_l
    do r_z = 1, LZ
      n = shift_z_SQ(n)
      do r_y = 1, LY
        n = shift_y_SQ(n)
        do r_x = 1, LX
          n = shift_x_SQ(n)
          call insertion_sort(n,NOD)
          !We obtain now |b>=T^r|a>_l

          !Look <S_site^z> in |b>
          lm_int=0
          if(any(n(1:NOD) == site))then !down
            lm_int=-1
          else!up
            lm_int=1
          end if

          !search T^r'|a>_l=|b>
          np = n0
          do rp_z = 1, LZ
            np = shift_z_SQ(np)
            p1 = pkz*dble(rp_z-r_z)
            do rp_y=1,LY
              np = shift_y_SQ(np)
              p2 = pky*dble(rp_y-r_y) + p1
              do rp_x=1,LX
                np = shift_x_SQ(np)
                call insertion_sort(np,NOD)
                !We obtain now T^r'|a>
                if( all(n==np) )then
                  lm_tmp=lm_tmp+c*dble(lm_int)*exp((0.0d0,1.0d0)*(pkx*dble(rp_x-r_x) + p2))
                end if
              end do
            end do
          end do
          !
        end do
      end do
    end do
  end subroutine calcu_lm_trans

  subroutine calcu_szsz_trans(NOD,psi_l,list_r_l,site1,site2,n0,pkx,pky,pkz,szsz_tmp) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)::NOD,site1,site2,n0(1:NOD)
    complex(8),intent(in)::psi_l
    real(8),intent(in)::list_r_l
    real(8),intent(in)::pkx,pky,pkz
    integer, allocatable :: n(:), np(:)
    integer ::r_x, r_y, r_z, rp_x, rp_y, rp_z
    integer::szsz_int
    complex(8),intent(out)::szsz_tmp
    logical ::f1, f2
    real(8) :: p1, p2, c

    allocate(n(1:NOD),np(1:NOD))

    szsz_int=0
    szsz_tmp=(0.0d0,0.0d0)
    c = dble( conjg(psi_l)*psi_l*0.25d0/(list_r_l**2) )
    n = n0 !n0 is |a>_l

    do r_z = 1, LZ
      n = shift_z_SQ(n)
      do r_y = 1, LY
        n = shift_y_SQ(n)
        do r_x = 1, LX
          n = shift_x_SQ(n)
          call insertion_sort(n,NOD)
          !We obtain now |b>=T^r|a>_l

          !Look <S_site1^zS_site2^z> in |b> state
          f1 = any(n(1:NOD)==site1)
          f2 = any(n(1:NOD)==site2)
          if(f1 .neqv. f2)then !up-down or down-up
            szsz_int =  - 1 
          else !up-up or down-down
            szsz_int =  1
          end if

          !search T^r'|a>_l=|b>
          np = n0
          do rp_z=1,LZ
            np = shift_z_SQ(np)
            p1 = pkz*dble(rp_z-r_z)
            do rp_y=1,LY
              np = shift_y_SQ(np)
              p2 = pky*dble(rp_y-r_y) + p1
              do rp_x=1,LX
                np = shift_x_SQ(np)
                call insertion_sort(np,NOD)
                !We obtain now T^r'|a>
                if( all(n==np) )then
                  szsz_tmp=szsz_tmp + c * dble(szsz_int) * &
                    exp((0.0d0,1.0d0)*(pkx*dble(rp_x-r_x) + p2))
                end if
              end do
            end do
          end do

        end do
      end do
    end do

  end subroutine calcu_szsz_trans

  subroutine calcu_spsm_trans(l,NOD,lv,psi,list_r,list_s,st_list,site1,site2,n0,pkx,pky,pkz,spsm_tmp) 
    !(calculate <Si^+Sj^->=<SxSx>+<SySy>)                                                         
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)::l,lv
    integer,intent(in)::NOD,site1,site2,n0(1:NOD) !lv=THS=dim
    complex(8),intent(in)::psi(1:lv)
    real(8),intent(in)::list_r(1:lv)
    integer,intent(in)::list_s(1:lv)
    integer, intent(in) :: st_list(1:NOD,1:lv)
    real(8),intent(in)::pkx,pky,pkz
    integer :: r_x, r_y, r_z, rp_x, rp_y, rp_z
    complex(8),intent(out)::spsm_tmp
    logical :: f1, f2
    integer :: s,id
    integer :: ell(3)
    integer, allocatable :: n(:), np(:), ni(:),ni2(:)
    real(8) :: p1, p2
    complex(8) :: c

    allocate(ni(NOD),ni2(NOD),n(NOD),np(NOD))

    spsm_tmp=(0.0d0,0.0d0)
    n = n0 !n0 is |a>_l

    do r_z = 1, LZ
      n = shift_z_SQ(n)

      do r_y = 1, LY
        n = shift_y_SQ(n)

        do r_x = 1, LX
          n = shift_x_SQ(n)
          call insertion_sort(n,NOD)
          id=0
          f1 = any(n(1:NOD)==site1) !S+
          f2 = any(n(1:NOD)==site2) !S-

          if(site1 /= site2)then  
            if(f1 .and. (.not. f2))then 
              !down-up state or up-down (we will calculate <Si^+Sj^->=<SxSx>+<SySy>)
              ni(1:NOD) = j_flip_ni(site1,site2,n(1:NOD),NOD) !flip site1, site2-spins
              ni2(1:NOD)=ni(1:NOD)
              !
              !We set ni2l here just for knowing "id". (ell will not used later)
              call representative_SQ(NOD,s,ell,ni2) 
              !Here, ni2 is changed to representative state |a>_id, that means ni2 is "not" |b>.
              call findstate(s,list_s,id,lv)
              if(id>0)then
                c = conjg(psi(id))*psi(l)/list_r(l)/list_r(id)
                ! We obtain now |b>=SpSmT^r|a>_l. Note that this |b> state belongs to 
                ! |a>_id (l/=id). |b>=T^[-ell(1),-ell(2)] |a>_id
                call insertion_sort(ni,NOD)

                !search T^r'|a>_id=|b>
                np = st_list(1:NOD,id)!n0(1:NOD)
                do rp_z=1, LZ
                  np = shift_z_SQ(np)
                  p1 = pkz*dble(rp_z-r_z)
                  do rp_y=1,LY
                    np = shift_y_SQ(np)
                    p2 = pky*dble(rp_y-r_y) + p1
                    do rp_x=1,LX
                      np = shift_x_SQ(np)
                      call insertion_sort(np,NOD)
                      !We obtain now T^r'|a>
                      if( all(ni==np) )then
                        spsm_tmp=spsm_tmp + c * exp((0.0d0,1.0d0)*(pkx*dble(rp_x-r_x) + p2))
                      end if
                    end do
                  end do
                end do

              end if
            end if

          else if(site1 == site2)then
            ! site1=site2 and this site has up spin !2020/7/30-memo(We will not calculate, being zero)
            !   ni(1:NOD)=n(1:NOD)
            ! call representative_SQ(NOD,s,ell,ni)
            ! call findstate(s,list_s,id,lv)
            ! if(id>0)then
            !  spsm_tmp=spsm_tmp+exp((0.0d0,-1.0d0)*(pkx*dble(ell(1)+i) + pky*dble(ell(2)+j)))* &
            ! & conjg(psi(id))*psi(i8)/list_r(id)/list_r(i8)
            ! end if
          end if

        end do
      end do
    end do

  end subroutine calcu_spsm_trans

  subroutine allocate_lists_omp_SQ 
    use input_param, only: NOD,rkx,rky,rkz,THS,NOS,list_s,list_r,st_list
    use omp_lib
    integer :: a, i, j
    real(8) :: r
    integer :: num = 1
    real(8), allocatable :: tmp_list_r(:)
    integer, allocatable :: tmp_st_list(:,:), order(:)
    num = omp_get_max_threads()
    write(*,*) "********************"
    write(*,*) "max_threads", num
    write(*,*) "********************"
    !
    write(*,'(" ### Count # of representative states. ")')
    a = 0
    !$omp parallel
    !$omp do private(j,i,r)
    do j = 1, num
      do i = j, THS, num
        call checkstate_SQ(i,r,NOD,rkx,rky,rkz) 
        if(r>1.0d-15) then
          !$omp atomic
          a = a + 1
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    !
    write(*,'(" ### Allocate work arrays for lists. ")')
    allocate(list_s(a),tmp_list_r(a),tmp_st_list(NOD,a))
    !
    write(*,'(" ### Store representative states. ")')
    a = 0
    !$omp parallel
    !$omp do private(j,i,r)
    do j = 1, num
      do i = j, THS, num
        call checkstate_SQ(i,r,NOD,rkx,rky,rkz) 
        if(r>1.0d-15) then
          !$omp critical
          a = a + 1
          list_s(a) = i
          tmp_list_r(a) = r
          tmp_st_list(:,a) = list_fly(i,NOD,NOS)
          !$omp end critical
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    !
    write(*,'(" ### Allocate arrays for lists. ")')
    THS = a
    allocate(order(THS),list_r(THS),st_list(NOD,THS))
    !
    write(*,'(" ### Store reordered lists. ")')
    !$omp parallel do
    do i = 1, THS
      order(i) = i
    end do
    call qsort_w_order(list_s, order, 1, THS)
    !$omp parallel do private(i)
    do i = 1, THS
      list_r(i) = sqrt( tmp_list_r(order(i)) )
      st_list(:,i) = tmp_st_list(:,order(i))
    end do
    !$omp end parallel do 
    return
  end subroutine allocate_lists_omp_SQ

  function list_fly(t,NOD,NOS) result(ni)
    integer, intent(in) :: t, NOD, NOS
    integer :: s, q
    integer :: i, j
    integer :: ni(NOD) 
    s = t
    j = NOS - 1
    do i = NOD, 2, -1
      do
        q = combination(j,i)
        if(q >= s)then
          j = j-1
        else
          ni(i) = j+1
          s = s - q
          exit
        end if
      end do
    end do
    ni(1) = s
  end function list_fly

  function j_flip_ni(i,j,n,NOD) result(nd)
    integer, intent(in) :: NOD, i, j, n(NOD)
    integer :: kr, kl
    integer :: nd(NOD)
    nd = 1
    do kr = NOD, 1, -1
      if(j < n(kr)) then
        cycle
      else if(j > n(kr)) then
        exit
      else
        nd = 0; exit
      end if
    end do
    if(nd(NOD) == 1)then ! S+_i S-_j
      do kl = 1, kr
        if(i == n(kl)) exit
      end do
      nd(kl:kr-1) = n(kl+1:kr)
      nd(kr) = j
    else ! S-_i S+_j
      do kl = 1, kr
        if(i < n(kl)) exit
      end do
      nd(kl) = i
      nd(kl+1:kr) = n(kl:kr-1)
    end if
    nd(1:kl-1) = n(1:kl-1)
    nd(kr+1:NOD) = n(kr+1:NOD)
  end function j_flip_ni

  subroutine findstate(s,list_s,b,m) 
    integer, intent(out) :: b
    integer, intent(in) :: s, m, list_s(m)
    integer :: bmin, bmax 
    bmin = 1; bmax = m
    do
      b = bmin + (bmax-bmin)/2
      if(s < list_s(b))then
        bmax = b-1
      else if(s > list_s(b))then
        bmin = b+1
      else
        return
      end if
      if(bmin>bmax)then
        b = -1; return
      end if
    end do
  end subroutine findstate

  subroutine representative_SQ(NOD,t,ell,n) 
    use input_param, only: LX,LY,LZ,shift_x_SQ,shift_y_SQ,shift_z_SQ
    integer,intent(in)::NOD 
    integer,intent(inout)::n(NOD) 
    integer,intent(out) :: t, ell(3)
    integer :: i,j,k,l
    integer, allocatable :: nd(:)
    allocate(nd(NOD))
    ell = 0
    nd = n
    do k = 1, LZ
      nd(:) = shift_z_SQ(nd(:))
      do j = 1, LY
        nd(:) = shift_y_SQ(nd(:))
        do i = 1, LX
          nd(:) = shift_x_SQ(nd(:))
          call insertion_sort(nd,NOD)
          !
          do l = NOD, 1, -1
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
    t = inv_list(n,NOD)
  end subroutine representative_SQ
  !
  subroutine allocate_lists_omp_SQ_new 
    use input_param, only: list_s_new,list_r_new,st_list_new,NOD_new,THS_new,NOS,&
      rkx_new,rky_new,rkz_new
    !$use omp_lib
    integer :: i, dim, a
    real(8) :: r
    integer :: num = 1, j
    real(8), allocatable :: tmp_list_r(:)
    integer, allocatable :: tmp_st_list(:,:), order(:)
    integer ::THS0 
    !$num = omp_get_max_threads()

    a = 0
    !$omp parallel
    !$omp do private(j,i,r)
    do j = 1, num
      do i = j, THS_new, num
        call checkstate_SQ(i,r,NOD_new,rkx_new,rky_new,rkz_new)
        if(r>1.0d-15) then
          !$omp critical
          a = a + 1
          !$omp end critical
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    dim = a           !update dim 
    allocate(list_s_new(dim),tmp_list_r(dim),tmp_st_list(NOD_new,dim))
    a = 0
    !$omp parallel
    !$omp do private(j,i,r)
    do j = 1, num
      do i = j, THS_new, num
        call checkstate_SQ(i,r,NOD_new,rkx_new,rky_new,rkz_new)
        if(r>1.0d-15) then
          !$omp critical
          a = a + 1
          if(a > dim) stop "dim is too small!"
          list_s_new(a) = i
          tmp_list_r(a) = r
          tmp_st_list(:,a) = list_fly(i,NOD_new,NOS)
          !$omp end critical
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    THS_new = a             !update THS_new
    allocate(order(THS_new),list_r_new(THS_new),st_list_new(NOD_new,THS_new))
    !$omp parallel do
    do i = 1, THS_new
      order(i) = i
    end do
    THS0=THS_new
    call qsort_w_order(list_s_new, order, 1, THS0)       
    !$omp parallel do
    do i = 1, THS_new
      list_r_new(i) = sqrt( tmp_list_r(order(i)) )
      st_list_new(:,i) = tmp_st_list(:,order(i))
    end do
    !
    return
  end subroutine allocate_lists_omp_SQ_new


  subroutine search_representative_Spi_state_SQ(i,NOD,NOD_new,dim_new,list_s_new,n,lp,ell) 
    integer,intent(in)::i,NOD,NOD_new 
    integer,intent(in)::dim_new
    integer,intent(in)::list_s_new(1:dim_new)
    integer,intent(in)::n(1:NOD)
    integer,intent(inout) :: ell(3)
    integer,intent(inout)::lp
    integer :: s
    integer::l,kr
    integer, allocatable :: bp(:)
    allocate(bp(1:NOD_new))

    !i-siteスピンがダウンかどうか、何番目のダウンスピン(l=kr)かを調べる。
    if(NOD>1)then !not FM state
      l = 1
      do kr = 1, NOD
        if(n(kr)==i)then
          l=kr
          exit  
        end if
      end do

      !新しいスピン状態|b'>
      bp(1:l-1)=n(1:l-1)
      !bp(l:NOD_new)=n(l+1:NOD)!when l=NOD then this is not operated
      bp(l:NOD-1)=n(l+1:NOD)

    else if (NOD==1.or.NOD==0)then
      print *, "You can not calculate...when sztot=NOS/2-1..." 
      stop
    end if

    call representative_SQ(NOD_new,s,ell,bp) !Here, bp is changed to representative state |b>_lp, that means ni2 is "not" |b>.
    call findstate(s,list_s_new,lp,dim_new)

  end subroutine search_representative_Spi_state_SQ

  subroutine search_representative_Smi_state_SQ(i,NOD,NOD_new,dim_new,list_s_new,n,lp,ell) 
    integer,intent(in)::i,NOD,NOD_new 
    integer,intent(in)::dim_new
    integer,intent(in)::list_s_new(1:dim_new)
    integer,intent(in)::n(1:NOD)
    integer,intent(inout) :: ell(3)
    integer,intent(inout)::lp
    integer :: s
    integer::l,kr
    integer, allocatable :: bp(:)
    allocate(bp(1:NOD_new))

    if(NOD/=0)then !not FM state
      if(i<n(1))then
        !新しいスピン状態|b'>
        bp(1)=i
        !bp(2:NOD_new)=n(1:NOD)
        bp(2:NOD+1)=n(1:NOD)
      else if(i>n(NOD))then
        !新しいスピン状態|b'>
        bp(1:NOD)=n(1:NOD)
        bp(NOD_new)=i
      else
        l = 1
        do kr = 1, NOD-1
          if(n(kr)<i)then
            l=kr  
          else if(n(kr)>i)then
            exit
          end if
        end do

        !新しいスピン状態|b'>
        bp(1:l)=n(1:l)
        bp(l+1)=i
        bp(l+2:NOD+1)=n(l+1:NOD)
      end if
    else if (NOD==0)then
      !新しいスピン状態|b'>
      bp(1)=i
    end if

    call representative_SQ(NOD_new,s,ell,bp) 
    !Here, bp is changed to representative state |b>_lp, that means ni2 is "not" |b>.
    call findstate(s,list_s_new,lp,dim_new)
    !
  end subroutine search_representative_Smi_state_SQ

  subroutine search_representative_Szi_state_SQ(NOD,NOD_new,dim_new,list_s_new,n,lp,ell) 
    integer,intent(in)::NOD,NOD_new 
    integer,intent(in)::dim_new
    integer,intent(in)::list_s_new(1:dim_new)
    integer,intent(in)::n(1:NOD)
    integer,intent(inout) :: ell(3)
    integer,intent(inout)::lp
    integer :: s
    integer, allocatable ::bp(:)
    allocate(bp(1:NOD_new))

    bp(1:NOD)=n(1:NOD)

    call representative_SQ(NOD_new,s,ell,bp) 
    !Here, bp is changed to representative state |b>_lp, that means ni2 is "not" |b>.
    call findstate(s,list_s_new,lp,dim_new)
  end subroutine search_representative_Szi_state_SQ

  !
end module state_lists

