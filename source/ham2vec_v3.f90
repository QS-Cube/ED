module ham2vec
  implicit none
contains

  subroutine make_full_hamiltonian(lv,Ham)
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly
    use input_param, only: NOD,NOxxz,p_xxz,list_s,list_r,explist,sJint,NOS
    integer, intent(in)::lv
    complex(8),intent(inout)::Ham(lv,lv)
    integer::s,id
    logical :: f1, f2
    integer :: i, j, ell(3)
    integer, allocatable :: ni(:), st_list(:)
    allocate(ni(NOD),st_list(NOD))
    !$omp parallel do private(i,j,f1,f2,s,id,ell,ni,st_list)
    do i=1, lv
      st_list = list_fly(list_s(i),NOD,NOS)
       do j=1, Noxxz
          f1 = any(st_list == p_xxz(1,j))
          f2 = any(st_list == p_xxz(2,j))
          if(f1 .neqv. f2)then
             Ham(i,i)=Ham(i,i)-sJint(j,3)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list,NOD)
             call representative_SQ(NOD,s,ell,ni)
             call findstate(s,list_s,id,lv)
             if(id > 0) Ham(i,id) = Ham(i,id) + sJint(j,2)*list_r(id)/list_r(i)*&
               explist(ell(1),ell(2),ell(3)) 
          else
             Ham(i,i)=Ham(i,i)+sJint(j,3)
          end if
       end do
    end do
    !
    return
  end subroutine make_full_hamiltonian
  !
  subroutine ham_to_vec_wave_vector(v0,v1,lv,NOD,list_s,list_r,explist) 
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly
    use input_param, only: NOxxz,p_xxz,sJint,LX,LY,LZ,NOS
    integer, intent(in) :: lv
    complex(8), intent(in) :: v1(lv)
    complex(8), intent(out) :: v0(lv)
    integer, intent(in) :: NOD
    real(8), intent(in) :: list_r(1:lv)           
    integer, intent(in)::list_s(1:lv)          
    complex(8),intent(in)::explist(0:LX,0:LY,0:LZ)
    integer :: j, i, ell(3)
    integer, allocatable ::ni(:), st_list(:)
    integer :: s,id
    logical :: f1, f2
    allocate(ni(NOD),st_list(NOD))
    !$omp parallel do private(f1,f2,s,ni,id,ell,i,j,st_list)
    do i = 1, lv
       v0(i) = (0.0d0, 0.0d0)
       st_list = list_fly(list_s(i),NOD,NOS)
       do j = 1, Noxxz
          f1 = any(st_list == p_xxz(1,j))
          f2 = any(st_list == p_xxz(2,j))
          if(f1 .neqv. f2)then
             v0(i) = v0(i) - sJint(j,3) * v1(i)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list,NOD)
             call representative_SQ(NOD,s,ell,ni)
             call findstate(s,list_s,id,lv)
             if(id > 0) v0(i) = v0(i) + sJint(j,2) * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3)) * v1(id)
          else
             v0(i) = v0(i) + sJint(j,3) * v1(i)
          end if
       end do
    end do
    return
  end subroutine ham_to_vec_wave_vector
  !
  subroutine calcu_lm_trans_2(v1,lv,NOD,expe_sz,NOS,list_n,non)
    use state_lists, only: list_fly
    use input_param, only: LX,LY,LZ,list_s
    integer, intent(in) :: lv, NOD, NOS, non
    complex(8), intent(in) :: v1(lv)
    integer, intent(in) :: list_n(non)
    real(8), intent(inout) :: expe_sz(NOS)
    integer :: j, i
    logical :: f1
    real(8) :: expval
    integer, allocatable :: st_list(:)
    allocate(st_list(NOD))
    !
    expval = 0.0d0
    !$omp parallel do private(f1,i,j,st_list) reduction(+:expval)
    do i = 1, lv
       st_list = list_fly(list_s(i),NOD,NOS)
       do j = 1, non
          f1 = any(st_list == list_n(j))
          if(f1)then
             expval = expval - dble( v1(i)*conjg(v1(i)) )
          else
             expval = expval + dble( v1(i)*conjg(v1(i)) )
          end if
       end do
    end do
    !
    expval = expval / non * 0.5
    !
    !$omp parallel do
    do i = 1, non
      expe_sz(list_n(i)) = expval
    end do
    !
    return
  end subroutine calcu_lm_trans_2
  !
  subroutine calcu_cf_trans_2(v1,lv,NOD,list_s,list_r,explist,NOS,list_ij,noij,&
    expe_szsz,expe_spsm)
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly
    use input_param, only: NOxxz,LX,LY,LZ
    integer, intent(in) :: lv, NOS, NOD, noij
    complex(8), intent(in) :: v1(lv)
    integer, intent(in) :: list_ij(2,noij)
    real(8), intent(in) :: list_r(1:lv)
    integer, intent(in)::list_s(1:lv)
    complex(8),intent(in)::explist(0:LX,0:LY,0:LZ)
    real(8), intent(inout) :: expe_szsz(NOS**2)
    complex(8), intent(inout) :: expe_spsm(NOS**2)
    integer :: j, i, ell(3)
    integer, allocatable ::ni(:), st_list(:)
    integer :: s,id
    logical :: f1, f2
    real(8) :: expval_szsz
    complex(8) :: expval_spsm
    allocate(ni(NOD),st_list(NOD))
    expval_szsz=0.0d0
    expval_spsm=(0.0d0,0.0d0)
    if(list_ij(1,1)==list_ij(2,1))then
      expval_szsz = 0.25d0
      !$omp parallel do private(f1,i,j,st_list) reduction(+:expval_spsm)
      do i = 1, lv
       st_list = list_fly(list_s(i),NOD,NOS)
        do j = 1, noij
          f1 = any(st_list == list_ij(1,j))
          if(f1) cycle
          expval_spsm = expval_spsm + dble( v1(i)*conjg(v1(i)) )
        end do
      end do
    else
      !$omp parallel do private(f1,f2,s,ni,id,ell,i,j,st_list) reduction(+:expval_szsz,expval_spsm)
      do i = 1, lv
        st_list = list_fly(list_s(i),NOD,NOS)
        do j = 1, noij
          f1 = any(st_list == list_ij(1,j))
          f2 = any(st_list == list_ij(2,j))
          if((.not. f1) .and. f2)then
            expval_szsz = expval_szsz - dble(conjg(v1(i))*v1(i))
          else if(f1 .and. (.not. f2))then
            expval_szsz = expval_szsz - dble(conjg(v1(i))*v1(i))
            ni = j_flip_ni(minval(list_ij(:,j)),maxval(list_ij(:,j)),st_list,NOD)
            call representative_SQ(NOD,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0) expval_spsm = expval_spsm + list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3)) * conjg(v1(i))*v1(id)
          else
            expval_szsz = expval_szsz + dble(conjg(v1(i))*v1(i))
          end if
        end do
      end do
      !
      expval_szsz = expval_szsz / noij *0.25
    end if
    expval_spsm = expval_spsm / noij
    !$omp parallel do
    do i = 1, noij
      expe_szsz(list_ij(1,i)+(list_ij(2,i)-1)*NOS) = expval_szsz
      expe_spsm(list_ij(1,i)+(list_ij(2,i)-1)*NOS) = expval_spsm
    end do
    return
  end subroutine calcu_cf_trans_2
  !
end module ham2vec
