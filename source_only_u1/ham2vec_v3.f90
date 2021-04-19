module ham2vec
  implicit none
contains

  subroutine make_full_hamiltonian(lv,Ham)
    use state_lists, only: j_flip_ni, list_fly, inv_list
    use input_param, only: NOD,NOxxz,p_xxz,sJint,NOS
    integer, intent(in)::lv
    complex(8),intent(inout)::Ham(lv,lv)
    integer :: id
    logical :: f1, f2
    integer :: i, j
    integer, allocatable :: ni(:), st_list(:)
    allocate(ni(NOD),st_list(NOD))
    !$omp parallel do private(i,j,f1,f2,id,ni,st_list)
    do i=1, lv
      st_list = list_fly(i,NOD,NOS)
       do j=1, Noxxz
          f1 = any(st_list == p_xxz(1,j))
          f2 = any(st_list == p_xxz(2,j))
          if(f1 .neqv. f2)then
             Ham(i,i)=Ham(i,i)-sJint(j,3)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list,NOD)
             id = inv_list(ni,NOD)
             Ham(i,id) = Ham(i,id) + sJint(j,2) 
          else
            Ham(i,i)=Ham(i,i)+sJint(j,3)
          end if
       end do
    end do
    !
    return
  end subroutine make_full_hamiltonian
  !
  subroutine ham_to_vec_wave_vector(v0,v1,lv,NOD) 
    use state_lists, only: j_flip_ni,list_fly,inv_list
    use input_param, only: NOxxz,p_xxz,sJint,NOS
    integer, intent(in) :: lv
    complex(8), intent(in) :: v1(lv)
    complex(8), intent(out) :: v0(lv)
    integer, intent(in) :: NOD
    integer :: j, i
    integer, allocatable ::ni(:), st_list(:)
    integer :: id
    logical :: f1, f2
    allocate(ni(NOD),st_list(NOD))
    !$omp parallel do private(f1,f2,ni,id,i,j,st_list)
    do i = 1, lv
       v0(i) = (0.0d0, 0.0d0)
       st_list = list_fly(i,NOD,NOS)
       do j = 1, Noxxz
          f1 = any(st_list == p_xxz(1,j))
          f2 = any(st_list == p_xxz(2,j))
          if(f1 .neqv. f2)then
             v0(i) = v0(i) - sJint(j,3) * v1(i)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list,NOD)
             id = inv_list(ni,NOD)
             v0(i) = v0(i) + sJint(j,2) * v1(id)
          else
             v0(i) = v0(i) + sJint(j,3) * v1(i)
          end if
       end do
    end do
    return
  end subroutine ham_to_vec_wave_vector
  !
  subroutine calcu_lm(v1,lv,NOD,NOS,expe_sz,list_n)
    use state_lists, only: list_fly
    integer, intent(in) :: lv, NOD, NOS
    complex(8), intent(in) :: v1(lv)
    integer, intent(in) :: list_n
    real(8), intent(inout) :: expe_sz
    integer :: i
    logical :: f1
    real(8) :: expval
    integer, allocatable :: st_list(:)
    allocate(st_list(NOD))
    !
    expval = 0.0d0
    !$omp parallel do private(f1,i,st_list) reduction(+:expval)
    do i = 1, lv
       st_list = list_fly(i,NOD,NOS)
       f1 = any(st_list == list_n)
       if(f1)then
         expval = expval - dble( v1(i)*conjg(v1(i)) )
       else
         expval = expval + dble( v1(i)*conjg(v1(i)) )
       end if
    end do
    !
    expe_sz = expval * 0.5
    !
    return
  end subroutine calcu_lm
  !
  subroutine calcu_cf(v1,lv,NOD,NOS,list_ij,expval_szsz,expval_spsm)
    use state_lists, only: j_flip_ni, list_fly, inv_list
    use input_param, only: NOxxz
    integer, intent(in) :: lv, NOS, NOD
    complex(8), intent(in) :: v1(lv)
    integer, intent(in) :: list_ij(2)
    real(8), intent(inout) :: expval_szsz
    complex(8), intent(inout) :: expval_spsm
    integer :: i
    integer, allocatable ::ni(:), st_list(:)
    integer :: id
    logical :: f1, f2
    allocate(ni(NOD),st_list(NOD))
    expval_szsz=0.0d0
    expval_spsm=(0.0d0,0.0d0)
    if(list_ij(1)==list_ij(2))then
      expval_szsz = 0.25d0
      !$omp parallel do private(f1,i,st_list) reduction(+:expval_spsm)
      do i = 1, lv
        st_list = list_fly(i,NOD,NOS)
        f1 = any(st_list == list_ij(1))
        if(f1) cycle
        expval_spsm = expval_spsm + dble( v1(i)*conjg(v1(i)) )
      end do
    else
      !$omp parallel do private(f1,f2,ni,id,i,st_list) reduction(+:expval_szsz,expval_spsm)
      do i = 1, lv
        st_list = list_fly(i,NOD,NOS)
        f1 = any(st_list == list_ij(1))
        f2 = any(st_list == list_ij(2))
        if((.not. f1) .and. f2)then
          expval_szsz = expval_szsz - dble(conjg(v1(i))*v1(i))
        else if(f1 .and. (.not. f2))then
          expval_szsz = expval_szsz - dble(conjg(v1(i))*v1(i))
          ni = j_flip_ni(minval(list_ij(:)),maxval(list_ij(:)),st_list,NOD)
          id = inv_list(ni,NOD)
          expval_spsm = expval_spsm + conjg(v1(i))*v1(id)
        else
          expval_szsz = expval_szsz + dble(conjg(v1(i))*v1(i))
        end if
      end do
      !
      expval_szsz = expval_szsz *0.25
    end if
    return
  end subroutine calcu_cf
  !
end module ham2vec
