module ham2vec
  implicit none
contains

  subroutine make_full_hamiltonian(lv,Ham)
    use state_lists, only: j_flip_ni, representative_SQ, findstate
    use input_param, only: NOD,NOxxz,st_list,p_xxz,list_s,list_r,explist,Jint
    integer, intent(in)::lv
    complex(8),intent(inout)::Ham(lv,lv)
    integer::s,id
    logical :: f1, f2
    integer :: i, j, ell(3)
    integer, allocatable :: ni(:)
    allocate(ni(NOD))
    Jint(:,2) = Jint(:,2) * 0.50d0
    Jint(:,3) = Jint(:,3) * 0.25d0

    !$omp parallel do private(i,j,f1,f2,s,id,ell,ni)
    do i=1, lv
       do j=1, Noxxz
          f1 = any(st_list(:,i) == p_xxz(1,j))
          f2 = any(st_list(:,i) == p_xxz(2,j))
          if(f1 .neqv. f2)then
             Ham(i,i)=Ham(i,i)-Jint(j,3)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list(:,i),NOD)
             call representative_SQ(NOD,s,ell,ni)
             call findstate(s,list_s,id,lv)
             if(id > 0) Ham(i,id) = Ham(i,id) + Jint(j,2)*list_r(id)/list_r(i)*&
               explist(ell(1),ell(2),ell(3)) 
          else
             Ham(i,i)=Ham(i,i)+Jint(j,3)
          end if
       end do
    end do
    !
    Jint(:,2) = Jint(:,2) * 2.0d0
    Jint(:,3) = Jint(:,3) * 4.0d0
    !
    return
  end subroutine make_full_hamiltonian
  !
  subroutine ham_to_vec_wave_vector(v0,v1,lv,NOD,list_s,list_r,st_list,explist) 
    use state_lists, only: j_flip_ni, representative_SQ, findstate
    use input_param, only: NOxxz,p_xxz,Jint,LX,LY,LZ
    integer, intent(in) :: lv
    complex(8), intent(in) :: v1(lv)
    complex(8), intent(out) :: v0(lv)
    integer, intent(in) :: NOD                    
    integer, intent(in) :: st_list(1:NOD,1:lv)    
    real(8), intent(in) :: list_r(1:lv)           
    integer, intent(in)::list_s(1:lv)          
    complex(8),intent(in)::explist(0:LX,0:LY,0:LZ)
    integer :: j, i, l, ell(3)
    integer, allocatable ::ni(:)
    integer :: s,id
    logical :: f1, f2
    allocate(ni(NOD))
    Jint(:,2) = Jint(:,2) * 0.50d0
    Jint(:,3) = Jint(:,3) * 0.25d0
    !$omp parallel do private(f1,f2,s,ni,id,ell,i,j,l)
    do i = 1, lv
       v0(i) = (0.0d0, 0.0d0)
       do j = 1, Noxxz
          f1 = any(st_list(:,i) == p_xxz(1,j))
          f2 = any(st_list(:,i) == p_xxz(2,j))
          if(f1 .neqv. f2)then
             v0(i) = v0(i) - Jint(j,3) * v1(i)
             ni = j_flip_ni(p_xxz(1,j),p_xxz(2,j),st_list(:,i),NOD)
             call representative_SQ(NOD,s,ell,ni)
             call findstate(s,list_s,id,lv)
             if(id > 0) v0(i) = v0(i) + Jint(j,2) * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3)) * v1(id)
          else
             v0(i) = v0(i) + Jint(j,3) * v1(i)
          end if
       end do
    end do
    Jint(:,2) = Jint(:,2) * 2.0d0
    Jint(:,3) = Jint(:,3) * 4.0d0
    return
  end subroutine ham_to_vec_wave_vector
  !
end module ham2vec
