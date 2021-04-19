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

  function f_combination(n,k) result(nCk)
    integer, intent(in) :: n, k
    integer :: i
    integer ::nCk

    if(n < k)then
      nCk = 0
      return
    end if

    nCk=1
    do i=1,k
      nCk=nCk*(n-k+i)
      nCk=nCk/i
    end do

  end function f_combination

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

  function list_fly(t,NOD,NOS) result(ni)
    integer, intent(in) :: t, NOD, NOS
    integer :: s
    integer :: i, j, b, j0
    integer :: ni(NOD) 
    s = t
    j = NOS - 1
    do i = NOD, 2, -1
       call binary_search(s,combination(:,i), i, j, b, j0)
       j = j0 - 1
       ni(i) = j0
       s = s - combination(j,i)
    end do
    ni(1) = s
  end function list_fly
  !                                                                                                                                                                                                                
  subroutine binary_search(s,list_s,ls,le,b,bmin)
    integer, intent(out) :: b, bmin
    integer, intent(in) :: s,ls,le
    integer, intent(in) :: list_s(le)
    integer :: bmax
    bmin = ls; bmax = le
    do
       b = bmin + (bmax-bmin)/2
       if(s<list_s(b))then
          bmax = b-1
       else if(list_s(b)<s)then
          bmin = b+1
       else
          bmin = b
          return
       end if
       if(bmin>bmax)then
          b = -1
          return
       end if
    end do
  end subroutine binary_search

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
  !
  subroutine search_Spi_state_SQ(i,NOD,NOD_new,n,lp) 
    integer,intent(in)::i,NOD,NOD_new 
    integer,intent(in)::n(1:NOD)
    integer,intent(inout)::lp
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

    lp = inv_list(bp,NOD_new)

  end subroutine search_Spi_state_SQ

  subroutine search_Smi_state_SQ(i,NOD,NOD_new,n,lp) 
    integer,intent(in)::i,NOD,NOD_new 
    integer,intent(in)::n(1:NOD)
    integer,intent(inout)::lp
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

    lp = inv_list(bp,NOD_new)
    !
  end subroutine search_Smi_state_SQ
  !
end module state_lists

