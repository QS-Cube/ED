program DSF_analysis
  implicit none
  integer::n,om
  integer::dummy
  integer::ex_i0, itr
  real(8) :: sztot, pi
  integer :: itr_num
  integer, parameter :: dim = 1000
  integer,parameter::omega_num=300000
  real(8),parameter::eta=0.04d0 !for DSF calculation
  real(8),parameter::eps=1.0d-8 !for degeneracy check
  real(8)::alpha(1:dim), beta(1:dim)
  !
  ! alpha(i) and beta(i+1) correspond to 
  ! \alpha_i and \beta_i for Eq.~(23) in Sec.3.9 of arXiv:****.****  
  !
  real(8)::omega, norm1, E0
  real(8)::kakko_real, kakko_image, tmp_real, tmp_image, bunbo
  real(8)::DSF(0:omega_num)
  character(len=5)::itr_name
  character(len=5)::eta_name

  pi = acos(-1.0d0)
  DSF=0.0d0
  write(eta_name,  '(I1)') INT(eta*100)

  if(itr_num<100)then
    write(itr_name, '(I2)') itr_num
  else if(100<=itr_num .and. itr_num<1000)then
    write(itr_name, '(I3)') itr_num
  else if(1000<=itr_num .and. itr_num<10000)then
    write(itr_name, '(I4)') itr_num
  end if

  !calculation of DSF

  alpha=0.0d0
  beta=0.0d0

  open(162, file='alpha_beta_values_E0-tmp.d',status='old')
 
  n = 1
  do
    read(162,*,end=999) ex_i0, dummy, itr, alpha(n), beta(n), norm1, E0, sztot
    n = n + 1
  end do
999 close(162)

  itr_num = n

  !continued-fraction-method
  do om=0, omega_num
    kakko_real=0.0d0
    kakko_image=0.0d0
    bunbo=0.0d0
    tmp_real=0.0d0
    tmp_image=0.0d0

    omega=dble(om)/10000d0

    n = itr_num - 1
    kakko_real=omega+E0-alpha(n)-beta(n+1)
    bunbo=kakko_real**2 + eta**2
    tmp_real=beta(n)*kakko_real/bunbo
    tmp_image=-beta(n)*eta

    do n=itr_num-2, 2,-1
      kakko_real=omega+E0-alpha(n)-tmp_real
      kakko_image=eta-tmp_image
      bunbo=kakko_real**2 + kakko_image**2
      tmp_real=beta(n)*kakko_real/bunbo
      tmp_image=-beta(n)*kakko_image/bunbo
    end do

    n=1
    kakko_real=omega+E0-alpha(n)-tmp_real
    kakko_image=eta-tmp_image
    bunbo=kakko_real**2+kakko_image**2
    tmp_real=norm1*kakko_real/bunbo
    tmp_image=-norm1*kakko_image/bunbo

    DSF(om)=DSF(om)+(-tmp_image/pi)

  end do

  !normalization
  DSF(:)=DSF(:)/dble(1) !modify as you need 2019/6/21

  !output dsf result for each samp
  open(40,file='DSF_itrnum'//trim(itr_name)//'_eta00'//trim(eta_name)//'.dat',position='append')
  do om=0, omega_num
    write(40,'(2e23.15)') dble(om)/10000d0, DSF(om)
  end do
  close(40)

end program DSF_analysis
