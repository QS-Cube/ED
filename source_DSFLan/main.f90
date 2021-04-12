 program DSF_analysis
 implicit none
 integer::i,j,k,l,m,n,om
 integer::dummy
 integer::seed, samp, tsz, tot_samp
 integer::ex_i,ex_i0, itr
 real(8) :: sztot, pi
 integer,parameter::Nsize=64
 integer,parameter::lattice_type=2
 integer,parameter::seed_num=1
 integer,parameter::samp_num=1
 integer,parameter::tot_samp_num=seed_num*samp_num
 integer,parameter::omega_num=300000
 integer,parameter::itr_num=400
 double precision,parameter::eta=0.04d0 !for DSF calculation
 double precision,parameter::eps=1.0d-8 !for degeneracy check
 double precision::dum1, dum2, dum3, dum4
 double precision::E_tsz(0:3)
 integer::dge_gs=0, tsz_gs=0
 integer::dge_tsz00, dge_tsz01 
 double precision::alpha(1:itr_num), beta(1:itr_num)
 double precision::omega, norm1, E0
 double precision::ren_real, ren_image
 double precision::kakko_real, kakko_image, tmp_real, tmp_image, bunbo
 double precision::DSF(0:omega_num,1:tot_samp_num)
 double precision::ave_DSF(0:omega_num)
 double precision::err_DSF(0:omega_num)
 character(len=5)::size_name
 character(len=5)::type_name
 character(len=5)::seed_name
 character(len=5)::samp_name
 character(len=5)::tsz_name
 character(len=5)::ex_index
 character(len=5)::itr_name
 character(len=5)::tot_samp_name
 character(len=5)::eta_name

 pi = acos(-1.0d0)
 DSF=0.0d0
 tot_samp=1
 write(eta_name,  '(I1)') INT(eta*100)
 write(size_name, '(I2)') Nsize
 write(type_name, '(I1)') lattice_type

 if(itr_num<100)then
 write(itr_name, '(I2)') itr_num
 else if(100<=itr_num .and. itr_num<1000)then
 write(itr_name, '(I3)') itr_num
 else if(1000<=itr_num .and. itr_num<10000)then
 write(itr_name, '(I4)') itr_num
 end if


 !degeneracy check
 !You should input by hand... 2019/6/21-memo
 !tsz_gs=0
 !dge_tsz00=2
 !dge_tsz01=0
 !dge_gs=2
!  if(samp<=9)then
!  write(samp_name, '(I1)') samp
!  open(161, file='seed0'//trim(seed_name)//'/ene_levelsdel00samp000'//trim(samp_name)//'full_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  open(161, file='ene_levelsfull_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  else if(10<=samp.and.samp<=99)then
!  write(samp_name, '(I2)') samp
!  !open(161, file='seed0'//trim(seed_name)//'/ene_levelsdel00samp00'//trim(samp_name)//'full_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  open(161, file='ene_levelsfull_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  else if(100<=samp.and.samp<=999)then
!  write(samp_name, '(I3)') samp
!  !open(161, file='seed0'//trim(seed_name)//'/ene_levelsdel00samp0'//trim(samp_name)//'full_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  open(161, file='ene_levelsfull_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'.d')
!  end if

!  do l=1, 4
!  read(161,*) dum1, E_tsz(l-1), dum2, dum3, dum4                                                    
!  end do
!  close(161)

!  if(abs(E_tsz(0)-E_tsz(1))>eps)then
!  print *, "seed=", seed, "samp=", samp, "gs is singlet"
!  dge_gs=1 
!  tsz_gs=0
!  else if(abs(E_tsz(0)-E_tsz(1))<eps .and. abs(E_tsz(0)-E_tsz(2))>eps )then
!  print *, "seed=", seed, "samp=", samp, "gs is triplet"
!  dge_gs=3
!  tsz_gs=1
!  else if(abs(E_tsz(0)-E_tsz(1))<eps .and. abs(E_tsz(0)-E_tsz(2))<eps .and. abs(E_tsz(0)-E_tsz(3))>eps )then
!  print *, "seed=", seed, "samp=", samp, "gs is quintet"
!  dge_gs=5
!  tsz_gs=2
!  else if(abs(E_tsz(0)-E_tsz(1))<eps .and. abs(E_tsz(0)-E_tsz(2))<eps .and. abs(E_tsz(0)-E_tsz(3))<eps)then
!  print *, "seed=", seed, "samp=", samp, "degerate of gs is larger than septet. You must recalc.."
!  dge_gs=7
!  tsz_gs=3
!  end if

 ! open(21, file='degeneracy-check.dat',position='append')
 ! write(21,'(3I4)') tot_samp, tsz_gs, dge_gs
 ! close(21)


 !calculation of DSF
  do 200 tsz=0, 0 !degeneracy between total Sz spaces !modify as you need!!! 2019/6/21
   E0=0.0d0
   norm1=0.0d0
   alpha=0.0d0
   beta=0.0d0
   write(tsz_name,'(I1)') tsz


   do ex_i=1, 1  !degeneracy in total Sz space !modify as you need!!! 2019/6/21 
   write(ex_index,'(I1)') ex_i-1
   E0=0.0d0
   norm1=0.0d0
   alpha=0.0d0
   beta=0.0d0

  ! if(samp<=9)then
   open(162, file='alpha_beta_values_E0-tmp.d',status='old')
  ! else if(10<=samp.and.samp<=99)then
  ! open(162, file='alpha_beta_values_sztot00'//trim(tsz_name)//'0_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'_E'//trim(ex_index)//'.d',status='old')
  ! else if(100<=samp.and.samp<=999)then
  ! open(162, file='alpha_beta_values_sztot00'//trim(tsz_name)//'0_rand_2DQHSAF_q_lar_N0'//trim(size_name)//'0'//trim(type_name)//'_E'//trim(ex_index)//'.d',status='old')
  ! end if

   do n=1, itr_num
   read(162,*,end=999) ex_i0, dummy, itr, alpha(n), beta(n), norm1, E0, sztot
   end do
999 close(162)

   !continued-fraction-method
   do 300 om=0, omega_num
    kakko_real=0.0d0
    kakko_image=0.0d0
    bunbo=0.0d0
    tmp_real=0.0d0
    tmp_image=0.0d0

    omega=dble(om)/10000d0

    n=itr_num-1
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

    !print *,tmp_real, tmp_image 
    n=1
    kakko_real=omega+E0-alpha(n)-tmp_real
    kakko_image=eta-tmp_image
    bunbo=kakko_real**2+kakko_image**2
    tmp_real=norm1*kakko_real/bunbo
    tmp_image=-norm1*kakko_image/bunbo

    !print *, tmp_image 
    !DSF(om,tot_samp)=-tmp_image/pi
    if(tsz==0)then
    DSF(om,tot_samp)=DSF(om,tot_samp)+(-tmp_image/pi)
    else if(tsz>0)then
    DSF(om,tot_samp)=DSF(om,tot_samp)+(-tmp_image/pi)*2.0d0 !+tsz and -tsz
    end if
 
  300 end do !om do
   end do !ex_i do

 200 end do  !tsz do

 !normalization
 DSF(:,tot_samp)=DSF(:,tot_samp)/dble(1) !modify as you need 2019/6/21

 !output dsf result for each samp
 open(40,file='DSF_itrnum'//trim(itr_name)//'_eta00'//trim(eta_name)//'.dat',position='append')
 do om=0, omega_num
 write(40,'(2e23.15,2I15)') dble(om)/10000d0, DSF(om,tot_samp), tsz_gs, dge_gs
 end do
 close(40)


 !sampling average and error
! ave_DSF=0.0d0
! err_DSF=0.0d0

! do samp=1, tot_samp_num
! ave_DSF(:)=ave_DSF(:)+DSF(:,samp)
! end do
! ave_DSF(:)=ave_DSF(:)/dble(tot_samp_num)

! do samp=1, tot_samp_num
! err_DSF(:)=err_DSF(:)+(DSF(:,samp)-ave_DSF(:))**2
! end do
! err_DSF(:)=sqrt(err_DSF(:)/(dble(tot_samp_num)*dble(tot_samp_num-1)))

 !output result
!  open(50,file='DSF_itrnum'//trim(itr_name)//'_seed0'//trim(seed_name)//'_all_sample_eta00'//trim(eta_name)//'.dat',position='append')
! do om=0, omega_num
! write(50,'(3e)') dble(om)/10000d0, ave_DSF(om), err_DSF(om)
! end do
! close(50)

end program
