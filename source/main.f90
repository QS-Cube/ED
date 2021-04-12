!******************************************************************************************************
! 2020/12/23 coded by Hiroshi Ueda (Osaka Univ.) & Tokuro Shimokawa (OIST)
!    This numerical program can be applied to S=1/2 XXZ spin systems 
! on regular lattices near saturation.
!
! 2021/01/02 major update of get_lm_2_wave_vector, get_cf_2_wave_vector, and their branches.
! 2021/01/03 major update of list_fly and do not use st_list & st_list_new
!
!******************************************************************************************************
program main
  !$ use omp_lib
  use state_lists, only: get_combination,combination,mk_shift_x_SQ, mk_shift_y_SQ, mk_shift_z_SQ &
    ,allocate_lists_omp_SQ
  use input_param, only: read_ip, THS, NOS, NOD, NOV, LX, LY, LZ, shift_x_SQ, shift_y_SQ, &
    shift_z_SQ,ene,psi,ALG,cal_lm,cal_cf,cal_dsf, &
    NOLM,NOCF,FILElm,FILECF,FILEpos,list_r,list_s,qx,qy,qz, &
    kvec_calc,rfield,itr_dsf,wr_wf, write_wf,re_wf,read_wf,OUTDIR
  use eigen_solver, only: my_trlanczos_routines_complex,Full_diag_routines,calcu_FM_energy,&
    i_vec_min, i_vec_max,read_input_TRLan,NOE
  use lanczos, only: read_lanczos_para, lanczos_routines
  use get_expectation_values, only: allocate_expe_mem,&
    calcu_DSF_wavevector,get_lm_2_wave_vector, get_cf_2_wave_vector
  implicit none
  real(8) :: eneFM
  complex(8) :: prdct
  integer :: i8
  integer :: k
  !$ real(8) :: st, en
  !$ st = omp_get_wtime()
  write(*,'("********************************* Start QS^3 ***************************************")')
  call read_ip
  !
  write(*,'(" ### Store combination numbers. ")')
  call get_combination(NOS,NOD+1)
  !
  write(*,'(" ### Set THS. ")')
  if(NOD>0)then
     THS = combination(NOS,NOD)
  else
     THS = 1
  end if
  write(*,*) "  THS   = ", THS
  if(THS < 0) stop "THS is bigger than the upper limit of 8-bite integer!"
  !
  if(THS.eq.1)then
    write(*,'(" ### Focus on the FM state. ")')
    eneFM=0.0d0
    call calcu_FM_energy(eneFM)
    open(10,file=trim(adjustl(OUTDIR))//'energy.dat',position='append')
    write(10,*) eneFM
    close(10)
    write(*,*) "We are calculating eneFM only"
    write(*,*) "eneFM=", eneFM
    stop
  end if
  !
  write(*,'(" ### Set translational operations. ")')
  shift_x_SQ = mk_shift_x_SQ(LX,LY,LZ) 
  shift_y_SQ = mk_shift_y_SQ(LX,LY,LZ) 
  shift_z_SQ = mk_shift_z_SQ(LX,LY,LZ)
  !
  write(*,'(" ### Allocate and Set arrays for state_list. ")')
  call allocate_lists_omp_SQ 
  write(*,*) "  THS(k)   = ", THS

  If(ALG.eq.1)then !Conventional Lanczos
    if(re_wf.ne.1)then
      write(*,'(" ### Start the Lanczos method. ")')
      call lanczos_routines(ene,THS,psi)
      !
      if(wr_wf.eq.1)then
        write(*,'(" ### Write wavevectors. ")')
        call write_wf(THS,1,1,1,ene)
      end if
    else
      write(*,'(" ### Read lanczos_para. ")')
      call read_lanczos_para
      !
      write(*,'(" ### Read wavefunctions. ")')
      call read_wf(THS,1,1,1)
    end if
    !
    if(cal_lm==1 .or. cal_cf==1)then
      write(*,'(" ### Allocate arrays for expectation values. ")')
      call allocate_expe_mem(NOS,1,NOLM,NOCF,FILElm,FILECF,FILEpos) 
    end if
    !
    if(cal_lm==1)then
      write(*,'(" ### Get local magnetizations. ")')
      call get_lm_2_wave_vector(psi,NOD,THS,1,NOS)
    end if
    !
    if(cal_cf==1)then
      write(*,'(" ### Get correlation functions. ")')
      call get_cf_2_wave_vector(psi,NOCF,NOD,THS,1,list_r,list_s,NOS)
    end if
    !
  else if(ALG.eq.2 )then 
    if(re_wf.ne.1)then
      write(*,'(" ### Start the thick-restarted Lanczos method. ")')
      call my_trlanczos_routines_complex(ene,psi,THS) 
      !
      if(wr_wf.eq.1)then
        write(*,'(" ### Write wavevectors. ")')
        call write_wf(THS,i_vec_min, i_vec_max,NOE,ene)
      end if
      !
    else
      write(*,'(" ### Read input_TRLan. ")')
      call read_input_TRLan
      !
      write(*,'(" ### Read wavefunctions. ")')
      call read_wf(THS,i_vec_min,i_vec_max,NOE)
    end if
    !
    write(*,'(" ### Check orthogonality. ")')
    prdct=(0.0d0,0.0d0)
    !$omp parallel do private(i8) reduction(+:prdct) 
    do i8 = 1, THS
      prdct = prdct + conjg(psi(i8,1))*psi(i8,2)
    end do
    print *, prdct
    !
    write(*,'(" ### Check normalizations. ")')
    do k=1, NOV
      prdct=(0.0d0,0.0d0)
      !$omp parallel do private(i8) reduction(+:prdct) 
      do i8=1, THS
        prdct=prdct+ conjg(psi(i8,k))*psi(i8,k)
      end do
      print *, k, prdct
      psi(1:THS,k)=psi(1:THS,k)/sqrt(abs(prdct))
    end do
    !
    write(*,'(" ### Check normalizations again. ")')
    do k=1, NOV
      prdct=(0.0d0,0.0d0)
      !$omp parallel do private(i8) reduction(+:prdct)
      do i8=1, THS
        prdct=prdct+ conjg(psi(i8,k))*psi(i8,k)
      end do
      print *, k, prdct
    end do
    !
    if(cal_lm==1 .or. cal_cf==1)then
      write(*,'(" ### Allocate arrays for expectation values. ")')
      call allocate_expe_mem(NOS,NOV,NOLM,NOCF,FILElm,FILECF,FILEpos) 
    end if
    !
    if(cal_lm==1)then
      write(*,'(" ### Get local magnetizations. ")')
      call get_lm_2_wave_vector(psi,NOD,THS,NOV,NOS)
    end if
    !
    if(cal_cf==1)then
      write(*,'(" ### Get correlation functions. ")')
      call get_cf_2_wave_vector(psi,NOCF,NOD,THS,NOV,list_r,list_s,NOS)
    end if
    !
  else if(ALG.eq.3)then
    write(*,'(" ### Start the full diagonalization. ")')
    call Full_diag_routines(ene,psi,THS)
    !
    if(cal_lm==1 .or. cal_cf==1)then
      write(*,'(" ### Allocate arrays for expectation values. ")')
      call allocate_expe_mem(NOS,THS,NOLM,NOCF,FILElm,FILECF,FILEpos) 
    end if
    !
    if(cal_lm==1)then
      write(*,'(" ### Get local magnetizations. ")')
      call get_lm_2_wave_vector(psi,NOD,THS,THS,NOS)
    end if
    !
    if(cal_cf==1)then
      write(*,'(" ### Get correlation functions. ")')
      call get_cf_2_wave_vector(psi,NOCF,NOD,THS,THS,list_r,list_s,NOS)
    end if
  end if
  !
  if(cal_dsf==1)then
    write(*,'(" ### Set kvec_calc. ")')
    kvec_calc(1)=qx
    kvec_calc(2)=qy
    kvec_calc(3)=qz
    !
    write(*,'(" ### Calculate DSF. ")')
    call calcu_DSF_wavevector(psi,ene(1:2),kvec_calc,rfield,NOS,NOD,itr_dsf,THS) 
  end if
  !
  !$ en = omp_get_wtime()
  !$ write(*,'("Time:",f10.3," [sec]")') en-st
  write(*,'("****************************************************************************************")')
  stop
end program main

