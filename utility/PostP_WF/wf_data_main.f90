PROGRAM wf_data_main

  use parallel_module
  use atom_module
  use wf_data_module

  implicit none

  integer :: n1,n2,i,j

  call start_mpi_parallel

  call read_parameters
  call read_atom(myrank,970)

  nlst=0
  do i=1,3,2
     n1 = lst_wf_tmp(i)
     n2 = lst_wf_tmp(i+1)
     if ( 0 < n1 .and. n1 <= n2 ) then
        nlst = nlst + n2-n1+1
     end if
  end do
  if ( nlst == 0 ) stop "stop@main(1)"
  if ( myrank == 0 ) then
     write(*,*) "nlst=",nlst
  end if

  allocate( lst_wf(nlst) ) ; lst_wf=0
  nlst=0
  do i=1,3,2
     n1 = lst_wf_tmp(i)
     n2 = lst_wf_tmp(i+1)
     if ( 0 < n1 .and. n1 <= n2 ) then
        do j=n1,n2
           nlst=nlst+1
           lst_wf(nlst)=j
        end do
     end if
  end do

  call init_parallel(myrank==0,Ngrid(1),Nband,Nbzsm,Nspin)

  call wf_data

  call end_mpi_parallel

END PROGRAM wf_data_main
