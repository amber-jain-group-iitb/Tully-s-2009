program tulli12
use atom_interaction 
implicit none
integer :: io,jo,ko,total_trajectories,number_of_cores,ntraj,TT,t,loop,q_c
real*8 :: x1,dx1,start,finish
real*8, dimension(:),allocatable :: population_matrix
real*8, dimension(:),allocatable :: Time_matrix
!call setup_initial_values1
open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)
allocate(Time_matrix(TT))
allocate(population_matrix(TT))
do jo=1,TT
population_matrix(TT)=0
enddo
do io=1,ntraj
   call setup_initial_values2

   call CPU_TIME(start)
   call classical_evolution 


  ! do ko=1,TT
   !   population_matrix(k)=population_matrix(k)+population_mat(k)
  ! enddo
   call CPU_time(finish)
   write(119,*)finish-start     
   call vibrational_energy
   
   open(34,file="quantum_count")
   do q_c=1,16
   write(34,*) quantum_count(q_c)
   enddo
   close(34)                                                                                                                                                                                                                              
enddo

do t=1,TT
   Time_matrix(t)=10*t
enddo

write(77,*)   0.0000000000,     1.000000000
do loop=1,TT
    write(77,*) Time_matrix(loop),1-population_matrix(loop)/ntraj
enddo

end program
!..............................................................................


