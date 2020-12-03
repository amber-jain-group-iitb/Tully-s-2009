module atom_interaction
        implicit none
        integer::i,j,k,l,Au_nei_list(528,12),step_dynamics,step_equilibrium,step,vib_state,Au_i(528),n_O_Au,n_N_Au,n_Au,traj
        real*8,dimension(396,3)::acc_Au,acc_Au_old,Au_vel
        real*8,dimension(528,3)::Au,Au_eq
        real*8,dimension(3)::N,O,N_vel,O_vel,COM,COM_v,acc_N,acc_N_old,acc_O,acc_O_old
        real*8::r_Au_O,r_Au_N
        integer,dimension(528)::N_Au_list,O_Au_list
        real*8::r_N_O

        real*8::a0,alpha0,b0,beta0,f0,gamma0,r0_N_O,a1,alpha1,b1,beta1,r1_Au_N,d,c1,zimage,f1,gamma1,r1_N_O,phi,ea
        real*8::a2,a3,gamma2,gamma3,b2,b3,r_cut

        real*8::mass_N,mass_O,red_mass,av,kb,temp,conv,wavenumber_2_j,amu2kg,c_light,wavenumber,mass_Au
        real*8::e_kin,pi,planck_const,total_mass

        real*8,dimension(528,528,3,3)::ten_mat
        real*8::l_x,r_x,l_y,r_y,l_z,r_z,box_len_x,box_len_y,box_len_z,force_const,dt1,dt,z_COM
        integer::iflag_write,iflag_writevmd,adsorbed,no_traj

        real*8::diab_pot(2,2),gr_vec(1,2),gr_vec_t(2,1),ex_vec(1,2),ex_vec_t(2,1),gr_eig,ex_eig
        real*8::der_pot_N(2,2,3),der_pot_O(2,2,3),der_pot_Au(528,3,2,2),pot_Au_Au(1,1),der_pot_Au_Au(528,3)
        !.............................................................................................................. 

        integer::nsize,INFO,total_dimensions,Hi,Ne,nold,number_of_atoms,movable_layer!(be careful with setting dt)
        real*8 :: tau,planks_constant,kT,Vr,Energy_of_diabat,dtc,total_time
        complex*16 :: population
        real*8 :: time,rnd,omega,gh,dG,Band_width,RAND
        complex*16,dimension(:,:), allocatable :: c,b,A
        real*8,dimension(:,:), allocatable :: Energy_hamil,g,Identity,H,Grd,vdij
        real*8,dimension(:,:,:), allocatable :: Gradient_hamil,acw
        real*8,dimension(:), allocatable ::Energy,acc1,E_levels,w_levels,knot_x
        real*8,dimension(:), allocatable :: pos,v,old_acc2,momentum,population_mat,mass
        real*8 :: m_acc1(1590,5)
        integer,dimension(:), allocatable :: lambda
        complex*16, dimension(:), allocatable :: cwork
        integer(kind=4) ::iseed
        integer, dimension(16) :: quantum_count
        real*8, dimension(3) :: O_pos,N_pos
 

contains

subroutine main_sub
        implicit none
        call parameter_value

        open(2,file='528atom.dat')
        do i=1,528
        read(2,*)Au_eq(i,:)
        enddo
        close(2)
       
        open(1,file='gold_neigh_pos.dat')
        do i=1,528
        read(1,*)Au_i(i),Au_nei_list(i,:)
        enddo
        close(1)

        call tensor_matrix1
end subroutine main_sub

subroutine parameter_value
        implicit none
        pi=4.d0*datan(1.d0)
        amu2kg=1.66053892d-27
        av=6.0221367d23
        conv=1000.d0/av
        kb=8.314d0/av
        mass_O=15.999d0*amu2kg
        mass_N=14.0067d0*amu2kg
        mass_Au=196.96d0*amu2kg
        red_mass=(mass_O*mass_N)/(mass_O+mass_N)
        total_mass=mass_O+mass_N
        planck_const=6.62607d-34
        c_light=2.99792458d10
        wavenumber_2_j=c_light*planck_const

        a0=457095d0*conv
        alpha0=3.7594*1.d10
        b0=30707*conv
        beta0=3.0082*1.d10
        f0=638.5*conv
        gamma0=2.743*1.d10
        r0_N_O=1.15077d-10

        a1=a0
        alpha1=alpha0
        b1=24.056*conv
        beta1=1.9649*1.d10
        r1_Au_N=2.3491*1.d-10
        d=347.22*conv*1.d-10
        c1=1.2423*1.d-10
        zimage=1.1536*1.d-10
        f1=495.98*conv
        gamma1=2.4890*1.d10
        r1_N_O=1.2904*1.d-10
        phi=511.37*conv
        ea=-0.67540*conv
        a2=11.842*conv
        a3=0.0061803
        gamma2=1.3693*1.d10
        b2=50*conv
        b3=0.0047213
        gamma3=2.0194*1.d10
        r_cut=10.d-10

        l_x=2.9499787804321573E-010
        r_x=3.3924755974969813E-009
        l_y=8.5158552149309472E-011
        r_y=3.2360249816737607E-009
        l_z=-7.230000000000003E-010
        r_z=0.d0
        box_len_x=r_x-l_x
        box_len_y=r_y-l_y
        box_len_z=r_z-l_z

        vib_state=0                     !vibrational state of the NO
        e_kin=conv*5.d0                 !it converts to joule per atom
end subroutine parameter_value


subroutine potential_force_calculate
        implicit none
        
        diab_pot=0.d0
        der_pot_N=0.d0
        der_pot_O=0.d0
        der_pot_Au=0.d0

        call interaction_list
        call nitrogen_oxygen_int
        call gold_nitrogen_int
        call gold_oxygen_int
        call image_potential
        call gold_gold_int
        call mat_diag(diab_pot,gr_eig,ex_eig,gr_vec,ex_vec)

        diab_pot(2,1)=diab_pot(1,2)
        diab_pot(2,2)=diab_pot(2,2)+phi-ea

        der_pot_N(2,1,:)=der_pot_N(1,2,:)
        der_pot_O(2,1,:)=der_pot_O(1,2,:)
        der_pot_Au(:,:,2,1)=der_pot_Au(:,:,1,2)

end subroutine potential_force_calculate

subroutine mat_diag(h,gr_eig,ex_eig,gr_vec,ex_vec)
        implicit none
        real*8,dimension(2,2),intent(in)::h
        real*8,intent(out)::gr_eig,ex_eig
        real*8,dimension(1,2),intent(out)::gr_vec,ex_vec
        real*8::tmp

        gr_eig=0.5d0*((h(1,1)+h(2,2))-dsqrt((h(1,1)-h(2,2))**2+4.d0*(h(1,2))**2))
        ex_eig=0.5d0*((h(1,1)+h(2,2))+dsqrt((h(1,1)-h(2,2))**2+4.d0*(h(1,2))**2))

        tmp=(h(1,2)**2)+(h(1,1)-gr_eig)**2
        gr_vec(1,1)=h(1,2)/dsqrt(tmp)
        gr_vec(1,2)=-(h(1,1)-gr_eig)/dsqrt(tmp)

        tmp=(h(1,2)**2)+(h(1,1)-ex_eig)**2
        ex_vec(1,1)=h(1,2)/dsqrt(tmp)
        ex_vec(1,2)=(ex_eig-h(1,1))/dsqrt(tmp)

end subroutine mat_diag


subroutine nitrogen_oxygen_int
        implicit none
        real*8::tmp1,dist(3),deriv(3)

        call pbc_distance(N(:),O(:),dist(:))
        r_N_O=dsqrt(sum(dist(:)*dist(:)))
        dist=dist/r_N_O

        tmp1=dexp(-gamma0*(r_N_O-r0_N_O))
        deriv(:)=2.d0*f0*(1.d0-tmp1)*tmp1*gamma0*dist(:)
        diab_pot(1,1)=diab_pot(1,1)+f0*((1.d0-tmp1)**2)
        der_pot_N(1,1,:)=der_pot_N(1,1,:)+deriv(:)
        der_pot_O(1,1,:)=der_pot_O(1,1,:)-deriv(:)

        tmp1=dexp(-gamma1*(r_N_O-r1_N_O))
        deriv(:)=2.d0*f1*(1.d0-tmp1)*gamma1*tmp1*dist(:)
        diab_pot(2,2)=diab_pot(2,2)+f1*((1.d0-tmp1)**2)
        der_pot_N(2,2,:)=der_pot_N(2,2,:)+deriv(:)
        der_pot_O(2,2,:)=der_pot_O(2,2,:)-deriv(:)

end subroutine nitrogen_oxygen_int

subroutine gold_oxygen_int
        implicit none
        real*8::tmp1,dist(3),deriv(3)

        do i=1,n_O_Au
           j=O_Au_list(i)
           call pbc_distance(O(:),Au(j,:),dist(:))
           r_Au_O=dsqrt(sum(dist(:)*dist(:)))
           dist(:)=dist(:)/r_Au_O

           tmp1=dexp(-alpha0*r_Au_O)
           deriv(:)=a0*tmp1*(-alpha0)*dist(:)
           diab_pot(1,1)=diab_pot(1,1)+a0*(tmp1-dexp(-alpha0*r_cut))
           der_pot_O(1,1,:)=der_pot_O(1,1,:)+deriv(:)
           der_pot_Au(j,:,1,1)=der_pot_Au(j,:,1,1)-deriv(:)

           tmp1=dexp(-alpha1*r_Au_O)
           deriv(:)=a1*tmp1*(-alpha1)*dist(:)
           diab_pot(2,2)=diab_pot(2,2)+a1*(tmp1-dexp(-alpha1*r_cut))
           der_pot_O(2,2,:)=der_pot_O(2,2,:)+deriv(:)
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)-deriv(:)

           tmp1=dexp(gamma2*r_Au_O)
           deriv(:)=a2*a3*gamma2*tmp1*dist(:)/((1.d0+a3*tmp1)**2)
           diab_pot(1,2)=diab_pot(1,2)-a2*((1.d0/(1.d0+a3*tmp1))-(1.d0/(1.d0+a3*dexp(gamma2*r_cut))))
           der_pot_O(1,2,:)=der_pot_O(1,2,:)+deriv(:)
           der_pot_Au(j,:,1,2)=der_pot_Au(j,:,1,2)-deriv(:)
         enddo
end subroutine gold_oxygen_int

subroutine gold_nitrogen_int
        implicit none
        real*8::tmp1,deriv(3),dist(3),costheta,dcostheta(3)

        do i=1,n_N_Au
           j=N_Au_list(i)
           call pbc_distance(N(:),Au(j,:),dist(:))
           r_Au_N=dsqrt(sum(dist(:)*dist(:)))
           dist=dist/r_Au_N

           tmp1=dexp(-beta0*r_Au_N)
           deriv(:)=b0*tmp1*(-beta0)*dist(:)
           diab_pot(1,1)=diab_pot(1,1)+b0*(tmp1-dexp(-beta0*r_cut))
           der_pot_N(1,1,:)=der_pot_N(1,1,:)+deriv(:)
           der_pot_Au(j,:,1,1)=der_pot_Au(j,:,1,1)-deriv(:)

           tmp1=dexp(-2.d0*beta1*(r_Au_N-r1_Au_N))
           deriv(:)=b1*tmp1*(-2.d0*beta1)*dist(:)
           diab_pot(2,2)=diab_pot(2,2)+b1*(tmp1-dexp(-2.d0*beta1*(r_cut-r1_Au_N)))
           der_pot_N(2,2,:)=der_pot_N(2,2,:)+deriv(:)
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)-deriv(:)

           tmp1=dexp(-beta1*(r_Au_N-r1_Au_N))
           deriv(:)=b1*tmp1*(-beta1)*dist(:)
           costheta=(O(3)-N(3))/r_N_O
           dcostheta(:)=-((costheta)*(N(:)-O(:)))/(r_N_O**2)
           dcostheta(3)=dcostheta(3)-1.d0/r_N_O
           diab_pot(2,2)=diab_pot(2,2)-2.d0*b1*(costheta**2)*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))
           der_pot_N(2,2,:)=der_pot_N(2,2,:)-4.d0*b1*costheta*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))*dcostheta(:)
           der_pot_O(2,2,:)=der_pot_O(2,2,:)-4.d0*b1*costheta*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))*(-dcostheta(:))
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)+2.d0*(costheta**2)*deriv(:)
           der_pot_N(2,2,:)=der_pot_N(2,2,:)-2.d0*(costheta**2)*deriv(:)
           tmp1=dexp(gamma3*r_Au_N)
           deriv(:)=b2*b3*tmp1*gamma3*dist(:)/((1.d0+b3*tmp1)**2)
           diab_pot(1,2)=diab_pot(1,2)-b2*((1.d0/(1.d0+b3*tmp1))-(1.d0/(1.d0+b3*dexp(gamma3*r_cut))))
           der_pot_N(1,2,:)=der_pot_N(1,2,:)+deriv(:)
           der_pot_Au(j,:,1,2)=der_pot_Au(j,:,1,2)-deriv(:)
           enddo

end subroutine gold_nitrogen_int

subroutine gold_gold_int
        implicit none
        real*8::pot(1,3),dist(1,3),dist_eq(1,3),tensor_mat(3,3),dist_trans(3,1)
        integer::l,m
        der_pot_Au_Au(:,:)=0.d0
        pot_Au_Au=0.d0
        do i=1,528
           res:do j=1,12
               tensor_mat(:,:)=0.d0
               if(Au_nei_list(i,j)==0) cycle
               k=Au_nei_list(i,j)
               call pbc_distance(Au_eq(i,:),Au_eq(k,:),dist_eq(1,:))
               call pbc_distance(Au(i,:),Au(k,:),dist(1,:))
               do l=1,3
               do m=1,3
               tensor_mat(l,m)=ten_mat(i,k,l,m)
               enddo
               enddo
               dist=dist-dist_eq
               dist_trans=transpose(dist)
               pot=matmul(dist,tensor_mat)
               pot_Au_Au=pot_Au_Au+0.5d0*(matmul(pot,dist_trans))
               der_pot_Au_Au(i,:)=der_pot_Au_Au(i,:)+pot(1,:)
               der_pot_Au_Au(k,:)=der_pot_Au_Au(k,:)-pot(1,:)
               enddo res
         enddo

end subroutine gold_gold_int

subroutine image_potential
        implicit none
        real*8::tmp1,tmp2

        z_COM=(mass_N*N(3)+mass_O*O(3))/(mass_N+mass_O)
        tmp1=dsqrt(c1*c1+(z_COM-zimage)**2)
        diab_pot(2,2)=diab_pot(2,2)-d/tmp1
        der_pot_N(2,2,3)=d*(z_COM-zimage)*(mass_N/total_mass)/(tmp1**3)+der_pot_N(2,2,3)
        der_pot_O(2,2,3)=d*(z_COM-zimage)*(mass_O/total_mass)/(tmp1**3)+der_pot_O(2,2,3)

end subroutine image_potential

subroutine tensor_matrix1
        implicit none
        real*8::rij(1,3),t(3,3),t_trans(3,3),t_mat(3,3),ll,a,thetay,thetaz,rij1(1,3)
        real*8::alpha=-4.94d0,beta=17.15d0,gamma4=19.4d0!,ten_mat(528,528,3,3)
        ten_mat(:,:,:,:)=0.d0
        do i=1,528
        res: do j=1,12
                if(Au_nei_list(i,j)==0) cycle
                k=Au_nei_list(i,j)
                rij(1,:)=(Au_eq(i,:)-Au_eq(k,:))
                call check_pbc1(rij,rij1)
                a=dsqrt(sum(rij1*rij1))

                thetaz=dasin(rij1(1,2)/a)
                !thetay=datan(rij1(1,3)/rij1(1,1))
                thetay=dacos(rij1(1,1)/dsqrt(rij1(1,1)**2+rij1(1,3)**2))
                t(1,1)=dcos(thetay)*dcos(thetaz)
                t(1,2)=dcos(thetay)*dsin(thetaz)
                t(1,3)=dsin(thetay)
                t(2,1)=-dsin(thetaz)
                t(2,2)=dcos(thetaz)
                t(2,3)=0.d0
                t(3,1)=-dsin(thetay)*dcos(thetaz)
                t(3,2)=-dsin(thetay)*dsin(thetaz)
                t(3,3)=dcos(thetay)
                t_mat=0.d0
                t_mat(1,1)=beta+gamma4
                t_mat(2,2)=beta-gamma4
                t_mat(3,3)=alpha
                t=transpose(t)
                t_trans=transpose(t)
                ten_mat(i,k,:,:)=matmul(t_trans,matmul(t_mat,t))


             enddo res
         enddo
end subroutine tensor_matrix1

subroutine check_pbc1(x,y)              !check pbc of r vector
        implicit none
        real*8,intent(in)::x(1,3)
        real*8,intent(out)::y(1,3)
        if(abs(x(1,1))>box_len_x*0.5d0)then
                if(x(1,1)>0.d0) y(1,1)=x(1,1)-box_len_x
                if(x(1,1)<0.d0) y(1,1)=x(1,1)+box_len_X
        else
                y(1,1)=x(1,1)
        endif
        if(abs(x(1,2))>box_len_y*0.5d0)then
                if(x(1,2)>0.d0) y(1,2)=x(1,2)-box_len_y
                if(x(1,2)<0.d0) y(1,2)=x(1,2)+box_len_y
        else
                y(1,2)=x(1,2)
        endif
        y(1,3)=x(1,3)
end subroutine check_pbc1

subroutine pbc_distance(x,y,rij)
  implicit none
  real*8,intent(in)::x(3),y(3)
  real*8,intent(out)::rij(3)

  rij=x-y
  if(rij(1)>box_len_x/2.d0) rij(1)=rij(1)-box_len_x
  if(rij(1)<-box_len_x/2.d0) rij(1)=rij(1)+box_len_x

  if(rij(2)>box_len_y/2.d0) rij(2)=rij(2)-box_len_y
  if(rij(2)<-box_len_y/2.d0) rij(2)=rij(2)+box_len_y

end subroutine pbc_distance


subroutine interaction_list
        implicit none
        real*8::dist(3),r
        N_Au_list=0
        n_N_Au=1
        do i=1,528
        call pbc_distance(N(:),Au(i,:),dist(:))
        r=dsqrt(sum(dist(:)*dist(:)))
        if(r<r_cut)then
                N_Au_list(n_N_Au)=i
                n_N_Au=n_N_Au+1
        endif
        enddo
        O_Au_list=0
        n_O_Au=1
        do i=1,528
        call pbc_distance(O(:),Au(i,:),dist(:))
        r=dsqrt(sum(dist(:)*dist(:)))
        if(r<r_cut)then
                O_Au_list(n_O_Au)=i
                n_O_Au=n_O_Au+1
        endif
        enddo

end subroutine interaction_list

subroutine set_NO2                      !any random NO orientation. random theta and random phi.
        implicit none
        real*8::rel_x,rel_v,theta,phi,rnd,omega0,vib_energy

        call random_number(rnd)
        COM(1)=l_x+box_len_x*rnd!*0.33d0+box_len_x*0.33d0
        call random_number(rnd)
        COM(2)=l_y+box_len_y*rnd!*0.33d0+box_len_y*0.33d0
        COM(3)=11.d-10

        omega0=dsqrt(2.d0*f0*gamma0**2/red_mass)
        vib_energy=(0.5+real(vib_state))*planck_const*omega0/(2.d0*pi)

        call random_number(rnd)
        theta=2.d0*pi*rnd
        rel_x=dsqrt(2.d0*vib_energy/(red_mass*omega0**2))*dcos(theta)
        rel_v=-dsqrt(2.d0*vib_energy/red_mass)*dsin(theta)

        call random_number(rnd)
        theta=pi*rnd
        call random_number(rnd)
        phi=2.d0*pi*rnd

        N(1)=COM(1)-(r0_N_O+rel_x)*dsin(theta)*dcos(phi)*mass_O/(mass_O+mass_N)
        N(2)=COM(2)-(r0_N_O+rel_x)*dsin(theta)*dsin(phi)*mass_O/(mass_O+mass_N)
        N(3)=COM(3)-(r0_N_O+rel_x)*dcos(theta)*mass_O/(mass_O+mass_N)

        O(1)=COM(1)+(r0_N_O+rel_x)*dsin(theta)*dcos(phi)*mass_N/(mass_O+mass_N)
        O(2)=COM(2)+(r0_N_O+rel_x)*dsin(theta)*dsin(phi)*mass_N/(mass_O+mass_N)
        O(3)=COM(3)+(r0_N_O+rel_x)*dcos(theta)*mass_N/(mass_O+mass_N)
        N_vel(1)=rel_v*dcos(phi)*dsin(theta)*0.5d0*mass_O/total_mass
        N_vel(2)=rel_v*dsin(phi)*dsin(theta)*0.5d0*mass_O/total_mass
        N_vel(3)=-dsqrt(2*e_kin/total_mass)+rel_v*dcos(theta)*0.5d0*mass_O/total_mass

        O_vel(1)=-rel_v*dcos(phi)*dsin(theta)*0.5d0*mass_N/total_mass
        O_vel(2)=-rel_v*dsin(phi)*dsin(theta)*0.5d0*mass_N/total_mass
        O_vel(3)=-dsqrt(2*e_kin/total_mass)-rel_v*dcos(theta)*0.5d0*mass_N/total_mass

end subroutine

!...........................................................................................................................................
subroutine reshaping_array(array,reshaped,ND_1d)
        real*8, intent(inout), dimension(530,3) :: array
        real*8, intent(inout), dimension(3*530) :: reshaped
        integer :: i8,i9,i10
        integer, intent(in) :: ND_1d

        i10=0

  do i8=1,530
     do i9=1,3
        if (i10 .le. 530*3) then
                i10=i10+1
                if (ND_1d.eq.1) then
                reshaped(i10)=array(i8,i9)
                else
                array(i8,i9)=reshaped(i10)
                endif
        endif
     enddo
  enddo




end subroutine reshaping_array

!...........................................................................................................................................

subroutine setting_initial_position(Au_N_O,cordinates_Au,ND_2)
real*8, intent(inout) :: cordinates_Au(528,3)
real*8, intent(inout) :: Au_N_O(530,3)
integer :: i5
integer, intent(in) :: ND_2


do i5=3,530
   if (ND_2.eq.1) then
     Au_N_O(i5,:)=cordinates_Au(i5-2,:)
   else
     cordinates_Au(i5-2,:)=Au_N_O(i5,:)
   end if
enddo



end subroutine setting_initial_position


!...............................................................................................................................................

subroutine combined_derivatives(com,d1,d2)
                        real*8, intent(out) :: com(1590)
                        integer :: i
                        integer, intent(in) :: d1,d2
                        real*8 :: vertical_pot_Au(528*3)

        do i=1,1590
           if (i.le.3) then
              com(i)=der_pot_N(d1,d2,i)
           else if (i.le.6) then
              com(i)=der_pot_O(d1,d2,i-3)
           else
              call reshaping_array2(der_pot_Au(:,:,d1,d2),vertical_pot_Au,1)
              com(i)=vertical_pot_Au(i-6)
           end if
        enddo


end subroutine combined_derivatives
!..........................................................................................................................................................
subroutine returning_N_O_Au(bool)
implicit none
integer, intent(in) :: bool
integer :: i,j,a,y,z
real*8, dimension(530,3) :: xr

if (bool.eq.1) then
   call reshaping_array(xr,pos,2)
   do i=1,3
      N(i)=xr(1,i)
      O(i)=xr(2,i)
   enddo
   call setting_initial_position(xr,Au,2)
else
    do j=1,3
    xr(1,j)=N(j)
    xr(2,j)=O(j)
    end do
    call setting_initial_position(xr,Au,1)
    call reshaping_array(xr,pos,1)
end if



end subroutine returning_N_O_Au
!..............................................................................................................................................................
  subroutine reshaping_array2(array_528_3,reshaped_1584,ND_1d)
        real*8, intent(inout), dimension(528,3) :: array_528_3
        real*8, intent(inout), dimension(3*528) :: reshaped_1584
        integer :: i8,i9,i10
        integer, intent(in) :: ND_1d

        i10=0

  do i8=1,528
     do i9=1,3
        if (i10 .le. 528*3) then
                i10=i10+1
                if (ND_1d.eq.1) then
                reshaped_1584(i10)=array_528_3(i8,i9)
                else
                array_528_3(i8,i9)=reshaped_1584(i10)
                endif
        endif
     enddo
  enddo




end subroutine reshaping_array2


!...........................................................................................................................................................
subroutine setup_initial_values2
implicit none
real*8 :: rnd2,rnd1
real*8 :: vel(530,3),x(530,3)
integer ::n1,loop_x,loop_v,i1,j1,k1,l1,t1,read_vel,g1,t12


call main_sub

call set_NO2


do j1=1,3
x(1,j1)=N(j1)
x(2,j1)=O(j1)
end do

do loop_x=3,530
        open(656,file="gold_positions_fort.20")
        read(656,*) x(loop_x,1), x(loop_x,2), x(loop_x,3)
end do





call setting_initial_position(x,Au,2)
call reshaping_array(x,pos,1)


!pos=pos/(5.2917E-011)



vel(1,1)=N_vel(1)
vel(1,2)=N_vel(2)
vel(1,3)=N_vel(3)

vel(2,1)=O_vel(1)
vel(2,2)=O_vel(2)
vel(2,3)=O_vel(3)


do read_vel=3,530
   if (read_vel.le.396) then
        open(658,file="gold_velocities_fort.21")
        read(658,*) vel(read_vel,1), vel(read_vel,2), vel(read_vel,3)
   else
           vel(read_vel,:)=0
   end if

end do



v=v/(2187691.26379)


call potential_force_calculate


time=0.00

c=0
do i=1,Ne
   c(i,i+1)=1
enddo 

do n1=1,Ne
lambda(n1)=n1
enddo

call potential
!write(87,*) sum(c(:,1)*conjg(c(:,1)))
c=matmul(c,Energy_hamil)
!write(88,*) sum(c(:,1)*conjg(c(:,1)))  


end subroutine 
!........................................................................

subroutine setup_initial_values1
implicit none
integer(kind=4) :: O,j
integer(kind=4), dimension(:),allocatable::x
real*8 :: rhos,wt
integer:: ut,xt,ip,kn,loop_mass1,loop_mass2,loop_mass3
real*8, dimension(14) :: inpot

 open(25,file='fort.23')
 do ip=1,15
 read(25,*) inpot(ip)
 enddo
close(25)



Hi=inpot(1)
Ne=int(Hi/2)
tau=inpot(4)
Band_width=inpot(5)
number_of_atoms=inpot(6)
omega=inpot(7)
gh=inpot(8)
dG=inpot(9)
KT=inpot(10)
dtc=inpot(11)
dt=inpot(12)
total_dimensions=inpot(13)
wt=inpot(14)
movable_layer=inpot(15)
total_time=wt*5000
rhos=real(Hi)/Band_width
!Vr=sqrt(tau/(2*3.1416))
  allocate(pos(total_dimensions))
  allocate(v(total_dimensions))
  allocate(momentum(total_dimensions))
  allocate(acc1(total_dimensions))
  allocate(old_acc2(total_dimensions))
  allocate(acw(Hi,Hi,total_dimensions))
  allocate(H(Hi,Hi))
  allocate(Energy_hamil(Hi,Hi)) 
  allocate(Energy(Hi))
  allocate(Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c(Ne,Hi))
  allocate(b(Ne,Hi))
  allocate(A(Ne,Hi))
  allocate(g(Ne,Hi))
  allocate(lambda(Ne))
  allocate(Identity(Hi,Hi))
  allocate(Grd(Hi,Hi))
  allocate(vdij(Hi,Hi))
  allocate(population_mat(int(total_time/dtc)))
  allocate(E_levels(int(Hi/2)))
  allocate(knot_x(int(Hi/2)))
  allocate(w_levels(int(Hi/2)))
  allocate(mass(number_of_atoms))
call random_seed(size=O)
  allocate(x(O))
   do j=1,O
   x(j)=j**6*iseed+2777772
   enddo
  call random_seed(put=x)
 open(2,file='rndom',status='new')
 write(2,*) iseed,x,O
 

open(167,file='raw_x.txt')
do ut=1,int(Hi/2)
read(167,*) knot_x(ut)
enddo


open(169,file='raw_w.txt')
do xt=1,int(Hi/2)
read(169,*)w_levels(xt)
enddo
close(169)

do loop_mass1=1,3
      mass(loop_mass1)=14.0067
end do
do loop_mass2=4,6
        mass(loop_mass2)=15.999
enddo

do loop_mass3=7,1590
mass(loop_mass3)=196.96
end do


  

end subroutine
!.........................................................
  
subroutine gaussian_random_number(rnd0)
!USE IFPORT
   !! generates gaussian distribution with center 0, sigma 1
   !! q0+sig*rnd gives center=q0, sigma=sig
   implicit none
   integer(kind=4) :: n,j,M,O,k
   real*8,intent(out)::rnd0
   real*8 rnd1,rnd2,pi
   pi=dacos(-1.d0)
   call random_number(rnd1)
   call random_number(rnd2)
   rnd0 = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!.............................................................
!subroutine setup_initial_values3
!implicit none
!integer :: i,n1

!time=0.00

!c=0
!do i=1,Ne
!   c(i,i+1)=1
!enddo 

!do n1=1,Ne
!lambda(n1)=n1
!enddo

!call potential
!c=matmul(c,Energy_hamil)



!end subroutine 
!........................................................................
subroutine diag_wrapper(matrix,nsize,eigen_values,eigen_vectors)
real*8, intent(inout) :: matrix(nsize,nsize)
real*8, dimension(nsize,nsize) :: mat
integer LWORK,nsize
real*8, allocatable :: WORK(:)
real*8, intent(out) :: eigen_vectors(nsize,nsize),eigen_values(nsize)
mat=matrix
LWORK=3*nsize-1
allocate(WORK(LWORK))
call dsyev('V','U',nsize,mat,nsize,eigen_values,WORK,LWORK,INFO)
eigen_vectors=mat

end subroutine

!..........................................................................
subroutine logm(mat,log_mat,n)
   !! http://arxiv.org/pdf/1203.6151v4.pdf
   implicit none
   integer,intent(in):: n
   real*8,intent(in):: mat(n,n)
   real*8,intent(out):: log_mat(n,n)
   integer i
   complex*16 T(n,n),en(n),vect(n,n)
   complex*16 dd(n,n)

   call schur(mat,T,n,en,vect,nold,cwork)

   dd=0.d0
   do i=1,n
     dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
   enddo

   log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm


!...................................................................................

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
   !! Diaganalizing matrix using dsyevr. First m_values eigen values and
!eigenvectors computed.
   !! The module's common variables should contain:

   !! Initialize nold=0

   !! nold makes sure that everytime value of n changes, work and iwork
!are re-allocated for optimal performance.
   !! mat is destroyed after use.

  implicit none
   integer,intent(in) :: n
   integer,intent(inout) :: nold
   complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
   real*8,intent(in) :: mat(n,n)
   complex*16,intent(out) :: T(n,n)
   complex*16,allocatable,intent(inout):: cwork(:)
   real*8 rwork(n)
   complex*16 mat_c(n,n)

   integer lwork
   logical:: select
   logical bwork(n)
   integer sdim,info,AllocateStatus

   T=mat

   info=0
   sdim=0

   if(nold.ne.n .or. .not.allocated(cwork)) then
   !if(nold.ne.n) then
     lwork=-1
     if(allocated(cwork))deallocate(cwork)
     allocate(cwork(n))
     call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
     lwork=int(cwork(1))
     deallocate(cwork)
     allocate(cwork(lwork),STAT=AllocateStatus)
     if(allocatestatus.ne.0) write(6,*)"problem in schur, allocation"
     nold=n
   endif
!   write(151,*) time,T
   lwork=size(cwork)
   call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
   if(info.ne.0) then
     write(6,*) "problem in scur",info
     stop
   endif
end subroutine schur


!...................................................................................................................
subroutine potential
implicit none
real*8 :: U0,U1,Vr
real*8 :: der_U0(total_dimensions), der_U1(total_dimensions), der_Vak(total_dimensions)
integer :: i,j,k,l
real*8 :: pot_SI_atomic,force_SI_atomic



call combined_derivatives(der_U0,1,1)
call combined_derivatives(der_U1,2,2)
call combined_derivatives(der_Vak,1,2)

U0=gr_eig
U1=ex_eig
Vr=Diab_pot(1,2)

pot_SI_atomic=1/4.359744E-18
force_SI_atomic=12137804.1108


U1=U1*pot_SI_atomic
U0=U0*pot_SI_atomic
Vr=Vr*pot_SI_atomic

der_U0=der_U0*force_SI_atomic
der_U1=der_U1*force_SI_atomic
der_Vak=der_Vak*force_SI_atomic

H(1,1)=U1-U0
    do i=2,Hi
        if (i.le.int(Hi/2)) then
                H(1,i)=sqrt(Band_width*w_levels(i)/2)*Vr
                H(i,1)=sqrt(Band_width*w_levels(i)/2)*Vr
        else
                H(1,i)=sqrt(Band_width*w_levels(i-int(Hi/2))/2)*Vr
                H(i,1)=sqrt(Band_width*w_levels(i-int(Hi/2))/2)*Vr
        end if
        do j=2,Hi
                if (i.eq.j) then
                        if (i.le.(int(Hi/2)+1)) then
                             H(i,j)=-(Band_width/2)*(0.5+0.5*knot_x(int(Hi/2)-i+2))
                        else
                             H(i,j)=(Band_width/2)*(0.5+0.5*knot_x(i-int(Hi/2)-1))
                                 
                        end if
                else
                        H(i,j)=0.0
                end if
        end do
   enddo




!do i=1,Hi
! do j=1,Hi
!   write(153,'(f12.5$)')H(i,j)
! enddo
! write(153,*)
!enddo

nsize=Hi
call diag_wrapper(H,nsize,Energy,Energy_hamil)

Gradient_hamil=0

Gradient_hamil(1,1,:)= der_U1 -der_U0
do k=2,Hi
   if (k.le.int(Hi/2)) then
       Gradient_hamil(1,k,:)=der_Vak*sqrt(Band_width*w_levels(k)/2)
       Gradient_hamil(k,1,:)=der_Vak*sqrt(Band_width*w_levels(k)/2)
   else
       Gradient_hamil(1,k,:)=der_Vak*sqrt(Band_width*w_levels(k-int(Hi/2))/2)
       Gradient_hamil(k,1,:)=der_Vak*sqrt(Band_width*w_levels(k-int(Hi/2))/2)      
   end if
 enddo

!write(211,*) der_pot_Au(:,:,1,2)
end subroutine
!........................................................................................
subroutine nonadiabaticvector(t,u)
implicit none
integer :: acwi
integer, intent(in) :: t,u
do acwi=1,total_dimensions

if (t.ne.u) then
acw(t,u,acwi)=(sum(Energy_hamil(:,t)*matmul(Gradient_hamil(:,:,acwi),Energy_hamil(:,u))))/(Energy(t)-Energy(u))
else
acw(t,u,acwi)=0
endif

end do
end subroutine
!..........................................................................
subroutine force
implicit none
integer :: acc1i,z,u,q
!real*8, dimension(total_dimensions,Ne) :: m_acc1
real*8, dimension(total_dimensions) :: grad_U0
do acc1i=1,total_dimensions

do z=1,Ne
m_acc1(acc1i,z)=-((sum((Energy_hamil(:,lambda(z)))*matmul(Gradient_hamil(:,:,acc1i),Energy_hamil(:,lambda(z)))))/mass(acc1i))

end do
end do


call combined_derivatives(grad_U0,1,1)


acc1=sum(m_acc1(:,1:Ne))-grad_U0/mass

!do u=2,Ne
!  m_acc1(:,u)=m_acc1(:,u-1)+m_acc1(:,u)
!enddo

!acc1=m_acc1(:,Ne)-grad_U0



end subroutine
!..........................................................................
subroutine Rungekutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi 
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-sum(v*acw(p,q,:))
enddo
enddo

k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine
!.........................................................................

subroutine Rungefutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-vdij(p,q)
enddo
enddo



k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine



!...................................................................................
subroutine gs
implicit none
integer :: i,d,r,s,l,m,o,t
integer :: findloc
integer, dimension(Ne) :: p
complex*16 :: det_S1,det_S2
complex*16, dimension(Ne,Ne) ::  S1,S2

nold=0
p=lambda
do i=1,Ne
do d=1,Hi

if (FINDLOC(p,d,1).eq.0) then

do r=1,Ne
S1(:,r)=c(:,p(r))
enddo

call modulus(S1,Ne,det_S1)

p(i)=d


do s=1,Ne
S2(:,s)=c(:,p(s))
enddo
call modulus(S2,Ne,det_S2)
A(d,i)=det_S1*conjg(det_S2)
b(d,i)=-2*real(A(d,i)*vdij(d,i))
A(i,i)=det_S1*conjg(det_S1)
g(i,d)=dt*real(b(d,i)/A(i,i))

else
A(i,d)=0
b(i,d)=0
g(i,d)=0
endif

enddo
p(i)=lambda(i)
enddo
end subroutine
!...........................................................................
subroutine quantum_evolution

implicit none
integer :: k,j,r,s,i,t,u,M
real*8 :: EE,TE,KE,rnd
integer(kind=4), dimension(:), allocatable :: z

call random_number(rnd)
call Rungefutta
call gs
jloop : do t=1,Ne
do u=1,Hi
if (g(t,u)>rnd) then
call nonadiabaticvector(t,u) 
call hop(t,u)

exit jloop
endif
enddo
enddo jloop

end subroutine
!.........................................................................
subroutine classical_evolution
integer ::no,ip,r,TT,yt,i,nan1,nan2
real*8, dimension(Hi,Hi) :: old_Energy_hamil
real*8, dimension(530,3) :: sq_pos
real*8,dimension(Hi) :: signature
real*8 :: KE,TE,EE
logical :: tr

!write(79,*) sum(c(:,2)*conjg(c(:,2)))

TT=int(total_time/dtc)
open(206, file='VMD_file.xyz')
do while(time.le.total_time)

  old_Energy_hamil=Energy_hamil
  call boundary_conditions
  call velocity_verlet2
  call boundary_conditions
  !do i=1,6
 !    write(204,*) time, pos(i)/5.2917E-011, i,'a', 32.5E-10/5.2917E-011
 ! enddo
  call Make_Video 



  call signt(old_Energy_hamil,signature)
  do r=1,Hi
     if (signature(r)<0) then
         Energy_hamil(:,r)=-Energy_hamil(:,r)
     endif
  enddo

!write(95,*) old_Energy_hamil
call vdotd(old_Energy_hamil)

no=int(dtc/dt)
  do while(no>0) 
     call quantum_evolution
   !  write(80,*) sum(c(:,1)*conjg(c(:,1))) , time, pos(1) 
     no=no-1
  enddo



call populations
yt=int(time/dtc)
population_mat(yt)=real(population)




 !open(35,file="poso")
 !    write(35,*)time,0.5*sum(mass*v*v)+sum(Energy_hamil(:,1:Ne))+gr_eig,sum(c(1,:)*conjg(c(1,:)))
 !close(35)

time=time+dtc 


enddo
!write(79,*) sum(c(:,1)*conjg(c(:,1))) 

close(206)

end subroutine
!................................................................................
!subroutine hop(t,u)
!implicit none
!integer, intent(in) :: t,u
!real*8 :: para_v
!real*8, dimension(:),allocatable :: perp_v


!allocate(perp_v(total_dimensions))
!if (((sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))) &  !frustrated hops are forbidden
!**2+(2*Energy(lambda(t))/mass)-(2*Energy(u)/mass))>0) then
!    para_v=sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))
!    perp_v=v-para_v*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))  
!    para_v=((para_v)/abs(para_v))*sqrt(para_v**2+(2*Energy(lambda(t))/mass)-(2*Energy(u)/mass))!initial_KE+initial_PE=final_KE+final_PE,
!    v=(para_v)*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))+perp_v
!    lambda(t)=u
!end if

!Note: (sum(v*acw(lambda(t),u,:))/norm(acw(lambda(t),u,:))) is the parallel component of velocity in the direction of nonadiabatic coupling vector
!end subroutine
!................................................................................

subroutine classical_evolution2
integer :: p,r,TT,yt,i,no,q_int,nt,q_t
real*8, dimension(Hi,Hi) :: old_Energy_hamil
real*8, dimension(530,3) :: sq_pos
real*8,dimension(Hi) :: signature
real*8, dimension(1590) :: poso
!real*, dimension(1590) :: velo
real*8 :: KE,TE,EE
integer :: ipos,iham1,iham2
!call potential(U0,U1,Vr,der_U0,der_U1,der_Vak)
!call force


!pos=pos/(5.2917E-011) 

!v=v/(2187691.26379)  


do nt=1,Ne
!     write(89,*) sum(c(:,nt)*conjg(c(:,nt))) , time, pos(1)
 enddo
TT=int(total_time/dtc)


!write(35,*)time,pos(ipos),0.5*sum(mass*v*v)+sum(Energy_hamil(:,1:Ne))+gr_eig

 !  open(40,file="poso1")
 !   do ipos=1,1590
 !    write(40,*)time,pos(ipos),0.5*sum(mass*v*v)+sum(Energy_hamil(:,1:Ne))+gr_eig
 !  enddo
 ! close(40)






do while(time.le.total_time)


old_Energy_hamil=Energy_hamil
call velocity_verlet2


!do q_t=1,Ne
!     write(90,*) sum(c(:,q_t)*conjg(c(:,q_t))) , time, pos(1)
!enddo


call signt(old_Energy_hamil,signature)
do r=1,Hi
 if (signature(r)<0) then
    Energy_hamil(:,r)=-Energy_hamil(:,r)
 endif
enddo



call vdotd(old_Energy_hamil)

no=int(dtc/dt)
  do while(no>0)
     call quantum_evolution
 !    do q_int=1,Ne
 !    write(85,*) sum(c(:,q_int)*conjg(c(:,q_int))) , time, pos(1)
 !    enddo
     no=no-1
  enddo










!write(35,*) time,0.5*sum(mass*v*v)+sum(Energy_hamil(:,1:Ne))+gr_eig,pos(1)
!write(35,*) pos 

time=time+dtc
enddo



end subroutine classical_evolution2

!............................................................................................................................

subroutine hop(t,u)
implicit none
integer, intent(in) :: t,u
integer :: reciprocal_mass_loop,v_loop
real*8, dimension(total_dimensions) :: reciprocal_mass
real*8 :: a,b,gama,frustrated_condition




do reciprocal_mass_loop=1,total_dimensions
     reciprocal_mass(reciprocal_mass_loop)=1/mass(reciprocal_mass_loop)
end do



       a=0.5*sum(reciprocal_mass*acw(lambda(t),u,:))
       b=sum(v*acw(lambda(t),u,:))
       frustrated_condition=b**2+4*a*(Energy(lambda(t))-Energy(u))
       if (frustrated_condition>0) then
           if (b<0) then
           gama=(b+sqrt(frustrated_condition))/2*a
           else
           gama=(b-sqrt(frustrated_condition))/2*a
           end if
           do v_loop=1,total_dimensions
             v(v_loop)=v(v_loop)-gama*acw(lambda(t),u,v_loop)/mass(v_loop)
           end do
           lambda(t)=u
        end if
 



end subroutine hop

!................................................................................................................
subroutine velocity_verlet
real*8 :: delr(total_dimensions),delv(total_dimensions),fix_pos(total_dimensions-movable_layer)
real*8 :: gama_dt,gamma_B,c0,c10,c20,pos_atm,v_atm
integer :: ir,iv,ips


pos_atm=5.2917E-011
pos=pos/pos_atm

gamma_B=2*omega
gama_dt=gamma_B*dtc
c0=dexp(-gama_dt)
c10=1.d0/gama_dt*(1.d0-c0)
c20=1.d0/gama_dt*(1.d0-c10)
call stochastic_force(delr,delv)

do ips=1,(total_dimensions-movable_layer)
    fix_pos(ips)=pos(movable_layer+ips)
 !   write(199,*) ips, movable_layer+ips
enddo


pos=pos+c10*dtc*v+c20*dtc*dtc*acc1+delr

do ips=movable_layer+1,total_dimensions
pos(ips)=fix_pos(ips-movable_layer)
!write(200,*) ips, ips-movable_layer
enddo


old_acc2=acc1
pos=pos*pos_atm

call main_sub
call returning_N_O_Au(1)                                                                                                                                                                                           
call potential_force_calculate        

                                                                                                                                                                            
call returning_N_O_Au(2)          

                                                                                                                                                                                 

call potential
call force
v=c0*v+(c10-c20)*dtc*old_acc2+c20*dtc*acc1+delv





end subroutine
!................................................................................
subroutine vdotd(old_Energy_hamil)
real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave
!real*8, intent(out),dimension(Hi,Hi) :: dij
real*8, dimension(Hi,Hi) :: Ut,logarithm_Ut
integer :: i

Ut=matmul(transpose(old_Energy_hamil),Energy_hamil)
!write(152,*)time,Energy_hamil
call logm(Ut,logarithm_Ut,Hi)
vdij=logarithm_Ut/dtc




end subroutine
!........................................................................................
subroutine signt(old_Energy_hamil,signature)

real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave

integer :: i
real*8, intent(out),dimension(Hi) :: signature
do i=1,Hi
signature(i)=sum(old_Energy_hamil(:,i)*Energy_hamil(:,i))
enddo



end subroutine
!....................................................................................
subroutine velocity_verlet2
real*8 :: pos_atm
integer :: q_int,ips
real*8 :: fix_pos(total_dimensions-movable_layer)



!write(91,*) pos, v
 !do q_int=1,Ne
 !    write(93,*) sum(c(:,q_int)*conjg(c(:,q_int))) , time, pos(1)
 !enddo

pos_atm=5.2917E-011
pos=pos/pos_atm

do ips=1,(total_dimensions-movable_layer)
    fix_pos(ips)=pos(movable_layer+ips)
  !  write(199,*) ips, movable_layer+ips
enddo


pos=pos+v*dtc+0.5*acc1*dtc*dtc

do ips=movable_layer+1,total_dimensions
pos(ips)=fix_pos(ips-movable_layer)
!write(200,*) ips, ips-movable_layer
enddo



old_acc2=acc1

pos=pos*pos_atm

call main_sub
call returning_N_O_Au(1)
call potential_force_calculate
call returning_N_O_Au(2)

call potential
call force
v=v+0.5*(acc1+old_acc2)*dtc

!write(92,*) pos, v

 !do q_int=1,Ne
 !    write(94,*) sum(c(:,q_int)*conjg(c(:,q_int))) , time, pos(1)
!enddo




end subroutine

!...................................................................................
subroutine stochastic_force(delr,delv)
real*8, intent(out) :: delr(total_dimensions),delv(total_dimensions)
integer :: i
real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv,gdt,gamma_B
gamma_B=2*omega
gdt=gamma_B*dtc

do i=1,total_dimensions

sig_r=dtc*dsqrt(KT/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
sig_v=dsqrt(KT/mass(i)*(1-dexp(-2*gdt)))
sig_rv=(dtc*KT/mass(i)*1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v) !correlation coefficient

call gaussian_random_number(rnd1)
call gaussian_random_number(rnd2)
delr(i)=sig_r*rnd1
delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
enddo

end subroutine stochastic_force
!....................................................................................
subroutine populations
integer :: a,i,j,k
real*8, dimension(Hi,Hi,Ne) :: rho_a
real*8, dimension(Hi,Hi) :: rho_d
real*8, dimension(Ne) :: LUMO_population
do a=1,Ne
do i=1,Hi
do j=1,Hi
rho_a(i,j,a)=c(a,i)*conjg(c(a,j))
enddo
enddo
enddo

do k=1,Ne
rho_d=matmul(Energy_hamil,(matmul(rho_a(:,:,k),transpose(Energy_hamil))))
LUMO_population(k)=rho_d(1,1)
enddo


population=sum(LUMO_population(1:Ne))

end subroutine
!...................................................................................................................

subroutine modulus(matrix,n,determinant)
 IMPLICIT NONE
     complex*16, DIMENSION(n,n) :: matrix
     INTEGER, INTENT(IN) :: n
     complex*16 :: m, temp
     INTEGER :: i, j, k, l
     LOGICAL :: DetExists = .TRUE.
     complex*16,intent(out) :: determinant
     l = 1
     !Convert to upper triangular form
     
     DO k = 1, n-1
         IF (matrix(k,k) == 0) THEN
             DetExists = .FALSE.
             DO i = k+1, n
                 IF (matrix(i,k) /= 0) THEN
                     DO j = 1, n
                         temp = matrix(i,j)
                         matrix(i,j)= matrix(k,j)
                         matrix(k,j) = temp
                     END DO
                     DetExists = .TRUE.
                     l=-l
                     EXIT
                 ENDIF
             END DO
             IF (DetExists .EQV. .FALSE.) THEN
                 determinant= 0
                 return
             END IF
         ENDIF
         DO j = k+1, n
             m = matrix(j,k)/matrix(k,k)
             DO i = k+1, n
                 matrix(j,i) = matrix(j,i) - m*matrix(k,i)
             END DO
         END DO
     END DO

     !Calculate determinant by finding product of diagonal elements
     determinant= l
     DO i = 1, n
         determinant= determinant* matrix(i,i)
     END DO

END subroutine modulus
!.........................................................................................................


subroutine vibrational_energy

implicit none
real*8 :: F0,gammma,r_NO,V00,muu,vib_KE,TE_vib
!real*8, dimension(3) :: O_pos,N_pos
real*8, dimension(16) :: q_number
real*8, dimension(16) :: closest
!integer, dimension(16) :: quantum_count
integer :: q_num,cls,op,np
F0=638.5
gammma=2.743
r_NO=1.1507*1.88973

do np=1,3
N_pos(np)=pos(np)
O_pos(np)=pos(np+3)
enddo






V00=F0*(1-exp(-gammma*(Norm2(N_pos-O_pos)-r_NO)))**2

!V001=F0*(1-exp(-gammma*(sqrt((pos(1)-pos(2))**2)-r_NO)))**2

muu=mass(1)*mass(4)/(mass(1)+mass(4))

vib_KE=0.5*muu*(v(1)**2+v(2)**2+v(3)**2+v(4)**2+v(5)**2+v(6)**2)


TE_vib=vib_KE+V00


do q_num=1,16
open(56,file="q_number.dat")
read(56,*) q_number(q_num)
close(56)
enddo

quantum_count=0

do cls=1,16
closest(cls)=abs(TE_vib-q_number(cls))
end do


quantum_count(minloc(closest))=quantum_count(minloc(closest))+1



 end subroutine vibrational_energy 




!.......................................................................................................................
subroutine Boundary_conditions
integer :: i,j
real*8 :: posa(530,3)
real*8 :: Lx,Ly

Lx=32.5E-10
Ly=30.7E-10


call reshaping_array(posa,pos,2)

do i=1,530
   if (posa(i,1)>Lx) then 
        posa(i,1)=posa(i,1)-Lx
   end if
   if (posa(i,2)>Ly) then
        posa(i,2)=posa(i,2)-Ly
   end if



   if (posa(i,1)<0) then
        posa(i,1)=posa(i,1)+Lx
   end if

   if (posa(i,2)<0) then
        posa(i,2)=posa(i,2)+Ly
   end if

enddo

call reshaping_array(posa,pos,1)
    

end subroutine Boundary_conditions
!.............................................................................................................................................................................
Subroutine Make_Video
real*8 :: pos_cordinate(530,3)
integer :: i
call reshaping_array(pos_cordinate,pos,2)

!open(206, file='VMD_file.xyz')
write(206,*) 'N', pos_cordinate(1,:)/1E-10
write(206,*) 'O', pos_cordinate(2,:)/1E-10
do i=3,530
   write(206,*) '3', pos_cordinate(i,:)/1E-10
enddo
!close(206)

call reshaping_array(pos_cordinate,pos,1)













end subroutine Make_Video
!.............................................................................................................................................................................
end module atom_interaction
