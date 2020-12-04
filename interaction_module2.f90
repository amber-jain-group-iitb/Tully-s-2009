module atom_interaction
        implicit none
        integer::i,j,k,l,Au_nei_list(528,12),step_dynamics,step_equilibrium,step,vib_state,Au_i(528),n_O_Au,n_N_Au,n_Au,traj
        real*8,dimension(396,3)::acc_Au,acc_Au_old,Au_vel
        real*8,dimension(528,3)::Au,Au_eq,gold_vel
        real*8,dimension(3)::N,O,N_vel,O_vel,COM,COM_v,acc_N,acc_N_old,acc_O,acc_O_old
        real*8::r_Au_O,r_Au_N
        integer,dimension(528)::N_Au_list,O_Au_list
        real*8::r_N_O

        real*8::a0,alpha0,b0,beta0,f0,gamma0,r0_N_O,a1,alpha1,b1,beta1,r1_Au_N,d,c,zimage,f1,gamma1,r1_N_O,phi,ea
        real*8::a2,a3,gamma2,gamma3,b2,b3,r_cut

        real*8::mass_N,mass_O,red_mass,av,kb,temp,conv,wavenumber_2_j,amu2kg,c_light,wavenumber,mass_Au
        real*8::e_kin,pi,planck_const,total_mass,time

        real*8,dimension(528,528,3,3)::ten_mat
        real*8::l_x,r_x,l_y,r_y,l_z,r_z,box_len_x,box_len_y,box_len_z,force_const,dt1,dt,z_COM
        integer::iflag_write,iflag_writevmd,adsorbed,no_traj

        real*8::diab_pot(2,2),gr_vec(1,2),gr_vec_t(2,1),ex_vec(1,2),ex_vec_t(2,1),gr_eig,ex_eig
        real*8::der_pot_N(2,2,3),der_pot_O(2,2,3),der_pot_Au(528,3,2,2),pot_Au_Au(1,1),der_pot_Au_Au(528,3)
        real*8 :: vel(530,3),x(530,3)
        real*8 :: pos(1590)
        real*8 :: time_atomic_SI



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
        c=1.2423*1.d-10
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



        time_atomic_SI= 2.418884326509E-17


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
        tmp1=dsqrt(c*c+(z_COM-zimage)**2)
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
!real*8 :: vel(530,3),x(530,3)
integer ::n1,loop_x,loop_v,i1,j1,k1,l1,t1,read_vel,g1,t12


call main_sub

call set_NO2



do j1=1,3
x(1,j1)=N(j1)
x(2,j1)=O(j1)
end do
open(656,file="gold_positions_fort.20")
do loop_x=3,530
!        open(656,file="gold_positions_fort.20")
        read(656,*) x(loop_x,1), x(loop_x,2), x(loop_x,3)
end do
close(656)




call setting_initial_position(x,Au,2)


vel(1,1)=N_vel(1)
vel(1,2)=N_vel(2)
vel(1,3)=N_vel(3)

vel(2,1)=O_vel(1)
vel(2,2)=O_vel(2)
vel(2,3)=O_vel(3)

open(658,file="gold_velocities_fort.21")
do read_vel=3,530
   if (read_vel.le.396) then
        read(658,*) vel(read_vel,1), vel(read_vel,2), vel(read_vel,3)
   else
           vel(read_vel,:)=0
   end if

end do
close(658)

call potential_force_calculate

!call potential





end subroutine 
!........................................................................
subroutine setup_initial_values3
!parameter*8 :: time_atomic_SI= 2.418884326509E-17
integer :: i,j

call main_sub



call set_NO2

dt=10*time_atomic_SI
open(13,file="gold_positions_fort.20")
do i=1,528
    read(13,*) Au(i,1), Au(i,2), Au(i,3)
end do
close(13)

write(16,*) Au
write(17,*) N
write(18,*) O


gold_vel=0

open(114,file="gold_velocities_fort.21")
do j=1,396
   read(114,*) gold_vel(j,1), gold_vel(j,2), gold_vel(j,3)
enddo
close(114)
    

write(19,*) O_vel
write(20,*) N_vel
write(21,*) gold_vel




call potential_force_calculate

time=0





end subroutine setup_initial_values3
!.............................................................................................................................................................................
subroutine velocity_verlet
real*8 :: old_O_acc(3), old_N_acc(3), old_Au_acc(528,3)
real*8 :: fix_Au(132,3)
integer :: i,j,k
O=O+O_vel*dt+0.5*(der_pot_O(1,1,:)/mass_O)*dt*dt
N=N+N_vel*dt+0.5*(der_pot_N(1,1,:)/mass_N)*dt*dt

do i=1,132
  fix_Au(i,:)=Au(396+i,:)
end do

Au=Au+gold_vel*dt +0.5*(der_pot_Au(:,:,1,1)/mass_Au)*dt*dt

do j=397,528
 Au(j,:)=fix_Au(j-396,:)
enddo


old_O_acc=(der_pot_O(1,1,:)/mass_O)
old_N_acc=(der_pot_N(1,1,:)/mass_N)
old_Au_acc=(der_pot_Au(:,:,1,1)/mass_Au)

call potential_force_calculate

O_vel=O_vel+0.5*((der_pot_O(1,1,:)/mass_O)+old_O_acc)*dt        
N_vel=N_vel+0.5*((der_pot_N(1,1,:)/mass_N)+old_N_acc)*dt
gold_vel=gold_vel+0.5*((der_pot_Au(:,:,1,1)/mass_Au)+old_Au_acc)*dt



write(15,*) der_pot_Au(:,:,2,2)         
write(25,*) der_pot_O(1,1,:)
write(26,*) der_pot_N(1,1,:)                                                                                                                                                                                                                   
write(215,*) mass_O, mass_N, mass_Au



end subroutine velocity_verlet
!............................................................................................................................................
subroutine classical_evolution
integer :: i
call setup_initial_values3
open(27,file='VMD_simple_MD.xyz')
do while(time<50000*time_atomic_SI)
   open(27,file='VMD_simple_MD.xyz') 
   write(27,*) 'N ', N/1E-10
   write(27,*) 'O ', O/1E-10
   do i=1,528
     write(27,*) 'Au', Au(i,:)/1E-10
   enddo
   call velocity_verlet
  
   time=time+dt           
end do
close(27)




end subroutine classical_evolution
!..............................................................................................................................
subroutine set_NO1              !NO is perpendicular to the surface. vibration only in z direction
        implicit none
        real*8::rel_x,rel_v,theta,phi,rnd,omega0,vib_energy,r_N_z,r_O_z,rand_ori
        
        call random_number(rnd)
        COM(1)=l_x+box_len_x*rnd*0.33d0+box_len_x*0.33d0
        call random_number(rnd)
        COM(2)=l_y+box_len_y*rnd*0.33d0+box_len_y*0.33d0
        COM(3)=10d-10
        
        omega0=dsqrt(2.d0*f0*gamma0**2/red_mass)
        !vib_energy=(real(vib_state)+0.5d0)*wavenumber*wavenumber_2_j
        vib_energy=0.5*planck_const*omega0/(2.d0*pi)
!        write(*,*)'vib_energy',vib_energy/conv
        call random_number(rnd)
        theta=pi*rnd*2.d0
        rel_x=dsqrt(2.d0*vib_energy/(red_mass*omega0**2))*dcos(theta)
        rel_v=dsqrt(2.d0*vib_energy/red_mass)*dsin(theta)
        write(9,*)'rel_x',rel_x,'rel_v',rel_v
        call random_number(rnd)
        if(rnd<0.5d0) rand_ori=1.d0
        if(rnd>0.5d0) rand_ori=-1.d0
        r_N_z=r0_N_O*mass_O/total_mass
        r_O_z=r0_N_O*mass_N/total_mass

        N(1)=COM(1)
        N(2)=COM(2)
        N(3)=COM(3)-(r0_N_O+rel_x)*mass_O/total_mass

        O(1)=COM(1)
        O(2)=COM(2)
        O(3)=COM(3)+(r0_N_O+rel_x)*mass_N/total_mass

        N_vel(1)=0.d0
        N_vel(2)=0.d0
        N_vel(3)=-rel_v*0.5d0
        O_vel(1)=0.d0
        O_vel(2)=0.d0
        O_vel(3)=+rel_v*0.5d0

        N_vel(3)=N_vel(3)-dsqrt(2.d0*e_kin/total_mass)
        O_vel(3)=O_vel(3)-dsqrt(2.d0*e_kin/total_mass)
end subroutine set_NO1






























!............................................................................................................................
end module atom_interaction
