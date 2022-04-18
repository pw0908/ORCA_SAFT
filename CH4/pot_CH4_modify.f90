! CH4 potential energy code, adapted from the original code written by Dr Robert Hellmann
! the adaptation facilitates the potential average calculations
module potpar_CH4

real(8),dimension(9,9) :: A,alpha,C6,C8,b
real(8)                :: dC6,dC8,b_iso,rfac(8)

real(8),parameter,dimension(1:33) :: p = &       
 (/  0.262373609870D+07, & ! A
     0.265413949306D+07, &
     0.241399202548D+06, &
    -0.271732286369D+06, &
    -0.749715217817D+05, &
     0.123654939014D+06, &
     0.168784214727D+01, & ! alpha
     0.288272189383D+01, &
     0.359175610364D+01, &
     0.164907472587D+01, &
     0.205930856206D+01, &
     0.214516406393D+01, &
     0.112317356393D+07, & ! C6
    -0.139633536709D+07, &
     0.294147229989D+06, &
     0.127844394118D+07, &
     0.169329264726D+06, &
    -0.590727146203D+06, &
    -0.120939118595D+09, & ! C8
     0.385078059812D+08, &
    -0.264781785792D+07, &
     0.174762763583D+07, &
    -0.810401688137D+07, &
     0.679543866185D+07, &
     0.168275674697D+01, & ! b
     0.288261053675D+01, &
     0.384703187791D+01, &
     0.155011959551D+01, &
     0.266424602826D+01, &
     0.304993944136D+01, &
     45346.8d0         , & ! dC6
     432463.3d0        , & ! dC8
     1.77d0             /) ! bcorr

real(8),parameter,dimension(1:9) :: x0 = &
 (/  0.000000000000d0, &
     0.000000000000d0, &
    -0.911809480295d0, &
     0.455904740147d0, &
     0.455904740147d0, &
     0.000000000000d0, &
     0.683857110221d0, &
    -0.341928555111d0, &
    -0.341928555111d0 /)

real(8),parameter,dimension(1:9) :: y0 = &
 (/  0.000000000000d0, &
     0.000000000000d0, &
     0.000000000000d0, &
    -0.789650173347d0, &
     0.789650173347d0, &
     0.000000000000d0, &
     0.000000000000d0, &
     0.592237630010d0, &
    -0.592237630010d0 /)

real(8),parameter,dimension(1:9) :: z0 = &
 (/  0.000000000000d0, &
     0.967120000000d0, &
    -0.322373333333d0, &
    -0.322373333333d0, &
    -0.322373333333d0, &
    -0.725340000000d0, &
     0.241780000000d0, &
     0.241780000000d0, &
     0.241780000000d0 /)

real(8),parameter,dimension(1:9) :: q = &
 (/ -0.379012d3, &
     0.947530d2, &
     0.947530d2, &
     0.947530d2, &
     0.947530d2, &
     0.000000d0, &
     0.000000d0, &
     0.000000d0, &
     0.000000d0 /)

end module potpar_CH4

!------------------------------------------------------------------------------

subroutine Test_FN(NDim,X,NumFun,ssExpV)
! use sin_exp_pot_CH4
implicit none

real(8),intent(in) :: NDim,NumFun
real(8),dimension(1:5),intent(in) :: X
real(8),intent(out) :: ssExpV

ssExpV = 1.0d0

end subroutine Test_FN

!------------------------------------------------------------------------------

subroutine init_CH4
use potpar_CH4

implicit none
integer :: i,j

A(1,1)     = p(1)
alpha(1,1) = p(7)
C6(1,1)    = p(13)
C8(1,1)    = p(19)
b(1,1)     = p(25)

do i = 2,5
  A(1,i)     = p(2)
  A(i,1)     = p(2)
  alpha(1,i) = p(8)
  alpha(i,1) = p(8)
  C6(1,i)    = p(14)
  C6(i,1)    = p(14)
  C8(1,i)    = p(20)
  C8(i,1)    = p(20)
  b(1,i)     = p(26)
  b(i,1)     = p(26)
enddo

do i = 2,5
  do j = 2,5
    A(i,j)     = p(3)
    alpha(i,j) = p(9)
    C6(i,j)    = p(15)
    C8(i,j)    = p(21)
    b(i,j)     = p(27)
  enddo
enddo

do i = 6,9
  A(1,i)     = p(4)
  A(i,1)     = p(4)
  alpha(1,i) = p(10)
  alpha(i,1) = p(10)
  C6(1,i)    = p(16)
  C6(i,1)    = p(16)
  C8(1,i)    = p(22)
  C8(i,1)    = p(22)
  b(1,i)     = p(28)
  b(i,1)     = p(28)
enddo

do i = 2,5
  do j = 6,9
    A(j,i)     = p(5)
    A(i,j)     = p(5)
    alpha(j,i) = p(11)
    alpha(i,j) = p(11)
    C6(j,i)    = p(17)
    C6(i,j)    = p(17)
    C8(j,i)    = p(23)
    C8(i,j)    = p(23)
    b(j,i)     = p(29)
    b(i,j)     = p(29)
  enddo
enddo

do i = 6,9
  do j = 6,9
    A(i,j)     = p(6)
    alpha(i,j) = p(12)
    C6(i,j)    = p(18)
    C8(i,j)    = p(24)
    b(i,j)     = p(30)
  enddo
enddo

dC6   = p(31)
dC8   = p(32)
b_iso = p(33)

rfac(1) = 1.0d0
do i = 2,8
  rfac(i) = rfac(i-1)/dble(i) 
enddo

end subroutine init_CH4

! modified based on Hellmann's code
subroutine sin_exp_pot_CH4(Rcom,T,ThetaA,ThetaB,PsiA,PsiB,Phi,ssExpV)
use potpar_CH4

implicit none

! Rcom (distance between centers of mass) in Angstrom, angles in radians
real(8),intent(in) :: Rcom,T,ThetaA,ThetaB,PsiA,PsiB,Phi

! V in Kelvin
real(8) :: V
real(8),intent(out) :: ssExpV

real(8) :: sinThetaA,sinThetaB,sinPsiA,sinPsiB,sinPhi
real(8) :: cosThetaA,cosThetaB,cosPsiA,cosPsiB,cosPhi

real(8),dimension(9) :: xa,ya,za,xb,yb,zb

real(8) :: EulerA(3,3),EulerB(3,3)
real(8) :: Dx,Dy,Dz
real(8) :: r,rk,rk2,rk3,rk6,rk8
real(8) :: br,dampexp,brpowk(0:8),dampsum(0:8)

integer :: i,j,k

sinThetaA = sin(ThetaA)
sinThetaB = sin(ThetaB)
sinPhi    = sin(Phi)
sinPsiA   = sin(PsiA)
sinPsiB   = sin(PsiB)

cosThetaA = cos(ThetaA) 
cosThetaB = cos(ThetaB) 
cosPhi    = cos(Phi)
cosPsiA   = cos(PsiA)
cosPsiB   = cos(PsiB)

EulerA(1,1) = cosThetaA*cosPsiA
EulerA(1,2) = -cosThetaA*sinPsiA
EulerA(1,3) = sinThetaA

EulerA(2,1) = sinPsiA
EulerA(2,2) = cosPsiA
EulerA(2,3) = 0.0d0

EulerA(3,1) = -sinThetaA*cosPsiA
EulerA(3,2) = sinThetaA*sinPsiA
EulerA(3,3) = cosThetaA

EulerB(1,1) = -sinPhi*sinPsiB + cosPhi*cosThetaB*cosPsiB
EulerB(1,2) = -sinPhi*cosPsiB - cosPhi*cosThetaB*sinPsiB
EulerB(1,3) = cosPhi*sinThetaB

EulerB(2,1) = cosPhi*sinPsiB + sinPhi*cosThetaB*cosPsiB
EulerB(2,2) = cosPhi*cosPsiB - sinPhi*cosThetaB*sinPsiB
EulerB(2,3) = sinPhi*sinThetaB

EulerB(3,1) = -sinThetaB*cosPsiB
EulerB(3,2) = sinThetaB*sinPsiB
EulerB(3,3) = cosThetaB

V = 0.0d0

do j = 1,9

  xb(j) = EulerB(1,1)*x0(j) + EulerB(1,2)*y0(j) + EulerB(1,3)*z0(j)
  yb(j) = EulerB(2,1)*x0(j) + EulerB(2,2)*y0(j) + EulerB(2,3)*z0(j)
  zb(j) = EulerB(3,1)*x0(j) + EulerB(3,2)*y0(j) + EulerB(3,3)*z0(j) + Rcom

  do i = 1,9

    xa(i) = EulerA(1,1)*x0(i) + EulerA(1,2)*y0(i) + EulerA(1,3)*z0(i)
    ya(i) = EulerA(2,1)*x0(i) + EulerA(2,2)*y0(i) + EulerA(2,3)*z0(i)
    za(i) = EulerA(3,1)*x0(i) + EulerA(3,2)*y0(i) + EulerA(3,3)*z0(i)

    Dx = xb(j) - xa(i)
    Dy = yb(j) - ya(i)
    Dz = zb(j) - za(i)
    
    r = sqrt(Dx*Dx + Dy*Dy + Dz*Dz)
    
    rk  = 1.0d0/r
    rk2 = rk*rk
    rk3 = rk2*rk
    rk6 = rk3*rk3
    rk8 = rk6*rk2

    br         = b(i,j) * r
    dampexp    = exp(-br)
    dampsum(0) = 1.0d0
    brpowk(0)  = 1.0d0
    do k = 1,8
      brpowk(k)  = brpowk(k-1) * br
      dampsum(k) = dampsum(k-1) + brpowk(k) * rfac(k)
    enddo
    
    V = V + A(i,j) * exp(-alpha(i,j)*r) &
          - C6(i,j) * rk6 * (1.0d0 - dampexp * dampsum(6)) &
          - C8(i,j) * rk8 * (1.0d0 - dampexp * dampsum(8)) &
          + q(i) * q(j) * rk

    if(i == 1 .and. j == 1) then
    
      br         = b_iso * r
      dampexp    = exp(-br)
      dampsum(0) = 1.0d0
      brpowk(0)  = 1.0d0
      do k = 1,8
        brpowk(k)  = brpowk(k-1) * br
        dampsum(k) = dampsum(k-1) + brpowk(k) * rfac(k)
      enddo

      V = V - dC6 * rk6 * (1.0d0 - dampexp * dampsum(6)) &
            - dC8 * rk8 * (1.0d0 - dampexp * dampsum(8))
      
    endif

  enddo
enddo

ssExpV = sinThetaA*sinThetaB*exp(-V/T)

end subroutine sin_exp_pot_CH4

subroutine pot_CH4(Rcom,ThetaA,ThetaB,PsiA,PsiB,Phi,V)
use potpar_CH4

implicit none

! Rcom (distance between centers of mass) in Angstrom, angles in radians
real(8),intent(in) :: Rcom,ThetaA,ThetaB,PsiA,PsiB,Phi

! V in Kelvin
real(8),intent(out) :: V

real(8) :: sinThetaA,sinThetaB,sinPsiA,sinPsiB,sinPhi
real(8) :: cosThetaA,cosThetaB,cosPsiA,cosPsiB,cosPhi

real(8),dimension(9) :: xa,ya,za,xb,yb,zb

real(8) :: EulerA(3,3),EulerB(3,3)
real(8) :: Dx,Dy,Dz
real(8) :: r,rk,rk2,rk3,rk6,rk8
real(8) :: br,dampexp,brpowk(0:8),dampsum(0:8)

integer :: i,j,k

sinThetaA = sin(ThetaA)
sinThetaB = sin(ThetaB)
sinPhi    = sin(Phi)
sinPsiA   = sin(PsiA)
sinPsiB   = sin(PsiB)

cosThetaA = cos(ThetaA) 
cosThetaB = cos(ThetaB) 
cosPhi    = cos(Phi)
cosPsiA   = cos(PsiA)
cosPsiB   = cos(PsiB)

EulerA(1,1) = cosThetaA*cosPsiA
EulerA(1,2) = -cosThetaA*sinPsiA
EulerA(1,3) = sinThetaA

EulerA(2,1) = sinPsiA
EulerA(2,2) = cosPsiA
EulerA(2,3) = 0.0d0

EulerA(3,1) = -sinThetaA*cosPsiA
EulerA(3,2) = sinThetaA*sinPsiA
EulerA(3,3) = cosThetaA

EulerB(1,1) = -sinPhi*sinPsiB + cosPhi*cosThetaB*cosPsiB
EulerB(1,2) = -sinPhi*cosPsiB - cosPhi*cosThetaB*sinPsiB
EulerB(1,3) = cosPhi*sinThetaB

EulerB(2,1) = cosPhi*sinPsiB + sinPhi*cosThetaB*cosPsiB
EulerB(2,2) = cosPhi*cosPsiB - sinPhi*cosThetaB*sinPsiB
EulerB(2,3) = sinPhi*sinThetaB

EulerB(3,1) = -sinThetaB*cosPsiB
EulerB(3,2) = sinThetaB*sinPsiB
EulerB(3,3) = cosThetaB

V = 0.0d0

do j = 1,9

  xb(j) = EulerB(1,1)*x0(j) + EulerB(1,2)*y0(j) + EulerB(1,3)*z0(j)
  yb(j) = EulerB(2,1)*x0(j) + EulerB(2,2)*y0(j) + EulerB(2,3)*z0(j)
  zb(j) = EulerB(3,1)*x0(j) + EulerB(3,2)*y0(j) + EulerB(3,3)*z0(j) + Rcom

  do i = 1,9

    xa(i) = EulerA(1,1)*x0(i) + EulerA(1,2)*y0(i) + EulerA(1,3)*z0(i)
    ya(i) = EulerA(2,1)*x0(i) + EulerA(2,2)*y0(i) + EulerA(2,3)*z0(i)
    za(i) = EulerA(3,1)*x0(i) + EulerA(3,2)*y0(i) + EulerA(3,3)*z0(i)

    Dx = xb(j) - xa(i)
    Dy = yb(j) - ya(i)
    Dz = zb(j) - za(i)
    
    r = sqrt(Dx*Dx + Dy*Dy + Dz*Dz)
    
    rk  = 1.0d0/r
    rk2 = rk*rk
    rk3 = rk2*rk
    rk6 = rk3*rk3
    rk8 = rk6*rk2

    br         = b(i,j) * r
    dampexp    = exp(-br)
    dampsum(0) = 1.0d0
    brpowk(0)  = 1.0d0
    do k = 1,8
      brpowk(k)  = brpowk(k-1) * br
      dampsum(k) = dampsum(k-1) + brpowk(k) * rfac(k)
    enddo
    
    V = V + A(i,j) * exp(-alpha(i,j)*r) &
          - C6(i,j) * rk6 * (1.0d0 - dampexp * dampsum(6)) &
          - C8(i,j) * rk8 * (1.0d0 - dampexp * dampsum(8)) &
          + q(i) * q(j) * rk

    if(i == 1 .and. j == 1) then
    
      br         = b_iso * r
      dampexp    = exp(-br)
      dampsum(0) = 1.0d0
      brpowk(0)  = 1.0d0
      do k = 1,8
        brpowk(k)  = brpowk(k-1) * br
        dampsum(k) = dampsum(k-1) + brpowk(k) * rfac(k)
      enddo

      V = V - dC6 * rk6 * (1.0d0 - dampexp * dampsum(6)) &
            - dC8 * rk8 * (1.0d0 - dampexp * dampsum(8))
      
    endif

  enddo
enddo

end subroutine pot_CH4