! the main program
program CH4_Int
implicit none
External FUNSUB_pot_CH4

logical           :: exist
real(8)           :: expV
integer             :: ir,iNr
character(len=3)    :: arg2
character(len=1024) :: filename
integer,parameter   :: Nr = 50
real(8),parameter   :: rInit = 3.15d0 
real(8),parameter   :: rFinal = 7.00d0 
real(8),allocatable :: rVal(:)
real(8),Dimension(1) :: result,abserr
integer           :: nsub,neval,ifail
real(8),parameter :: pi = 3.141592653589793238d0
real(8),parameter :: pid180 = pi / 180.0d0
integer           :: NMAX = 500000000
integer           :: NW 
real(8),Dimension(1:5) :: LowLim,UpLim
real(8),allocatable :: WorkArray(:)

allocate(rVal(Nr))
rVal = (/(dble(iNr-1)*(rFinal-rInit)/dble(Nr-1)+rInit, iNr= 1,Nr)/)

NW = 10*NMAX
allocate(WorkArray(NW))

call get_command_argument(1,filename)
call get_command_argument(2,arg2)
read(arg2,*)  ir

LowLim = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)
UpLim  = (/pi,pi,2.0d0*pi,2.0d0*pi,2.0d0*pi/)

call init_CH4

call DCUHRE(5, 1, LowLim, UpLim, 1, NMAX, FUNSUB_pot_CH4, 0.0d0,  &
            0.001d0, 0, NW, 0, result, abserr, neval, ifail, WorkArray)

inquire(file=filename, exist=exist)
if (exist) then
  open(1, file=filename, status="old", position="append", action="write")
else
  open(1, file=filename, status="new", action="write")
  write(1, '(A)') '# r/A   V_direct   abserr    NEval   Exit_statues'
end if

Write(1, '(F15.5,F15.4, F20.4, I16, I4)') &
     rVal(ir), result/(32.0d0*pi*pi*pi), abserr/(32.0d0*pi*pi*pi), NEval, iFail
close(1) 

end program CH4_Int

subroutine FUNSUB_pot_CH4(NDim,X,NumFun,Vr)
! use sin_exp_pot_CH4
implicit none

integer             :: ir,iNr, iT
character(len=3)    :: arg2
integer,parameter   :: Nr = 50
real(8),parameter   :: rInit = 3.15d0 
real(8),parameter   :: rFinal = 7.00d0 
real(8),allocatable :: rVal(:)
real(8),intent(in)                :: NDim,NumFun
real(8),intent(in),dimension(1:5) :: X
real(8),intent(out)               :: Vr


allocate(rVal(Nr))
rVal = (/(dble(iNr-1)*(rFinal-rInit)/dble(Nr-1)+rInit, iNr= 1,Nr)/)

call get_command_argument(2,arg2)
read(arg2,*)  ir
call pot_CH4(rVal(ir),X(1),X(2),X(3),X(4),X(5),Vr)
Vr = Vr * sin(X(1))*sin(X(2))

end subroutine FUNSUB_pot_CH4