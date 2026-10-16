! Reference-scaled two-parameter tumour--immune--DC model.
! PAR(1)=sigma=eta*sT/r0^2, PAR(9)=rho=r/r0.
! r0 is the internal reference clock; runtime parameters come from the inputs.
! Thus changing rho does not silently alter the other reference parameters.
subroutine FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
  implicit none
  integer NDIM,IJAC,ICP(*),i
  double precision U(NDIM),PAR(*),F(NDIM),DFDU(*),DFDP(*)
  double precision x,y,z,sigma,rho,alpha,delta,beta,a,gamma,kappa,sat
  double precision xn,kn,h,hx,q,qa
  x=U(1)
  y=U(2)
  z=U(3)
  sigma=PAR(1)
  alpha=PAR(2)
  delta=PAR(3)
  beta=PAR(4)
  a=PAR(5)
  gamma=PAR(6)
  kappa=PAR(7)
  rho=PAR(9)
  sat=PAR(10)
  xn=x**4
  kn=kappa**4
  h=xn/(kn+xn)
  q=x*y
  F(1)=rho*x*(1d0-x)-x*y
  F(2)=sigma*h+alpha*z*y/(1d0+z)-delta*y
  F(3)=beta*q/(a+sat*q)-gamma*z
  if (IJAC.ge.1) then
    hx=4d0*kn*x**3/(kn+xn)**2
    qa=beta*a/(a+sat*q)**2
    DFDU(1)=rho*(1d0-2d0*x)-y
    DFDU(2)=sigma*hx
    DFDU(3)=qa*y
    DFDU(4)=-x
    DFDU(5)=alpha*z/(1d0+z)-delta
    DFDU(6)=qa*x
    DFDU(7)=0d0
    DFDU(8)=alpha*y/(1d0+z)**2
    DFDU(9)=-gamma
  endif
  if (IJAC.ge.2) then
    do i=1,33
      DFDP(i)=0d0
    enddo
    DFDP(2)=h
    DFDP(5)=z*y/(1d0+z)
    DFDP(8)=-y
    DFDP(12)=q/(a+sat*q)
    DFDP(15)=-beta*q/(a+sat*q)**2
    DFDP(18)=-z
    DFDP(20)=-4d0*sigma*kappa**3*xn/(kn+xn)**2
    DFDP(25)=x*(1d0-x)
    DFDP(30)=-beta*q*q/(a+sat*q)**2
  endif
end subroutine FUNC

subroutine STPNT(NDIM,U,PAR,T)
  implicit none
  integer NDIM
  double precision U(NDIM),PAR(*),T
  PAR(1)=1.60000000000000003d-02
  PAR(2)=3.75000000000000000d+00
  PAR(3)=5.99999999999999978d-01
  PAR(4)=3.05967261904761978d-01
  PAR(5)=1.00000000000000000d+00
  PAR(6)=5.00000000000000000d-01
  PAR(7)=1.50000000000000022d+00
  PAR(9)=1.94360750690608763d+00
  PAR(10)=1d0
  U(1)=1d0
  U(2)=0d0
  U(3)=0d0
end subroutine STPNT

subroutine BCND
end subroutine BCND
subroutine ICND
end subroutine ICND
subroutine FOPT
end subroutine FOPT
subroutine PVLS(NDIM,U,PAR)
  implicit none
  integer NDIM
  double precision U(*),PAR(*),GETP
  PAR(8)=GETP('MIN',1,U)
end subroutine PVLS
