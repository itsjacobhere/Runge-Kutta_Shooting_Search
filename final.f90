module var
    implicit none
    integer, parameter :: wp = selected_real_kind(15)
    real(wp) :: tstart, tend ,h, hmin, hmax, epsi1, epsi2
    integer :: i, k, l , j, exitflag, neqn
    real(wp) :: h1, snew
    !function variables:
    real(wp) :: mu
    real(wp), parameter :: ka1=0.003, kp1=0.05,kp2=1.0,kp3=0.03, mumax=0.43
    real(wp) :: timereport, dreportinteger, velocity, distance
    !variable declaration:
    !ka1 = the algal death rate
    !kp1 = the constant settling rate of phosphorus
    !kp2 = the yield coefficient
    !kp3 = saturation concentration of P
    !mumax = max algal growth rate
    !h = step size controlled with truncation error in rfk sub
    !hmin = minimum step size within rung-kutta-fehlberg
    !hmax = maximum step size within rkf algo changes with h control
    !tstart = starting time initialized at zero
    !tend = time at 44.5 km
    !epsi1-2 = values to compare error and make solution more accurate
    !neqn = number of equations in the system 
    !c values & a values = coeffs descirbed in chapra and kanale pg 747
    !e = trunctation error found from: 5th order - 4th order estimate , 
    !t = time
    !
end module
!_______________________________________________________________________________
Module Rung_Kutta
use var
implicit none
contains
    subroutine rungk(tstart,tend,y, epsi1,epsi2, h, hmin, hmax, D, exitflag, neqn)
        implicit none
        integer, intent(in) :: neqn
        integer, intent(out) :: exitflag
        real(wp), intent(in) :: hmin, hmax, epsi1, epsi2, tstart, tend
        real(wp), intent(inout) :: h
        !local vars
        real(wp) ::  t,e, hsave
        real(wp), dimension(neqn) :: k1,k2,k3,k4,k5,k6,y4,y5,ysave
        !see page 747 Chapra Runge-Kutta Fehlberg
        real(wp), parameter :: c1=37._wp/378._wp, c2 = 250._wp/621._wp,c3=125._wp/594._wp,c4=512._wp/1771._wp
        real(wp), parameter :: c5 = 2825._wp/27648._wp
        real(wp), parameter :: c6 = 18575._wp/48384._wp,c7=13525._wp/55296._wp,c8=277._wp/14336._wp,c9=1._wp/4._wp
        real(wp), parameter :: a1=1._wp/5._wp, a2 =3._wp/10._wp,a22 = 3._wp/40._wp,a3=9._wp/40._wp,a4=3._wp/5._wp
        real(wp), parameter :: a5=-9._wp/10._wp, a6 = 6._wp/5._wp,a7=-11._wp/54._wp,a8=5._wp/2._wp
        real(wp), parameter :: a9=-70._wp/27._wp, a10 = 35._wp/27._wp,a11=7._wp/8._wp,a12=1631._wp/55296._wp
        real(wp), parameter :: a13=175._wp/512._wp, a14 = 575._wp/13824._wp
        real(wp), parameter :: a15=44275._wp/110592._wp,a16=253._wp/4096._wp
        real(wp), dimension(:), intent(inout) :: y
        !interface to function here
        interface
            Function D(t,y,neqns)
                use var
                implicit none
                integer, intent(in) :: neqns
                real(wp), intent(in) :: t
                real(wp), dimension(:), intent(in) :: y
                real(wp), dimension(neqns) :: D
            end function D
        end interface
        t = 0
        !checking to make sure answer makes sense
        !if (tstart>tend) .or. (epsi1>epsi2) .or. (hmin>hmax) exit
        if (tstart>tend) then
            write(*,*) "tstart cannot be greater than tend"
            return
        else if (epsi1>epsi2) then
            write(*,*) "epsi1 must be smaller than epsi2"
            return
        else if (hmin>hmax) then
            write(*,*) "hmin should be less than hmax"
            return
        end if 
        if (h>(tend-tstart)) then
            h = tend-tstart
        end if
        hsave = h
        t = tstart !time begins at zero
        do 
            j = j+1 !tracking iterations
            !write(*,*) t, y
            ysave=y(:neqn)!!
            !RFK fehlberg eqns
            k1 = D(t,y,neqn)
            k2 = D((t+(a1*h)), (y+(a1*k1*h)), neqn)
            k3 = D((t+(a2*h)), (y+(a22*k1*h)+(a3*k2*h)), neqn)
            k4 = D((t+(a4*h)), (y+(a2*k1*h)+(a5*k2*h)+(a6*k3*h)), neqn)
            k5 = D((t+h), (y+(a7*k1*h)+(a8*k2*h)+(a9*k3*h)+(a10*k4*h)),neqn)
            k6 = D(t+(a11*h),y+(a12*k1*h)+(a13*k2*h)+(a14*k3*h)+(a15*k4*h)+(a16*k5*h),neqn)
            y4 = y+(c1*k1*h)+((c2*k3*h)+(c3+k4*h)+(c4*k6*h))
            y5 = y+((c5*k1*h)+(c6*k3*h)+(c7*k4*h)+(c8*k5*h)+(c9*k6*h))
            !compute error estimate
            e = maxval((abs(y5-y4))/y5)
            !controlling error
                if (e>epsi2) then
                    if (abs(hmin-h)<0.000001) then
                        write(*,*) "step size (h) is too small"
                        return
                    !end if
                    h = h/2._wp
                    if (h<hmin) then 
                        h=hmin
                    end if
                    y = ysave
                else
                    t = t + h
                    y = y5
                    if(e<epsi1) then
                        h = h*2 !double step size
                        if (h>hmax) h = hmax
                    end if
                    !within this range calculations are no longer necessary
                    if(abs(t-tend)<1.e-6) then
                        exitflag = 10
                        return
                    end if
                    !t cannot exceed tend
                    if (t+h>tend) then
                        hsave = h
                        h = tend-t
                        exit
                    end if
                end if
            end if
        end do
    end subroutine rungk
end module Rung_Kutta
!__________________________________________________________________________________________________
module shooting
    use var
    use Rung_Kutta
    implicit none
    contains
    subroutine shootingmeth(solder,esolder,sold,esold,lim)
        real(wp), intent(inout) :: solder, esolder
        real(wp), intent(inout) :: sold, esold
        real(wp), intent(in) :: lim
        !local vars
        real(wp)::  esnew
        integer :: iteration
        real(wp), dimension(neqn) :: y
        character(len=20) :: afmt
        !interface function
        interface
            Function D(t,y,neqns) !function 
                use var
                implicit none
                integer, intent(in) :: neqns
                real(wp), intent(in) :: t
                real(wp), dimension(:), intent(in) :: y
                real(wp), dimension(neqns) :: D
            end function D
        end interface
        !------------------------------------------------------------------------------------------
        iteration = 3 !start at third iteration
        
        do

            h = .001
            snew = sold - esold*((sold - solder) / (esold-esolder))
            
            !write(*,'(I0,",",f10.4)') iteration, snew
            write(20,'(I0,",",f10.4)') iteration, snew
            
            y(1) = snew
            y(2) = 0.002
            call rungk(tstart,tend,y, epsi1,epsi2, h, hmin, hmax, D, exitflag, neqn)
            esnew = y(1) - lim !new error
            if ((abs(esnew)) <= 0.00001) then !eps=0.0001 precision of error
                write(*,*) "shooting method complete."
                return
                exit
            end if
            if (iteration > 200) then !itermax = 500 iteration limit
                write(*,*) "Something is wrong, too many iterations."
                STOP
            else !updating snew
                solder = sold
                esolder = esold
                sold = snew
                esold = esnew
                iteration = iteration +1
            end if
        end do
    end subroutine shootingmeth
end module shooting
!___________________________________________________________________________________________________
program wastewater
    use var
    use Rung_Kutta
    use shooting
    implicit none
    interface
        Function D(t,y,neqns)
            use var
            implicit none
            integer, intent(in) :: neqns
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(neqns) :: D
        end function D
    end interface
    real(wp) :: solder, esolder
    real(wp) :: sold, esold
    real(wp), dimension(:), allocatable :: y
    !real(wp), dimension(:,:), allocatable :: f 
    
    real(wp) :: p_in, lim,conc_algal
    integer :: nsteps, iterations, iteration
    character(len=40) :: output="plotsvals9.csv",afmt
!___________________________________________________________________
    !create output files to write results to.
    open(unit=20,file=output,status="replace")
    !initializations:
    iterations = 89
    !target limit
    lim = 0.05_wp
    neqn = 2

    y = [0.9 , 0.002] !initial phosphorus concentration, and algal uptake
    p_in = y(1)
    conc_algal = y(2)
    !precision
    epsi1 = 0.000001
    epsi2 = 0.00001
    nsteps = iterations
    !allocate(f(nsteps+1,3)) !allocate final array
    !parameters
    velocity = 0.6_wp !(m/s)
    distance = 44500._wp !(meters)
    tstart = 0._wp !start at time zero
    tend = (distance/velocity)/86400._wp !calc time and convert to days
    
    timereport = (tend-tstart)/nsteps !time interval
    dreportinteger = distance/nsteps !distance interval
    
    hmax = tend
    hmin = hmax/100 !try diff values
    h = .01
    !**************************************************************************
    
    
    
    !begin shooting method
    solder = p_in
    call rungk(tstart,tend,y, epsi1,epsi2, h, hmin, hmax, D, exitflag, neqn)

    esolder = y(1) - lim !initial difference 
    
    if (esolder>0) then
        sold = 0.8*solder
    else
        sold = 1.2*solder
    end if
    
    y(1) = sold
    y(2) = .002
    
    !write(*,*) "before finding esold, y=",y
    call rungk(tstart,tend,y, epsi1,epsi2, h, hmin, hmax, D, exitflag, neqn)

    esold = y(1) - lim
    !write first two iterations here
    afmt = '(I0,",",f10.4)'
    write(20,*) "Iteration,S_value"
    write(20,'(A,",",f10.4)') "1", solder
    write(20,'(A,",",f10.4)') "2", sold
!    write(*,*) "Iteration,S_value"
!    write(*,'(A,",",f10.4)') "1", solder
!    write(*,'(A,",",f10.4)') "2", sold
    
    call shootingmeth(solder,esolder,sold,esold,lim)!
    write(*,*) "snew=", snew



end program wastewater
!______________________________________________________________________
Function D(t,y,neqns) !function with the problem specific equations
    use var
    implicit none
    integer, intent(in) :: neqns
    real(wp), intent(in) :: t
    real(wp), dimension(:), intent(in) :: y
    real(wp), dimension(neqns) :: D
    mu = mumax*(y(1)/(kp3+y(1))) !algal growth rate
    !System of differential equations
    D(1) = (-kp1*y(1)) - (kp2*mu*y(2))
    D(2) = (-ka1*y(2)) + (mu*y(2))
end function D