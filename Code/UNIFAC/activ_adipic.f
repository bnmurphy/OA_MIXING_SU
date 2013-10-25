      program activ

c     This this code gives the activity coefficients for adipic acid/water
c     mixture using UNIFAC Dortmund

c     Anca Gaman 
c     January 2005


c     a for adipic acid
c     w for water

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!ADIPIC_ACID!!!!!!!!!!!!!!!!!!!
      !!!!!!!COOH-(CH2)_4-COOH!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!water!!!!!!!!!!!!!!!!!
      !!!!!!!H2O!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
        open(16,file='activ2.res',status='unknown') 
c      xn=1
C     1 stands for COOH
c     2 stands for CH2
c     3 stands for H2O
      !!!!!!!!MOLE FRACTIONS!!!!!!!!!!!!!!!!!!!
      xlow=0.
      xup=1.0
      step=0.01
      xa=xlow
      do
         if (xa.gt.xup) EXIT
         write(*,*) xa, aa, aw
         write(16,*) xa, gammaa*xa, gammaw*xw
         xa=xa+step
      xw=1.-xa    
      T=277.0
      ua1=2
      ua2=4
      uw1=1
      R1=0.8000
      R2=0.6325
      R3=1.7334
      Q1=0.9215
      Q2=0.7081
      Q3=2.4561
      a13=624.97
      b13=-4.6878
      c13=0.005237
      a31=-1795.2
      b31=12.708
      c31=-0.15146*1d-1
      a12=2017.7
      b12=-9.0933
      c12=0.1024*1d-1
      a21=1182.2
      b21=-3.2647
      c21=0.9198*1d-2
      a23=1391.3
      b23=-3.6156
      c23=0.1144*1d-2
      a32=-17.253
      b32=0.8389
      c32=0.9021*1d-3
      
      xi11=1
      xi12=exp(-(a12+b12*T+c12*T**2)/T)
      xi21=exp(-(a21+b21*T+c21*T**2)/T)
      xi22=1
      xi13=exp(-(a13+b13*T+c13*T**2)/T)
      xi23=exp(-(a13+b13*T+c13*T**2)/T)
      xi31=exp(-(a31+b31*T+c31*T**2)/T)
      xi32=exp(-(a31+b31*T+c31*T**2)/T)
      xi33=1
      !!!!!!!!!! COMBINATORIAL PART!!!!!!!!!!!!
 
      qa=ua1*Q1+ua2*Q2
      qw=uw1*Q1


      fa=qa/(xa*qa+xw*qw)
      fw=qw/(xa*qa+xw*qw)

      ra=ua1*R1+ua2*R2
      rw=uw1*R3

      va=ra/(xa*ra+xw*rw)
      vw=rw/(xa*ra+xw*rw)
      !!!parameter w (V' in the article)can be calculated using relative van der Waals
      !!! volumes Rk of different groups
      wa=ra**(0.75)/(xa*ra**(0.75)+xw*rw**(0.75))
      ww=rw**(0.75)/(xa*ra**(0.75)+xw*rw**(0.75))

      gammaac=exp(1-wa+log(wa)-5*qa*(1-va/fa+log(va/fa)))
      gammawc=exp(1-ww+log(ww)-5*qw*(1-vw/fw+log(vw/fw)))

      !!!!!!!!!! REZIDUAL PART !!!!!!!!!!!!!!!!!!

      X1ref=ua1/(ua1+ua2)   ! mole fraction of COOH in pure adipic acid
      X2ref=ua2/(ua1+ua2)   ! mole fraction of CH2 in pure adipic acid
      X3ref=uw1/uw1 

      X1=(xa*ua1)/(xa*(ua1+ua2)+xw*uw1)
      X2=(xa*ua2)/(xa*(ua1+ua2)+xw*uw1)
      X3=(xw*uw1)/(xa*(ua1+ua2)+xw*uw1)

      teta1=Q1*X1/(Q1*X1+Q2*X2+Q3*X3)
      teta2=Q2*X2/(Q1*X1+Q2*X2+Q3*X3)
      teta3=Q3*X3/(Q1*X1+Q2*X2+Q3*X3)
      teta1ref=Q1*X1ref/(Q1*X1ref+Q2*X2ref)!+Q3*X3ref)  
      teta2ref=Q2*X2ref/(Q1*X1ref+Q2*X2ref)!+Q3*X3ref)
      teta3ref=Q3*X3ref/Q3*X3ref

      !!!!!!gamma!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      T1=exp(Q1*(1-log(teta1*xi11+teta2*xi21+teta3*xi31)- 
     1    (teta1*xi11/(teta1*xi11+teta2*xi21+teta3*xi31)+ 
     2     teta2*xi12/(teta1*xi12+teta2*xi22+teta3*xi32)+ 
     3     teta3*xi13/(teta1*xi13+teta2*xi23+teta3*xi33))))


      T2=exp(Q2*(1-log(teta1*xi12+teta2*xi22+teta3*xi32)-
     1  (teta1*xi21/(teta1*xi11+teta2*xi21+teta3*xi31)+
     2   teta2*xi22/(teta1*xi12+teta2*xi22+teta3*xi32)+
     3   teta3*xi23/(teta1*xi13+teta2*xi23+teta3*xi33))))


      T3=exp(Q3*(1-log(teta1*xi13+teta2*xi23+teta3*xi33)-
     1  (teta1*xi31/(teta1*xi11+teta2*xi21+teta3*xi31)+
     2   teta2*xi32/(teta1*xi12+teta2*xi22+teta3*xi32)+
     3   teta3*xi33/(teta1*xi13+teta2*xi23+teta3*xi33))))

      !!!!!!the reference case
      T1w_ref=1  !!! because there is only one main group
 
      T1a_ref=exp(Q1*(1-
     1     log(teta1ref*xi11+teta2ref*xi21)-
     2   (teta1ref*xi11/(teta1ref*xi11+teta2ref*xi21)+
     3    teta2ref*xi12/(teta1ref*xi12+teta2ref*xi22))))


      T2a_ref=exp(Q2*(1-
     1   log(teta1ref*xi12+teta2ref*xi22)-
     2  (teta1ref*xi21/(teta1ref*xi11+teta2ref*xi21)+
     3   teta2ref*xi22/(teta1ref*xi12+teta2ref*xi22))))




      gammaar=exp(ua1*(log(T1)-log(T1a_ref))+ua2*(log(T2)-log(T2a_ref)))
      gammawr=exp(uw1*(log(T3)-log(T1w_ref)))


c     Activity coefficients : sum of combinatorial -c- and rezidual -r- part
      gammaa=exp(log(gammaac)+log(gammaar))
      gammaw=exp(log(gammawc)+log(gammawr))

      aa=gammaa*xa
      aw=gammaw*xw
      ge=8.314*T*(xa*log(gammaa)+xw*log(gammaw))
      end do
c       write(16,*) xw, gammaa*xa, gammaw*xw
    

      end program activ
