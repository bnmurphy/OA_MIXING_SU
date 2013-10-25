      program activ

c     This this code gives the activity coefficients for succinic acid/water
c     mixture using UNIFAC Dortmund

c     Anca Gaman 
c     August 2005
c	Modified for succinic acid
c	Ilona Riipinen
c	September 2005


c     a for acid
c     w for water

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!SUCCINIC ACID!!!!!!!!!!!!!!!
      !!!!!!!COOH-CH2-CH2-COOH!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!WATER!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!H2O!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
        open(16,file='activ_succinic.res',status='unknown') 
c      xn=1
C     1 stands for CH2
c     2 stands for COOH
c     3 stands for H2O
      !!!!!!!!MOLE FRACTIONS!!!!!!!!!!!!!!!!!!!
      xlow=0.
      xup=1.0
      step=0.0001
      xa=xlow
      do
         if (xa.gt.xup) EXIT
         write(*,*) xa,aa,aw
         write(16,*) xa, gammaa*xa, gammaw*xw
         xa=xa+step
      xw=1.-xa    
      T=293.15
      ua1=2
      ua2=2
      uw=1
      
      
      R1=0.6325
      R2=0.8
      R3=1.7334
      Q2=0.9215
      Q1=0.7081
      Q3=2.4561

      a12=1182.2
      b12=-3.2647
      c12=0.009198
      
      a21=2017.7
      b21=-9.0933
      c21=0.01024

      a13=1391.3
      b13=-3.6156
      c13=0.001144

      a31=-17.253
      b31=0.8389
      c31=0.0009021

      a23=624.97
      b23=-4.6878
      c23=0.005237

      a32=-1795.2
      b32=12.708
      c32=-0.01546

      xi11=1
      xi22=1
      xi33=1

      xi12=exp(-(a12+b12*T+c12*T**2)/T)
      xi21=exp(-(a21+b21*T+c21*T**2)/T)

      xi13=exp(-(a13+b13*T+c13*T**2)/T)
      xi31=exp(-(a31+b31*T+c31*T**2)/T)

      xi23=exp(-(a23+b23*T+c23*T**2)/T)
      xi32=exp(-(a32+b32*T+c32*T**2)/T)





    !!!!!!!!!! COMBINATORIAL PART!!!!!!!!!!!!
 
      qa=ua1*Q1+ua2*Q2
      qw=uw*Q3


      fa=qa/(xa*qa+xw*qw)
      fw=qw/(xa*qa+xw*qw)

      ra=ua1*R1+ua2*R2
      rw=uw*R3

      va=ra/(xa*ra+xw*rw)
      vw=rw/(xa*ra+xw*rw)
      !!!parameter w can be calculated using relative van der Waals
      !!! volumes Rk of different groups
      wa=ra**(3/4)/(xa*ra**(3/4)+xw*rw**(3/4))
      ww=rw**(3/4)/(xa*ra**(3/4)+xw*rw**(3/4))

      gammaac=exp(1-wa+log(wa)-5*qa*(1-va/fa+log(va/fa)))
      gammawc=exp(1-ww+log(ww)-5*qw*(1-vw/fw+log(vw/fw)))

       !!!!!!!!!! REZIDUAL PART !!!!!!!!!!!!!!!!!!


      !!!!! the reference state!!!!
      X1ref=ua1/(ua1+ua2)
      X2ref=ua2/(ua1+ua2)
      X3ref=1

      teta1ref=Q1*X1ref/(Q1*X1ref+Q2*X2ref)
      teta2ref=Q2*X2ref/(Q1*X1ref+Q2*X2ref)
      teta3ref=1

      T1ref=exp(Q1*(1-
     1     log(teta1ref*xi11+teta2ref*xi21)-
     2     (teta1ref*xi11/(teta1ref*xi11+teta2ref*xi21)+
     3      teta2ref*xi12/(teta1ref*xi12+teta2ref*xi22))))

      T2ref=exp(Q2*(1-
     1     log(teta1ref*xi12+teta2ref*xi22)-
     2     (teta1ref*xi21/(teta1ref*xi11+teta2ref*xi21)+
     3      teta2ref*xi22/(teta1ref*xi12+teta2ref*xi22))))

      T3ref=1

      !!!!!Same thing for the solution!!!!

      X1=xa*ua1/(xa*(ua1+ua2)+xw*uw)
      X2=xa*ua2/(xa*(ua1+ua2)+xw*uw)
      X3=xw*uw/(xa*(ua1+ua2)+xw*uw)
      
      teta1=Q1*X1/(Q1*X1+Q2*X2+Q3*X3)
      teta2=Q2*X2/(Q1*X1+Q2*X2+Q3*X3)
      teta3=Q3*X3/(Q1*X1+Q2*X2+Q3*X3)

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

      gammaar=exp(ua1*(log(T1)-log(T1ref))+ua2*(log(T2)-log(T2ref)))
      gammawr=exp(uw*(log(T3)-log(T3ref)))




c     Activity coefficients : sum of combinatorial -c- and rezidual -r- part
      gammaa=exp(log(gammaac)+log(gammaar))
      gammaw=exp(log(gammawc)+log(gammawr))
      aa=gammaa*xa
      aw=gammaw*xw     
      end do
c       write(16,*) xp, gamman*xn, gammap*xp
    

      end program activ
