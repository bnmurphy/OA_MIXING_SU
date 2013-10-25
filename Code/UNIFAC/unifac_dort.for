	PROGRAM UNIFAC_DORT

	! Activity coefficient predictions with the UNIFAC Dortmund model

	IMPLICIT NONE

	! The presicion parameters just for the main program version
	INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
	INTEGER, PARAMETER :: dk = SELECTED_INT_KIND(9)
	INTEGER(KIND = dk) :: i, j, l ! The running indices 

	! The total number of functional groups, given as input
	INTEGER(KIND = dk) :: k
	! The number of the groups k in molecules 1 and 2
	REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: nk

	! The parameters
	REAL(KIND = dp), PARAMETER :: RGAS = 8.314510E3

	! The rest of input data
	! The group-specific parameters and mole fractions
	REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: Qk, Rk, Xg, Theta
	REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: help71, help72
	REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: loggammak  
	! The matrix containing the interaction parameters
	REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: ak, psi
	! The mole fraction of compound 2
	REAL(KIND = dp) :: Xmole, summa
	! The temperature
	REAL(KIND = dp) :: Temp
	! The molecule-specific parameters and mole fractions
	REAL(KIND = dp), DIMENSION(2) :: q, r, V, Vdot, F, Xx, help9, Xapu
	! The helping variables needed to calculate derivatives for
	! mixing enthalpies
	REAL(KIND = dp) :: deltaT, T
	REAL(KIND = dp), DIMENSION(2,3) :: gamma_diff
	REAL(KIND = dp), DIMENSION(2) :: deltagamma, M
	! The output data
	! The activity coefficients
	REAL(KIND = dp), DIMENSION(2) :: gamma, loggamma
	REAL(KIND = dp), DIMENSION(2) :: loggammac, loggammar
	! The mixing enthalpies
	REAL(KIND = dp), DIMENSION(2) :: mix_ent

	! Parameters for this program version only
	REAL(KIND = dp) :: X_low, X_up, step


	k = 3 ! For the propanol-nonane case
		
	! Allocating for the sizes of the matrices
	ALLOCATE(Qk(k),Rk(k), Xg(k), Theta(k))
	ALLOCATE(help71(k), help72(k))
	ALLOCATE(loggammak(k,3))
	ALLOCATE(ak(k,k), psi(k,k))
	ALLOCATE(nk(k,2))
	! Molar masses of propanol and nonane
	M(1) = 60.11
	M(2) = 128.2

	OPEN(16,FILE = 'activ_nonprop_ud.res', STATUS = 'unknown') 

      X_low=0.0
      step=0.0001
	X_up=1.0 - step
	Xmole = 0.0 + step

	

      DO

	T = 277.11
	deltaT = 0.01 
	Temp = T - deltaT
       
	Xx(2) = Xmole
	Xx(1) = 1.0 - Xmole
	

	DO l = 1,3

	nk(1,1) = 1 ! CH3 groups for propanol
	nk(1,2) = 2 ! CH3 groups for nonane
	nk(2,1) = 2 ! CH2 groups for propanol
	nk(2,2) = 7 ! CH2 groups for nonane
	nk(3,1) = 1 ! OH group for propanol

	Rk(1) = 0.6325 
	Rk(2) = 0.6325 
	Rk(3) = 1.2302 !From Gmehling et al.,Ind. Eng. Chem. Res., 32, 178, 1993

	Qk(1)=1.0608
      Qk(2)=0.7081
      Qk(3)=0.8927 !From Gmehling et al.,Ind. Eng. Chem. Res., 32, 178, 1993

	ak(1,1) = 1.0 !CH3 & CH3
	ak(1,2) = 1.0 !CH3 & CH2
	ak(1,3) = EXP(-(2777.0 - 4.674*Temp + 0.1551E-2*Temp**2.0)/Temp) !CH3 & OH
	ak(2,1) = 1.0 !CH2 & CH3
	ak(2,2) = 1.0 !CH2 & CH2
	ak(2,3) = ak(1,3) ! CH2 & OH, as CH3 and CH2 are the same main group
	ak(3,1) = EXP(-(1606.0 - 4.746*Temp + 0.9181E-3*Temp**2.0)/Temp) !OH & CH3
	ak(3,2) = ak(3,1) ! OH & CH2, as CH3 and CH2 are the same main group	
	ak(3,3) = 1.0 !OH & OH

	! Here the main part that should be implemented into the subroutine starts
	! Eqs. in Gmehling et al., 1993:

	psi = ak

	DO i = 1,2
	! Eq. 5a:
	   q(i) = SUM(nk(:,i)*Qk)
	! Eq. 4a:
	   r(i) = SUM(nk(:,i)*Rk)
	END DO

	DO i = 1,2
	! Eq. 4:
		V(i) = r(i)/SUM(Xx*r)
	! Eq. 5:
		F(i) = q(i)/SUM(Xx*q)   	
	! Eq. 3:
		Vdot(i) = r(i)**(3./4.)/SUM(Xx*(r**(3./4.)))
	! Eq. 2: The combinatorial term
		loggammac(i) = 1.0 - Vdot(i) + LOG(Vdot(i)) - 
	1	5.*q(i)*(1.0 - V(i)/F(i) + LOG(V(i)/F(i))) 
	! Helps for Eq. 9:
		help9(i) = SUM(nk(:,i)*Xx(i))
	END DO

	DO j = 1,3 ! The pure reference cases are calculated simultaneously
	
		IF (j == 1) THEN
			Xapu(1) = 1.0
			Xapu(2) = 0.0 ! Pure propanol ref. case
		ELSE IF (j == 2) THEN
			Xapu(1) = 0.0
			Xapu(2) = 1.0 ! Pure nonane ref. case 
		ELSE IF (j == 3) THEN
			Xapu = Xx
		END IF

		DO i = 1,k
		! Eq. 9:
			Xg(i) = SUM(nk(i,:)*Xapu)/SUM(help9)
		END DO

		DO i = 1,k
		! Eq. 8:
			Theta(i) = Qk(i)*Xg(i)/SUM(Qk*Xg)
		END DO
	

		DO i = 1,k
		! Helps for Eq. 7:
			help71(i) = SUM(Theta*psi(:,i))
		END DO

		DO i = 1,k
		!Helps for Eq. 7:
			help72(i) = LOG(help71(i)) + SUM((Theta*psi(i,:))/help71)
		! Eq. 7:
			loggammak(i,j) = Qk(i)*(1.0 - help72(i))	 
		END DO

	END DO

	! Eq. The calculation of the residual part:
	DO i = 2,3
		loggammar(i-1) = SUM(nk(:,i-1)*
	1	(loggammak(:,3) - loggammak(:,i-1)))
		loggamma(i-1) = loggammac(i-1) + loggammar(i-1)
	END DO

	IF (l<3) THEN
	gamma_diff(:,l) = loggamma
		DO i = 1,2
		deltagamma(i) = (gamma_diff(i,2) - gamma_diff(i,1))/deltaT
		END DO
		Temp = T + deltaT
	END IF
	IF (l==2) THEN
	T = Temp
	END IF
		gamma = EXP(loggamma)
	END DO

	mix_ent = -RGAS*(Temp**2)*deltagamma/M
	
	IF (Xmole > X_up) EXIT
          WRITE(*,*) Xmole, gamma*Xx
          WRITE(16,*) Xmole, gamma*Xx
          Xmole = Xmole + step


      
	END DO



	END PROGRAM UNIFAC_DORT