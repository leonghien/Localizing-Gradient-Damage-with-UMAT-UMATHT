C	************************************************************************																				 C
C	*	FULLY COUPLED THERMO-MECHANICAL LOCALIZING GRADIENT DAMAGE MODEL   *
C	************************************************************************
C	von-Mises equivalent strain
C	For 2D and 3D analyses
C	POH L.H. (ceeplh@nus.edu.sg)
C 	© 2021 National University of Singapore
C	************************************************************************		
		SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

		INCLUDE 'ABA_PARAM.INC'

		CHARACTER*80 CMNAME
		DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
	 
		DOUBLE PRECISION CC(NTENS,NTENS),delta(NTENS),
     $ II1,JJ2,II1_d(NTENS),JJ2_d(NTENS),strstr(3),tr_strstr,
     $ H44(NTENS,NTENS),HSTRAN(NTENS),H66(NTENS,NTENS),AA,
     $ e_tilde1,kappa,omega,ee,gg,coeff_omega,coeff_gg,coeff_ee(NTENS)
			

C   VARIABLES FOR PROPS
		DOUBLE PRECISION E,NU,k,kappa_o
        
		INTEGER i, j, a, b       
		
		PARAMETER(Zero=0d0, One=1.0d0, Two=2.0d0, Three=3.0d0, 
     $ Six=6.0d0, TOL=1d-25, Alpha=1.0d0, Beta=9.0d0, RR=5.0d-3, 
     $ Eta=5.0d0) 
		
C	PROPS
        E=PROPS(1)          ! Young's modulus
        NU=PROPS(2)         ! Poisson's ratio
		k=PROPS(3)          ! Ratio f_t/f_c 
		kappa_o=props(4)    ! Damage initiation threshold
		
C	INITIATION
		CC=0d0
		delta=0d0
		DDSDDE=0d0
		DDSDDT=0d0
		DRPLDE=0d0
		DRPLDT=0d0

C	TOTAL STRAIN 
		DO a=1,NTENS		
			STRAN(a)=STRAN(a)+DSTRAN(a)	
		END DO
					
	! FOR PLANE STRESS (CONVENTIONS:11,22,12)
		IF (NDI==2 .and. NSHR==1) THEN 
		! Elastic stiffness CC      
			DO i=1,2
            	DO j=1,2
               		CC(i,j)=NU
            	END DO
        	END DO

        	DO i=1,2
				CC(i,i)=One
        	END DO 

			CC(3,3)=(One-NU)/Two 
			CC=CC*E/((One-NU)*(One+NU))
		
		! Kronecker delta 
			delta(1)=One
        	    delta(2)=One
        	    delta(3)=Zero
		
		! Invariants II1&JJ2
			II1=(STRAN(1)+STRAN(2))*((One-Two*NU)/(One-NU))

			II1_d=((One-Two*NU)/(One-NU))*delta

			JJ2=Two*(STRAN(1)**Two+STRAN(2)**Two-STRAN(1)*STRAN(2))
     $ +Three*((STRAN(3)/Two)**Two+(STRAN(3)/Two)**Two)
     $ +(Two*NU/((One-NU)**Two))*((STRAN(1)+STRAN(2))**Two)
		 	
		!Deviatoriate part of JJ2
        	JJ2_d(1)=Two*((Two+Two*NU/((One-NU)**Two))*STRAN(1) 
     $ +(-One+Two*NU/((One-NU)**Two))*STRAN(2)) 	
              
        	JJ2_d(2)=Two*((Two+Two*NU/((One-NU)**Two))*STRAN(2) 
     $ +(-One+Two*NU/((One-NU)**Two))*STRAN(1))
              
        	JJ2_d(3)=Six*(STRAN(3)/Two)
			
		ELSE
	! FOR PLANE STRAIN (CONVENTIONS:11,22,33,12)
		IF (NDI==3 .and. NSHR==1) THEN
		! Elastic stiffness CC       
			DO i=1,NDI
            	DO j=1,NDI
               		CC(i,j)=NU
            	END DO
        	END DO

        	DO i=1,NDI
				CC(i,i)=1.d0-NU
        	END DO

			CC(4,4)=(One-Two*NU)/Two                   		
			CC=CC*E/((One+NU)*(One-Two*NU))

		! Kronecker delta
			delta(1)=One
        	delta(2)=One
        	delta(3)=One
        	delta(4)=Zero

		! Invariants II1&JJ2
			II1=0d0
        	DO a=1,NTENS
				II1=II1+STRAN(a)*delta(a)
			END DO
                    
          ! strstr: ???
        	strstr(1)=STRAN(1)*STRAN(1)+(STRAN(4)/Two)*(STRAN(4)/Two)	 	
        	strstr(2)=STRAN(2)*STRAN(2)+(STRAN(4)/Two)*(STRAN(4)/Two)			
        	strstr(3)=STRAN(3)*STRAN(3) 		

			tr_strstr=0d0
    		DO a=1,3
        		tr_strstr=tr_strstr+strstr(a)
    		END DO

			JJ2=Three*tr_strstr-II1*II1

			H44=0d0
			DO a=1,NDI
            	DO b=1,NDI
                	H44(a,b)=-One
            	END DO
        	END DO

        	DO a=1,NDI
		    	H44(a,a)=Two
        	END DO

			H44(4,4)=Three/Two

			Hstran=0d0
        	DO a=1,NTENS
				DO b=1,NTENS	 
					Hstran(a)=Hstran(a)+H44(a,b)*STRAN(b)
            	END DO
        	END DO
		
		ELSE
	! FOR 3D (CONVENTIONS:11,22,33,12,13,23)
		IF (NDI==3 .and. NSHR==3) THEN
		! Elastic stiffness CC      
			CC(1,1)=1.d0-NU
			CC(1,2)=NU
			CC(1,3)=NU
			CC(2,1)=NU
			CC(2,2)=1.d0-NU
			CC(2,3)=NU
			CC(3,1)=NU
			CC(3,2)=NU
			CC(3,3)=1.d0-NU
			CC(4,4)=(1.d0-2.d0*NU)/Two
			CC(5,5)=(1.d0-2.d0*NU)/Two
			CC(6,6)=(1.d0-2.d0*NU)/Two
			CC=CC*E/((1.d0+NU)*(1.d0-2.d0*NU))	

		! Kronecker delta
			delta(1)=One
			delta(2)=One
			delta(3)=One
			delta(4)=Zero		
			delta(5)=Zero
			delta(6)=Zero		
			
		! Invariants II1&JJ2		
			II1=0d0
			DO a=1,NTENS
				II1 = II1+STRAN(a)*delta(a)
			END DO

			strstr(1)=stran(1)*stran(1)+(stran(4)/Two)*(stran(4)/Two)+(stran(5)/Two)*(stran(5)/Two)		
			strstr(2)=(stran(4)/Two)*(stran(4)/Two)+stran(2)*stran(2)+(stran(6)/Two)*(stran(6)/Two)				
			strstr(3)=(stran(5)/Two)*(stran(5)/Two)+(stran(6)/Two)*(stran(6)/Two)+stran(3)*stran(3)		

			tr_strstr=0d0
			DO a=1,3
				tr_strstr = tr_strstr + strstr(a)
			END DO

			JJ2=Three*tr_strstr-II1*II1	

C   H66(6,6) = 3*delta_ik*delta_jl - delta_kl*delta_ij   
			H66=0d0 
			DO a=1,3
				DO b=1,3
					H66(a,b)=-One
				END DO
			END DO 			
		
			DO a=1,3
				H66(a,a)=Two
			END DO
      
			DO a=4,6
				H66(a,a)=Three/Two
			END DO		

			Hstran=0d0
			DO a=1,NTENS
				DO b=1,NTENS	 
					Hstran(a) = Hstran(a)+H66(a,b)*stran(b)
				END DO
			END DO 			
				
		ENDIF		
		ENDIF
	    ENDIF

C	EQUIV STRAIN
        AA=((((k-One)*II1)/(One-Two*NU))**Two)+((Two*k*JJ2)/((One+NU)**Two))
			 
	!Make sure no zero in the denominator   
        IF(AA .le. Zero) AA=TOL 

		ee=((k-One)*II1)/(Two*k*(One-Two*NU))+(sqrt(AA))/(Two*k)

      !Update e_tilde(TEMP)	
		TEMP=TEMP+DTEMP

	!Make sure no zero in the denominator
		IF (TEMP .eq. Zero) TEMP=TOL
			        
	!kappa (histrory parameter)
		e_tilde1=STATEV(1)

        IF (e_tilde1 .lt. TEMP) THEN
			kappa=TEMP        
        ELSE
			kappa=e_tilde1
	    ENDIF
									
	!omega (damage variable)
        IF  (kappa .lt. kappa_o) THEN
            omega=Zero							
        ELSE 
            omega=One-(kappa_o/kappa)*((One-Alpha)+Alpha*(exp(-Beta*(kappa-kappa_o))))				
        ENDIF

	!gg (interction parameter)
        gg=((One-RR)*(exp(-Eta*omega))+RR-(exp(-Eta)))/(One-exp(-Eta))	

C	TOTAL STRESS				
		STRESS=0d0
		DO i=1, NTENS
			DO j=1, NTENS
				STRESS(i)=STRESS(i)+(One-omega)*CC(i,j)*STRAN(j)
			END DO
		END DO

C	************************************************************************
C   COEFFICIENTS 							
C	DDSDDE
		DDSDDE=(One-omega)*CC

C	DDSDDT				
	! coeff_omega - domega/dkappa
        IF (kappa .lt. kappa_o) THEN
                coeff_omega=Zero			
        ELSE							
		IF (e_tilde1 .lt. TEMP) THEN
                coeff_omega=(kappa_o/(kappa**Two))*(One-Alpha
     $ 	+ Alpha*(exp(-Beta*(kappa-kappa_o))))+(kappa_o/kappa)*(Alpha*
     $ (exp(-Beta*(kappa-kappa_o)))*Beta)					
        ELSE
                coeff_omega=Zero
        ENDIF
        ENDIF

	! coeff_gg - dgg/domega
		coeff_gg=((One-RR)*(exp(-Eta*omega))*(-Eta))/(One-exp(-Eta))			
		
		DO i=1, NTENS
			DO j=1, NTENS
				DDSDDT(j)=-CC(j,i)*STRAN(i)*coeff_omega
			END DO
		END DO	

C	RPL
		RPL=ee-TEMP

C	DRPLDE
	! coeff_ee - dee/dstran
	! FOR PLANE STRESS
		IF (NDI==2 .and. NSHR==1) THEN
			DO a=1,NTENS
        		coeff_ee(a)=((k-One)/(Two*k*(One-Two*NU)))*
     $ (One+((k-One)/(One-Two*NU))*(One/sqrt(AA))*II1)*II1_d(a)			
     $ +((One/sqrt(AA))/(Two*((One+NU)**Two)))*JJ2_d(a)		
        	END DO			
		ELSE		
	! FOR PLANE STRESS & 3D 
			DO a=1,NTENS
				coeff_ee(a)=((k-One)/(Two*k*(One-Two*NU)))*(One+	
     $ ((k-One)/(One-Two*NU))*(One/sqrt(AA))*II1)*delta(a) 
     $ 	+(One/sqrt(AA))*Hstran(a)/((One+NU)**Two)
        	END DO
        ENDIF		
				
		DRPLDE=coeff_ee
		
C	DRPLDT
		DRPLDT=-One

C	STORE
		STATEV(1)=kappa
		STATEV(2)=omega
		STATEV(3)=gg
		STATEV(4)=coeff_gg
		STATEV(5)=coeff_omega
			
      
		RETURN
    	END
		
C	*********************************************************************	
		SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

		INCLUDE 'ABA_PARAM.INC'

		CHARACTER*80 CMNAME
		DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)

		COND=PROPS(1)
		SPECHT=PROPS(2)

		DUDT=SPECHT
		DUDG=0d0
		DU=DUDT*DTEMP
		U=U+DU

		DO I=1, NTGRD
		FLUX(I)=-STATEV(3)*(COND*COND)*DTEMDX(I)
		DFDG(I,I)=-STATEV(3)*(COND*COND)
		DFDT(I)=-STATEV(4)*STATEV(5)*(COND*COND)*DTEMDX(I)
		END DO					

     	RETURN
     	END