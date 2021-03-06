BEMM.4PL=function(data, 				          #A matrix of response [n.examinees * n.items]
                  PriorA=c(0,0.25), 	    #The log normal prior for a parameters with default mean 0 and variance 0.25
                  PriorB=c(0,4),          #The normal prior for b parameters with default mean 0 and variance 4
                  PriorC=c(4,16),         #The beta prior for c parameters with default hyper-parameter 4 and 16
                  PriorS=c(4,16),         #The beta prior for s parameters with default hyper-parameter 4 and 16
                  InitialA=NA,			      #Initial values for a parameters, default is NA 
                  InitialB=NA,			      #Initial values for b parameters, default is NA 
                  InitialC=NA,			      #Initial values for c parameters, default is NA 
                  InitialS=NA,			      #Initial values for s parameters, default is NA 
                  Tol=0.0001,			        #The tolerate threshold for convergence, default is 0.0001
                  max.ECycle=2000L,		    #The max of E-step iteration, default is 2000L
                  max.MCycle=100L,		    #The max of M-step iteration, default is 100L
                  n.decimal=3L,           #The decimal length of outputs parameters, default is 3L
                  n.Quadpts =31L,		      #The number of quadrature, default is 31L
                  Theta.lim=c(-6,6),      #The range the Theta, default is [-6,6]
                  Missing=-9,             #A number to indicate missing value, default is -9
                  ParConstraint=FALSE,    #A logical value to determine whether imposing a range limitation for parameters, default is FALSE
                  BiasSE=FALSE){          #A logical value to determine whether directly estimating SEs from inversed Hession matrix, default is FALSE
  
  Time.Begin=Sys.time()  #Recording starting time
  Model='4PL'            #Set the model to be 4PLM
  D=1.702                #Set the constant D to 1.702
  
  ###Check Input variables and return processed results###        
  Check.results=Input.Checking(Model=Model, data=data, PriorA=PriorA, PriorB=PriorB, PriorC=PriorC, PriorS=PriorS,
                               InitialA=InitialA, InitialB=InitialB, InitialC=InitialC, InitialS=InitialS,
                               Tol=Tol, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                               Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint, BiasSE=BiasSE)
  
  data=Check.results$data                     #Return the processed data
  data.simple=Check.results$data.simple       #Return the simplified data
  CountNum=Check.results$CountNum             #Return the counts of simplified data
  I=Check.results$I                           #Return the nrow of data
  J=Check.results$J                           #Return the ncol of data
  n.class=Check.results$n.class               #the ncol of data.simple
  
  PriorA=Check.results$PriorA                 #Return the prior setting for A parameter
  PriorB=Check.results$PriorB                 #Return the prior setting for B parameter
  PriorC=Check.results$PriorC                 #Return the prior setting for C parameter
  PriorS=Check.results$PriorS                 #Return the prior setting for S parameter
  Prior=list(PriorA=PriorA, PriorB=PriorB, PriorC=PriorC, PriorS=PriorS)   #create a list to save Priors
  
  InitialA=Check.results$InitialA             #Return the starting values for A parameter
  InitialB=Check.results$InitialB             #Return the starting values for B parameter
  InitialC=Check.results$InitialC             #Return the starting values for C parameter
  InitialS=Check.results$InitialS             #Return the starting values for S parameter
  
  max.ECycle=Check.results$max.ECycle         #Return the max of Expectation cycles
  max.MCycle=Check.results$max.MCycle         #Return the max of Maximization cycles
  n.Quadpts=Check.results$n.Quadpts           #Return the number of Quadpts
  n.decimal=Check.results$n.decimal           #Return the number of decimal
  
  ParConstraint=Check.results$ParConstraint   #Return the value of ParConstraint
  BiasSE=Check.results$BiasSE                 #Return the value of BiasSE
  
  Par.est0=list(A=InitialA, B=InitialB, C=InitialC, S=InitialS) 			#Assemble parameters estimates to a list
  Par.SE0=list(SEA=InitialA*0, SEB=InitialB*0, SEC=InitialC*0, SES=InitialS*0)  #Assemble parameters SEs to a list
  np=J*4		#Obtain the number of estimated parameters
  
  ###Run the BEMM algorithms###
  Est.results=BEMM.4PL.est(Model=Model, data=data, data.simple=data.simple, CountNum=CountNum, n.class=n.class,
                           Prior=Prior, Par.est0=Par.est0, Par.SE0=Par.SE0, D=D, np, 
                           Tol=Tol, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                           Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint, BiasSE=BiasSE, I=I, J=J, Time.Begin=Time.Begin)
  
  #Print important information into screen after estimation						
  if (Est.results$StopNormal==1){
    message('PROCEDURE TERMINATED NORMALLY')
  }else{
    message('PROCEDURE TERMINATED WITH ISSUES')
  }
  message('IRTEMM version: 1.0.7') 
  message('Item Parameter Calibration for the 4PLM.','\n')
  message('Quadrature: ', n.Quadpts, ' nodes from ', Theta.lim[1], ' to ', Theta.lim[2], ' were used to approximate Gaussian distribution.') 
  message('Method for Items: Ability-based Bayesian Expectation-Maximization-Maximization (BEMM) Algorithm.')
  if (BiasSE){
    message('Method for Item SEs: directly estimating SEs from inversed Hession matrix.')
    warning('Warning: The SEs maybe not trustworthy!', sep = '')
  }else{
    message('Method for Item SEs: Updated SEM algorithm.')
  }
  message('Method for Theta: Expected A Posteriori (EAP).')
  if (Est.results$StopNormal==1){
    message('Converged at LL-Change < ', round(Est.results$cr, 6), ' after ', Est.results$Iteration, ' EMM iterations.', sep = '')
  }else{
    warning('Warning: Estimation cannot converged under current max.ECycle and Tol!', sep = '')
    warning('Warning: The reults maybe not trustworthy!', sep = '')
    message('Terminated at LL-Change = ', round(Est.results$cr, 6), ' after ', Est.results$Iteration, ' EMM iterations.', sep = '')
  }
  message('Running time:', Est.results$Elapsed.time, '\n')
  message('Log-likelihood (LL):', as.character(round(Est.results$Loglikelihood, n.decimal)))
  message('Estimated Parameters:', as.character(np))
  message('AIB: ', round(Est.results$fits.test$AIC, n.decimal), ', BIC: ', round(Est.results$fits.test$BIC, n.decimal), ', RMSEA = ', round(Est.results$fits.test$RMSEA, n.decimal))
  message('G2 (', round(Est.results$fits.test$G2.df, n.decimal), ') = ', round(Est.results$fits.test$G2, n.decimal), ', p = ', round(Est.results$fits.test$G2.P, n.decimal), ', G2/df = ', round(Est.results$fits.test$G2.ratio, n.decimal), sep='')		
  
  return(Est.results)  #Return results
  
}




BEMM.4PL.est=function(Model=Model, data=data, data.simple=data.simple, CountNum=CountNum, n.class=n.class,
                      Prior=Prior, Par.est0=Par.est0, Par.SE0=Par.SE0, D=D, np, 
                      Tol=Tol, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                      Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint, BiasSE=BiasSE, I=I, J=J, Time.Begin=Time.Begin){
  
  #Generating the nodes and their weights for approximating the integration
  node.Quadpts=seq(Theta.lim[1],Theta.lim[2],length.out = n.Quadpts)	#Generating the nodes
  weight.Quadpts=dnorm(node.Quadpts,0,1)                                #Generating the weight of nodes
  weight.Quadpts=weight.Quadpts/sum(weight.Quadpts)
  
  #Initializing the program settings
  InitialValues=Par.est0    #InitialValues
  LH=rep(0,max.ECycle)		#Log-likelihood
  IA=rep(0,J)  #EM iteration history for inverse second derivative of a parameters
  IB=rep(0,J)  #EM iteration history for inverse second derivative of b parameters
  IC=rep(0,J)  #EM iteration history for inverse second derivative of c parameters
  IS=rep(0,J)  #EM iteration history for inverse second derivative of s parameters
  IAB=rep(0,J) #EM iteration history for inverse second derivative of ab parameters
  TA=matrix(0,max.ECycle,J)  #EM iteration history for a parameters
  TB=matrix(0,max.ECycle,J)  #EM iteration history for b parameters
  TC=matrix(0,max.ECycle,J)  #EM iteration history for c parameters
  TS=matrix(0,max.ECycle,J)  #EM iteration history for s parameters
  deltahat.A=rep(0,J)     #SEM iteration history for a parameters
  deltahat.B=rep(0,J)     #SEM iteration history for b parameters
  deltahat.C=rep(0,J)     #SEM iteration history for c parameters
  deltahat.S=rep(0,J)     #SEM iteration history for s parameters
  n.ECycle=1L				#The first E-step iteration
  StopNormal=0L			#Whether program terminate normally
  E.exit=0L         #Whether estimation should be ended 
  
  #Calculate Log-likelihood & Update artificial data
  LLinfo=LikelihoodInfo(data.simple, CountNum, Model, Par.est0, n.Quadpts, node.Quadpts, weight.Quadpts, D)
  LH0=LLinfo$LH
  f=LLinfo$f
  r=LLinfo$r
  fz=LLinfo$fz
  rz=LLinfo$rz
  
  while(E.exit==0L && (n.ECycle <= max.ECycle)){#E step iteration
    for (j in 1:J){
      #estimate c parameters parameters
      c.part1 = sum(r[j,] - rz[j,])
      c.part2 = sum(f - fz[j,] - r[j,] + rz[j,])
      c.part3 = sum(f - fz[j,])
      c.part4 = 1.0 - Par.est0$C[j]
      c.part5 = c.part4 * c.part4
      c.part6 = Par.est0$C[j] * Par.est0$C[j]
      if (Prior$PriorC[j]==-9 || Prior$PriorC[j + J]==-9) {
        lcc = - (c.part1 / c.part6) - (c.part2) / c.part5	#lcc=-(r-rz)/c^2-(f-fz-r+rz)/(1-c)^2
        ct0 = c.part1 / c.part3							#c=(r-rz)/(f-fz)
      } else {
        lcc = - ((c.part1 + Prior$PriorC[j] -1) / c.part6) - ((c.part2 + Prior$PriorC[j + J] - 1.0) / c.part5)
        ct0 = (Prior$PriorC[j] -1 + c.part1) / (((Prior$PriorC[j] + Prior$PriorC[j + J]) - 2.0) + c.part3)
      }
      #estimate s parameters parameters
      s.part1 = sum(fz[j,] - rz[j,])
      s.part2 = sum(rz[j,])
      s.part3 = sum(fz[j,])
      s.part4 = 1.0 - Par.est0$S[j]
      s.part5 = c.part4 * c.part4
      s.part6 = Par.est0$S[j] * Par.est0$S[j]
      if (Prior$PriorS[j]==-9 || Prior$PriorS[j + J]==-9) {
        lss = - ( s.part1 / s.part6) - ( s.part2 / s.part5)	#lss=-(fz-rz)/s^2-rz/(1-s)^2
        st0 = s.part1 / s.part3							                #s=(fz-rz)/fz
      } else {
        lss = - (( s.part1 + Prior$PriorS[j] -1) / s.part6) - (( s.part2+ Prior$PriorS[j + J] - 1.0) / s.part5)
        st0 = (Prior$PriorS[j] -1 + s.part1) / (((Prior$PriorS[j] + Prior$PriorS[j + J]) - 2.0) + s.part3)
      }
      
      #estimate a and b parameter
      at0 = Par.est0$A[j]
      bt0 = Par.est0$B[j]
      n.MCycle = 1
      M.exit = 0L
      #M step iteration
      while (M.exit==0L && (n.MCycle <= max.MCycle)){
        Da = D * at0
        D2a2 = Da * Da
        x.bt = node.Quadpts - bt0   	#x-bt
        x.bt2 = x.bt * x.bt  				#(x-bt)^2
        pstar = 1 / (1 + exp(-Da * x.bt))	#ps
        psf = pstar * f					#f * ps
        wsf = psf * (1-pstar)    			#f * ws
        fz.psf = fz[j,] - psf  	#rz - f *ps
        la1 = Da * sum(fz.psf * x.bt)			#la1: (fz - f*ps)* (x-bt)
        lb1 = -Da * sum(fz.psf)       			#lb1: fz - f* ps
        laa = -D2a2 * sum(wsf * x.bt2)			#laa: (f*ws)* (x-bt)^2 
        lbb = -D2a2 * sum(wsf)          			#lbb: (f*ws)
        lab = D2a2 * sum(wsf * x.bt)  			#lab: wsf*(x-bt) 
        # Maximize a and b
        if (Prior$PriorA[j]!=-9 && Prior$PriorA[j + J]!=-9) {
          la1 = la1 -((log(at0)-Prior$PriorA[j])/Prior$PriorA[j + J])
          laa = laa -1/Prior$PriorA[j + J]
        }
        if (Prior$PriorB[j]!=-9 && Prior$PriorB[j + J]!=-9) {
          lb1 = lb1 -((bt0-Prior$PriorB[j])/Prior$PriorB[j + J])
          lbb = lbb -1/Prior$PriorB[j + J]
        }
        Weight = laa * lbb - lab * lab
        Iaa = - lbb / Weight
        Ibb = - laa / Weight
        Iab = lab / Weight
        Icc = - 1 / lcc
        Iss = - 1 / lss
        at1 = log(at0) + (Iaa * la1 + Iab * lb1)
        bt1 = bt0 + (Ibb * lb1 + Iab * la1)
        if (sqrt((at1-log(at0))^2 + (bt1-bt0)^2) < 0.01){
          M.exit = 1L
          at0 = exp(at1)
          bt0 = bt1
        }else{
          at0 = exp(at1)
          bt0 = bt1
          n.MCycle = n.MCycle+1
        }
      }
      if (is.finite(at0) && is.finite(bt0) && is.finite(ct0) && is.finite(st0)){
        if (ParConstraint){
          if (at0>=0.001 && at0<=6 && bt0>=-6 && bt0<=6){
            Par.est0$A[j] = at0
            Par.est0$B[j] = bt0
            TA[n.ECycle,j] = at0
            TB[n.ECycle,j] = bt0
            IA[j] = Iaa
            IB[j] = Ibb
            IAB[j] = Iab
          }
          if (ct0>=0.0001 && ct0<=0.5){
            Par.est0$C[j] = ct0
            TC[n.ECycle,j] = ct0
            IC[j] = Icc
          }
          if (st0>=0.0001 && st0<=0.5){
            Par.est0$S[j] = st0
            TS[n.ECycle,j] = st0
            IS[j] = Iss
          }
        }else{
          if (at0>=0.001){
            Par.est0$A[j] = at0
            TA[n.ECycle,j] = at0
          }else{
            TA[n.ECycle,j] = TA[n.ECycle-1,j]
          }
          Par.est0$B[j] = bt0
          Par.est0$C[j] = ct0
          Par.est0$S[j] = st0
          TB[n.ECycle,j] = bt0
          TC[n.ECycle,j] = ct0
          TS[n.ECycle,j] = st0
          IA[j] = Iaa
          IB[j] = Ibb
          IAB[j] = Iab
          IC[j] = Icc
          IS[j] = Iss
        }
      }else{
        if (n.ECycle!=1){
          TA[n.ECycle,j] = TA[n.ECycle-1,j]
          TB[n.ECycle,j] = TB[n.ECycle-1,j]
          TC[n.ECycle,j] = TC[n.ECycle-1,j]
          TS[n.ECycle,j] = TS[n.ECycle-1,j]
        }else{
          TA[n.ECycle,j] = Par.est0$A[j]
          TB[n.ECycle,j] = Par.est0$B[j]
          TC[n.ECycle,j] = Par.est0$C[j]
          TS[n.ECycle,j] = Par.est0$S[j]
        }
      }
    }
    #Calculate Log-likelihood & Update artificial data
    LLinfo=LikelihoodInfo(data.simple, CountNum, Model, Par.est0, n.Quadpts, node.Quadpts, weight.Quadpts, D)
    LH[n.ECycle]=LLinfo$LH
    f=LLinfo$f
    r=LLinfo$r
    fz=LLinfo$fz
    rz=LLinfo$rz  
    cr=LH[n.ECycle]-LH0
    LH0=LH[n.ECycle]   
    if (abs(cr)<Tol){
      n.ECycle = n.ECycle + 1
      E.exit=1
      StopNormal=1L
    }else{
      n.ECycle = n.ECycle + 1
    }
  }
  n.ECycle = n.ECycle - 1  
  if (BiasSE==FALSE){
    # Estimating Standard Errors via Supplement EM algorithm
    start.SEM=0
    end.SEM=n.ECycle
    delta=rep(0,6)
    delta0=rep(0,6)
    delta1=rep(0,6)
    cr.SEM0=1
    cr.SEM1=1
    cr.SEM2=1
    cr.SEM3=1
    cr.SEM4=1
    for (i in 1:(n.ECycle-1)){
      deltatemp=exp(-(LH[i+1]-LH[i]))
      if (deltatemp>=0.9 && deltatemp<=0.999){
        if (cr.SEM0==0){
          end.SEM=i
        }else{
          start.SEM=i
          cr.SEM0=0
        }
      }
    }   
    Time.Mid=Sys.time()		#End time
    message(paste('Estimating SEs via USEM algorithm (Requires about ',as.character(round(difftime(Time.Mid, Time.Begin, units="mins"), 2)), ' mins).', sep=''),'\n')
    for (j in 1:J){
      z=start.SEM
      SEM.exit=0
      cr.SEM1=1
      cr.SEM2=1
      cr.SEM3=1
      cr.SEM4=1
      while (SEM.exit==0 && z<=end.SEM){
        for (ParClass in 1:4){
          if (ParClass==1 && z>=2 && cr.SEM1<sqrt(Tol)){
            next
          }
          if (ParClass==2 && z>=2 && cr.SEM2<sqrt(Tol)){
            next
          }
          if (ParClass==3 && z>=2 && cr.SEM3<sqrt(Tol)){
            next
          }
          if (ParClass==4 && z>=2 && cr.SEM4<sqrt(Tol)){
            next
          }
          deltahat=Par.est0
          if (ParClass==1){deltahat$A[j]=TA[z,j]}
          if (ParClass==2){deltahat$B[j]=TB[z,j]}
          if (ParClass==3){deltahat$C[j]=TC[z,j]}
          if (ParClass==4){deltahat$S[j]=TS[z,j]}
          LLinfo=LikelihoodInfo(data.simple, CountNum, Model, deltahat, n.Quadpts, node.Quadpts, weight.Quadpts, D)
          f=LLinfo$f
          r=LLinfo$r
          fz=LLinfo$fz
          rz=LLinfo$rz  
          #estimate c parameters
          if (ParClass==3){
            c.part1 = sum(r[j,] - rz[j,])
            c.part3 = sum(f - fz[j,])
            if (Prior$PriorC[j]==-9 || Prior$PriorC[j + J]==-9) {
              ct0 = c.part1 / c.part3							#c=(r-rz)/(f-fz)
            } else {
              ct0 = (Prior$PriorC[j] -1 + c.part1) / (((Prior$PriorC[j] + Prior$PriorC[j + J]) - 2.0) + c.part3)
            }
            if (is.finite(ct0)){
              if (ParConstraint){
                if (ct0>=0.0001 && ct0<=0.5){
                  deltahat$C[j] = ct0
                }
              }else{
                deltahat$C[j] = ct0
              }
            }
          }	
          #estimate s parameters
          if (ParClass==4){
            s.part1 = sum(fz[j,] - rz[j,])
            s.part3 = sum(fz[j,])
            if (Prior$PriorS[j]==-9 || Prior$PriorS[j + J]==-9) {
              st0 = s.part1 / s.part3							                #s=(fz-rz)/fz
            } else {
              st0 = (Prior$PriorS[j] -1 + s.part1) / (((Prior$PriorS[j] + Prior$PriorS[j + J]) - 2.0) + s.part3)
            }
            if (is.finite(st0)){
              if (ParConstraint){
                if (st0>=0.0001 && st0<=0.5){
                  deltahat$S[j] = st0
                }
              }else{
                deltahat$S[j] = st0
              }
            }
          }
          
          if (ParClass==1 || ParClass==2){
            #estimate a and b parameter
            at0 = deltahat$A[j]
            bt0 = deltahat$B[j]
            n.MCycle = 1
            M.exit = 0
            #M step iteration
            while (M.exit==0 && (n.MCycle <= max.MCycle)) {
              Da = D * at0
              D2a2 = Da * Da
              x.bt = node.Quadpts - bt0   	#x-bt
              x.bt2 = x.bt * x.bt  				#(x-bt)^2
              pstar = 1 / (1 + exp(-Da * x.bt))	#ps
              psf = pstar * f					#f * ps
              wsf = psf * (1-pstar)    			#f * ws
              fz.psf = fz[j,] - psf  	#rz - f *ps
              la1 = Da * sum(fz.psf * x.bt)			#la1: (fz - f*ps)* (x-bt)
              lb1 = -Da * sum(fz.psf)       			#lb1: fz - f* ps
              laa = -D2a2 * sum(wsf * x.bt2)			#laa: (f*ws)* (x-bt)^2 
              lbb = -D2a2 * sum(wsf)          			#lbb: (f*ws)
              lab = D2a2 * sum(wsf * x.bt)  			#lab: wsf*(x-bt) 
              # Maximize a and b
              if (Prior$PriorA[j]!=-9 && Prior$PriorA[j + J]!=-9) {
                la1 = la1 -((log(at0)-Prior$PriorA[j])/Prior$PriorA[j + J])
                laa = laa -1/Prior$PriorA[j + J]
              }
              if (Prior$PriorB[j]!=-9 && Prior$PriorB[j + J]!=-9) {
                lb1 = lb1 -((bt0-Prior$PriorB[j])/Prior$PriorB[j + J])
                lbb = lbb -1/Prior$PriorB[j + J]
              }
              Weight = laa * lbb - lab * lab
              Iaa = - lbb / Weight
              Ibb = - laa / Weight
              Iab = lab / Weight
              Icc = - 1 / lcc
              Iss = - 1 / lss
              at1 = log(at0) + (Iaa * la1 + Iab * lb1)
              bt1 = bt0 + (Ibb * lb1 + Iab * la1)
              if (sqrt((at1-log(at0))^2 + (bt1-bt0)^2) < 0.01){
                M.exit = 1L
                at0 = exp(at1)
                bt0 = bt1
              }else{
                at0 = exp(at1)
                bt0 = bt1
                n.MCycle = n.MCycle+1
              }
            }
            if (is.finite(at0) && is.finite(bt0)){
              if (ParConstraint){
                if (at0>=0.001 && at0<=6 && bt0>=-6 && bt0<=6){
                  deltahat$A[j] = at0
                  deltahat$B[j] = bt0
                }
              }else{
                if (at0>=0.001){
                  deltahat$A[j] = at0
                }else{
                }
                deltahat$B[j] = bt0
              }
            }
          }
          if (ParClass==1){
            delta1[1]=(log(deltahat$A[j])-log(Par.est0$A[j]))/(log(TA[z,j])-log(Par.est0$A[j])+0.0001)
            delta1[2]=(deltahat$B[j]-Par.est0$B[j])/(log(TA[z,j])-log(Par.est0$A[j])+0.0001)
            break
          }
          if (ParClass==2){
            delta1[3]=(log(deltahat$A[j])-log(Par.est0$A[j]))/(TB[z,j]-Par.est0$B[j]+0.0001)
            delta1[4]=(deltahat$B[j]-Par.est0$B[j])/(TB[z,j]-Par.est0$B[j]+0.0001)
            break
          }
          if (ParClass==3){
            delta1[5]=(deltahat$C[j]-Par.est0$C[j])/(TC[z,j]-Par.est0$C[j]+0.0001)
            break
          }
          if (ParClass==4){
            delta1[6]=(deltahat$S[j]-Par.est0$S[j])/(TS[z,j]-Par.est0$S[j]+0.0001)
            break
          }
        }
        cr.SEM1=sqrt((delta1[1]-delta0[1])^2+(delta1[2]-delta0[2])^2)
        cr.SEM2=sqrt((delta1[3]-delta0[3])^2+(delta1[4]-delta0[4])^2)
        cr.SEM3=abs(delta1[5]-delta0[5])
        cr.SEM4=abs(delta1[6]-delta0[6])
        if (cr.SEM1<sqrt(Tol) && cr.SEM2<sqrt(Tol) && cr.SEM3<sqrt(Tol)  && cr.SEM4<sqrt(Tol) && z>=2){
          SEM.exit=1
        }else{
          z=z+1
        }
        delta0[is.finite(delta1)] = delta1[is.finite(delta1)]
      }
      delta[1]=1-delta0[1]
      delta[2]= -delta0[2]
      delta[3]= -delta0[3]
      delta[4]=1-delta0[4]
      delta[5]=1-delta0[5]
      delta[6]=1-delta0[6]
      Weight = delta[1] * delta[4] - delta[2] * delta[3]
      delta1[1] =  delta[4] / Weight
      delta1[2] = -delta[2] / Weight
      delta1[3] = -delta[3] / Weight
      delta1[4] =  delta[1] / Weight
      delta1[5] =  1 / delta[5]
      delta1[6] =  1 / delta[6]
      if (is.finite(delta1[1])==F || delta1[1]<=0){delta1[1] = 1}
      if (is.finite(delta1[2])==F || delta1[2]<=0){delta1[2] = 0}
      if (is.finite(delta1[3])==F || delta1[3]<=0){delta1[3] = 0}
      if (is.finite(delta1[4])==F || delta1[4]<=0){delta1[4] = 1}
      if (is.finite(delta1[5])==F || delta1[5]<=0){delta1[5] = 1}
      if (is.finite(delta1[6])==F || delta1[6]<=0){delta1[6] = 1}
      Par.SE0$SEA[j]= sqrt(Par.est0$A[j] * Par.est0$A[j] * IA[j] * delta1[1] + IAB[j] * delta1[3])
      Par.SE0$SEB[j]= sqrt(IB[j] * delta1[4] + IAB[j] * delta1[2])
      Par.SE0$SEC[j]= sqrt(IC[j] * delta1[5])
      Par.SE0$SES[j]= sqrt(IS[j] * delta1[6])
      if (is.finite(Par.SE0$SEA[j])){if (Par.SE0$SEA[j]>1){Par.SE0$SEA[j]= sqrt(Par.est0$A[j] * Par.est0$A[j] * IA[j])}}
      if (is.finite(Par.SE0$SEB[j])){if (Par.SE0$SEB[j]>1){Par.SE0$SEB[j]= sqrt(IB[j])}}
	  if (is.finite(Par.SE0$SEC[j])){if (Par.SE0$SEC[j]>1){Par.SE0$SEC[j]= sqrt(IC[j])}}
      if (is.finite(Par.SE0$SES[j])){if (Par.SE0$SES[j]>1){Par.SE0$SES[j]= sqrt(IS[j])}}
    }
  }else{
    message('Directly estimating SEs from inversed Hession matrix.', '\n')
    for (j in 1:J){
      Par.SE0$SEA[j]= sqrt(Par.est0$A[j] * Par.est0$A[j] * IA[j])
      Par.SE0$SEB[j]= sqrt(IB[j])
      Par.SE0$SEC[j]= sqrt(IC[j])
      Par.SE0$SES[j]= sqrt(IS[j])
    }
  }
  
  Par.est0$A=round(Par.est0$A, n.decimal)		#Save the estimated a parameters
  Par.est0$B=round(Par.est0$B, n.decimal)		#Save the estimated b parameters
  Par.est0$C=round(Par.est0$C, n.decimal)		#Save the estimated c parameters
  Par.est0$S=round(Par.est0$S, n.decimal)		#Save the estimated s parameters
  
  Par.SE0$SEA=round(Par.SE0$SEA, n.decimal)		#Save the estimated SEs of a parameters 
  Par.SE0$SEB=round(Par.SE0$SEB, n.decimal)		#Save the estimated SEs of b parameters
  Par.SE0$SEC=round(Par.SE0$SEC, n.decimal)		#Save the estimated SEs of c parameters
  Par.SE0$SES=round(Par.SE0$SES, n.decimal)		#Save the estimated SEs of s parameters
  
  EM.Map=list(Map.A=TA[1:n.ECycle,],Map.B=TB[1:n.ECycle,], Map.C=TC[1:n.ECycle,], Map.S=TS[1:n.ECycle,]) #Transfer the variable EM.map to a list
  
  Est.ItemPars=as.data.frame(list(est.a=Par.est0$A, est.b=Par.est0$B, est.c=Par.est0$C, est.s=Par.est0$S, se.a=Par.SE0$SEA, se.b=Par.SE0$SEB, se.c=Par.SE0$SEC, se.s=Par.SE0$SES))
  #Transfer the variable estimated item parameters to a dataframe
  
  #Calculate examinees ability via EAP methods.
  P.Quadpts=lapply(as.list(node.Quadpts), Prob.model, Model=Model, Par.est0=Par.est0, D=D)	
  Joint.prob=mapply('*',lapply(P.Quadpts, function(P,data){apply(data*P+(1-data)*(1-P),2,prod,na.rm = T)}, data=t(data)),
                    as.list(weight.Quadpts), SIMPLIFY = FALSE)		
  Whole.prob=Reduce("+", Joint.prob)		
  LogL=sum(log(Whole.prob))      #Obtain the final Log-likelihood
  Posterior.prob=lapply(Joint.prob, '/', Whole.prob)
  #calculate the posterior probability	
  EAP.JP=simplify2array(Joint.prob)   #Save the joint proability.	
  EAP.Theta=rowSums(matrix(1,I,1)%*%node.Quadpts*EAP.JP)/rowSums(EAP.JP)	#Save the examinees ability.	
  EAP.WP=EAP.JP*simplify2array(lapply(as.list(node.Quadpts), function(node.Quadpts, Est.Theta){(node.Quadpts-Est.Theta)^2}, Est.Theta=EAP.Theta)) 
  #Calculate the weighted joint proability.	
  hauteur=node.Quadpts[2:n.Quadpts]-node.Quadpts[1:(n.Quadpts-1)]	
  base.JP=colSums(t((EAP.JP[,1:(n.Quadpts-1)]+EAP.JP[,2:n.Quadpts])/2)*hauteur)	
  base.WP=colSums(t((EAP.WP[,1:(n.Quadpts-1)]+EAP.WP[,2:n.Quadpts])/2)*hauteur)
  EAP.Theta.SE=sqrt(base.WP/base.JP)	#Calculate the SEs of estimated EAP theta
  Est.Theta=as.data.frame(list(Theta=EAP.Theta, Theta.SE=EAP.Theta.SE))	#Transfer the variable Theta to a dataframe
  
  #Compute model fit information
  N2loglike=-2*LogL			#-2Log-likelihood
  AIC=2*np+N2loglike		#AIC
  BIC=N2loglike+log(I)*np	#BIC
  Theta.uni=sort(unique(EAP.Theta))
  Theta.uni.len=length(Theta.uni)
  G2=NA
  df=NA
  G2.P=NA
  G2.ratio=NA
  G2.size=NA
  RMSEA=NA
  if (Theta.uni.len>=11){
    n.group=10
    cutpoint=rep(NA,n.group)
    cutpoint[1]=min(Theta.uni)-0.001
    cutpoint[11]=max(Theta.uni)+0.001
    for (i in 2:n.group){
      cutpoint[i]=Theta.uni[(i-1)*Theta.uni.len/n.group]
    }
    Index=cut(EAP.Theta, cutpoint, labels = FALSE)
  }
  if (Theta.uni.len>=3 & Theta.uni.len<11){
    n.group=Theta.uni.len-1
    cutpoint=rep(NA,n.group)
    cutpoint[1]=min(Theta.uni)-0.001
    cutpoint[n.group]=max(Theta.uni)+0.001
    if (Theta.uni.len>=4){
      for (i in 2:(Theta.uni.len-2)){
        cutpoint[i]=Theta.uni[i]
      }
    }
    Index=cut(EAP.Theta, cutpoint, labels = FALSE)
  }
  if (Theta.uni.len>=3){
    X2.item=matrix(NA, n.group, J)
    G2.item=matrix(NA, n.group, J)
    Index.Uni=unique(Index)
    for (k in 1:n.group){
      data.group=data[Index==Index.Uni[k],]
      Theta.group=EAP.Theta[Index==Index.Uni[k]]
      Obs.P=colMeans(data.group)
      Exp.P=Reduce('+',lapply(Theta.group, Prob.model, Model=Model, Par.est0=Par.est0, D=D))/nrow(data.group)
      Obs.P[Obs.P>=1]=0.9999
      Obs.P[Obs.P<=0]=0.0001
      Exp.P[Exp.P>=1]=0.9999
      Exp.P[Exp.P<=0]=0.0001
      X2.item[k,]=nrow(data.group)*(Obs.P-Exp.P)^2/(Exp.P*(1-Exp.P))
      Odds1=log(Obs.P/Exp.P)
      Odds2=log((1-Obs.P)/(1-Exp.P))
      G2.item[k,]=nrow(data.group)*(Obs.P*Odds1+(1-Obs.P)*Odds2)
    }
    X2=sum(colSums(X2.item, na.rm = T))
    G2=sum(2*colSums(G2.item, na.rm = T))
    df=J*(n.group-3)
    G2.P=1-pchisq(G2,df)
    G2.ratio=G2/df
    RMSEA=sqrt(((X2-df)/(nrow(data)-1))/X2)
  }else{
    warning('The frequence table is too small to do fit tests.')
  }
  fits.test=list(G2=G2, G2.df=df, G2.P=G2.P, G2.ratio=G2.ratio, RMSEA=RMSEA, AIC=AIC, BIC=BIC)
  
  Time.End=Sys.time()		#End time
  Elapsed.time=paste('Elapsed time:', as.character(round(difftime(Time.End, Time.Begin, units="mins"), digits = 4)), 'mins')
  #Running time
  return(list(Est.ItemPars=Est.ItemPars, Est.Theta=Est.Theta, Loglikelihood=LogL, Iteration=n.ECycle, EM.Map=EM.Map,
              fits.test=fits.test, Elapsed.time=Elapsed.time, StopNormal=StopNormal, InitialValues=InitialValues, cr=cr))
}
 