BEMM.1PLAG=function(data, 					#A matrix of response [n.examinees * n.items]
                 PriorAlpha=c(-1.9,1), 	    #The log normal prior for alpha parameters with default mean -1.9 and variance 1
                 PriorBeta=c(0,4),          #The normal prior for beta parameters with default mean 0 and variance 4
                 PriorGamma=c(-1.39,0.25),  #The normal prior for gamma parameters with default mean -1.39 and variance 0.25
                 InitialAlpha=NA,			#Initial values for alpha parameters, default is NA 
                 InitialBeta=NA,			#Initial values for beta parameters, default is NA 
                 InitialGamma=NA,			#Initial values for gamma parameters, default is NA 
                 Tol=0.0001,				#The tolerate threshold for convergence, default is 0.0001
                 max.ECycle=2000L,		    #The max of E-step iteration, default is 2000L
                 max.MCycle=100L,		    #The max of M-step iteration, default is 100L
                 n.decimal=3L,              #The decimal length of outputs parameters, default is 3L
                 n.Quadpts =31L,		    #The number of quadrature, default is 31L
                 Theta.lim=c(-6,6),         #The range the Theta, default is [-6,6]
                 Missing=-9,                #A number to indicate missing value, default is -9
                 ParConstraint=FALSE,   #A logical value to determine whether impose a range limitation for parameters, default is FALSE
                 BiasSE=FALSE){         #A logical value to determine whether directly estimating SEs from inversed Hession matrix, default is FALSE     
  
   
  
  Time.Begin=Sys.time()       #Recording starting time
  Model='1PLAG'               #Set the model to be 1PLAG
  D=1                         #Set Constant D to be 1
  
  ###Check Input variables and return processed results###        
  Check.results=Input.Checking(Model=Model, data=data, PriorAlpha=PriorAlpha, PriorBeta=PriorBeta, PriorGamma=PriorGamma, 
                               InitialAlpha=InitialAlpha, InitialBeta=InitialBeta, InitialGamma=InitialGamma, Tol=Tol, 
                               max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                               Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint, BiasSE=BiasSE)
  
  data=Check.results$data                    #Return the processed data
  data.simple=Check.results$data.simple      #Return the simplified data
  CountNum=Check.results$CountNum            #Return the counts of simplified data
  I=Check.results$I                          #Return the nrow of data
  J=Check.results$J                          #Return the ncol of data
  n.class=Check.results$n.class              #the ncol of data.simple
  
  PriorAlpha=Check.results$PriorAlpha        #Return the prior setting for Alpha parameter
  PriorBeta=Check.results$PriorBeta          #Return the prior setting for Beta parameter
  PriorGamma=Check.results$PriorGamma        #Return the prior setting for Gamma parameter
  Prior=list(PriorAlpha=PriorAlpha, PriorBeta=PriorBeta, PriorGamma=PriorGamma)   #create a list to save Priors
  
  InitialAlpha=Check.results$InitialAlpha    #Return the starting values for Alpha parameter
  InitialBeta=Check.results$InitialBeta      #Return the starting values for Beta parameter
  InitialGamma=Check.results$InitialGamma    #Return the starting values for Gamma parameter
  
  max.ECycle=Check.results$max.ECycle        #Return the max of Expectation cycles
  max.MCycle=Check.results$max.MCycle        #Return the max of Maximization cycles
  n.Quadpts=Check.results$n.Quadpts          #Return the number of Quadpts
  n.decimal=Check.results$n.decimal          #Return the number of decimal
  
  ParConstraint=Check.results$ParConstraint  #Return the value of ParConstraint
  BiasSE=Check.results$BiasSE                 #Return the value of BiasSE
  
  Par.est0=list(Alpha=InitialAlpha, Beta=InitialBeta, Gamma=InitialGamma) 			#Assemble parameters estimates to a list
  Par.SE0=list(SEAlpha=InitialAlpha*0, SEBeta=InitialBeta*0, SEGamma=InitialGamma*0)    #Assemble parameters SEs to a list
  np=J*2+1		#Obtain the number of estimated parameters
  
  ###Run the BEMM algorithms###
  Est.results=BEMM.1PLAG.est(Model=Model, data=data, data.simple=data.simple, CountNum=CountNum, n.class=n.class,
                            Prior=Prior, Par.est0=Par.est0, Par.SE0=Par.SE0, D=D, np, 
                            Tol=Tol, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                            Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint, BiasSE=BiasSE, I=I, J=J, Time.Begin=Time.Begin)
  
  #Print important information into screen after estimation						
  if (Est.results$StopNormal==1){
    message('PROCEDURE TERMINATED NORMALLY')
  }else{
    message('PROCEDURE TERMINATED WITH ISSUES')
  }
  message('IRTEMM version: 1.0.5') 
  message('Item Parameter Calibration for the 1PL-AG Model.','\n')
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



BEMM.1PLAG.est=function(Model=Model, data=data, data.simple=data.simple, CountNum=CountNum, n.class=n.class,
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
  IAlpha=0                #EM iteration history for inverse second derivative of alpha parameters
  IBeta=rep(0,J)              #EM iteration history for inverse second derivative of beta parameters
  IGamma=rep(0,J)             #EM iteration history for inverse second derivative of gamma parameters
  TAlpha=rep(0,max.ECycle)    #EM iteration history for alpha parameters
  TBeta=matrix(0,max.ECycle,J)   #EM iteration history for beta parameters
  TGamma=matrix(0,max.ECycle,J)  #EM iteration history for gamma parameters
  deltahat.Alpha=0            #SEM iteration history for alpha parameters
  deltahat.Beta=rep(0,J)      #SEM iteration history for beta parameters
  deltahat.Gamma=rep(0,J)     #SEM iteration history for gamma parameters
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
      #estimate beta and gamma parameter
      gt0 = Par.est0$Gamma[j]
      bt0 = Par.est0$Beta[j]
      n.MCycle = 1
      M.exit = 0L
      #M step iteration
      while (M.exit==0L && (n.MCycle <= max.MCycle)){
        pstar = 1 / (1 + exp(-(node.Quadpts - bt0)))
        ag = 1 / (1 + exp(-(Par.est0$Alpha * node.Quadpts + gt0)))
        lg1 = sum(r[j,] - rz[j,] - (f - fz[j,]) * ag) #lg1: r - rz - (f - fz)*ps
        lgg = sum(-(f - fz[j,]) * ag * (1 - ag)) 			#lgg: -(f-fz)*(ag*(1-ag))
        lb1 = sum(-(rz[j,] - f * pstar))        			#lb1: - (rz - f * ps)
        lbb = sum(-(f * pstar * (1 - pstar)))  		    #lbb: - (f * ps * (1-ps))
        # Maximize beta and gamma
        if (Prior$PriorBeta[j]!=-9 && Prior$PriorBeta[j + J]!=-9) {
          lb1 = lb1 -((bt0-Prior$PriorBeta[j])/Prior$PriorBeta[j + J])
          lbb = lbb -1/Prior$PriorBeta[j + J]
        }
        Ibb = - 1 / lbb
        bt1 = bt0 + (Ibb * lb1)
        if (abs(bt1-bt0) >= 0.01){
          bt0 = bt1
        }
        
        if (Prior$PriorGamma[j]!=-9 && Prior$PriorGamma[j + J]!=-9) {
          lg1 = lg1 -((gt0-Prior$PriorGamma[j])/Prior$PriorGamma[j + J])
          lgg = lgg -1/Prior$PriorGamma[j + J]
        }
        Igg = - 1 / lgg
        gt1 = gt0 + (Igg * lg1)
        if (abs(gt1-gt0) >= 0.01){
          gt0 = gt1
        }
        if (abs(bt1-bt0) < 0.01 && abs(gt1-gt0) < 0.01){
          bt0 = bt1
          gt0 = gt1
          M.exit = 1L
        } else{
          bt0 = bt1
          gt0 = gt1
          n.MCycle = n.MCycle+1
        }
      }
      if (is.finite(gt0) && is.finite(bt0)){
        if (ParConstraint){
          if (bt0>=-6 && bt0<=6){
            Par.est0$Beta[j] = bt0
            TBeta[n.ECycle,j] = bt0
            IBeta[j] = Ibb
          }
          if (gt0>=-7 && gt0<=0){
            Par.est0$Gamma[j] = gt0
            TGamma[n.ECycle,j] = gt0
            IGamma[j] = Igg
          }
        }else{
          Par.est0$Beta[j] = bt0
          Par.est0$Gamma[j] = gt0
          TBeta[n.ECycle,j] = bt0
          TGamma[n.ECycle,j] = gt0
          IBeta[j] = Ibb
          IGamma[j] = Igg
        }
      }else{
        if (n.ECycle!=1){
          TBeta[n.ECycle,j] = TBeta[n.ECycle-1,j]
          TGamma[n.ECycle,j] = TGamma[n.ECycle-1,j]
        }else{
          TBeta[n.ECycle,j] = Par.est0$Beta[j]
          TGamma[n.ECycle,j] = Par.est0$Gamma[j]
        }
      }
    }
    #M step iteration for alpha
    at0 = Par.est0$Alpha
    n.MCycle = 1
    M.exit = 0
    while (M.exit==0L && (n.MCycle <= max.MCycle)){
      la1=0
      laa=0
      for (j in 1:J){
        ag = 1 / (1 + exp(-(at0 * node.Quadpts + Par.est0$Gamma[j])))
        la1 = la1 + at0 * node.Quadpts * (r[j,] - rz[j,] - (f - fz[j,]) * ag)
        #la1: exp(log(at))*X*(r - rz - (f - fz)*ps)
        laa = laa - at0 * at0 * node.Quadpts * node.Quadpts * (f - fz[j,]) * ag * (1 - ag)				  
        #laa: -exp(2*log(at))*X*X*(f-fz)*(ag*(1-ag))
      }
      la1=sum(la1)
      laa=sum(laa)
      # Maximize alpha
      if (Prior$PriorAlpha[1]!=-9 && Prior$PriorAlpha[2]!=-9) {
        la1 = la1 -((log(at0)-Prior$PriorAlpha[1])/Prior$PriorAlpha[2])
        laa = laa -1/Prior$PriorAlpha[2]
      }
      Iaa = - 1 / laa
      at1 = log(at0) + (Iaa * la1)
      if (abs(at1-log(at0)) < 0.01){
        at0 = exp(at1)
        M.exit = 1
      } else{
        at0 = exp(at1)
        n.MCycle = n.MCycle+1
      }
    }  
    if (is.finite(at0)){
      if (at0>=0.001){
        Par.est0$Alpha = at0
        TAlpha[n.ECycle] = at0
        IAlpha = Iaa
      }
    }else{
      if (n.ECycle!=0){
        TAlpha[n.ECycle] = TAlpha[n.ECycle-1]
      }else{
        TAlpha[n.ECycle] = Par.est0$Alpha
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
    delta=rep(0,3)
    delta0=rep(0,3)
    delta1=rep(0,3)
    cr.SEM0=1
    cr.SEM1=1
    cr.SEM2=1
    cr.SEM3=1
    
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

    # Estimating SEs of alpha
    z=start.SEM
    SEM.exit=0
    cr.SEM1=1
    deltahat=Par.est0
    while (SEM.exit==0 && z<=end.SEM){
      deltahat$Alpha=TAlpha[z]
      # Update Log-likelihood
      LLinfo=LikelihoodInfo(data.simple, CountNum, Model, deltahat, n.Quadpts, node.Quadpts, weight.Quadpts, D)
      f=LLinfo$f
      r=LLinfo$r
      fz=LLinfo$fz
      rz=LLinfo$rz 
      
      #Estimating alpha	
      at0=deltahat$Alpha
      n.MCycle = 1
      M.exit = 0
      while (M.exit==0L && (n.MCycle <= max.MCycle)){
        la1=0
        laa=0
        for (j in 1:J){
          ag = 1 / (1 + exp(-(at0 * node.Quadpts + deltahat$Gamma[j])))
          la1 = la1 + at0 * node.Quadpts * (r[j,] - rz[j,] - (f - fz[j,]) * ag)
          #la1: exp(log(at))*X*(r - rz - (f - fz)*ps)
          laa = laa - at0 * at0 * node.Quadpts * node.Quadpts * (f - fz[j,]) * ag * (1 - ag)				  
          #laa: -exp(2*log(at))*X*X*(f-fz)*(ag*(1-ag))
        }
        la1=sum(la1)
        laa=sum(laa)
        # Maximize alpha
        if (Prior$PriorAlpha[1]!=-9 && Prior$PriorAlpha[2]!=-9) {
          la1 = la1 -((log(at0)-Prior$PriorAlpha[1])/Prior$PriorAlpha[2])
          laa = laa -1/Prior$PriorAlpha[2]
        }
        Iaa = - 1 / laa
        at1 = log(at0) + (Iaa * la1)
        if (abs(at1-log(at0)) < 0.01){
          at0 = exp(at1)
          M.exit = 1
        } else{
          at0 = exp(at1)
          n.MCycle = n.MCycle+1
        }
      }  
      if (is.finite(at0)){
        if (at0>=0.001){
          deltahat$Alpha = at0
        }
      }
      delta1[1]=(log(deltahat$Alpha)-log(Par.est0$Alpha))/(log(TAlpha[z])-log(Par.est0$Alpha)+0.0001)
      cr.SEM1=abs(delta1[1]-delta0[1])	
      if (cr.SEM1<sqrt(Tol) && z>=2){
        SEM.exit=1
      }else{
        z=z+1
      }
      if (is.finite(delta1[1])){
        delta0[1] = delta1[1]
      }
    }
    delta[1]=1-delta0[1]
    delta1[1] =  1 / delta[1]
    if (is.finite(delta1[1])==0 || delta1[1]<=0){delta1[1] = 1}
    Par.SE0$SEAlpha= sqrt(Par.est0$Alpha * Par.est0$Alpha * IAlpha * delta1[1])
    if (Par.SE0$SEAlpha>1){Par.SE0$SEAlpha= sqrt(Par.est0$Alpha * Par.est0$Alpha * IAlpha)}
    
    #Estimating SEs of beta & gamma
    for (j in 1:J){
      z=start.SEM
      SEM.exit=0
      cr.SEM2=1
      cr.SEM3=1
      while (SEM.exit==0 && z<=end.SEM){
        for (ParClass in 1:2){
          if (ParClass==1 && z>=2 && cr.SEM2<sqrt(Tol)){
            next
          }
          if (ParClass==2 && z>=2 && cr.SEM3<sqrt(Tol)){
            next
          }
          deltahat=Par.est0
          if (ParClass==1){deltahat$Beta[j]=TBeta[z,j]}
          if (ParClass==2){deltahat$Gamma[j]=TGamma[z,j]}
          LLinfo=LikelihoodInfo(data.simple, CountNum, Model, deltahat, n.Quadpts, node.Quadpts, weight.Quadpts, D)
          f=LLinfo$f
          r=LLinfo$r
          fz=LLinfo$fz
          rz=LLinfo$rz  
          #estimate beta & gamma parameters
          bt0 = deltahat$Beta[j]
          gt0 = deltahat$Gamma[j]
          n.MCycle = 1
          M.exit = 0
          #M step iteration
          while (M.exit==0 && (n.MCycle <= max.MCycle)) {
            if (ParClass==1){
              #estimate beta parameter
              pstar = 1 / (1 + exp(-(node.Quadpts - bt0)))
              lb1 = sum(-(rz[j,] - f * pstar))        					            #lb1: - (rz - f * ps)
              lbb = sum(-(f * pstar * (1 - pstar)))  		                    #lbb: - (f * ps * (1-ps))
              if (Prior$PriorBeta[j]!=-9 && Prior$PriorBeta[j + J]!=-9) {
                lb1 = lb1 -((bt0-Prior$PriorBeta[j])/Prior$PriorBeta[j + J])
                lbb = lbb -1/Prior$PriorBeta[j + J]
              }
              Ibb = - 1 / lbb
              bt1 = bt0 + (Ibb * lb1)
              if (abs(bt1-bt0) < 0.01){
                bt0 = bt1
                M.exit = 1
              } else{
                bt0 = bt1
                n.MCycle = n.MCycle + 1
              }
            }
            if (ParClass==2){
              #estimate gamma parameter
              ag = 1 / (1 + exp(-(Par.est0$Alpha * node.Quadpts + gt0)))
              lg1 = sum(r[j,] - rz[j,] - (f - fz[j,]) * ag)   #lg1: r - rz - (f - fz)*ps
              lgg = sum(-(f - fz[j,]) * ag * (1 - ag)) 			  #lgg: -(f-fz)*(ag*(1-ag))
              if (Prior$PriorGamma[j]!=-9 && Prior$PriorGamma[j + J]!=-9) {
                lg1 = lg1 -((gt0-Prior$PriorGamma[j])/Prior$PriorGamma[j + J])
                lgg = lgg -1/Prior$PriorGamma[j + J]
              }
              Igg = - 1 / lgg
              gt1 = gt0 + (Igg * lg1)
              if (abs(gt1-gt0) < 0.01){
                gt0 = gt1
                M.exit =1
              } else{
                gt0 = gt1
                n.MCycle = n.MCycle + 1
              }
            }
          }     
          if (is.finite(gt0) && is.finite(bt0)){
            if (ParConstraint){
              if (bt0>=-6 && bt0<=6){
                deltahat$Beta[j] = bt0
              }
              if (gt0>=-7 && gt0<=0){
                deltahat$Gamma[j] = gt0
              }
            }else{
              deltahat$Beta[j] = bt0
              deltahat$Gamma[j] = gt0
            }
          }
          if (ParClass==1){
            delta1[2]=(deltahat$Beta[j]-Par.est0$Beta[j])/(TBeta[z,j]-Par.est0$Beta[j]+0.0001)
            break
          }
          if (ParClass==2){
            delta1[3]=(deltahat$Gamma[j]-Par.est0$Gamma[j])/(TGamma[z,j]-Par.est0$Gamma[j]+0.0001)
            break
          }
        }
        cr.SEM2=abs(delta1[2]-delta0[2])
        cr.SEM3=abs(delta1[3]-delta0[3])
        if (cr.SEM2<sqrt(Tol) && cr.SEM3<sqrt(Tol) && z>=2){
          SEM.exit=1
        }else{
          z=z+1
        }
        delta0[is.finite(delta1)] = delta1[is.finite(delta1)]
      }
      delta[2]=1-delta0[2]
      delta[3]=1-delta0[3]
      delta1[2] =  1 / delta[2]
      delta1[3] =  1 / delta[3]
      if (is.finite(delta1[2])==F || delta1[2]<=0){delta1[2] = 1}
      if (is.finite(delta1[3])==F || delta1[3]<=0){delta1[3] = 0}
      Par.SE0$SEBeta[j]= sqrt(IBeta[j] * delta1[2])
      Par.SE0$SEGamma[j]= sqrt(IGamma[j] * delta1[3])
      if (Par.SE0$SEBeta[j]>1){Par.SE0$SEBeta[j]= sqrt(IBeta[j])}
      if (Par.SE0$SEGamma[j]>1){Par.SE0$SEGamma[j]= sqrt(IGamma[j])}
    }
  }else{
    message('Directly estimating SEs from inversed Hession matrix.', '\n')
    Par.SE0$SEAlpha= sqrt(Par.est0$Alpha * Par.est0$Alpha * IAlpha)
    for (j in 1:J){
      Par.SE0$SEBeta[j]= sqrt(IBeta[j])
      Par.SE0$SEGamma[j]= sqrt(IGamma[j])
    }
  }
  
  Par.est0$Alpha=round(Par.est0$Alpha, n.decimal)		#Save the estimated beta parameters
  Par.est0$Beta=round(Par.est0$Beta, n.decimal)		#Save the estimated beta parameters
  Par.est0$Gamma=round(Par.est0$Gamma, n.decimal)		#Save the estimated gamma parameters
  
  Par.SE0$SEAlpha=round(Par.SE0$SEAlpha, n.decimal)		#Save the estimated beta parameters
  Par.SE0$SEBeta=round(Par.SE0$SEBeta, n.decimal)		#Save the estimated SEs of b parameters
  Par.SE0$SEGamma=round(Par.SE0$SEGamma, n.decimal)		#Save the estimated SEs of g parameters
  
  EM.Map=list(Map.Alpha=TAlpha[1:n.ECycle], Map.Beta=TBeta[1:n.ECycle,], Map.Gamma=TGamma[1:n.ECycle,]) #Transfer the variable EM.map to a list
  
  Est.ItemPars=as.data.frame(list(est.alpha=Par.est0$Alpha, est.beta=Par.est0$Beta, est.gamma=Par.est0$Gamma, se.alpha=Par.SE0$SEAlpha, se.beta=Par.SE0$SEBeta, se.gamma=Par.SE0$SEGamma))
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
