BE3M.4PL=function(data, 				#A matrix of respoonse [n.examinees * n.items]
                 PriorA=c(0,0.25), 	    #The log normal prior for a parameters with default mean 0 and variance 0.25
                 PriorB=c(0,4),         #The normal prior for b parameters with default mean 0 and variance 4
                 PriorC=c(4,16),        #The beta prior for c parameters with default hyper-parameter 4 and 16
                 PriorS=c(4,16),        #The beta prior for s parameters with default hyper-parameter 4 and 16
                 InitialA=NA,			#Initial values for a parameters, default is NA 
                 InitialB=NA,			#Initial values for b parameters, default is NA 
                 InitialC=NA,			#Initial values for c parameters, default is NA 
                 InitialS=NA,			#Initial values for s parameters, default is NA 
                 Tol=0.001,				#The tolerate threshold for convergnece, default is 0.001
                 max.ECycle=2000L,		#The max of Estem interation, default is 2000L
                 max.MCycle=30L,		#The max of Mstep interation, default is 30L
                 n.decimal=3L,          #The decimal length of outputs parameters, default is 3L
                 n.Quadpts =31L,		#The number of quadratures, default is 31L
                 Theta.lim=c(-6,6),     #The range the Theta, default is [-6,6]
                 Missing=-9,            #A number to indicate missing value, default is -9
                 ParConstraint=FALSE){  #A logical value to determine whether impose a range limitation for parameters, default is FALSE
  
  Time.Begin=Sys.time()  #Recoding starting time
  Model='4PL'            #Set the model to be 3PLM
  D=1.702                #Set the constant D to 1.702
  
  ###Check Input variables and return processed results###        
  Check.results=Input.Checking(Model=Model, data=data, PriorA=PriorA, PriorB=PriorB, PriorC=PriorC, PriorS=PriorS, 
                 InitialA=InitialA, InitialB=InitialB, InitialC=InitialC, InitialS=InitialS, 
                 Tol=Tol, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.Quadpts =n.Quadpts, n.decimal=n.decimal,
                 Theta.lim=Theta.lim, Missing=Missing, ParConstraint=ParConstraint)
  
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
  max.MCycle=Check.results$max.MCycle         #Return the max of Mxpectation cycles
  n.Quadpts=Check.results$n.Quadpts           #Return the number of Quadpts
  n.decimal=Check.results$n.decimal           #Return the number of decimal
  
  ParConstraint=Check.results$ParConstraint   #Return the value of ParConstraint
  
  Par.est0=list(A=InitialA, B=InitialB, C=InitialC, S=InitialS) 			    #Assemble parameters estimates to a list
  Par.SE0=list(SEA=InitialA*0, SEB=InitialB*0, SEC=InitialC*0, SES=InitialS*0)  #Assemble parameters SEs to a list
  
  
  #Generating the nodes and their weights for approximating the integration
  node.Quadpts=seq(Theta.lim[1],Theta.lim[2],length.out = n.Quadpts)	#Generating the nodes
  weight.Quadpts=dnorm(node.Quadpts,0,1)                                #Generating the weight of nodes
  weight.Quadpts=weight.Quadpts/sum(weight.Quadpts)
  node.Quadpts.list=as.list(node.Quadpts)                               #Create a list for nodes
  weight.Quadpts.list=as.list(weight.Quadpts)                           #Create a list for weights of nodes
  
  #Initializing the program settings
  LH=rep(0,max.ECycle)		#Log-likelihood
  IA=rep(0,J)  #EM iteration history for invered second derivative of a parameters
  IB=rep(0,J)  #EM iteration history for invered second derivative of b parameters
  IC=rep(0,J)  #EM iteration history for invered second derivative of c parameters
  IS=rep(0,J)  #EM iteration history for invered second derivative of s parameters
  IAB=rep(0,J) #EM iteration history for invered second derivative of ab parameters
  TA=rep(0,max.ECycle*J)  #EM iteration history for a parameters
  TB=rep(0,max.ECycle*J)  #EM iteration history for b parameters
  TC=rep(0,max.ECycle*J)  #EM iteration history for c parameters
  TS=rep(0,max.ECycle*J)  #EM iteration history for s parameters  
  deltahat.A=rep(0,J)     #SEM iteration history for a parameters
  deltahat.B=rep(0,J)     #SEM iteration history for b parameters
  deltahat.C=rep(0,J)     #SEM iteration history for c parameters
  deltahat.S=rep(0,J)     #SEM iteration history for S parameters
  
  

  f=rep(0,n.Quadpts)			#Aritifical data f
  r=rep(0,n.Quadpts*J)          #Aritifical data r
  fz=rep(0,n.Quadpts*J)         #Aritifical data fz
  rz=rep(0,n.Quadpts*J)			#Aritifical data rz
  P=rep(0,n.Quadpts*J)          #Probability for 3PLM
  Pstar=rep(0,n.Quadpts*J)      #Probability for 2PLM
  LL=rep(0,n.Quadpts*n.class)   #Joint probability
  LL0=rep(0,n.class)            #Sum of joint probability
  Posterior.prob=rep(0,n.Quadpts*n.class)   #Posterior probability
  cr=0							#The difference between last and current log-likelihood
  n.ECycle=0L					#The first E-step iteration
  StopNormal=1L					#Whether program terminate normally
  
  #Call BE3M estimation program from C
  
  res=.C("BE3M4PL", data.simple=data.simple, CountNum=CountNum, n.class=n.class, J=J,  LH=LH, 
         IA=IA, IB=IB, IC=IC, IS=IS, IAB=IAB, TA=TA, TB=TB, TC=TC, TS=TS,
         deltahat.A=deltahat.A, deltahat.B=deltahat.B, deltahat.C=deltahat.C, deltahat.S=deltahat.S,
         est.A=Par.est0$A, est.B=Par.est0$B, est.C=Par.est0$C, est.S=Par.est0$S, Tol=Tol, cr=cr, StopNormal=StopNormal,
         se.A=Par.SE0$SEA, se.B=Par.SE0$SEB, se.C=Par.SE0$SEC, se.S=Par.SE0$SES, ParConstraint=ParConstraint,  
         f=f,  r=r, fz=fz,  rz=rz,  P=P, Pstar=Pstar,  LL=LL, LL0=LL0,  Posterior.prob=Posterior.prob,
         PriorA=PriorA, PriorB=PriorB, PriorC=PriorC, PriorS=PriorS, n.Quadpts=n.Quadpts, node.Quadpts=node.Quadpts, 
         weight.Quadpts=weight.Quadpts, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.ECycle=n.ECycle, D=D)   #return the results  
  
  cr=res$cr		#Save the convergent point
  LH=res$LH		#Save the final log-likelihood
  StopNormal=1-res$StopNormal		#Save the results whether program terminate normally
  
  
  n.ECycle=res$n.ECycle                         #Save the number of final E-step iteration
  Par.est0$A=round(res$est.A, n.decimal)		#Save the estimated a parameters
  Par.est0$B=round(res$est.B, n.decimal)		#Save the estimated b parameters
  Par.est0$C=round(res$est.C, n.decimal)		#Save the estimated c parameters
  Par.est0$S=round(res$est.S, n.decimal)		#Save the estimated c parameters
  
  Par.SE0$SEA=round(res$se.A, n.decimal)		#Save the estimated SEs of a parameters 
  Par.SE0$SEB=round(res$se.B, n.decimal)		#Save the estimated SEs of b parameters
  Par.SE0$SEC=round(res$se.C, n.decimal)		#Save the estimated SEs of c parameters
  Par.SE0$SES=round(res$se.S, n.decimal)		#Save the estimated SEs of c parameters
  
  EM.Map.A=matrix(res$TA[1:(n.ECycle*J)], ncol=J, byrow=T)  #Save the EM Map for a parameters
  EM.Map.B=matrix(res$TB[1:(n.ECycle*J)], ncol=J, byrow=T)  #Save the EM Map for b parameters
  EM.Map.C=matrix(res$TC[1:(n.ECycle*J)], ncol=J, byrow=T)  #Save the EM Map for c parameters
  EM.Map.S=matrix(res$TS[1:(n.ECycle*J)], ncol=J, byrow=T)  #Save the EM Map for s parameters
  EM.Map=list(Map.A=EM.Map.A,Map.B=EM.Map.B, Map.C=EM.Map.C, Map.S=EM.Map.S) #Transfer the variable EM.map to a list
  
  Est.ItemPars=as.data.frame(list(est.a=Par.est0$A, est.b=Par.est0$B, est.c=Par.est0$C, est.s=Par.est0$S, se.a=Par.SE0$SEA, se.b=Par.SE0$SEB, se.c=Par.SE0$SEC, se.s=Par.SE0$SES))
		#Transfer the variable estimated item parameters to a dataframe
		
 #Calulate examinees ability via EAP methods.
		
  P.Quadpts=lapply(node.Quadpts.list, Prob.model, Model=Model, Par.est0=Par.est0, D=D)	
  Joint.prob=mapply('*',lapply(P.Quadpts, function(P,data){apply(data*P+(1-data)*(1-P),2,prod,na.rm = T)}, data=t(data)),
                    weight.Quadpts.list, SIMPLIFY = FALSE)		
  Whole.prob=Reduce("+", Joint.prob)		
  LogL=sum(log(Whole.prob))      #Obtain the final Log-likelihood
  Posterior.prob=lapply(Joint.prob,function(P,Whole.prob){P/Whole.prob}, Whole.prob=Whole.prob)
			#calculate the posterior probability	
  EAP.JP=simplify2array(Joint.prob)   #Save the joint proability.	
  EAP.Theta=rowSums(matrix(1,I,1)%*%node.Quadpts*EAP.JP)/rowSums(EAP.JP)	#Save the examinees ability.	
  EAP.WP=EAP.JP*simplify2array(lapply(node.Quadpts.list, function(node.Quadpts, Est.Theta){(node.Quadpts-Est.Theta)^2}, Est.Theta=EAP.Theta)) 
		#Calculate the weighted joint proability.	
  hauteur=node.Quadpts[2:n.Quadpts]-node.Quadpts[1:(n.Quadpts-1)]	
  base.JP=colSums(t((EAP.JP[,1:(n.Quadpts-1)]+EAP.JP[,2:n.Quadpts])/2)*hauteur)	
  base.WP=colSums(t((EAP.WP[,1:(n.Quadpts-1)]+EAP.WP[,2:n.Quadpts])/2)*hauteur)
  EAP.Theta.SE=sqrt(base.WP/base.JP)	#Calculate the SEs of estimated EAP theta
  Est.Theta=as.data.frame(list(Theta=EAP.Theta, Theta.SE=EAP.Theta.SE))	#Transfer the variable Theta to a dataframe
  
  
  
  #Compute model fit information
  np=J*3		#Obtain the number of estimated parameters
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
    df=J*(n.group-4)
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
  							
  #Print important information into screen after estimation						

  if (StopNormal==1){
    message('PROCEDURE TERMINATED NORMALLY')
  }else{
    message('PROCEDURE TERMINATED WITH ISSUES')
  }
  message('IRTEMM version: 1.0.2') 
  message('Item Parameter Calibration for the 4PLM.','\n')
  message('Quadrature: ', n.Quadpts, ' nodes from ', Theta.lim[1], ' to ', Theta.lim[2], ' were used to approximate Gaussian distribution.') 
  message('Method for Items: Ability-based Bayesian Expectation-Maximization-Maximization-Maximization (BE3M) Algorithm.')
  message('Method for Item SEs: Updated Supplemented EM.')
  message('Method for Theta: Expected A Posteriori (EAP).')
  if (StopNormal==1){
    message('Converged at LL-Change=', cr, ' after ', n.ECycle, ' EMM iterations.', sep = '')
  }else{
    warning('Warning: Estimation cannot converged under current max.ECycle and Tol!', sep = '')
    warning('Warning: The reults maybe not trustworthy!', sep = '')
    message('Terminated at LL-Change=', cr, ' after ', n.ECycle, ' EMM iterations.', sep = '')
  }
  message('Running time:', Elapsed.time, '\n')
  message('Log-likelihood (LL):', as.character(round(LogL, n.decimal)))
  message('Estimated Parameters:', as.character(np))
  message('AIB: ', round(fits.test$AIC, n.decimal), ', BIC: ', round(fits.test$BIC, n.decimal), ', RMSEA = ', round(fits.test$RMSEA, n.decimal))
  message('G2 (', round(fits.test$G2.df, n.decimal), ') = ', round(fits.test$G2, n.decimal), ', p = ', round(fits.test$G2.P, n.decimal), ', G2/df = ', round(fits.test$G2.ratio, n.decimal), sep='')		
											
							
  return(list(Est.ItemPars=Est.ItemPars, Est.Theta=Est.Theta, Loglikelihood=LogL, Iteration=n.ECycle, EM.Map=EM.Map,
              fits.test=fits.test, Elapsed.time=Elapsed.time))
							#Return results
}

