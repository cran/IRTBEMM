Prob.model<- function(X, Model, Par.est0, D=1.702){
  if (Model=='3PL' || Model =='4PL' || Model =='1PLAG' || Model =='1PLG'){
	  if (Model=='3PL'){
		Prob=Par.est0$C+(1-Par.est0$C)/(1+exp(-D*Par.est0$A*(X-Par.est0$B)))    #the prob of 3PLM
	  }
	  if (Model=='4PL'){
		Prob=Par.est0$C+(1-Par.est0$S-Par.est0$C)/(1+exp(-D*Par.est0$A*(X-Par.est0$B)))    #the prob of 4PLM
	  }
	  if (Model=='1PLAG'){
		P.1pl=1/(1+exp(-(X-Par.est0$Beta)))    #the prob of 1PLM
		P.ag=1/(1+exp(-(Par.est0$Alpha*X+Par.est0$Gamma)))    #the prob of ability-based guessing
		Prob=P.1pl+(1-P.1pl)*P.ag    #the prob of 1PLAG
	  }
	  if (Model=='1PLG'){
		P.1pl=1/(1+exp(-(X-Par.est0$Beta)))    #the prob of 1PLM
		P.g=1/(1+exp(-(Par.est0$Gamma)))    #the prob of ability-based guessing
		Prob=P.1pl+(1-P.1pl)*P.g    #the prob of 1PLAG
	  }
	  
	  Prob[Prob>0.9999]=0.9999   #constrain the max value to 0.9999
	  Prob[Prob<0.0001]=0.0001   #constrain the min value to 0.9999
	  
	  return(Prob)  
  }else{
	stop('The Model user specified does not exist!')
  }
}

###Checking whether the input variables satisfy the requirements###
Input.Checking <- function(Model, 
                            data, 
							PriorA=c(0,0.25), 
							PriorB=c(0,4), 
							PriorC=c(4,16), 
							PriorS=c(4,16), 
							PriorAlpha=c(-1.9,1), 
							PriorBeta=c(0,4), 
							PriorGamma=c(-1.39,0.25), 
							InitialA=NA, 
							InitialB=NA, 
							InitialC=NA,
							InitialS=NA, 
							InitialAlpha=NA, 
							InitialBeta=NA, 
							InitialGamma=NA, 
							Tol=0.0001, 
							max.ECycle=1000L, 
							max.MCycle=100L, 
							n.Quadpts=31L, 
							n.decimal=3L,
							Theta.lim=c(-6,6), 
							Missing=-9, 
							ParConstraint=FALSE){
  
  ###Checking data###
  if (is.data.frame(data)){
    data=data.matrix(data)
  }   
  if (is.matrix(data)){
    I=as.integer(nrow(data))
    J=as.integer(ncol(data))
    if (I==1 | J==1){
      stop('Error: The ncol and nrow of data must bigger than 1!')
    }else{
      if (sum(is.na(data))!=0){
        stop('Error: Some elements in data are not 1, 0 or Missing!')
      }
      if (sum(data!=0 & data!=1 & data!=Missing)!=0){
        stop('Error: Some elements in data are not 1, 0 or Missing!')
      }else{
        if (sum(data==Missing)!=0){
			Index.miss=which(data==Missing, arr.ind = T)
			data[data==Missing]=0
			PI=rowMeans(data)
			PJ=colMeans(data)
			for (i in 1:nrow(Index.miss)){
				PI.tmp=PI[Index.miss[i,1]]
				PJ.tmp=PJ[Index.miss[i,2]]
				P.correct=0.5+PI.tmp-PJ.tmp
				if (is.na(P.correct)){
					data[Index.miss[i,1],Index.miss[i,2]]=0
				}else{
					if (P.correct>=1){
						data[Index.miss[i,1],Index.miss[i,2]]=1
					}else{
						if (P.correct<=0){
							data[Index.miss[i,1],Index.miss[i,2]]=0
						}else{
							data[Index.miss[i,1],Index.miss[i,2]]=rbinom(1,1,P.correct)
						}
					}
				}
			}
        }
		datafull=as.data.frame(data)
		data.class=data.matrix(aggregate(list(Num=rep(1,I)), datafull, length))
		data.simple=as.numeric(data.class[,1:J])
		n.class=as.integer(nrow(data.class))
		CountNum=as.integer(data.class[,J+1])
		data.list=list(data=data, data.simple=data.simple, CountNum=CountNum, n.class=n.class, I=I, J=J)
      }
    }
  }else{
    stop('Error: The type of data must be a matrix or a data.framework!')
  }
  
  
  ###Checking Priors for each parameters###
  if (Model=='3PL' | Model=='4PL'){
    ###Checking Variable PriorA###
    if (is.numeric(PriorA) & length(PriorA)==2 & is.na(sum(PriorA))==FALSE){  #Whether PriorA=c(0,0.25)
      if (PriorA[2]<=0){
        stop('Error: PriorA[2] is the variance, and it must bigger than 0!')
      }else{
        PriorA=rep(PriorA, each = J)
      }
    }else{
      if (length(PriorA)==1){ #Whether PriorA=NA
        if (is.na(PriorA)){
          PriorA=rep(-9, 2*J)
        }else{
          stop('Error: PriorA must have two input values unless PriorA=NA!')
        }
      }else{
        if (is.matrix(PriorA) & length(PriorA)==J*2){                          #Whether PriorA is a matrix
          if (ncol(PriorA)==2 | nrow(PriorA)==J){
            if (sum(rowSums(is.na(PriorA))==1)!=0){
              stop('Error: The type of PriorA[1] and PriorA[2] are different!')
            }
            if (sum(PriorA[,2]<=0,na.rm = T)!=0){
              stop('Error: PriorA[2] is the variance, and it must bigger than 0!')
            }
            PriorA[is.na(PriorA)]=-9
            PriorA=as.numeric(PriorA)
          }else{
            stop('Error: The dim of matrix PriorA must be c(n.item, 2)!')
          }
        }else{
          stop('Error: The class of PriorA must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
        }
      }
    }  
    
    
    ###Checking Variable PriorB###
    if (is.numeric(PriorB) & length(PriorB)==2 & is.na(sum(PriorB))==FALSE){  #Whether PriorB=c(0,0.25)
      if (PriorB[2]<=0){
        stop('Error: PriorB[2] is the variance, and it must bigger than 0!')
      }else{
        PriorB=rep(PriorB, each = J)
      }
    }else{
      if (length(PriorB)==1){ #Whether PriorB=NA
        if (is.na(PriorB)){
          PriorB=rep(-9, 2*J)
        }else{
          stop('Error: PriorB must have two input values unless PriorB=NA!')
        }
      }else{
        if (is.matrix(PriorB) & length(PriorB)==J*2){                          #Whether PriorB is a matrix
          if (ncol(PriorB)==2 | nrow(PriorB)==J){
            if (sum(rowSums(is.na(PriorB))==1)!=0){
              stop('Error: The type of PriorB[1] and PriorB[2] are different!')
            }
            if (sum(PriorB[,2]<=0,na.rm = T)!=0){
              stop('Error: PriorB[2] is the variance, and it must bigger than 0!')
            }
            PriorB[is.na(PriorB)]=-9
            PriorB=as.numeric(PriorB)
          }else{
            stop('Error: The dim of matrix PriorB must be c(n.item, 2)!')
          }
        }else{
          stop('Error: The class of PriorB must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
        }
      }
    }
    
    ###Checking Variable PriorC###
    if (is.numeric(PriorC) & length(PriorC)==2 & is.na(sum(PriorC))==FALSE){  #Whether PriorC=c(4,16)
      if (PriorC[1]<=0 | PriorC[2]<=0){
        stop('Error: The elements in Beta Prior must bigger than 0!')
      }else{
        PriorC=rep(PriorC, each = J)
      }
    }else{
      if (length(PriorC)==1){ #Whether PriorC=NA
        if (is.na(PriorC)){
          PriorC=rep(-9, 2*J)
        }else{
          stop('Error: PriorC must have two input values unless PriorC=NA!')
        }
      }else{
        if (is.matrix(PriorC) & length(PriorC)==J*2){                          #Whether PriorC is a matrix
          if (ncol(PriorC)==2 | nrow(PriorC)==J){
            if (sum(rowSums(is.na(PriorC))==1)!=0){
              stop('Error: The type of PriorC[1] and PriorC[2] are different!')
            }
            if (sum(PriorC[,1]<=0,na.rm = T)!=0 | sum(PriorC[,2]<=0,na.rm = T)!=0){
              stop('Error: The elements in Beta Prior must bigger than 0!')
            }
            PriorC[is.na(PriorC)]=-9
            PriorC=as.numeric(PriorC)
          }else{
            stop('Error: Error: The dim of matrix PriorC must be c(n.item, 2)!')
          }
        }else{
          stop('Error: The class of PriorC must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
        }
      }
    }
    if (Model=='4PL'){
      ###Checking Variable PriorS###
      if (is.numeric(PriorS) & length(PriorS)==2 & is.na(sum(PriorS))==FALSE){  #Whether PriorS=c(4,16)
        if (PriorS[1]<1 | PriorS[2]<1){
          stop('Error: The elements in Beta Prior must bigger than 1!')
        }else{
          PriorS=rep(PriorS, each = J)
        }
      }else{
        if (length(PriorS)==1){ #Whether PriorS=NA
          if (is.na(PriorS)){
            PriorS=rep(-9, 2*J)
          }else{
            stop('Error: PriorS must have two input values unless PriorS=NA!')
          }
          PriorS[is.na(PriorS)]=-9
          PriorS=as.numeric(PriorS)
        }else{
          if (is.matrix(PriorS) & length(PriorS)==J*2){                          #Whether PriorS is a matrix
            if (ncol(PriorS)==2 | nrow(PriorS)==J){
              if (sum(rowSums(is.na(PriorS))==1)!=0){
                stop('Error: The type of PriorS[1] and PriorS[2] are different!')
              }
              if (sum(PriorS[,1]<1,na.rm = T)!=0 | sum(PriorS[,2]<1,na.rm = T)!=0){
                stop('Error: The elements in Beta Prior must bigger than 1!')
              }
            }else{
              stop('Error: The dim of matrix PriorS must be c(n.item, 2)!')
            }
          }else{
            stop('Error: The class of PriorS must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
          }
        }
      }
	  Prior.list=list(PriorA=PriorA, PriorB=PriorB, PriorC=PriorC, PriorS=PriorS)
    }else{
      Prior.list=list(PriorA=PriorA, PriorB=PriorB, PriorC=PriorC)
    }
  }
  
  #If Model is 1PL-AG
  if (Model=='1PLAG' | Model=='1PLG'){
    ###Checking Variable PriorAlpha###
    if (Model=='1PLAG'){
      if (is.numeric(PriorAlpha) & length(PriorAlpha)==2 & is.na(sum(PriorAlpha))==FALSE){  #Whether PriorAlpha=c(-1.9,1)
        if (PriorAlpha[2]<=0){
          stop('Error: The variance in Normal Prior must bigger than 0!')
        }else{
          PriorAlpha=PriorAlpha
        }
      }else{
        if (length(PriorAlpha)==1){ #Whether PriorAlpha=NA
          if (is.na(PriorAlpha)){
            PriorAlpha=rep(-9, 2)
          }else{
            stop('Error: PriorAlpha must have two input values unless PriorAlpha=NA!')
          }
          PriorAlpha[is.na(PriorAlpha)]=-9
          PriorAlpha=as.numeric(PriorAlpha)
        }else{
            stop('Error: The class of PriorAlpha must be NA or a numeric with length=2!')
        }
      }
    }
    
    ###Checking Variable PriorBeta###
    if (is.numeric(PriorBeta) & length(PriorBeta)==2 & is.na(sum(PriorBeta))==FALSE){  #Whether PriorBeta=c(0,0.25)
      if (PriorBeta[2]<=0){
        stop('Error: PriorBeta[2] is the variance, and it must bigger than 0!')
      }else{
        PriorBeta=rep(PriorBeta, each = J)
      }
    }else{
      if (length(PriorBeta)==1){ #Whether PriorBeta=NA
        if (is.na(PriorBeta)){
          PriorBeta=rep(-9, 2*J)
        }else{
          stop('Error: PriorBeta must have two input values unless PriorBeta=NA!')
        }
      }else{
        if (is.matrix(PriorBeta) & length(PriorBeta)==J*2){                          #Whether PriorBeta is a matrix
          if (ncol(PriorBeta)==2 | nrow(PriorBeta)==J){
            if (sum(rowSums(is.na(PriorBeta))==1)!=0){
              stop('Error: The type of PriorBeta[1] and PriorBeta[2] are different!')
            }
            if (sum(PriorBeta[,2]<=0,na.rm = T)!=0){
              stop('Error: PriorBeta[2] is the variance, and it must bigger than 0!')
            }
            PriorBeta[is.na(PriorBeta)]=-9
            PriorBeta=as.numeric(PriorBeta)
          }else{
            stop('Error: The dim of matrix PriorBeta must be c(n.item, 2)!')
          }
        }else{
          stop('Error: The class of PriorBeta must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
        }
      }
    }
  

    ###Checking Variable PriorGamma###
    if (is.numeric(PriorGamma) & length(PriorGamma)==2 & is.na(sum(PriorGamma))==FALSE){  #Whether PriorGamma=c(-1.39,0.25)
      if (PriorGamma[2]<=0){
        stop('Error: PriorGamma[2] is the variance, and it must bigger than 0!')
      }else{
        PriorGamma=rep(PriorGamma, each = J)
      }
    }else{
      if (length(PriorGamma)==1){ #Whether PriorGamma=NA
        if (is.na(PriorGamma)){
          PriorGamma=rep(-9, 2*J)
        }else{
          stop('Error: PriorGamma must have two input values unless PriorGamma=NA!')
        }
        }else{
        if (is.matrix(PriorGamma) & length(PriorGamma)==J*2){                          #Whether PriorGamma is a matrix
          if (ncol(PriorGamma)==2 | nrow(PriorGamma)==J){
            if (sum(rowSums(is.na(PriorGamma))==1)!=0){
              stop('Error: The type of PriorGamma[1] and PriorGamma[2] are different!')
            }
            if (sum(PriorGamma[,2]<=0,na.rm = T)!=0){
              stop('Error: PriorGamma[2] is the variance, and it must bigger than 0!')
            }
            PriorGamma[is.na(PriorGamma)]=-9
            PriorGamma=as.numeric(PriorGamma)
          }else{
            stop('Error: The dim of matrix PriorGamma must be c(n.item, 2)!')
          }
        }else{
          stop('Error: The class of PriorGamma must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!')
        }
      }
    }
    if (Model=='1PLAG'){
      Prior.list=list(PriorAlpha=PriorAlpha, PriorBeta=PriorBeta, PriorGamma=PriorGamma)
    }
    if (Model=='1PLG'){
      Prior.list=list(PriorBeta=PriorBeta, PriorGamma=PriorGamma)
    }
  }
  

  
  
  ###Checking Variable Initial###
  total_score=matrix(rowSums(data),ncol=1)     #calculate the total score for each individual
  corr0=t(cor(data,total_score,use = "complete.obs"))   #calculate the correlation between the each item and the total score
  PassRate=matrix(colSums(data)/I,nrow=1) #calculate the correct ratio for each item
  PassRate[PassRate>0.9999]=0.9999    #constrain the max value to 0.9999
  PassRate[PassRate<0.0001]=0.0001    #constrain the min value to 0.0001
  Zscore=-qnorm(PassRate,0,1)  # transform the correct ratio to the normal scale for each item
  if (Model=='3PL' | Model=='4PL'){
    Initial.B=as.numeric(Zscore/corr0)                 # estimating initial item difficult parameter
    Initial.A=as.numeric(corr0/sqrt(1-corr0^2))        # estimating initial item discirmation parameter
    Initial.B[Initial.B>2]=1.8    #constrain the max value to 2
    Initial.B[Initial.B< -2]=-1.8 #constrain the min value to 2
    Initial.A[Initial.A>2]=1.8    #constrain the max value to 1.8
    Initial.A[Initial.A<0.3]=0.5  #constrain the min value to 0.5
    Initial.C=rep(0.2,1,J)         # giving initial item guessing parameter
    if (Model=='3PL'){
      InitialValue=list(InitialA=Initial.A,InitialB=Initial.B,InitialC=Initial.C)
    }
    if (Model=='4PL'){
      Initial.S=rep(0.2,1,J)                 # giving initial item slipping parameter
      InitialValue=list(InitialA=Initial.A,InitialB=Initial.B,InitialC=Initial.C,InitialS=Initial.S)
    }
  }
  if (Model=='1PLAG'){
    Initial.Alpha=0.15                          # giving initial item alpha parameter
    Initial.Beta=as.numeric(Zscore/corr0)       # giving initial item beta parameter
    Initial.Gamma=rep(-1.39,1,J)                # giving initial item gamma parameter
    Initial.Beta[Initial.Beta>2]=1.8      #constrain the max value to 2
    Initial.Beta[Initial.Beta<-2]=-1.8    #constrain the min value to 2
    InitialValue=list(InitialAlpha=Initial.Alpha,InitialBeta=Initial.Beta,InitialGamma=Initial.Gamma)
  }
  if (Model=='1PLG'){
    Initial.Beta=as.numeric(Zscore/corr0)       # giving initial item beta parameter
    Initial.Gamma=rep(-1.39,1,J)                # giving initial item gamma parameter
    Initial.Beta[Initial.Beta>2]=1.8      #constrain the max value to 2
    Initial.Beta[Initial.Beta<-2]=-1.8    #constrain the min value to 2
    InitialValue=list(InitialBeta=Initial.Beta,InitialGamma=Initial.Gamma)
  }

  if (Model=='3PL' | Model =='4PL'){
    ###Checking Variable InitialA###
    if (is.numeric(InitialA) & length(InitialA)==1){  #Whether InitialA=1
      if (InitialA<=0){
        stop('Error: InitialA must bigger than 0!')
      }else{
        InitialValue$A=rep(InitialA,J)
      }
    }else{
      if (is.numeric(InitialA) & length(InitialA)==J){       #Whether InitialA=rep(1,J)
        for (j in 1:J){
          if (is.na(InitialA[j])==FALSE){InitialValue$A[j]=InitialA[j]}
        }
      }else{
        if (length(InitialA)==1){
          if (is.na(InitialA)==F){
            stop('Error: The class of InitialA must be NA or numeric with length= 1 or n.item!')
          }
        }else{
          stop('Error: The class of InitialA must be NA or numeric with length= 1 or n.item!')
        }
      }
    }
    
    ###Checking Variable InitialB###
    if (is.numeric(InitialB) & length(InitialB)==1){  #Whether InitialB=1
      if (InitialB<=0){
        stop('Error: InitialB must bigger than 0!')
      }else{
        InitialValue$B=rep(InitialB,J)
      }
    }else{
      if (is.numeric(InitialB) & length(InitialB)==J){       #Whether InitialB=rep(1,J)
        for (j in 1:J){
          if (is.na(InitialB[j])==FALSE){InitialValue$B[j]=InitialB[j]}
        }
      }else{
        if (length(InitialB)==1){
          if (is.na(InitialB)==F){
            stop('Error: The class of InitialB must be NA or numeric with length= 1 or n.item!')
          }
        }else{
          stop('Error: The class of InitialB must be NA or numeric with length= 1 or n.item!')
        }
      }
    }
    
    ###Checking Variable InitialC###
    if (is.numeric(InitialC) & length(InitialC)==1){  #Whether InitialC=1
      if (InitialC<=0 | InitialC>=0.5){
        stop('Error: InitialC must bigger than 0 and less than 0.5!')
      }else{
        InitialValue$C=rep(InitialC,J)
      }
    }else{
      if (is.numeric(InitialC) & length(InitialC)==J){       #Whether InitialC=rep(1,J)
        for (j in 1:J){
          if (is.na(InitialC[j])==FALSE & InitialC[j]<0.5){
            InitialValue$C[j]=InitialC[j]
          }else{
            stop('Error: InitialC must bigger than 0 and less than 0.5!')
          }
        }
      }else{
        if (length(InitialC)==1){
          if (is.na(InitialC)==F){
            stop('Error: The class of InitialC must be NA or numeric with length= 1 or n.item!')
          }
        }else{
          stop('Error: The class of InitialC must be NA or numeric with length= 1 or n.item!')
        }
      }
    }
    
    
    if (Model =='4PL'){
      ###Checking Variable InitialS###
      if (is.numeric(InitialS) & length(InitialS)==1){  #Whether InitialS=1
        if (InitialS<=0 | InitialS>0.5){
          stop('Error: InitialS must bigger than 0 and less than 0.5!')
        }else{
          InitialValue$S=rep(InitialS,J)
        }
      }else{
        if (is.numeric(InitialS) & length(InitialS)==J){       #Whether InitialS=rep(1,J)
          for (j in 1:J){
            if (is.na(InitialS[j])==FALSE & InitialS[j]<0.5){
              InitialValue$S[j]=InitialS[j]
            }else{
              stop('Error: InitialS must bigger than 0 and less than 0.5!')
            }
          }
        }else{
          if (length(InitialS)==1){
            if (is.na(InitialS)==F){
              stop('Error: The class of InitialS must be NA or numeric with length= 1 or n.item!')
            }
          }else{
            stop('Error: The class of InitialS must be NA or numeric with length= 1 or n.item!')
          }
        }
      }
    }
  }
  
  
  
  if (Model=='1PLAG' | Model=='1PLG'){
    if (Model=='1PLAG'){
    ###Checking Variable InitialAlpha###
      if (is.numeric(InitialAlpha) & length(InitialAlpha)==1){  #Whether InitialAlpha=0.2
        if (InitialAlpha<=0 | InitialAlpha>0.4){
          stop('Error: InitialAlpha must bigger than 0 and less than 0.4!')
        }else{
        InitialValue$Alpha=rep(InitialAlpha,J)
        }
      }else{
        if (length(InitialAlpha)==1){
          if (is.na(InitialAlpha)==F){
            stop('Error: The class of InitialAlpha must be NA or numeric with length= 1!')
          }
        }else{
          stop('Error: The class of InitialAlpha must be NA or numeric with length= 1!')
        }
      }
    }
    
    ###Checking Variable InitialBeta###
    if (is.numeric(InitialBeta) & length(InitialBeta)==1){  #Whether InitialBeta=1
      if (InitialBeta<=0){
        stop('Error: InitialBeta must bigger than 0!')
      }else{
        InitialValue$Beta=rep(InitialBeta,J)
      }
    }else{
      if (is.numeric(InitialBeta) & length(InitialBeta)==J){       #Whether InitialBeta=rep(1,J)
        for (j in 1:J){
          if (is.na(InitialBeta[j])==FALSE){InitialValue$Beta[j]=InitialBeta[j]}
        }
      }else{
        if (length(InitialBeta)==1){
          if (is.na(InitialBeta)==F){
            stop('Error: The class of InitialBeta must be NA or numeric with length= 1 or n.item!')
          }
        }else{
          stop('Error: The class of InitialBeta must be NA or numeric with length= 1 or n.item!')
        }
      }
    }
    
    ###Checking Variable InitialGamma###
    if (is.numeric(InitialGamma) & length(InitialGamma)==1){  #Whether InitialGamma=1
      if (InitialGamma<=0){
        stop('Error: InitialGamma must bigger than 0!')
      }else{
        InitialValue$Gamma=rep(InitialGamma,J)
      }
    }else{
      if (is.numeric(InitialGamma) & length(InitialGamma)==J){       #Whether InitialGamma=rep(1,J)
        for (j in 1:J){
          if (is.na(InitialGamma[j])==FALSE){InitialValue$Gamma[j]=InitialGamma[j]}
        }
      }else{
        if (length(InitialGamma)==1){
          if (is.na(InitialGamma)==F){
            stop('Error: The class of InitialGamma must be NA or numeric with length= 1 or n.item!')
          }
        }else{
          stop('Error: The class of InitialGamma must be NA or numeric with length= 1 or n.item!')
        }
      }
    }
  }
  
  
  
  ###Checking Variable Tol###
  if (is.numeric(Tol)){
    if(length(Tol)==1){
      if (Tol<=0){
        stop('Error: The min of Tol must bigger than 0!')
      }
    }else{
      stop('Error: The length of Tol must be 1!')
    }
  }else{
    stop('Error: The type of Tol must be numeric!')
  }
  
  
  ###Checking Variable max.ECycle###
  if (length(max.ECycle)==1){
    if (is.numeric(length(max.ECycle)) | is.integer(length(max.ECycle))){
      max.ECycle=as.integer(max.ECycle)
      if (max.ECycle<1){
        stop('Error: The min of max.ECycle is 1!')
      }
    }else{
      stop('Error: The type of max.ECycle must be integer!')
    }
  }else{
    stop('Error: The length of max.ECycle must be 1!')
  }
  
  ###Checking Variable max.MCycle###
  if (length(max.MCycle)==1){
    if (is.numeric(length(max.MCycle)) | is.integer(length(max.MCycle))){
      max.MCycle=as.integer(max.MCycle)
      if (max.MCycle<1){
        stop('Error: The min of max.MCycle is 1!')
      }
    }else{
      stop('Error: The type of max.MCycle must be integer!')
    }
  }else{
    stop('Error: The length of max.MCycle must be 1!')
  }
  
  
  ###Checking Variable n.Quadpts###
  if (length(n.Quadpts)==1){
    if (is.numeric(length(n.Quadpts)) | is.integer(length(n.Quadpts))){
      n.Quadpts=as.integer(n.Quadpts)
      if (n.Quadpts<5){
        stop('Error: The min of n.Quadpts is 5!')
      }
      if (n.Quadpts>256){
        warning('Too large n.Quadpts will cause produre to be extremely time-consuming!')
      }
    }else{
      stop('Error: The type of n.Quadpts must be integer!')
    }
  }else{
    stop('Error: The length of n.Quadpts must be 1!')
  }
  
  ###Checking Variable n.decimal###
  if (length(n.decimal)==1){
    if (is.numeric(length(n.decimal)) | is.integer(length(n.decimal))){
      n.decimal=as.integer(n.decimal)
      if (n.decimal<0){
        stop('Error: The min of n.decimal is 0!')
      }
      if (n.decimal>8){
        stop('Error: The max of n.decimal is 8!')
      }
    }else{
      stop('Error: The type of n.decimal must be integer!')
    }
  }else{
    stop('Error: The length of n.decimal must be 1!')
  }
  
  Integer.tot=list(n.Quadpts=n.Quadpts, max.ECycle=max.ECycle, max.MCycle=max.MCycle, n.decimal=n.decimal)
  
  
  ###Checking Variable Theta.lim###
  if (is.numeric(Theta.lim)){
    if (length(Theta.lim)==2){
      if (Theta.lim[1]>=Theta.lim[2]){
        stop('Error: Theta.lim[1] must bigger than Theta.lim[2]!')
      }
    }else{
      stop('Error: The length of Theta.lim must be 2!')
    }
  }else{
    stop('Error: The type of Theta.lim must be numeric!')
  }
  
  ###Checking Variable Missing###
  if (length(Missing)==1){
    if (is.numeric(Missing) | is.na(Missing)){
      if (Missing==1 | Missing==0){
        stop('Error: The value of Missing cannot be 1 or 0!')
      }
    }else{
      stop('Error: The type of Missing must be numeric or NA!')
    }
  }else{
    stop('Error: The length of Missing must be 1!')
  }
  
  
  ###Checking Variable ParConstraint###
  if (length(ParConstraint)==1){
    if (is.logical(ParConstraint)==F){
      stop('Error: The type of ParConstraint must be logical!')
    }else{
      if (ParConstraint){ParConstraint=list(ParConstraint=1L)}else{ParConstraint=list(ParConstraint=0L)}
    }
  }else{
    stop('Error: The length of ParConstraint must be 1!')
  }  
  
  
  return(c(data.list, Prior.list, InitialValue, Integer.tot, ParConstraint))
}


