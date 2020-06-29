#include <math.h>
#include <string.h>
#include <stdio.h>
//BEMM function for 1PLAG
void BEMM1PLAG(double *data, int *CountNum, int *n_class, int *n_item, double *LH,
			  double *IAlpha, double *IBeta, double *IGamma, 
			  double *TAlpha, double *TBeta, double *TGamma, 
			  double *deltahat_Alpha, double *deltahat_Beta, double *deltahat_Gamma,
			  double *Alpha, double *Beta, double *Gamma, double *Tol, double *cr, int *E_exit,
			  double *SEAlpha, double *SEBeta, double *SEGamma, int *ParConstraint,   
			  double *f, double *r, double *fz, double *rz, double *P, double *Pstar, double *AG, 
			  double *LL, double *LL0, double *Posterior_prob, const double *PriorAlpha, 
			  const double *PriorBeta, const double *PriorGamma, const int *n_Quadpts, 
			  const double *node_Quadpts, const double *weight_Quadpts,		  
			  const int *MECycle, const int *MMCycle, int *n_ECycle){
	
	int j;
	int i;
	int k;
	int v;
	int m;
	int n;
    int n_MCycle=1;
	int M_exit=1;
	int I = n_class[0];
	int J = n_item[0];
	int max_ECycle = MECycle[0];
	int max_MCycle = MMCycle[0];
	int nq = n_Quadpts[0];
	int Jnq = nq * J;
	
	
	double TOL = Tol[0];
	double LH0 = 0;
	double Iaa=0;
	double Ibb=0;
	double Igg=0;
	double laa = 0;
	double lbb = 0;
	double lgg = 0;
	double EZ;
	double L;	
	double la1;
	double lb1;
	double lg1;
	double pstar;
	double ag;

	double at0;
	double bt0;
	double gt0;
	double at1;
	double bt1;
	double gt1;

	
	memset(f, 0, sizeof(double) * nq);
	memset(r, 0, sizeof(double) * Jnq);
	memset(rz, 0, sizeof(double) * Jnq);
	memset(fz, 0, sizeof(double) * Jnq);
	memset(LL0, 0, sizeof(double) * I);
	
	
	//Calculate Log-likelihood
	for (k = 0; k < nq; k++) {
		for (i = 0; i<I; i++){
			L = 1.0;
			for (j = 0; j < J; j++) { 
				n = i + I * j;
				v = k + nq * j;
				Pstar[v] = 1 / (1 + exp(-(node_Quadpts[k] - Beta[j])));
				AG[v] = 1 / (1 + exp(-(Alpha[0] * node_Quadpts[k] + Gamma[j])));
				if (Pstar[v] >=1){Pstar[v] = 0.9999;}
				if (Pstar[v] <=0){Pstar[v] = 0.0001;}
				if (AG[v] >=1){AG[v] = 0.9999;}
				if (AG[v] <=0){AG[v] = 0.0001;}
				P[v] = Pstar[v] + (1-Pstar[v]) * AG[v];
				L *= P[v] * data[n] + (1 - P[v]) * (1 - data[n]);
			}
			v = i + I * k;
			LL[v] = L * weight_Quadpts[k];
			LL0[i] += LL[v];
		}
	}
	//Calculate artificial data f r fz rz
	for (i = 0; i<I; i++){
		LH0 += CountNum[i] * log(LL0[i]);
		for (k = 0; k < nq; k++) {
			v = i + I * k;
			Posterior_prob[v] = LL[v] / LL0[i];
			f[k] += Posterior_prob[v] * CountNum[i];
			for (j = 0; j < J; j++) { 
				m = k + nq * j;
				n = i + I * j;
				r[m] += Posterior_prob[v] * data[n] * CountNum[i];
				EZ = Pstar[m] / P[m] * data[n];
				fz[m] += Posterior_prob[v] * EZ * CountNum[i];
				rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
			}
		}
	}
	//E step iteration
	while (E_exit[0] && (n_ECycle[0] < max_ECycle)) {
		for (j = 0; j < J; j++) {
			gt0 = Gamma[j];
			bt0 = Beta[j];
			n_MCycle = 1;
			M_exit = 1;
			//M step iteration for beta & gamma respectively
			while (M_exit && (n_MCycle <= max_MCycle)) {
				for (k = 0; k < nq; k++) {
					pstar = 1 / (1 + exp(-(node_Quadpts[k] - bt0)));
					ag = 1 / (1 + exp(-(Alpha[0] * node_Quadpts[k] + gt0)));
					lg1 += r[k + nq * j] - rz[k + nq * j] - (f[k] - fz[k + nq * j]) * ag; //lg1: r - rz - (f - fz)*ps
					lgg += -(f[k] - fz[k + nq * j]) * ag * (1 - ag);					  //lgg: -(f-fz)*(ag*(1-ag))
					lb1 += -(rz[k + nq * j] - f[k] * pstar);       					      //lb1: - (rz - f * ps)
					lbb += -(f[k] * pstar * (1 - pstar)); 								  //lbb: - (f * ps*(1-ps))
				}
				// Maximize beta & gamma
				if (PriorBeta[j]!=-9 && PriorBeta[j + J]!=-9) {
					lb1 += -((bt0-PriorBeta[j])/PriorBeta[j + J]);
					lbb += -1/PriorBeta[j + J];
				}
				Ibb = - 1 / lbb;
				bt1 = bt0 + (Ibb * lb1);
				lb1 = 0;
				lbb = 0;
				if (fabs(bt1-bt0) >= 0.01){
					bt0 = bt1;
				}
				if (PriorGamma[j]!=-9 && PriorGamma[j + J]!=-9) {
					lg1 += -((gt0-PriorGamma[j])/PriorGamma[j + J]);
					lgg += -1/PriorGamma[j + J];
				}
				Igg = - 1 / lgg;
				gt1 = gt0 + (Igg * lg1);
				lg1 = 0;
				lgg = 0;
				if (fabs(gt1-gt0) >= 0.01){
					gt0 = gt1;
				}
				if (fabs(bt1-bt0) < 0.01 && fabs(gt1-gt0) < 0.01){
					bt0 = bt1;
					gt0 = gt1;
					M_exit = 0;
				} else{
					bt0 = bt1;
					gt0 = gt1;
					n_MCycle += 1;
				}
			}
			m = j + J * n_ECycle[0];
			if (isnormal(bt0) && isnormal(gt0)){
				if (ParConstraint[0]){
					if (bt0>=-6 && bt0<=6){
						Beta[j] = bt0;
						TBeta[m] = bt0;
						IBeta[j] = Ibb;
					}
					if (gt0>=-7 && gt0<=0){
						Gamma[j] = gt0;
						TGamma[m] = gt0;
						IGamma[j] = Igg;
					}
				}else{
					Beta[j] = bt0;
					Gamma[j] = gt0;
					TBeta[m] = bt0;
					TGamma[m] = gt0;
					IBeta[j] = Ibb;
					IGamma[j] = Igg;
				}
			}else{
				if (n_ECycle[0]!=0){
					n=j + J * (n_ECycle[0]-1);
					TBeta[m] = TBeta[n];
					TGamma[m] = TGamma[n];
				}else{
					TBeta[m] = Beta[j];
					TGamma[m] = Gamma[j];
				}
			}
		}
			
		//M step iteration for alpha
		at0 = Alpha[0];
		n_MCycle = 1;
		M_exit = 1;
		while (M_exit && (n_MCycle <= max_MCycle)) {
			for (k = 0; k < nq; k++) {
				for (j = 0; j < J; j++){
					ag = 1 / (1 + exp(-(at0 * node_Quadpts[k] + Gamma[j])));
					la1 += at0 * node_Quadpts[k] * (r[k + nq * j] - rz[k + nq * j] - (f[k] - fz[k + nq * j]) * ag); 
																					  //la1: exp(log(at))*X*(r - rz - (f - fz)*ps)
					laa += - at0 * at0 * node_Quadpts[k] * node_Quadpts[k] * (f[k] - fz[k + nq * j]) * ag * (1 - ag);					  
				                                                                      //laa: -exp(2*log(at))*X*X*(f-fz)*(ag*(1-ag))
				}
			}
			// Maximize alpha
			if (PriorAlpha[0]!=-9 && PriorAlpha[1]!=-9) {
				la1 += -((log(at0)-PriorAlpha[0])/PriorAlpha[1]);
				laa += -1/PriorAlpha[1];
			}
			Iaa = - 1 / laa;
			at1 = log(at0) + (Iaa * la1);
			la1 = 0;
			laa = 0;
			if (fabs(at1-log(at0)) < 0.01){
				at0 = exp(at1);
				M_exit = 0;
			} else{
				at0 = exp(at1);
				n_MCycle += 1;
			}
		}
		if (isnormal(at0)){
			if (at0>=0.001){
				Alpha[0] = at0;
				TAlpha[n_ECycle[0]] = at0;
				IAlpha[0] = Iaa;
			}
		}else{
			if (n_ECycle[0]!=0){
				TAlpha[n_ECycle[0]] = TAlpha[n_ECycle[0]-1];
			}else{
				TAlpha[n_ECycle[0]] = Alpha[0];
			}
		}
		memset(f, 0, sizeof(double) * nq);
		memset(r, 0, sizeof(double) * Jnq);
		memset(rz, 0, sizeof(double) * Jnq);
		memset(fz, 0, sizeof(double) * Jnq);
		memset(LL0, 0, sizeof(double) * I);
		// Update Log-likelihood
		for (k = 0; k < nq; k++) {
			for (i = 0; i<I; i++){
				L = 1.0;
				for (j = 0; j < J; j++) { 
					n = i + I * j;
					v = k + nq * j;
					Pstar[v] = 1 / (1 + exp(-(node_Quadpts[k] - Beta[j])));
					AG[v] = 1 / (1 + exp(-(Alpha[0] * node_Quadpts[k] + Gamma[j])));
					if (Pstar[v] >=1){Pstar[v] = 0.9999;}
					if (Pstar[v] <=0){Pstar[v] = 0.0001;}
					if (AG[v] >=1){AG[v] = 0.9999;}
					if (AG[v] <=0){AG[v] = 0.0001;}
					P[v] = Pstar[v] + (1-Pstar[v]) * AG[v];
					L *= P[v] * data[n] + (1 - P[v]) * (1 - data[n]);
				}
				v = i + I * k;
				LL[v] = L * weight_Quadpts[k];
				LL0[i] += LL[v];
			}
		}
		// Update artificial data
		for (i = 0; i<I; i++){
			LH[n_ECycle[0]] += CountNum[i] * log(LL0[i]);
			for (k = 0; k < nq; k++) {
				v = i + I * k;
				Posterior_prob[v] = LL[v] / LL0[i];
				f[k] += Posterior_prob[v] * CountNum[i];
				for (j = 0; j < J; j++) { 
					m = k + nq * j;
					n = i + I * j;
					r[m] += Posterior_prob[v] * data[n] * CountNum[i];
					EZ = Pstar[m] / P[m] * data[n];
					fz[m] += Posterior_prob[v] * EZ * CountNum[i];
					rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
				}
			}
		}
		cr[0]=LH[n_ECycle[0]]-LH0;
		LH0=LH[n_ECycle[0]];
		if (fabs(cr[0])<TOL){
			n_ECycle[0] = n_ECycle[0] + 1;
			E_exit[0]=0;
		}else{
			n_ECycle[0] = n_ECycle[0] + 1;
		}
	}
	n_ECycle[0] = n_ECycle[0] - 1;
	
	
	// Estimating Standard Errors via Supplement EM algorithm
	int z;
	int mm;
	int jj;
	int ParClass;
	int SEM_exit;
	int start_SEM=0;
	int end_SEM=n_ECycle[0];
	double delta[3]={0};
	double delta0[3]={0};
	double delta1[3]={0};
	int cr_SEM=n_ECycle[0];
	double cr_SEM0=1;
	double cr_SEM1=1;
	double cr_SEM2=1;
	double cr_SEM3=1;
	double deltatemp;
	
	for (i = 0; i<=cr_SEM; i++){
		deltatemp=exp(-(LH[i+1]-LH[i]));
		if (deltatemp>=0.9 && deltatemp<=0.999){
			if (cr_SEM0==0){
				end_SEM=i;
			}else{
				start_SEM=i;
				cr_SEM0=0;
			}
		}
	}
	// Estimating SEs of alpha
	z=start_SEM;
	SEM_exit=1;
	cr_SEM1=1;
	while (SEM_exit && z<=end_SEM){
		deltahat_Alpha[0]=TAlpha[z];
		memset(f, 0, sizeof(double) * nq);
		memset(r, 0, sizeof(double) * Jnq);
		memset(rz, 0, sizeof(double) * Jnq);
		memset(fz, 0, sizeof(double) * Jnq);
		memset(LL0, 0, sizeof(double) * I);
		// Update Log-likelihood
		for (k = 0; k < nq; k++) {
			for (i = 0; i<I; i++){
				L = 1.0;
				for (j = 0; j < J; j++) { 
					n = i + I * j;
					v = k + nq * j;
					Pstar[v] = 1 / (1 + exp(-(node_Quadpts[k] - Beta[j])));
					AG[v] = 1 / (1 + exp(-(deltahat_Alpha[0] * node_Quadpts[k] + Gamma[j])));
					if (Pstar[v] >=1){Pstar[v] = 0.9999;}
					if (Pstar[v] <=0){Pstar[v] = 0.0001;}
					if (AG[v] >=1){AG[v] = 0.9999;}
					if (AG[v] <=0){AG[v] = 0.0001;}
					P[v] = Pstar[v] + (1-Pstar[v]) * AG[v];
					L *= P[v] * data[n] + (1 - P[v]) * (1 - data[n]);
				}
				v = i + I * k;
				LL[v] = L * weight_Quadpts[k];
				LL0[i] += LL[v];
			}
		}
		// Update artificial data
    	for (i = 0; i<I; i++){
			for (k = 0; k < nq; k++) {
				v = i + I * k;
				Posterior_prob[v] = LL[v] / LL0[i];
				f[k] += Posterior_prob[v] * CountNum[i];
				for (j = 0; j < J; j++) { 
					m = k + nq * j;
					n = i + I * j;
					r[m] += Posterior_prob[v] * data[n] * CountNum[i];
					EZ = Pstar[m] / P[m] * data[n];
					fz[m] += Posterior_prob[v] * EZ * CountNum[i];
					rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
				}
			}
		}
		//Estimating alpha	
		at0=deltahat_Alpha[0];
		while (M_exit && (n_MCycle <= max_MCycle)) {
			for (k = 0; k < nq; k++) {
				for (j = 0; j < J; j++){
					ag = 1 / (1 + exp(-(at0 * node_Quadpts[k] + Gamma[j])));
					la1 += at0 * node_Quadpts[k] * (r[k + nq * j] - rz[k + nq * j] - (f[k] - fz[k + nq * j]) * ag); 
																					  //la1: exp(log(a))*X*(r - rz - (f - fz)*ps)
					laa += - at0 * at0 * node_Quadpts[k] * node_Quadpts[k] * (f[k] - fz[k + nq * j]) * ag * (1 - ag);					  
																						 //laa: -exp(2log(a))*X*X*(f-fz)*(ag*(1-ag))
				}
			}
			// Maximize alpha
			if (PriorAlpha[0]!=-9 && PriorAlpha[1]!=-9) {
				la1 += -((log(at0)-PriorAlpha[0])/PriorAlpha[1]);
				laa += -1/PriorAlpha[1];
			}
			Iaa = - 1 / laa;
			at1 = log(at0) + (Iaa * la1);
			la1 = 0;
			laa = 0;
			if (fabs(at1-log(at0)) < 0.01){
				at0 = exp(at1);
				M_exit = 0;
			} else{
				at0 = exp(at1);
				n_MCycle += 1;
			}
		}
		if (isnormal(at0)){
			if (at0>=0.001){
				deltahat_Alpha[0] = at0;
			}
		}
		delta1[0]=(log(deltahat_Alpha[0])-log(Alpha[0]))/(log(TAlpha[z])-log(Alpha[0])+0.0001);
		cr_SEM1=fabs(delta1[0]-delta0[0]);	
		if (cr_SEM1<0.0001 && z>=2){
			SEM_exit=0;
		}else{
			z=z+1;
		}
		if (isnormal(delta1[0])){
			delta0[0] = delta1[0];
		}
	}
	delta[0]=1-delta0[0];
	delta1[0] =  1 / delta[0];
	if (isnormal(delta1[0])==0 || delta1[0]<=0){delta1[0] = 1;}
	SEAlpha[0]= sqrt(Alpha[0] * Alpha[0] * IAlpha[0] * delta1[0]);
	if (SEAlpha[0]>1){SEAlpha[0]= sqrt(IAlpha[0]);}
	
	// Estimating SEs of beta & gamma
	for (jj = 0; jj < J; jj++) {
		z=start_SEM;
		SEM_exit=1;
		cr_SEM2=1;
		cr_SEM3=1;
		while (SEM_exit && z<=end_SEM){
			mm = jj + J * z;
			for (ParClass = 2; ParClass <= 3; ParClass++){
				if (ParClass==2 && z>=2 && cr_SEM1<0.0001){
					continue;
				}
				if (ParClass==3 && z>=2 && cr_SEM2<0.0001){
					continue;
				}
				deltahat_Alpha[0]=Alpha[0];
				for (j = 0; j < J; j++){
					deltahat_Beta[j]=Beta[j];
                    deltahat_Gamma[j]=Gamma[j];
				}
				switch (ParClass){
					case 2:
                        deltahat_Beta[jj]=TBeta[mm];
						break;
					case 3:
                        deltahat_Gamma[jj]=TGamma[mm];
						break;	
				}
				memset(f, 0, sizeof(double) * nq);
				memset(r, 0, sizeof(double) * Jnq);
				memset(rz, 0, sizeof(double) * Jnq);
				memset(fz, 0, sizeof(double) * Jnq);
				memset(LL0, 0, sizeof(double) * I);
				// Update Log-likelihood
				for (k = 0; k < nq; k++) {
					for (i = 0; i<I; i++){
						L = 1.0;
						for (j = 0; j < J; j++) { 
							n = i + I * j;
							v = k + nq * j;
							Pstar[v] = 1 / (1 + exp(-(node_Quadpts[k] - deltahat_Beta[j])));
							AG[v] = 1 / (1 + exp(-(deltahat_Alpha[0] * node_Quadpts[k] + deltahat_Gamma[j])));
							if (Pstar[v] >=1){Pstar[v] = 0.9999;}
							if (Pstar[v] <=0){Pstar[v] = 0.0001;}
							if (AG[v] >=1){AG[v] = 0.9999;}
							if (AG[v] <=0){AG[v] = 0.0001;}
							P[v] = Pstar[v] + (1-Pstar[v]) * AG[v];
							L *= P[v] * data[n] + (1 - P[v]) * (1 - data[n]);
						}
						v = i + I * k;
						LL[v] = L * weight_Quadpts[k];
						LL0[i] += LL[v];
					}
				}
				// Update artificial data
				for (i = 0; i<I; i++){
					for (k = 0; k < nq; k++) {
						v = i + I * k;
						m = k + nq * jj;
						n = i + I * jj;
						Posterior_prob[v] = LL[v] / LL0[i];
						f[k] += Posterior_prob[v] * CountNum[i];
						r[m] += Posterior_prob[v] * data[n] * CountNum[i];
						EZ = Pstar[m] / P[m] * data[n];
						fz[m] += Posterior_prob[v] * EZ * CountNum[i];
						rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
					}
				}		
				//estimate beta & gamma parameters
				gt0 = deltahat_Gamma[jj];
				bt0 = deltahat_Beta[jj];
				n_MCycle = 1;
				M_exit = 1;
				//M step iteration for beta & gamma respectively
				while (M_exit && (n_MCycle <= max_MCycle)) {
					if (ParClass==2){
						for (k = 0; k < nq; k++) {
							pstar = 1 / (1 + exp(-(node_Quadpts[k] - bt0)));
							lb1 += -(rz[k + nq * jj] - f[k] * pstar);       					      //lb1: - (rz - f * ps)
							lbb += -(f[k] * pstar * (1 - pstar)); 									  //lbb: - (f * ps)
						}
						// Maximize beta
						if (PriorBeta[jj]!=-9 && PriorBeta[jj + J]!=-9) {
							lb1 += -((bt0-PriorBeta[jj])/PriorBeta[jj + J]);
							lbb += -1/PriorBeta[jj + J];
						}
						Ibb = - 1 / lbb;
						bt1 = bt0 + (Ibb * lb1);
						lb1 = 0;
						lbb = 0;
						if (fabs(bt1-bt0) < 0.01){
							bt0 = bt1;
							M_exit = 0;
						} else{
							bt0 = bt1;
							n_MCycle += 1;
						}
					}
					if (ParClass==3){
						for (k = 0; k < nq; k++) {
							pstar = 1 / (1 + exp(-(node_Quadpts[k] - bt0)));
							ag = 1 / (1 + exp(-(deltahat_Alpha[0] * node_Quadpts[k] + gt0)));
							lg1 += r[k + nq * jj] - rz[k + nq * jj] - (f[k] - fz[k + nq * jj]) * ag; //lg1: r - rz - (f - fz)*ps
							lgg += -(f[k] - fz[k + nq * jj]) * ag * (1 - ag);				         //lgg: -(f-fz)*(ag*(1-ag))
						}
						// Maximize gamma
						if (PriorGamma[jj]!=-9 && PriorGamma[jj + J]!=-9) {
							lg1 += -((gt0-PriorGamma[jj])/PriorGamma[jj + J]);
							lgg += -1/PriorGamma[jj + J];
						}
						Igg = - 1 / lgg;
						gt1 = gt0 + (Igg * lg1);
						lg1 = 0;
						lgg = 0;
						if (fabs(gt1-gt0) < 0.01){
							gt0 = gt1;
							M_exit = 0;
						} else{
							gt0 = gt1;
							n_MCycle += 1;
						}
					}
				}
				if (ParClass==2){
					if (isnormal(bt0)){
						if (ParConstraint[0]){
							if (bt0>=-6 && bt0<=6){
								deltahat_Beta[jj] = bt0;
							}
						}else{
							deltahat_Beta[jj] = bt0;
						}
					}
				}
				if (ParClass==3){
					if (isnormal(gt0)){
						if (ParConstraint[0]){
							if (gt0>=-6 && gt0<=6){
								deltahat_Gamma[jj] = gt0;
							}
						}else{
							deltahat_Gamma[jj] = gt0;
						}
					}
				}
				switch (ParClass){
					case 2:
                        delta1[1]=(deltahat_Beta[jj]-Beta[jj])/(TBeta[mm]-Beta[jj]+0.0001);
						break;
					case 3:
                        delta1[2]=(deltahat_Gamma[jj]-Gamma[jj])/(TGamma[mm]-Gamma[jj]+0.0001);   
						break;	
				}
			}	
			cr_SEM2=fabs(delta1[1]-delta0[1]);
			cr_SEM3=fabs(delta1[2]-delta0[2]);
			if (cr_SEM2<0.0001 && cr_SEM3<0.0001 && z>=2){
				SEM_exit=0;
			}else{
				z=z+1;
			}
			for (m = 1; m<=2; m++){
				if (isnormal(delta1[m])){
					delta0[m] = delta1[m];
				}
			}
		}
		delta[1]=1-delta0[1];
		delta[2]=1-delta0[2];
		delta1[1] =  1 / delta[1];
		delta1[2] =  1 / delta[2];
		if (isnormal(delta1[1])==0 || delta1[1]<=0){delta1[1] = 1;}
		if (isnormal(delta1[2])==0 || delta1[2]<=0){delta1[2] = 1;}
		SEBeta[jj]= sqrt(IBeta[jj] * delta1[1]);
		SEGamma[jj]= sqrt(IGamma[jj] * delta1[2]);
		if (SEBeta[jj]>1){SEBeta[jj]= sqrt(IBeta[jj]);}
		if (SEGamma[jj]>1){SEGamma[jj]= sqrt(IGamma[jj]);}	
	}
}
