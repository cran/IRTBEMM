#include <math.h>
#include <string.h>
#include <stdio.h>
//BEMM function for 1PLG
void BEMM1PLG(double *data, int *CountNum, int *n_class, int *n_item, double *LH,
			  double *IBeta, double *IGamma, double *TBeta, double *TGamma, 
			  double *deltahat_Beta, double *deltahat_Gamma,
			  double *Beta, double *Gamma, double *Tol, double *cr, int *E_exit,
			  double *SEBeta, double *SEGamma, int *ParConstraint,   
			  double *f, double *r, double *fz, double *rz, double *P, double *Pstar, double *AG, 
			  double *LL, double *LL0, double *Posterior_prob, 
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
	double Ibb=0;
	double Igg=0;
	double lbb = 0;
	double lgg = 0;
	double EZ;
	double L;	
	double lb1;
	double lg1;
	double pstar;
	double ag;

	double bt0;
	double gt0;
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
				AG[v] = 1 / (1 + exp(-(Gamma[j])));
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
					ag = 1 / (1 + exp(-(gt0)));
					lg1 += r[k + nq * j] - rz[k + nq * j] - (f[k] - fz[k + nq * j]) * ag; //lg1: r - rz - (f - fz)*ps
					lgg += -(f[k] - fz[k + nq * j]) * ag * (1 - ag);					  //lgg: -(f-fz)*(ag*(1-ag))
					lb1 += -(rz[k + nq * j] - f[k] * pstar);       					      //lb1: - (rz - f * ps)
					lbb += -(f[k] * pstar * (1 - pstar)); 								  //lbb: - (f * ps * (1-ps))
				}
				// Maximize beta
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
				// Maximize gamma
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
					AG[v] = 1 / (1 + exp(-(Gamma[j])));
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
	double delta[2]={0};
	double delta0[2]={0};
	double delta1[2]={0};
	double cr_SEM1=1;
	double cr_SEM2=1;
	if (n_ECycle[0]>=10){
		start_SEM=floor(n_ECycle[0] * 0.2);
		end_SEM=floor(n_ECycle[0] * 0.8);
	}
	
	// Estimating SEs of beta & gamma
	for (jj = 0; jj < J; jj++) {
		z=start_SEM;
		SEM_exit=1;
		cr_SEM1=1;
		cr_SEM2=1;
		while (SEM_exit && z<=end_SEM){
			mm = jj + J * z;
			for (ParClass = 1; ParClass <= 2; ParClass++){
				if (ParClass==1 && z>=2 && cr_SEM1<0.01){
					continue;
				}
				if (ParClass==2 && z>=2 && cr_SEM2<0.01){
					continue;
				}
				for (j = 0; j < J; j++){
					deltahat_Beta[j]=Beta[j];
                    deltahat_Gamma[j]=Gamma[j];
				}
				switch (ParClass){
					case 1:
                        deltahat_Beta[jj]=TBeta[mm];
						break;
					case 2:
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
							AG[v] = 1 / (1 + exp(-(deltahat_Gamma[j])));
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
					if (ParClass==1){
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
					if (ParClass==2){
						for (k = 0; k < nq; k++) {
							pstar = 1 / (1 + exp(-(node_Quadpts[k] - bt0)));
							ag = 1 / (1 + exp(-(gt0)));
							lg1 += r[k + nq * jj] - rz[k + nq * jj] - (f[k] - fz[k + nq * jj]) * ag; //lg1: r - rz - (f - fz)*ps
							lgg += -(f[k] - fz[k + nq * jj]) * ag * (1 - ag);					     //lgg: -(f-fz)*(ag*(1-ag))
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
				if (ParClass==1){
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
				if (ParClass==2){
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
					case 1:
                        delta1[0]=(deltahat_Beta[jj]-Beta[jj])/(TBeta[mm]-Beta[jj]);
						break;
					case 2:
                        delta1[1]=(deltahat_Gamma[jj]-Gamma[jj])/(TGamma[mm]-Gamma[jj]);   
						break;	
				}
			}	
			cr_SEM1=fabs(delta1[0]-delta0[0]);
			cr_SEM2=fabs(delta1[1]-delta0[1]);
			if (cr_SEM1<0.01 && cr_SEM2<0.01 && z>=2){
				SEM_exit=0;
			}else{
				z=z+1;
			}
			for (m = 0; m<=1; m++){
				if (isnormal(delta1[m])){
					delta0[m] = delta1[m];
				}
			}
		}
		delta[0]=1-delta0[0];
		delta[1]=1-delta0[1];
		delta1[0] =  1 / delta[0];
		delta1[1] =  1 / delta[1];
		if (isnormal(delta1[0])==0){delta1[0] = 1;}
		if (isnormal(delta1[1])==0){delta1[1] = 1;}
		SEBeta[jj]= sqrt(fabs(IBeta[jj] * delta1[0]));
		SEGamma[jj]= sqrt(fabs(IGamma[jj] * delta1[1]));
		
	}
}
