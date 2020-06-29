#include <math.h>
#include <string.h>
#include <stdio.h>
//BEMM function for 3PLM
void BEMM3PL(double *data, int *CountNum, int *n_class, int *n_item, double *LH, 
			  double *IA, double *IB, double *IC, double *IAB, 
			  double *TA, double *TB, double *TC, 
			  double *deltahat_A, double *deltahat_B, double *deltahat_C, 
			  double *A, double *B, double *C, double *Tol, double *cr, int *E_exit,
			  double *SEA, double *SEB, double *SEC, int *ParConstraint,   
			  double *f, double *r, double *fz, double *rz, double *P, double *Pstar, 
			  double *LL, double *LL0, double *Posterior_prob, const double *PriorA, 
			  const double *PriorB, const double *PriorC, const int *n_Quadpts, 
			  const double *node_Quadpts, const double *weight_Quadpts,		  
			  const int *MECycle, const int *MMCycle, int *n_ECycle, const double *ConstD){
	
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
	
	
	double D = ConstD[0];
	double LH0 = 0;
	double TOL = Tol[0];
	double Iaa=0;
	double Ibb=0;
	double Icc=0;
	double Iab=0;
	double c_part1=0;
	double c_part2=0;
	double c_part3=0;
	double c_part4;
	double c_part5;
	double c_part6;
	double laa;
	double lbb;
	double lcc;
	double EZ;
	double L;	
	double Da;
	double D2a2;
	double Weight;
	double la1;
	double lb1;
	double lab;
	double x_bt;
	double x_bt2;
	double pstar;
	double psf;
	double wsf;
	double fz_psf;
	double la1_core=0;
	double laa_core=0;
	double lb1_core=0;
	double lbb_core=0;
	double lab_core=0;
	double at0;
	double bt0;
	double ct0;
	double at1;
	double bt1;

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
				Pstar[v] = 1 / (1 + exp(-D * A[j] * (node_Quadpts[k] - B[j])));
				if (Pstar[v] >=1){Pstar[v] = 0.9999;}
				if (Pstar[v] <=0){Pstar[v] = 0.0001;}
				P[v] = C[j] + (1-C[j]) * Pstar[v];
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
				EZ = (Pstar[m] / P[m]) * data[n];
				fz[m] += Posterior_prob[v] * EZ * CountNum[i];
				rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
			}
		}
	}
	//E step iteration
	while (E_exit[0] && (n_ECycle[0] < max_ECycle)) {
		for (j = 0; j < J; j++) {
			//estimate c parameters & s parameters
			for (k = 0; k < nq; k++) {
				m = k + nq * j;
				c_part1 += r[m] - rz[m];
				c_part2 += f[k] - fz[m] -r[m] + rz[m];
				c_part3 += f[k] - fz[m];
			}
			c_part4 = 1.0 - C[j];
			c_part5 = c_part4 * c_part4;
			c_part6 = C[j] * C[j];
			//estimate c parameters
			if (PriorC[j]==-9 || PriorC[j + J]==-9) {
				lcc = - (c_part1 / c_part6) - (c_part2) / c_part5;	//lcc=-(r-rz)/c^2-(f-fz-r+rz)/(1-c)^2
				ct0 = c_part1 / c_part3;							//c=(r-rz)/(f-fz)
			} else {
				lcc = - ((c_part1 + PriorC[j] -1) / c_part6) - ((c_part2 + PriorC[j + J] - 1.0) / c_part5);
				ct0 = (PriorC[j] -1 + c_part1) / (((PriorC[j] + PriorC[j + J]) - 2.0) + c_part3);
			}
			//estimate a and b parameter
			at0 = A[j];
			bt0 = B[j];
			n_MCycle = 1;
			c_part1 = 0;
			c_part2 = 0;
			c_part3 = 0;
			M_exit = 1;
			//M step iteration
			while (M_exit && (n_MCycle <= max_MCycle)) {
				Da = D * at0;
				D2a2 = Da * Da;
				for (k = 0; k < nq; k++) {
					x_bt = node_Quadpts[k] - bt0;   	//x-bt
					x_bt2 = x_bt * x_bt;  				//(x-bt)^2
					pstar = 1 / (1 + exp(-Da * x_bt));	//ps
					psf = pstar * f[k];					//f * ps
					wsf = psf * (1-pstar);    			//f * ws
					fz_psf = fz[k + nq * j] - psf;  	//fz - f *ps
					la1_core += fz_psf * x_bt;			//la1: (fz - f*ps)* (x-bt)
					lb1_core += fz_psf;       			//lb1: fz - f* ps
					laa_core += wsf * x_bt2;			//laa: (f*ws)* (x-bt)^2 
					lbb_core += wsf;          			//lbb: (f*ws)
					lab_core += wsf * x_bt;   			//lab: wsf*(x-bt)
				}
				la1 = Da * la1_core;
				lb1 = -Da * lb1_core;
				laa = -D2a2 * laa_core;
				lbb = -D2a2 * lbb_core;
				lab = D2a2 * lab_core;
				// Maximize a and b
				if (PriorA[j]!=-9 && PriorA[j + J]!=-9) {
					la1 += -((log(at0)-PriorA[j])/PriorA[j + J]);
					laa += -1/PriorA[j + J];
				}
				if (PriorB[j]!=-9 && PriorB[j + J]!=-9) {
					lb1 += -((bt0-PriorB[j])/PriorB[j + J]);
					lbb += -1/PriorB[j + J];
				}
				Weight = laa * lbb - lab * lab;
				Iaa = - lbb / Weight;
				Ibb = - laa / Weight;
				Iab = lab / Weight;
				Icc = - 1 / lcc;
				at1 = log(at0) + (Iaa * la1 + Iab * lb1);
				bt1 = bt0 + (Ibb * lb1 + Iab * la1);
				la1_core = 0;
				lb1_core = 0;
				laa_core = 0;
				lbb_core = 0;
				lab_core = 0;
				if (sqrt(pow(at1-log(at0), 2) + pow(bt1-bt0, 2)) < 0.01){
					M_exit = 0;
					at0 = exp(at1);
					bt0 = bt1;
				} else{
					at0 = exp(at1);
					bt0 = bt1;
					n_MCycle += 1;
				}
			}
			m = j + J * n_ECycle[0];
			if (isnormal(at0) && isnormal(bt0) && isnormal(ct0)){
				if (ParConstraint[0]){
					if (at0>=0.001 && at0<=6 && bt0>=-6 && bt0<=6){
						A[j] = at0;
						B[j] = bt0;
						TA[m] = at0;
						TB[m] = bt0;
						IA[j] = Iaa;
						IB[j] = Ibb;
						IAB[j] = Iab;
					}
					if (ct0>=0.0001 && ct0<=0.5){
						C[j] = ct0;
						TC[m] = ct0;
						IC[j] = Icc;
					}
				}else{
					if (at0>=0.001){
						A[j] = at0;
					}
					A[j] = at0;
					B[j] = bt0;
					C[j] = ct0;
					TA[m] = at0;
					TB[m] = bt0;
					TC[m] = ct0;
					IA[j] = Iaa;
					IB[j] = Ibb;
					IAB[j] = Iab;
					IC[j] = Icc;
				}
			}else{
				if (n_ECycle[0]!=0){
					n=j + J * (n_ECycle[0]-1);
					TA[m] = TA[n];
					TB[m] = TB[n];
					TC[m] = TC[n];
				}else{
					TA[m] = A[j];
					TB[m] = B[j];
					TC[m] = C[j];
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
					Pstar[v] = 1 / (1 + exp(-D * A[j] * (node_Quadpts[k] - B[j])));
					if (Pstar[v] >=1){Pstar[v] = 0.9999;}
					if (Pstar[v] <=0){Pstar[v] = 0.0001;}
					P[v] = C[j] + (1-C[j]) * Pstar[v];
					L *= P[v] * data[n] + (1 - P[v]) * (1 - data[n]);
				}
				v = i + I * k;
				LL[v] = L * weight_Quadpts[k];
				LL0[i] += LL[v];
			}
		}
		//Calculate artificial data f r fz rz
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
					EZ = (Pstar[m] / P[m]) * data[n];
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
	double delta[5]={0};
	double delta0[5]={0};
	double delta1[5]={0};
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
	for (jj = 0; jj < J; jj++) {
		z=start_SEM;
		SEM_exit=1;
		cr_SEM1=1;
		cr_SEM2=1;
		cr_SEM3=1;
		while (SEM_exit && z<=end_SEM){
			mm = jj + J * z;
			for (ParClass = 1; ParClass <= 3; ParClass++){
				if (ParClass==1 && z>=2 && cr_SEM1<0.0001){
					continue;
				}
				if (ParClass==2 && z>=2 && cr_SEM2<0.0001){
					continue;
				}
				if (ParClass==3 && z>=2 && cr_SEM3<0.0001){
					continue;
				}
				for (j = 0; j < J; j++){
					deltahat_A[j]=A[j];
					deltahat_B[j]=B[j];
                    deltahat_C[j]=C[j];
				}
				switch (ParClass){
					case 1:
                        deltahat_A[jj]=TA[mm];
						break;
					case 2:
                        deltahat_B[jj]=TB[mm];
						break;
					case 3:
                        deltahat_C[jj]=TC[mm];
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
							Pstar[v] = 1 / (1 + exp(-D * deltahat_A[j] * (node_Quadpts[k] - deltahat_B[j])));
							if (Pstar[v] >=1){Pstar[v] = 0.9999;}
							if (Pstar[v] <=0){Pstar[v] = 0.0001;}
							P[v] = deltahat_C[j] + (1-deltahat_C[j]) * Pstar[v];
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
						EZ = (Pstar[m] / P[m]) * data[n];
						fz[m] += Posterior_prob[v] * EZ * CountNum[i];
						rz[m] += Posterior_prob[v] * EZ * data[n] * CountNum[i];
					}
				}
				//estimate c parameters
				if (ParClass==3){
					//estimate c parameters
					for (k = 0; k < nq; k++) {
						m = k + nq * jj;
						c_part1 += r[m] - rz[m];
						c_part3 += f[k] - fz[m];
					}
					if (PriorC[jj]==-9 || PriorC[jj + J]==-9) {
						ct0 = c_part1 / c_part3;							    //c=(r-rz)/(f-fz)
					} else {
						ct0 = (PriorC[jj] -1 + c_part1) / (((PriorC[jj] + PriorC[jj + J]) - 2.0) + c_part3);
					}
					if (isnormal(ct0)){
						if (ParConstraint[0]){
							if (ct0>=0.0001 && ct0<=0.5){
								deltahat_C[jj] = ct0;
							}
						}else{
							deltahat_C[jj] = ct0;
						}
					}
					c_part1 = 0;
					c_part3 = 0;
				}	
				if (ParClass==1 || ParClass==2){
					//estimate a and b parameter
					at0 = deltahat_A[jj];
					bt0 = deltahat_B[jj];
					n_MCycle = 1;
					M_exit = 1;
					//M step iteration
					while (M_exit && (n_MCycle <= max_MCycle)) {
						Da = D * at0;
						D2a2 = Da * Da;
						for (k = 0; k < nq; k++) {
							x_bt = node_Quadpts[k] - bt0;   	//x-bt
							x_bt2 = x_bt * x_bt;  				//(x-bt)^2
							pstar = 1 / (1 + exp(-Da * x_bt));	//ps
							psf = pstar * f[k];					//f * ps
							wsf = psf * (1-pstar);    			//f * ws
							fz_psf = fz[k + nq * jj] - psf;  	//fz - f *ps
							la1_core += fz_psf * x_bt;			//la1: (fz - f*ps)* (x-bt)
							lb1_core += fz_psf;       			//lb1: fz - f* ps
							laa_core += wsf * x_bt2;			//laa: (f*ws)* (x-bt)^2 
							lbb_core += wsf;          			//lbb: (f*ws)
							lab_core += wsf * x_bt;   			//lab: wsf*(x-bt)
						}
						la1 = Da * la1_core;
						lb1 = -Da * lb1_core;
						laa = -D2a2 * laa_core;
						lbb = -D2a2 * lbb_core;
						lab = D2a2 * lab_core;
						// Maximize a and b
						if (PriorA[j]!=-9 && PriorA[j + J]!=-9) {
							la1 += -((log(at0)-PriorA[j])/PriorA[j + J]);
							laa += -1/PriorA[j + J];
						}
						if (PriorB[j]!=-9 && PriorB[j + J]!=-9) {
							lb1 += -((bt0-PriorB[j])/PriorB[j + J]);
							lbb += -1/PriorB[j + J];
						}
						Weight = laa * lbb - lab * lab;
						Iaa = - lbb / Weight;
						Ibb = - laa / Weight;
						Iab = lab / Weight;
						Icc = - 1 / lcc;
						at1 = log(at0) + (Iaa * la1 + Iab * lb1);
						bt1 = bt0 + (Ibb * lb1 + Iab * la1);
						la1_core = 0;
						lb1_core = 0;
						laa_core = 0;
						lbb_core = 0;
						lab_core = 0;
						if (sqrt(pow(at1-log(at0), 2) + pow(bt1-bt0, 2)) < 0.01){
							M_exit = 0;
							at0 = exp(at1);
							bt0 = bt1;
						} else{
							at0 = exp(at1);
							bt0 = bt1;
							n_MCycle += 1;
						}
					}
					if (isnormal(at0) && isnormal(bt0)){
						if (ParConstraint[0]){
							if (at0>=0.001 && at0<=6 && bt0>=-6 && bt0<=6){
								deltahat_A[jj] = at0;
								deltahat_B[jj] = bt0;
							}
						}else{
							if (at0>=0.001){
								deltahat_A[jj] = at0;
							}
							deltahat_B[jj] = bt0;
						}
					}
				}
				switch (ParClass){
					case 1:
						delta1[0]=(log(deltahat_A[jj])-log(A[jj]))/(log(TA[mm])-log(A[jj])+0.0001);
						delta1[1]=(deltahat_B[jj]-B[jj])/(log(TA[mm])-log(A[jj])+0.0001);
						break;
					case 2:
						delta1[2]=(log(deltahat_A[jj])-log(A[jj]))/(TB[mm]-B[jj]+0.0001);
						delta1[3]=(deltahat_B[jj]-B[jj])/(TB[mm]-B[jj]+0.0001);
						break;
					case 3:
						delta1[4]=(deltahat_C[jj]-C[jj])/(TC[mm]-C[jj]+0.0001);   
						break;						
				}
			}
			cr_SEM1=sqrt(pow(delta1[0]-delta0[0],2)+pow(delta1[1]-delta0[1],2));
			cr_SEM2=sqrt(pow(delta1[2]-delta0[2],2)+pow(delta1[3]-delta0[3],2));
			cr_SEM3=fabs(delta1[4]-delta0[4]);
			if (cr_SEM1<0.0001 && cr_SEM2<0.0001 && cr_SEM3<0.0001 && z>=2){
				SEM_exit=0;
			}else{
				z=z+1;
			}
			for (m = 0; m<=4; m++){
				if (isnormal(delta1[m])){
					delta0[m] = delta1[m];
				}
			}
		}
		delta[0]=1-delta0[0];
		delta[1]= -delta0[1];
		delta[2]= -delta0[2];
		delta[3]=1-delta0[3];
		delta[4]=1-delta0[4];
		Weight = delta[0] * delta[3] - delta[1] * delta[2];
		delta1[0] =  delta[3] / Weight;
		delta1[1] = -delta[1] / Weight;
		delta1[2] = -delta[2] / Weight;
		delta1[3] =  delta[0] / Weight;
		delta1[4] =  1 / delta[4];
		if (isnormal(delta1[0])==0 || delta1[0]<=0){delta1[0] = 1;}
		if (isnormal(delta1[1])==0 || delta1[1]<=0){delta1[1] = 0;}
		if (isnormal(delta1[2])==0 || delta1[2]<=0){delta1[2] = 0;}
		if (isnormal(delta1[3])==0 || delta1[3]<=0){delta1[3] = 1;}
		if (isnormal(delta1[4])==0 || delta1[4]<=0){delta1[4] = 1;}
		SEA[jj]= sqrt(A[jj] * A[jj] * IA[jj] * delta1[0] + IAB[jj] * delta1[2]);
        SEB[jj]= sqrt(IB[jj] * delta1[3] + IAB[jj] * delta1[1]);
        SEC[jj]= sqrt(IC[jj] * delta1[4]);
		if (SEA[jj]>1){SEA[jj]= sqrt(A[jj] * A[jj] * IA[jj]);}
		if (SEB[jj]>1){SEB[jj]= sqrt(IB[jj]);}
		if (SEC[jj]>1){SEC[jj]= sqrt(IC[jj]);}
	}
	n_ECycle[0] = n_ECycle[0] + 1;
}
