#pragma once
#include <stdint.h>
#include <inttypes.h>
#include <cassert>
#include "DFT.h"
#include "NumberTheory.h"
#include "BigInt.h"


//
template<int N, int L, int DNUM, int K>
void mul(const uint64_t q[L], const uint64_t p[K], 
	     const uint64_t evk_hat[DNUM][2][DNUM*K+K][N], 
	     const uint64_t A_hat[2][L][N], const uint64_t B_hat[2][L][N], uint64_t C_hat[2][L][N]){

	uint64_t d0_hat[L][N];
	uint64_t d1_hat[L][N];
	uint64_t d2_hat[L][N];

	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++){
		d0_hat[i][j] = mul_mod(A_hat[0][i][j], B_hat[0][i][j], q[i]);
		d1_hat[i][j] = ( mul_mod(A_hat[0][i][j], B_hat[1][i][j], q[i])
			            +mul_mod(A_hat[1][i][j], B_hat[0][i][j], q[i]) )%q[i];
		d2_hat[i][j] = mul_mod(A_hat[1][i][j], B_hat[1][i][j], q[i]);
	}

	uint64_t temp[2][L][N];
	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++){
		temp[0][i][j] = 0;
		temp[1][i][j] = d2_hat[i][j];
	}

	ks<N,L, DNUM, K>(q,p,evk_hat,temp,C_hat);
	
	for(int i=0; i<L; i++)
	for(int j=0; j<N; j++){
		C_hat[0][i][j] = (C_hat[0][i][j] + d0_hat[i][j]) %q[i];
		C_hat[1][i][j] = (C_hat[1][i][j] + d1_hat[i][j]) %q[i];
	}
}

//
void convert_poly_to_binarytreeform(double *c, int N) {
	if (N == 2) return;

	for (int i = 1; i < N / 2; i++) {
		c[N/2-i] -= c[N/2+i];
		c[N/2+i] *= 2;
	}

	convert_poly_to_binarytreeform(c    ,N/2);
	convert_poly_to_binarytreeform(c+N/2,N/2);
}

template<int N, int L, int DNUM, int K>
void mul_rs( const uint64_t q[L], const uint64_t p[K],
	         const uint64_t evk_hat[DNUM][2][DNUM*K+K][N],
	         const uint64_t Ahat[2][L  ][N],
	         const uint64_t Bhat[2][L  ][N],
	               uint64_t Chat[2][L-1][N]){
	uint64_t temp[2][L][N];
	mul<N,L,DNUM,K>(q,p,evk_hat,Ahat,Bhat,temp);
	RS_hat<N,L>(q,temp[0],Chat[0]);
	RS_hat<N,L>(q,temp[1],Chat[1]);
}

//
template<int N, int L, int DNUM, int K, int logd>
void EvalPoly(const uint64_t q[L],
			  const uint64_t p[K],
			  const uint64_t evk_hat[DNUM][2][K*DNUM+K][N],
			  const uint64_t ct_t_hat[2][L][N],
				    uint64_t Delta_L,
					uint64_t Delta_L_minus_logd,
			  const double* c,
				    uint64_t ct_p_hat      [2][L-logd  ][N],
					uint64_t ct_T_d_over_2 [2][L-logd+1][N],
					int      case_){
	int d = 1<< logd;
	if constexpr(logd==1){
		double k = ((((double)Delta_L_minus_logd)/Delta_L)*q[L-1])/Delta_L;
		double kc[2]; kc[0]=k*c[0]; kc[1]=k*c[1];
		uint64_t pt[2][L][N];
		for(int i=0;i<2;i++)
		for(int j=0;j<L;j++){
			for(int k=1;k<N;k++)
				pt[i][j][k]=0;
			if(kc[i]>=0) pt[i][j][0]=       uint64_t(  kc[i] *Delta_L+0.5)%q[j]       ;
			else         pt[i][j][0]=(q[j]-(uint64_t((-kc[i])*Delta_L+0.5)%q[j]))%q[j];
		}
		ntt<N,L>(q,pt[0]);
		ntt<N,L>(q,pt[1]);
		uint64_t temp[2][L][N];
		for(int j=0;j<L;j++)
		for(int k=0;k<N;k++){
			temp[0][j][k]=(mul_mod(pt[0][j][k],     Delta_L%q[j],q[j])
				          +mul_mod(pt[1][j][k],ct_t_hat[0][j][k],q[j]))%q[j];
			temp[1][j][k]= mul_mod(pt[1][j][k],ct_t_hat[1][j][k],q[j]);
		}
		RS_hat<N,L>(q,temp[0],ct_p_hat[0]);
		RS_hat<N,L>(q,temp[1],ct_p_hat[1]);

		if(case_==1){
			for(int i=0;i<2;i++)
			for(int j=0;j<L;j++)
			for(int k=0;k<N;k++)
				ct_T_d_over_2[i][j][k]=ct_t_hat[i][j][k];
		}
	}
	else{
		uint64_t ct_r_hat     [2][L-logd+1][N];
		uint64_t ct_T_d_over_4[2][L-logd+2][N];
		
		uint64_t Delta_tilde_L_minus_logd_plus_1 = Delta_L;
		for(int j=1; j<logd; j++)
			Delta_tilde_L_minus_logd_plus_1 = round_ab_over_c(Delta_tilde_L_minus_logd_plus_1,Delta_tilde_L_minus_logd_plus_1,q[L-j]);

		uint64_t Delta_L_minus_logd_plus_1 = uint64_t((((double)Delta_L_minus_logd)/Delta_tilde_L_minus_logd_plus_1)*q[L-logd]+0.5);
		EvalPoly<N,L,DNUM,K,logd-1>(q,p,evk_hat,ct_t_hat,Delta_L,Delta_L_minus_logd_plus_1,c    ,ct_r_hat,ct_T_d_over_4,1);
		uint64_t ct_q_hat [2][L-logd+1][N];
		EvalPoly<N,L,DNUM,K,logd-1>(q,p,evk_hat,ct_t_hat,Delta_L,Delta_L_minus_logd_plus_1,c+d/2,ct_q_hat,ct_T_d_over_4,2);

		if(case_==1){
			uint64_t temp[2][L-logd+2][N];
			mul<N,L-logd+2,DNUM,K>(q,p,evk_hat,ct_T_d_over_4,ct_T_d_over_4,temp);
			for(int i=0;i<L-logd+2;i++)
			for(int j=0;j <N;j++){
				temp[0][i][j] = (2*temp[0][i][j])%q[i];
				temp[1][i][j] = (2*temp[1][i][j])%q[i];
			}
			RS_hat<N,L-logd+2>(q,temp[0],ct_T_d_over_2[0]);
			RS_hat<N,L-logd+2>(q,temp[1],ct_T_d_over_2[1]);
			for(int i=0;i<L-logd+1;i++)
			for(int j=0;j<N;j++)
				ct_T_d_over_2[0][i][j] = (ct_T_d_over_2[0][i][j]+q[i]- (Delta_tilde_L_minus_logd_plus_1% q[i])) % q[i];
		}												
		uint64_t temp[2][L-logd+1][N];
		mul<N,L-logd+1,DNUM,K>(q,p,evk_hat,ct_T_d_over_2,ct_q_hat,temp);

		for(int i=0;i<L-logd+1;i++)
		for(int j=0;j<N;j++){
			temp[0][i][j] = (temp[0][i][j]+mul_mod(ct_r_hat[0][i][j],Delta_tilde_L_minus_logd_plus_1%q[i],q[i]))%q[i];
			temp[1][i][j] = (temp[1][i][j]+mul_mod(ct_r_hat[1][i][j],Delta_tilde_L_minus_logd_plus_1%q[i],q[i]))%q[i];
		}
		RS_hat<N,L-logd+1>(q,temp[0],ct_p_hat[0]);
		RS_hat<N,L-logd+1>(q,temp[1],ct_p_hat[1]);
	}
}
