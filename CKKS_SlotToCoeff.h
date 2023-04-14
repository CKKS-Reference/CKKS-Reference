#pragma once

#include "CKKS_ks_wih_dnum.h"
#include "CKKS_lineartransform.h"

//conventional STC
template< int L, int DNUM, int K>
void SlotToCoeff_logN_10( const uint64_t q[L],
	                      const uint64_t p[K], uint64_t Delta,
	                      const SparseComplexMatrix< 1<<9, 27> A[3],
						  const uint64_t rkey_hat[3][27][DNUM][2][DNUM*K+K][1<<10],
					      const uint64_t ct_evalmod_hat[2][2][L  ][1<<10],
	                            uint64_t ct_stc_hat       [2][L-3][1<<10]){
	const int N=1<<10;
	//
	uint64_t ct1_hat[2][2][L][N];

	for(int i=0;i<2;i++)
	for(int j=0;j<2;j++)
	for(int k=0;k<L;k++)
	for(int w=0;w<N;w++) ct1_hat[i][j][k][w] = ct_evalmod_hat[i][j][k][w];

	for(int n=2;n>=0;n--){
		uint64_t temp[2][2][L][N];
		for(int i=0;i<2;i++)
			lineartransform<N,10,L,DNUM,K,27>(A[n],Delta,q,p,rkey_hat[n],ct1_hat[i],temp[i]);
		for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
		for(int k=0;k<L;k++)
		for(int w=0;w<N;w++) ct1_hat[i][j][k][w] = temp[i][j][k][w];
	}

	//
	uint64_t ct2_hat[2][L-3][N];
	uint64_t ct3_hat[2][L-3][N];
	RS_hat<1<<10,L,L-3>(q,ct1_hat[0],ct2_hat);
	RS_hat<1<<10,L,L-3>(q,ct1_hat[1],ct3_hat);

	//
	uint64_t pti_hat[L-3][N];
	for(int i=0;i<L-3;i++){
		for(int j=0;j<N;j++)
			pti_hat[i][j]=0;
		pti_hat[i][N/2]=1;
	}
	ntt<N,L-3>(q,pti_hat);
	//
	for(int i=0;i<2;i++)
	for(int j=0;j<L-3;j++)
	for(int k=0;k<N;k++)
		ct_stc_hat[i][j][k] = (ct2_hat[i][j][k]+
			               mul_mod(pti_hat[j][k],ct3_hat[i][j][k],q[j]))%q[j];
}

// reordered STC
template< int L, int DNUM, int K>
void SlotToCoeff_logN_10(const uint64_t q[L],
						 const uint64_t p[K], uint64_t Delta,
						 const SparseComplexMatrix< 1 << 9, 27> A[3],
						 const uint64_t rkey_hat[3][27][DNUM][2][DNUM * K + K][1 << 10],
						 const uint64_t ct_evalmod_hat[2][L    ][1 << 10],
							   uint64_t ct_stc_hat    [2][L - 3][1 << 10]) {
	const int N = 1 << 10;
	//
	uint64_t ct1_hat[2][L][N];

	for (int j = 0; j < 2; j++)
	for (int k = 0; k < L; k++)
	for (int w = 0; w < N; w++) ct1_hat[j][k][w] = ct_evalmod_hat[j][k][w];

	for (int n = 2; n >= 0; n--) {
		uint64_t temp[2][L][N];
		lineartransform<N, 10, L, DNUM, K, 27>(A[n], Delta, q, p, rkey_hat[n], ct1_hat, temp);
		
		for (int j = 0; j < 2; j++)
		for (int k = 0; k < L; k++)
		for (int w = 0; w < N; w++) ct1_hat[j][k][w] = temp[j][k][w];
	}

	//
	RS_hat<1 << 10, L, L - 3>(q, ct1_hat, ct_stc_hat);
	
}

// assumption : logN-1 is multiple of MDNUM
template<int N, int logN, int L, int DNUM, int K, int MDNUM>
void SlotToCoeff(const uint64_t q[L],
				 const uint64_t p[K], uint64_t Delta,
				 const SparseComplexMatrix<N/2, 3> E[logN-1],
				 const uint64_t rkey_hat[logN - 1][3][DNUM][2][DNUM * K + K][N],
				 const uint64_t ct_hat      [2][L        ][N],
				       uint64_t ct_stc_hat  [2][L - MDNUM][N]) {
	
	const int NumMats = (logN - 1) / MDNUM;
	uint64_t ct1_hat[2][L][N];

	for (int j = 0; j < 2; j++)
	for (int k = 0; k < L; k++)
	for (int w = 0; w < N; w++) ct1_hat[j][k][w] = ct_hat[j][k][w];

	//-----------------------------------------------
	// errorless scaling
	//-----------------------------------------------
	double k = 1.0;
	for (int i = 1; i <= MDNUM; i++)
		k *= ((double)q[L - i]) / Delta;

	for (int n = MDNUM - 1; n >= 0; n--) 
	{
		uint64_t temp[2][L][N];
		lineartransform<N, logN, L, DNUM, K, 3, NumMats>(E + n * NumMats, Delta, q, p, rkey_hat + n * NumMats, ct1_hat, temp, n==0?k:1.);
		
		for (int j = 0; j < 2; j++)
		for (int k = 0; k < L; k++)
		for (int w = 0; w < N; w++) ct1_hat[j][k][w] = temp[j][k][w];

	}

	//
	RS_hat<N, L, L - MDNUM>(q, ct1_hat, ct_stc_hat);
	
}