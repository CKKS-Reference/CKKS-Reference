#pragma once

#include "CKKS_basic.h"
/////////////////////////////////////////////////////////////////////////////
// Key Switch with DNUM
/////////////////////////////////////////////////////////////////////////////

template<int N, int LplusK>
void modUP(const uint64_t qp[LplusK], int L, 
	             uint64_t  a[LplusK][N]) {
	int K = LplusK-L; 
	const uint64_t* q=qp;
	const uint64_t* p=qp+L;
	//
	uint64_t table1[LplusK][LplusK];
	for(int j=0; j<L; j++)
	for(int k=0; k<K; k++) {
		table1[j][k] = 1;
		for(int i = 0; i < L; i++)
			if (i != j) 
				table1[j][k] = mul_mod(table1[j][k], q[i] % p[k], p[k]);
	}
	//
	uint64_t table2[LplusK];
	for (int j = 0; j < L; j++) {
		table2[j] = 1;
		for (int i = 0; i < L; i++)
			if (i != j)
				table2[j] = mul_mod(table2[j], inv_mod(q[i]%q[j], q[j]), q[j]);
	}
	//
	uint64_t table3[LplusK];
	for (int k = 0; k < K; k++) {
		table3[k] = 1;
		for (int j = 0; j < L; j++)
			table3[k] = mul_mod(table3[k], q[j] % p[k], p[k]);
	}
	//
	for (int i = 0; i < N; i++) {
		uint64_t b[LplusK]; int count = 0;
		for (int j = 0; j < L; j++) {
			b[j] = mul_mod(a[j][i], table2[j], q[j]);
			if (2 * b[j] >= q[j]) count++;
		}
		for (int k = 0; k < K; k++) {
			a[L+k][i] = 0;
			for (int j = 0; j < L; j++)
				a[L+k][i] =(a[L+k][i]+ mul_mod(b[j] % p[k], table1[j][k], p[k]))%p[k];

			if (count > 0)
				a[L+k][i]= (a[L+k][i]+ mul_mod(table3[k], p[k] - count, p[k]))%p[k];
		}
	}
}

//
template<int L, int DNUM, int K>
void gadget_g( const uint64_t q[L],
			   const uint64_t p[K],
					 uint64_t g[DNUM][L+K]){
	for(int d = 0; d < DNUM ; d++)
	for(int j = 0; j < L + K; j++){
		if((d*K <= j) && (j<(d+1)*K) && (j<L)) g[d][j] = 1;
		else                                   g[d][j] = 0;
	}
}


//
template<int N, int L, int DNUM, int K>
void gadget_ginv(const uint64_t q[L],
				 const uint64_t p[K],
				 const uint64_t a[L][N],
					   uint64_t ginva[DNUM][L+K][N]) {
	for(int d=0; d<DNUM; d++){
		uint64_t QP[L+K];
		for(int i=0; i<L; i++) QP[i] = q[i];
		for(int i=0; i<K; i++) QP[L+i] = p[i];
		//
		int off1 = ( d   *K<L)?( d   *K) : L;
		int off2 = ((d+1)*K<L)?((d+1)*K) : L;
		uint64_t temp[K];
		for(int i=off1;i<off2;i++) temp[i-off1]=QP[i];
		for(int i=off1-1;i>=0;i--) QP[i+off2-off1]=QP[i];
		for(int i=off1;i<off2;i++) QP[i-off1]=temp[i-off1];
		//
		for(int i=off1; i<off2; i++)
		for(int j=0; j<N; j++)
			ginva[d][i-off1][j] = a[i][j];
		if(off2>off1)
			modUP<N,L+K>(QP, off2-off1,ginva[d]);
		//
		for(int j=0; j<N; j++){
			for(int i=0;i<off2-off1;i++) temp[i]=ginva[d][i][j];
			for(int i=0;i<     off1;i++) ginva[d][i][j]= ginva[d][i+off2-off1][j];
			for(int i=0;i<off2-off1;i++) ginva[d][i+off1][j]=temp[i];
		}
		if (off1 == off2)
			for (int i = 0; i < L + K; i++)
				for (int j = 0; j < N; j++) ginva[d][i][j] = 0;
	}
}

//
template<int N, int L, int DNUM>
void swkgen(const int sfr[N],
			const int sto[N],
			const uint64_t q[L],
			const uint64_t p[L / DNUM],
				  uint64_t swk_hat[DNUM][2][L + (L / DNUM)][N]) {
	const int K = L / DNUM; assert(L==K*DNUM);
	uint64_t g[DNUM][L + K];
	gadget_g<L,DNUM,K>(q, p, g);

	for (int n = 0; n < DNUM; n++) {
		uint64_t pt[L + K][N];
		for (int j = 0; j < L; j++) {
			uint64_t P = 1;
			for (int k = 0; k < K; k++)
				P = mul_mod(P, p[k] % q[j], q[j]);
			uint64_t Pg = mul_mod(P, g[n][j], q[j]);
			for (int i = 0; i < N; i++)
				pt[j][i] = mul_mod((q[j] + sfr[i]) % q[j], Pg, q[j]);
		}
		for (int j = 0; j < K; j++)
		for (int i = 0; i < N; i++)
			pt[j + L][i] = 0;
		uint64_t qp[L + K];
		for (int j = 0; j < L; j++) qp[    j] = q[j];
		for (int j = 0; j < K; j++) qp[L + j] = p[j];

		enc<N, L + K>(pt, sto, qp, swk_hat[n]);
	}
}

//
template<int N,int L, int DNUM, int K>
void ks( const uint64_t q[L],
		 const uint64_t p[K],
	     const uint64_t swk_hat[DNUM][2][DNUM*K+K][N],
		 const uint64_t  ct_hat[2][L][N],
	           uint64_t out_hat[2][L][N]){
	uint64_t a[L][N];
	for(int i=0;i<L;i++)
	for(int j=0;j<N;j++) a[i][j]=ct_hat[1][i][j];
	intt<N,L>(q,a);
	//
	uint64_t ginva[DNUM][L+K][N];
	gadget_ginv<N,L,DNUM,K>(q,p,a,ginva);
	//
	uint64_t qp[L+K];
	for(int i=0;i<L;i++) qp[i]=q[i];
	for(int i=0;i<K;i++) qp[L+i]=p[i];
	//
	uint64_t sum[2][L+K][N];
	for(int i=0;i<L+K;i++)
	for(int j=0;j<N  ;j++){ sum[0][i][j]=0;
	                        sum[1][i][j]=0;}
	for(int d=0;d<DNUM;d++){
		ntt<N,L+K>(qp,ginva[d]);
		for(int i=0;i<L;i++)
		for(int j=0;j<N;j++){
			sum[0][i][j]=(sum[0][i][j]+mul_mod(ginva[d][i][j],swk_hat[d][0][i][j],q[i]))%q[i];
			sum[1][i][j]=(sum[1][i][j]+mul_mod(ginva[d][i][j],swk_hat[d][1][i][j],q[i]))%q[i];
		}
		for(int i=0;i<K;i++)
		for(int j=0;j<N;j++){
			sum[0][L+i][j]=(sum[0][L+i][j]+mul_mod(ginva[d][L+i][j],swk_hat[d][0][DNUM*K+i][j],p[i]))%p[i];
			sum[1][L+i][j]=(sum[1][L+i][j]+mul_mod(ginva[d][L+i][j],swk_hat[d][1][DNUM*K+i][j],p[i]))%p[i];
		}
	}
	//
	RS_hat<N,L+K,L>(qp,sum,out_hat);
	//
	for(int i=0;i<L;i++)
	for(int j=0;j<N;j++)
		out_hat[0][i][j]=(out_hat[0][i][j]+ct_hat[0][i][j])%q[i];
}