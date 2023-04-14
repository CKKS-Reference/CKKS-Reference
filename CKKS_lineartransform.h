#pragma once

#include "CKKS_basic.h"
#include "CKKS_ks_wih_dnum.h"

template< int N, int L>
void rot( const uint64_t q[L], const uint64_t m   [L][N], int r,
	                                 uint64_t mrot[L][N]){
	if(r==0){
		for(int i=0;i<L;i++)
		for(int j=0;j<N;j++)
			mrot[i][j]=m[i][j];
	}
	else{
		uint64_t temp[L][N];
		for(int i=0;i<L;i++)
		for(int j=0;j<N;j++){
			uint64_t Aj=(5*j)/N;
			uint64_t rj=(5*j)%N;
			temp[i][rj]=(Aj%2==0)?m[i][j]:((q[i]-m[i][j])%q[i]);
		}
		rot<N,L>(q, temp,r-1,mrot);
	}
}
//
template<int N>
void rot( const int s[N], int r, int srot[N] ){
	if(r==0){
		for(int j=0;j<N;j++)
			srot[j]=s[j];
	}
	else{
		int temp[N];
		for(int j=0;j<N;j++){
			int Aj=(5*j)/N;
			int rj=(5*j)%N;
			temp[rj]=(Aj%2==0)?s[j]:-s[j];
		}
		rot<N>(temp,r-1,srot);
	}
}
//
template< int N, int L, int DNUM, int K, int NumDiags >
void rkey_gen( const SparseComplexMatrix<N/2,NumDiags>& A,
	           const uint64_t q[L],
	           const uint64_t p[K],
	           const  int     s[N],
	                 uint64_t rkey_hat[NumDiags][DNUM][2][DNUM*K+K][N]){
	for(int i=0;i<NumDiags;i++){
		int srot[N];
		rot<N>(s,A.shift[i],srot);
		swkgen<N,L,DNUM>(srot,s,q,p,rkey_hat[i]);
	}
}

//
template<int N, int logN, int L, int DNUM, int K, int NumDiags>
void lineartransform(const SparseComplexMatrix<N/2, NumDiags> & A,
				     const uint64_t Delta,
					 const uint64_t q[L],
					 const uint64_t p[K],
					 const uint64_t rkey_hat[NumDiags][DNUM][2][DNUM*K+K][N],
					 const uint64_t ct_hat[2][L][N],
						   uint64_t res_hat[2][L][N]){
	uint64_t ct[2][L][N];
	for(int i=0; i<2; i++)
	for(int j=0; j<L; j++)
	for(int k=0; k<N; k++){
		ct[i][j][k] = ct_hat[i][j][k];
		res_hat[i][j][k] = 0;
	}
	intt<N, L>(q, ct[0]);
	intt<N, L>(q, ct[1]);

	for(int d=0; d<NumDiags; d++){
		uint64_t pt[L][N];
		encode<N, logN, L>(A.diagr[d], A.diagi[d], Delta, q, pt);

		uint64_t temp[2][L][N];
		rot<N, L>(q, ct[0], A.shift[d], temp[0]);
		ntt<N, L>(q, temp[0]);

		rot<N, L>(q, ct[1], A.shift[d], temp[1]);
		ntt<N, L>(q, temp[1]);

		uint64_t ct_rot_hat[2][L][N];
		ks<N, L, DNUM, K>(q, p, rkey_hat[d], temp, ct_rot_hat);

		ntt<N, L>(q, pt);
		for(int i=0; i<2; i++)
		for(int j=0; j<L; j++)
		for(int k=0; k<N; k++)
			res_hat[i][j][k] = (res_hat[i][j][k] + mul_mod(pt[j][k], ct_rot_hat[i][j][k], q[j]))%q[j];
	}
}
	


template<int N>
void rot(int shift, double zr[N], double zi[N]){
	double zr_copy[N];
	double zi_copy[N];
	for(int i=0; i<N; i++){
		zr_copy[i] = zr[i];
		zi_copy[i] = zi[i];
	}

	for(int i=0; i<N; i++){
		zr[i] = zr_copy[(i+shift)%N];
		zi[i] = zi_copy[(i+shift)%N];
	}
}

template<int N, int logN, int L, int DNUM, int K, int NumDiags, int NumMats>
void lineartransform(const SparseComplexMatrix<N/2, NumDiags> E[NumMats],
					 const uint64_t Delta,
					 const uint64_t q[L],
					 const uint64_t p[K],
					 const uint64_t rkey_hat[NumMats][NumDiags][DNUM][2][DNUM*K+K][N],
					 const uint64_t ct_hat [2][L][N],
						   uint64_t res_hat [2][L][N], double k=1.0){
	uint64_t ct[2][L][N];
	for(int i=0; i<2; i++)
	for(int j=0; j<L; j++)
	for(int k=0; k<N; k++){
		ct[i][j][k] = ct_hat[i][j][k];
		res_hat[i][j][k] = 0;
	}
	intt<N, L>(q, ct[0]);
	intt<N, L>(q, ct[1]);

	int count = 1;
	for(int i=0; i<NumMats; i++) count*= NumDiags;

	for(int i=0; i<count; i++){
		int isave = i;

		double temp_vr[N/2], temp_vi[N/2];
		uint64_t temp_r[2][L][N];
		for(int j=0; j<N/2; j++){
			temp_vr[j] = k;
			temp_vi[j] = 0;
		}

		for(int j=0; j<2; j++)
		for(int k=0; k<L; k++)
		for(int l=0; l<N; l++)
			temp_r[j][k][l] = ct[j][k][l];

		for(int j=NumMats-1; j>=0; j--){
			int i_j = isave % NumDiags;
			isave = isave / NumDiags;

			int shift = E[j].shift[i_j];
			rot<N/2>(shift, temp_vr, temp_vi);

			for(int k=0; k<N/2; k++){
				double save = temp_vr[k];
				temp_vr[k] = temp_vr[k]*E[j].diagr[i_j][k] - temp_vi[k]*E[j].diagi[i_j][k];
				temp_vi[k] = save * E[j].diagi[i_j][k] + temp_vi[k]*E[j].diagr[i_j][k];
			}

			uint64_t temp_r_rot[2][L][N];
			rot<N, L>(q, temp_r[0], shift, temp_r_rot[0]); ntt<N, L>(q, temp_r_rot[0]);
			rot<N, L>(q, temp_r[1], shift, temp_r_rot[1]); ntt<N, L>(q, temp_r_rot[1]);
			ks<N, L, DNUM, K>(q, p, rkey_hat[j][i_j], temp_r_rot, temp_r);
			intt<N, L>(q, temp_r[0]);
			intt<N, L>(q, temp_r[1]);
		}

		uint64_t pt[L][N];
		encode<N, logN, L>(temp_vr, temp_vi, Delta, q, pt);

		ntt<N, L>(q, temp_r[0]);
		ntt<N, L>(q, temp_r[1]);
		ntt<N, L>(q, pt);

		for(int j=0; j<2; j++)
		for(int k=0; k<L; k++)
		for(int l=0; l<N; l++)
			res_hat[j][k][l] = (res_hat[j][k][l] + mul_mod(pt[k][l], temp_r[j][k][l], q[k]))%q[k];

		bool test = true;
	}
}

template<int N, int logN, int L, int DNUM, int K, int NumDiags, int NumMats>
void lineartransform_reverse(const SparseComplexMatrix<N / 2, NumDiags> E[NumMats],
	const uint64_t Delta,
	const uint64_t q[L],
	const uint64_t p[K],
	const uint64_t rkey_hat[NumMats][NumDiags][DNUM][2][DNUM * K + K][N],
	const uint64_t ct_hat[2][L][N],
	uint64_t res_hat[2][L][N], double k=1.0) {
	uint64_t ct[2][L][N];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < L; j++)
			for (int k = 0; k < N; k++) {
				ct[i][j][k] = ct_hat[i][j][k];
				res_hat[i][j][k] = 0;
			}
	intt<N, L>(q, ct[0]);
	intt<N, L>(q, ct[1]);

	int count = 1;
	for (int i = 0; i < NumMats; i++) count *= NumDiags;

	for (int i = 0; i < count; i++) {
		int isave = i;

		double temp_vr[N / 2], temp_vi[N / 2];
		uint64_t temp_r[2][L][N];
		for (int j = 0; j < N / 2; j++) {
			temp_vr[j] = k;
			temp_vi[j] = 0;
		}

		for (int j = 0; j < 2; j++)
			for (int k = 0; k < L; k++)
				for (int l = 0; l < N; l++)
					temp_r[j][k][l] = ct[j][k][l];

		for (int j = 0; j <NumMats; j++) {
			int i_j = isave % NumDiags;
			isave = isave / NumDiags;

			int shift = E[j].shift[i_j];
			rot<N / 2>(shift, temp_vr, temp_vi);

			for (int k = 0; k < N / 2; k++) {
				double save = temp_vr[k];
				temp_vr[k] = temp_vr[k] * E[j].diagr[i_j][k] - temp_vi[k] * E[j].diagi[i_j][k];
				temp_vi[k] = save * E[j].diagi[i_j][k] + temp_vi[k] * E[j].diagr[i_j][k];
			}

			uint64_t temp_r_rot[2][L][N];
			rot<N, L>(q, temp_r[0], shift, temp_r_rot[0]); ntt<N, L>(q, temp_r_rot[0]);
			rot<N, L>(q, temp_r[1], shift, temp_r_rot[1]); ntt<N, L>(q, temp_r_rot[1]);
			ks<N, L, DNUM, K>(q, p, rkey_hat[j][i_j], temp_r_rot, temp_r);
			intt<N, L>(q, temp_r[0]);
			intt<N, L>(q, temp_r[1]);
		}

		uint64_t pt[L][N];
		encode<N, logN, L>(temp_vr, temp_vi, Delta, q, pt);

		ntt<N, L>(q, temp_r[0]);
		ntt<N, L>(q, temp_r[1]);
		ntt<N, L>(q, pt);

		for (int j = 0; j < 2; j++)
			for (int k = 0; k < L; k++)
				for (int l = 0; l < N; l++)
					res_hat[j][k][l] = (res_hat[j][k][l] + mul_mod(pt[k][l], temp_r[j][k][l], q[k])) % q[k];

		bool test = true;
	}
}