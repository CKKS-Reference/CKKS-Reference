#pragma once

#include "linear_transform.h"
// L : level at the start of cts
// L_evalmod : level at the start of stc (= at the end of evalmod)
template <int N, int logN, int L, int L_evalmod>
void bootstrap_key_gen(
    const uint64_t q[L],
	const uint64_t p[L],
	const int s[N],
	uint64_t ctskey[logN-1][3][2][2*L][N],
	uint64_t stckey[logN-1][3][2][2*L_evalmod][N],
    uint64_t ckey[2][2*L][N]) {
	SparseComplexMatrix<N/2,3> E[logN-1];	
	splitU0R<logN>(E);
	{
		SparseComplexMatrix<N/2,3> U0R[logN-1];
		double cts_cnst = std::pow(1.0/N, 1.0 / (logN - 1));
		for(int i=0;i<logN-1;i++){
			U0R[i] = E[i];
			U0R[i].transpose();
			U0R[i].conjugate();
			U0R[i] *= cts_cnst;
		}
		rkey_gen<N,L,logN-1,3>(U0R,q,p,s,ctskey);		
	}
	{
		SparseComplexMatrix<N/2,3> U0[logN-1];
		for(int i = 0; i < logN-1; ++i) {
			U0[i] = E[logN-2-i];
		}
		rkey_gen<N,L_evalmod,logN-1,3>(U0,q,p,s,stckey);
	}

	int s_conj[N];
	conj<N>(s, s_conj);
	swkgen<N,L>(s_conj, s, q, p, ckey);
}

// Delta : scalefactor to encode iDFT matrix
template<int N, int logN, int L>
void CoeffToSlot ( const uint64_t ct_hat[2][L][N],
                   const uint64_t rkey[logN-1][3][2][2*L][N],
				   const uint64_t ckey[2][2*L][N],
                   const uint64_t Delta,
                   const uint64_t q[L],
	               const uint64_t p[L],
                   uint64_t ct_res_hat[2][2][L][N]) {
	SparseComplexMatrix<(N/2),3> E[logN-1]; 
	splitU0R<logN>(E);
    // multiply 1/N through the lin transform
    // *0.5 is done together with mat mul.
    double cts_cnst = std::pow(1.0/N, 1.0 / (logN - 1));
    for(int i=0;i<logN-1;i++){
		E[i].transpose();
		E[i].conjugate();
        E[i] *= cts_cnst;
	}

	uint64_t U0Rct_hat[2][L][N];
    serial_linear_transform<N, logN, L, logN-1, 3>(
		E, Delta, q, p, rkey, ct_hat, U0Rct_hat
	);

	uint64_t U0Rct_conj_hat[2][L][N], U0Rct_conj_hat_ks[2][L][N];
	conj_hat<N,L>(q, U0Rct_hat, U0Rct_conj_hat);
	ks<N,L>(q, p, ckey, U0Rct_conj_hat, U0Rct_conj_hat_ks);

	uint64_t pti_hat[L][N];
	for(int i = 0; i < L; ++i) {
		for(int j = 0; j < N; ++j) {
			pti_hat[i][j] = (j == N/2) ? q[i] - 1 : 0;
		}
	}
	ntt<N,L>(q, pti_hat);

	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < L; ++j) {
			for(int k = 0; k < N; ++k) {
				ct_res_hat[0][i][j][k] = (U0Rct_hat[i][j][k] + U0Rct_conj_hat_ks[i][j][k]) % q[j];
			}
		}
	}

	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < L; ++j) {
			for(int k = 0; k < N; ++k) {
				ct_res_hat[1][i][j][k] = (U0Rct_hat[i][j][k] + q[j] - U0Rct_conj_hat_ks[i][j][k]) % q[j];
				ct_res_hat[1][i][j][k] = mul_mod(ct_res_hat[1][i][j][k], pti_hat[j][k], q[j]);
			}
		}
	}
}

// Delta : scalefactor to encode DFT matrix
template<int N, int logN, int L>
void SlotToCoeff ( const uint64_t ct_r_hat[2][L][N],
				   const uint64_t ct_i_hat[2][L][N],
                   const uint64_t rkey[logN-1][3][2][2*L][N],
                   const uint64_t Delta,
                   const uint64_t q[L],
	               const uint64_t p[L],
                   uint64_t ct_res_hat[2][L][N],
				   const int s[N]) {	
	SparseComplexMatrix<(N/2),3> E[logN-1];
	SparseComplexMatrix<(N/2),3> U0[logN-1];
	SparseComplexMatrix<(N/2),3> iU0[logN-1];
	splitU0R<logN>(E);
	for(int i = 0; i < logN-1; ++i) {
		U0[i] = E[logN-2-i];
	}

	for(int i = 0; i < logN-1; ++i) {
		iU0[i] = U0[i];
		if(i == 0) {
			std::complex<double> imag_unit{0.0, 1.0};
			iU0[i] *= imag_unit;
		}
	}

	uint64_t U0ct_r_hat[2][L][N];
	uint64_t iU0ct_i_hat[2][L][N];
    serial_linear_transform<N, logN, L, logN-1, 3>(
		U0, Delta, q, p, rkey, ct_r_hat, U0ct_r_hat
	);
	serial_linear_transform<N, logN, L, logN-1, 3>(
		iU0, Delta, q, p, rkey, ct_i_hat, iU0ct_i_hat
	);

	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < L; ++j) {
			for(int k = 0; k < N; ++k) {
				ct_res_hat[i][j][k] = (U0ct_r_hat[i][j][k] + iU0ct_i_hat[i][j][k]) % q[j];
			}
		}
	}
}