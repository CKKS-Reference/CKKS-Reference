#include <stdio.h>
#include "CKKS_CoeffToSlot.h"
#include "CKKS_poly.h"
#include "CKKS_EvalMod.h"
#include "CKKS_SlotToCoeff.h"

#define logN 10
#define MDNUM 3
#define N (1<<logN)
#define L 25
#define DNUM 5
#define K (L/DNUM)
#define Delta (1ULL<<42)



double norm_square_exp(const double zr[N/2], const double zi[N/2]){
	double sum=0;
	for(int i=0; i<N/2; i++) sum+=zr[i]*zr[i]+zi[i]*zi[i];
	return (sum/(N/2));
}

double norm_square(const double zr[N / 2], const double zi[N / 2]) {
	double sum = 0;
	for (int i = 0; i < N / 2; i++)
		sum += zr[i] * zr[i] + zi[i] * zi[i];
	return (sum );
}

double norm_max(const double zr[N / 2], const double zi[N / 2]) {
	double max = 0;
	for (int i = 0; i < N / 2; i++) {
		double abs = sqrt(zr[i] * zr[i] + zi[i] * zi[i]);
		if (abs > max) max = abs;
	}
	return max;
}

//
void main(){
	//---------------------------------------------------------------------
	// Initial setting
	//---------------------------------------------------------------------
	uint64_t q[L];	uint64_t p[K];
	{
		uint64_t q58[18] = {288230376151760897ULL, 288230376152137729ULL, 288230376152154113ULL, 288230376152350721ULL, 288230376152727553ULL, 288230376152973313ULL, 288230376153726977ULL, 288230376153858049ULL, 288230376154267649ULL,
			288230376154300417ULL, 288230376154316801ULL, 288230376154562561ULL, 288230376154841089ULL, 288230376155185153ULL, 288230376155201537ULL, 288230376155250689ULL, 288230376155545601ULL, 288230376156610561ULL };
		uint64_t q42[12] = {4398047051777ULL,4398047232001ULL,4398047543297ULL,4398047772673ULL,4398048034817ULL,4398048526337ULL,4398048575489ULL,4398048706561ULL,4398048903169ULL,4398048919553ULL,4398048952321ULL,4398049116161ULL};
		// Base
		q[0] = q58[0];
		// STC
		for (int i = 0; i < 3; i++) q[1 + i] = q42[i];
		//CTS
		for (int i = 0; i < 3; i++) q[22 + i] = q58[1 + i];
		// EvalMod
		for (int i = 0; i < 9; i++) q[21 - i] = q58[4+i];
		//general
		for (int i = 0; i < 9; i++) q[4 + i] = q42[3 + i];
		//P
		for (int i = 0; i < 5; i++) p[i] = q58[13 + i];
	}
	int h=192; int s[N];
	keygen<N>(h,s);
	double zr[N/2], zi[N/2];
	for(int i=0;i<N/2;i++){
		zr[i]=((double)rand())/RAND_MAX/sqrt(2);
		zi[i]=((double)rand())/RAND_MAX/sqrt(2);
	}
	uint64_t pt[4][N];
	encode<N,logN,4>(zr,zi,Delta,q,pt);
	uint64_t ct_bot_hat[2][4][N];
	enc<N,4>(pt,s,q,ct_bot_hat);
	printf("Encryption complete\n");
	//---------------------------------------------------------------------
	// swkgen
	//---------------------------------------------------------------------
	SparseComplexMatrix<N / 2, 3> E[logN - 1];
	splitU0R<logN>(E);
	uint64_t rkey_hat[logN - 1][3][DNUM][2][DNUM * K + K][N];
	for (int n = 0; n < logN - 1; n++)
		rkey_gen<N, L, DNUM, K, 3>(E[n], q, p, s, rkey_hat[n]); 
	
	uint64_t  evk_hat[DNUM][2][DNUM * K + K][N];
	uint64_t ckey_hat[DNUM][2][DNUM * K + K][N];
	int ss[N]; conv<N>(s, s, ss);
	swkgen< N, L, DNUM>(ss, s, q, p, evk_hat);
	int sconj[N]; conj<N>(s, sconj);
	swkgen<N, L, DNUM>(sconj, s, q, p, ckey_hat);
	
	printf("KeyGen complete\n");

	
	//---------------------------------------------------------------------
	// SlotToCoeff
	//---------------------------------------------------------------------
	
	uint64_t ct_hat_stc[2][1][N];
	SlotToCoeff<N, logN, 4, DNUM, K, MDNUM>(q, p, Delta, E, rkey_hat, ct_bot_hat, ct_hat_stc);
	printf("SlotToCoeff complete\n");
	//-----------------------------------------------
	// debug
	//-----------------------------------------------
	double zDeltaBR[N];
	{
		uint64_t pt[1][N];
		double wr[N / 2], wi[N / 2];
		dec<N, 1>(ct_hat_stc, s, q, pt);
		decode<N, logN, 1>(pt, Delta, q, wr, wi);
		ifft<N, logN>(wr, wi, zDeltaBR);
		for (int i = 0; i < N; i++) zDeltaBR[i] *= Delta;

		double zDeltaBRmin = 1e8, zDeltaBRmax = -1e8;
		for (int i = 0; i < N; i++) {
			if (zDeltaBR[i] < zDeltaBRmin) zDeltaBRmin = zDeltaBR[i];
			if (zDeltaBR[i] > zDeltaBRmax) zDeltaBRmax = zDeltaBR[i];
		}
	}
	//---------------------------------------------------------------------
	// ModUp
	//---------------------------------------------------------------------
	uint64_t ct_hat[2][L][N];
	for(int i=0;i<2;i++)
	for(int j=0;j<N;j++)
		ct_hat[i][0][j] = ct_hat_stc[i][0][j];
	intt<N,1>(q,ct_hat[0]);
	intt<N,1>(q,ct_hat[1]);
	modUP<N, L>(q, 1, ct_hat[0]);
	modUP<N, L>(q, 1, ct_hat[1]);
	ntt<N,L>(q, ct_hat[0]);
	ntt<N,L>(q, ct_hat[1]);
	printf("ModUp complete\n");
	//----------------------------------------------------
	// debug : pt = ifft(z), pt+qI = ifft(decode(dec(ct))
	//----------------------------------------------------
	double I[N], I_BR[N];
	{
		double zDeltaBR_plus_qI_BR[N], qI_BR[N];
		uint64_t pt[L][N]; dec<N, L>(ct_hat, s, q, pt);
		double wr[N / 2], wi[N / 2]; decode<N, logN, L>(pt, Delta, q, wr, wi);
		ifft<N, logN>(wr, wi, zDeltaBR_plus_qI_BR);
		for (int i = 0; i < N; i++) {
			zDeltaBR_plus_qI_BR[i] *= Delta;
			qI_BR[i] = zDeltaBR_plus_qI_BR[i] - zDeltaBR[i];
			I_BR[i] = qI_BR[i] / q[0];
		}

		bitReverse<N / 2>(I_BR, I); 
		bitReverse<N / 2>(I_BR + N/2, I + N/2);

		double Imin = 1e8, Imax = -1e8;
		for (int i = 0; i < N; i++) {
			if (I[i] < Imin) Imin = I[i];
			if (I[i] > Imax) Imax = I[i];
		}
	}

	double t[2][N / 2];
	{
		double zDelta[N];
		bitReverse<N / 2>(zDeltaBR, zDelta);
		bitReverse<N / 2>(zDeltaBR + N / 2, zDelta + N / 2);
		
		for (int i = 0; i < N / 2; i++) {
			t[0][i] = (zDelta[i] / q[0] + I[i]) / 32;
			t[1][i] = (zDelta[i + N / 2] / q[0] + I[i + N / 2]) / 32;
		}

	}
	//---------------------------------------------------------------------
	// CoeffToSlot
	//---------------------------------------------------------------------
	uint64_t ct_cts_hat[2][2][L-MDNUM][N];
	CoeffToSlot<N, logN, L, DNUM, K, MDNUM>(q,p,q[0], E, ckey_hat, rkey_hat, ct_hat, ct_cts_hat,32);
	printf("CoeffToSlot complete\n");
	//-----------------------------------------------
	// debug
	//-----------------------------------------------
	double error_cts_r[2][N / 2];
	double error_cts_i[2][N / 2];
	{
		// dec and decode
		uint64_t pt[2][L - MDNUM][N];
		dec<N, L - MDNUM>(ct_cts_hat[0], s, q, pt[0]);
		dec<N, L - MDNUM>(ct_cts_hat[1], s, q, pt[1]);
		
		double wr[2][N / 2], wi[2][N / 2];
		decode<N, logN, L - MDNUM>(pt[0], q[0], q, wr[0], wi[0]);
		decode<N, logN, L - MDNUM>(pt[1], q[0], q, wr[1], wi[1]);
		
		double zDelta[N];
		bitReverse<N / 2>(zDeltaBR, zDelta);
		bitReverse<N / 2>(zDeltaBR+N/2, zDelta+N/2);
		for (int i = 0; i < N / 2; i++) {
			error_cts_r[0][i] = wr[0][i] * 32 - zDelta[i      ] / q[0] - I[i];
			error_cts_r[1][i] = wr[1][i] * 32 - zDelta[N/2 + i] / q[0] - I[i+N/2];
			error_cts_i[0][i] = wi[0][i];
			error_cts_i[1][i] = wi[1][i];
		}

		printf("error_cts in L2 norm: %e, %e\n", sqrt(norm_square(error_cts_r[0], error_cts_i[0])),
												 sqrt(norm_square(error_cts_r[1], error_cts_i[1])));
		printf("error_cts in L2_exp norm: %e\n", norm_square_exp(error_cts_r[0], error_cts_i[0]));
		printf("error_cts in L2_exp norm: %e\n", norm_square_exp(error_cts_r[1], error_cts_i[1]));
		printf("error_cts in Linf norm: %e, %e\n", norm_max(error_cts_r[0], error_cts_i[0]),
												    norm_max(error_cts_r[1], error_cts_i[1]));
	}
	
	//---------------------------------------------------------------------
	// EvalMod
	//---------------------------------------------------------------------
	uint64_t ct_hat_evalmod[2][2][L-12][N];
	
	EvalMod_h_192_mul_q_over_Delta<N,L-3,DNUM,K>(q,p,q[0],evk_hat,ct_cts_hat[0],ct_hat_evalmod[0],s,t[0]);
	EvalMod_h_192_mul_q_over_Delta<N,L-3,DNUM,K>(q,p,q[0],evk_hat,ct_cts_hat[1],ct_hat_evalmod[1],s,t[1]);
	printf("EvalMod complete\n");
	//-----------------------------------------------
	// debug
	//-----------------------------------------------
	double error_evalmod_r[2][N / 2];
	double error_evalmod_i[2][N / 2];
	{
		uint64_t pt[2][L - 12][N];
		dec<N, L - 12>(ct_hat_evalmod[0], s, q, pt[0]);
		dec<N, L - 12>(ct_hat_evalmod[1], s, q, pt[1]);

		double wr[2][N / 2], wi[2][N / 2];
		decode<N, logN, L - 12>(pt[0], Delta, q, wr[0], wi[0]);
		decode<N, logN, L - 12>(pt[1], Delta, q, wr[1], wi[1]);

		for (int i = 0; i < N / 2; i++) {
			error_evalmod_r[0][i] = wr[0][i] - zr[i];
			error_evalmod_r[1][i] = wr[1][i] - zi[i];
			error_evalmod_i[0][i] = wi[0][i];
			error_evalmod_i[1][i] = wi[1][i];
		}
		printf("error_evalmod in L2 norm: %e, %e\n", sqrt(norm_square(error_evalmod_r[0], error_evalmod_i[0])),
												 sqrt(norm_square(error_evalmod_r[1], error_evalmod_i[1])));
		printf("error_evalmod in L2_exp norm: %e\n", norm_square_exp(error_evalmod_r[0], error_evalmod_i[0]));
		printf("error_evalmod in L2_exp norm: %e\n", norm_square_exp(error_evalmod_r[1], error_evalmod_i[1]));
		//printf("error_evalmod in Linf norm: %e, %e\n", norm_max(error_evalmod_r[0], error_evalmod_r[0]),
		//										   norm_max(error_evalmod_r[1], error_evalmod_r[1]));
	}
	
	//---------------------------------------------------------------------
	// boot
	//---------------------------------------------------------------------
	uint64_t ct_hat_boot[2][L - 12][N];
	{
		uint64_t pti_hat[L - 12][N];
		for (int i = 0; i < L - 12; i++) {
			for (int j = 0; j < N; j++)
				pti_hat[i][j] = 0;
			pti_hat[i][N / 2] = 1;
		}
		ntt<N, L - 12>(q, pti_hat);

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < L - 12; j++)
				for (int k = 0; k < N; k++)
					ct_hat_boot[i][j][k] = (ct_hat_evalmod[0][i][j][k] + mul_mod(pti_hat[j][k], ct_hat_evalmod[1][i][j][k], q[j])) % q[j];
	}

	double er[N / 2], ei[N / 2];
	{
		uint64_t pt[L - 12][N];
		dec<N, L-12>(ct_hat_boot, s, q, pt);

		double wr[N / 2], wi[N / 2];
		decode<N, logN, L - 12>(pt, Delta, q, wr, wi);

		for (int i = 0; i < N / 2; i++) {
			er[i] = wr[i] - zr[i];
			ei[i] = wi[i] - zi[i];
		}

		printf("error_bootstrapping in Linf norm: %e\n", norm_max(er, ei));
		printf("error_bootstrapping in L2 norm: %e\n", sqrt(norm_square(er, ei)));
		printf("error_bootstrapping in L2_exp norm: %e\n", norm_square_exp(er, ei));
	}
	
}
