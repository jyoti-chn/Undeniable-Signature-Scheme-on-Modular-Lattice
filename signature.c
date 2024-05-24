#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"
#include "random.h"
#include "hash.h"

#include "cpucycles.h"


extern unsigned long long timing_sampleZ_KG;
extern unsigned long long timing_precomp_KG;
extern unsigned long long timing_arith_KG;

/*
	Generates a signing key (T, cplx_T, sch_comp) and an associated verification key A
*/
void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp    )
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	
	//scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];
	scalar *A_hat_coeffs = malloc(PARAM_D * PARAM_D * PARAM_N * sizeof(scalar)), *AprimeT_coeffs = malloc(PARAM_D * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the NTT_op domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D - 1);
	
	
	// T <- D_{R^{2d,dk},sigma}
	begin_timing = cpucycles_start();
	
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	end_timing = cpucycles_stop();
	timing_sampleZ_KG += (end_timing - begin_timing);

	// Compute the Schur complements
	begin_timing = cpucycles_start();
	
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	end_timing = cpucycles_stop();
	timing_precomp_KG += (end_timing - begin_timing);
	
	
	begin_timing = cpucycles_start();
	
	// Add q to each component of T and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		T[i] += PARAM_Q;
		}

	matrix_ntt(T, 2*PARAM_D, PARAM_D * PARAM_K);
	freeze_poly(T, PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K - 1);

	// T1 and T2 are the upper and lower half of T
	poly_matrix T1 = T, T2 = poly_matrix_element(T, PARAM_D * PARAM_K, PARAM_D, 0);
	
	// AprimeT <- A_hat * T2
	mul_crt_poly_matrix(AprimeT, A_hat, T2, PARAM_D, PARAM_D, PARAM_D * PARAM_K);
	
	// Convert AprimeT from the NTT_res domain into the NTT_op domain
	multiply_by_2pow32(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	// AprimeT <- AprimeT + T1
	add_to_poly_matrix(AprimeT, T1, PARAM_D, PARAM_D * PARAM_K);
	
	freeze_poly(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	
	// AprimeT <- - AprimeT
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		AprimeT[i] = 2 * PARAM_Q - AprimeT[i];
		}
	
	
	// AprimeT <- AprimeT + G
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix AprimeT_iik = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, i * PARAM_K);
		
		add_ring_gadget_vector(AprimeT_iik);
		}
	
	
	// A = (A_hat | -A'T) ( = (I | A_hat | -A'T) implicitly)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		poly_matrix A_hat_i = poly_matrix_element(A_hat, PARAM_D, i, 0);
		
		memcpy(A_i0, A_hat_i, PARAM_D * PARAM_N * sizeof(scalar));
		}
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_id = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D);
		poly_matrix AprimeT_i = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, 0);
		
		memcpy(A_id, AprimeT_i, PARAM_N * PARAM_D * PARAM_K * sizeof(scalar));
		}
	
	// Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	
	
	end_timing = cpucycles_stop();
	timing_arith_KG += (end_timing - begin_timing);
	
	free(A_hat);
	free(AprimeT);

	scalar* L = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	uint8_t sd[SALT_BYTES] ; 
	// for(int lp = 0 ; lp<10 ; lp++){
	// 	salt(sd) ;
	// for(uint8_t* j = sd  ; j<sd + SALT_BYTES ; j++   ){
	// 	printf( "%d" , *j ) ;

	// }
	// printf("%c",'\n') ;
	// }
	salt(sd) ;
	// for(uint8_t* j = sd  ; j<sd + SALT_BYTES ; j++   ){
	// 	printf( "%d" , *j ) ;

	// }
	// printf("%c",'\n') ;
	hash_com(L ,sd) ;
	
	scalar* v = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;
	SampleR_matrix_centered((signed_poly_matrix) v, (PARAM_M-PARAM_D), 1,PARAM_SIGMA);

	scalar* h = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;

	multiply_by_A(h , L , v );
	//  for(unsigned int* j = h ; j< h+PARAM_N*PARAM_D; j++  ){
    //     printf( "%d" , *j ) ;
    //   }
    //   printf("%c" , '\n') ;
	

	}
	// void mul(poly_matrix A , poly_matrix B, poly_matrix C){
	// 	int n = 
	// }
	void mul_by_A(poly_matrix y, poly_matrix A, poly_matrix x)
	{
	
	for(int i=0 ; i< PARAM_N ;i++){
		for(int j=0 ; j<PARAM_D ; j++){
			for(int k=0 ; k< PARAM_M-PARAM_D; k++){
				int64_t t = x[i+k*PARAM_N] ;
				t = t*(A[ i+j*PARAM_N*(PARAM_M-PARAM_D) + k*(PARAM_N)   ])  ;
				t = t%PARAM_Q ;
				y[i+j*PARAM_N] +=    t  ;
				y[i+j*PARAM_N]%=PARAM_Q ;

			}

		}
	}
	
	}
void KeyGen2(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp ,poly_matrix h , uint8_t * sd, poly_matrix v ,poly_matrix L   )
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	
	//scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];
	scalar *A_hat_coeffs = malloc(PARAM_D * PARAM_D * PARAM_N * sizeof(scalar)), *AprimeT_coeffs = malloc(PARAM_D * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the NTT_op domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D - 1);
	
	
	// T <- D_{R^{2d,dk},sigma}
	begin_timing = cpucycles_start();
	
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	end_timing = cpucycles_stop();
	timing_sampleZ_KG += (end_timing - begin_timing);

	// Compute the Schur complements
	begin_timing = cpucycles_start();
	
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	end_timing = cpucycles_stop();
	timing_precomp_KG += (end_timing - begin_timing);
	
	
	begin_timing = cpucycles_start();
	
	// Add q to each component of T and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		T[i] += PARAM_Q;
		}

	matrix_ntt(T, 2*PARAM_D, PARAM_D * PARAM_K);
	freeze_poly(T, PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K - 1);

	// T1 and T2 are the upper and lower half of T
	poly_matrix T1 = T, T2 = poly_matrix_element(T, PARAM_D * PARAM_K, PARAM_D, 0);
	
	// AprimeT <- A_hat * T2
	mul_crt_poly_matrix(AprimeT, A_hat, T2, PARAM_D, PARAM_D, PARAM_D * PARAM_K);
	
	// Convert AprimeT from the NTT_res domain into the NTT_op domain
	multiply_by_2pow32(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	// AprimeT <- AprimeT + T1
	add_to_poly_matrix(AprimeT, T1, PARAM_D, PARAM_D * PARAM_K);
	
	freeze_poly(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	
	// AprimeT <- - AprimeT
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		AprimeT[i] = 2 * PARAM_Q - AprimeT[i];
		}
	
	
	// AprimeT <- AprimeT + G
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix AprimeT_iik = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, i * PARAM_K);
		
		add_ring_gadget_vector(AprimeT_iik);
		}
	
	
	// A = (A_hat | -A'T) ( = (I | A_hat | -A'T) implicitly)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		poly_matrix A_hat_i = poly_matrix_element(A_hat, PARAM_D, i, 0);
		
		memcpy(A_i0, A_hat_i, PARAM_D * PARAM_N * sizeof(scalar));
		}

	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_id = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D);
		poly_matrix AprimeT_i = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, 0);
		
		memcpy(A_id, AprimeT_i, PARAM_N * PARAM_D * PARAM_K * sizeof(scalar));
		}
	

	freeze_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	
	
	
	printf("%llu",timing_arith_KG);
	printf("\n");
	free(A_hat);
	free(AprimeT);

	salt(sd) ;

	hash_com(L ,sd) ;
	
	SampleR_matrix_centered((signed_poly_matrix) v, (PARAM_M-PARAM_D), 1,PARAM_SIGMA);
	for(int i=0 ; i<PARAM_N*(PARAM_M-PARAM_D) ; i++){
		v[i] = (v[i]+PARAM_Q)%PARAM_Q ;
		
	}
	
	mul_crt_poly_matrix(h , L , v, PARAM_D,PARAM_M-PARAM_D,1 );

	end_timing = cpucycles_stop();
	timing_arith_KG += (end_timing - begin_timing);
	unsigned long long KeyGen2_cycles = timing_arith_KG+timing_precomp_KG+timing_sampleZ_KG;
	printf("KeyGen cycles: ");
	printf("%llu",KeyGen2_cycles);
	printf("\n");
	
	//printf("Printing H...\n");
	//print_poly_matrix(h, PARAM_D , 1 ) ;
	}
/*
	Signs a message m of length m_len using the signing key (T, cplx_T, sch_comp) and the verification key A
*/
void Sign(poly_matrix nu, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len)
	{
	// Generate a random salt

	salt(r);
	
	// Construct a target u using the message and the salt
	scalar u_coeffs[PARAM_D * PARAM_N];
	poly_matrix u = u_coeffs;
	
	H(u, m, m_len, r);

	// Sample nu
	sample_pre(nu, A, T, cplx_T, sch_comp, u);

	}

void Sign2(poly_matrix sig_2, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len, poly_matrix h , poly_matrix v , poly_matrix sig_1, poly_matrix sig_3, poly_matrix M)
{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	
	begin_timing = cpucycles_start();
	salt(r) ;

	Hash_com_message(M,m,m_len,r);

	salt(r) ;
	
	H(sig_1, (uint8_t*)h , PARAM_N*PARAM_D*sizeof(scalar)  ,r);
	
	sample_pre(sig_2, A, T, cplx_T, sch_comp, sig_1);
	
	for(int i=0 ; i<PARAM_N * PARAM_M ; i++){
		sig_2[i] += PARAM_Q  ;
		sig_2[i] %= PARAM_Q ;
	}
	
	mul_crt_poly_matrix(sig_3 , M , v, PARAM_D,PARAM_M-PARAM_D,1);
	
	end_timing = cpucycles_stop();
	
	unsigned long long Sign_cycles = end_timing-begin_timing;
	printf("Signcycles: ");
	
	printf("%llu",Sign_cycles);
	
	printf("\n");
	

	// printf("Printing Sig_1...\n");
	// print_poly_matrix(sig_1, PARAM_D , 1 ) ;
	// printf("Printing Sig_2...\n");
	// print_poly_matrix(sig_2, 	PARAM_M , 1 ) ;
	// printf("Printing Sig_3...\n");
	// print_poly_matrix(sig_3, PARAM_D , 1 ) ;

}
bool verifier0(poly_matrix permute_e,poly_matrix permute_v,poly_matrix comm_2,poly_matrix comm_3){
		scalar* tmcomm_2 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );
		scalar* tmcomm_3 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );
		H_2(tmcomm_2 , (uint8_t*)permute_e , PARAM_N*(PARAM_M-PARAM_D) *sizeof(scalar) );
		scalar* permute_v_e = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;
		for(int i= 0; i<PARAM_N*(PARAM_M-PARAM_D) ; i++  ){
			permute_v_e[ i ] = permute_e[i]+ permute_v[i] ;
			permute_v_e[ i ]%=PARAM_Q ;
		}

		H_2(tmcomm_3 , (uint8_t*)permute_v_e , PARAM_N*(PARAM_M-PARAM_D)  ) ;

		for(int i=0; i<PARAM_N*PARAM_D; i++){
			if(tmcomm_2[i]!=comm_2[i]) return false;
			if(tmcomm_3[i]!=comm_3[i]) return false;
		}

		return true;

	}
	bool verifier1(uint8_t* permute, poly_matrix add_v_e,poly_matrix comm_1, poly_matrix add_L_M, poly_matrix h, poly_matrix sig_3){
		scalar* tmcomm_1 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );
		scalar* mul_L_M_v_e = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;
		scalar* tmp1 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );
	
		mul_crt_poly_matrix(mul_L_M_v_e, add_L_M , add_v_e, PARAM_D,PARAM_M-PARAM_D,1);
		for(int i=0;i<PARAM_N*PARAM_D;i++){
			tmp1[i] = (mul_L_M_v_e[i] -h[i] + PARAM_Q )%PARAM_Q;
			tmp1[i] +=  PARAM_Q-sig_3[i] ;
			tmp1[i]%= PARAM_Q;
		}

		int mr_len = ((PARAM_M-PARAM_D)*sizeof(int)) + ( PARAM_N*PARAM_D*sizeof(scalar) );
		uint8_t mr[mr_len];
		memcpy(mr,permute,((PARAM_M-PARAM_D)*sizeof(int)));
		memcpy(&mr[((PARAM_M-PARAM_D)*sizeof(int))],tmp1,PARAM_N*PARAM_D*sizeof(scalar) );

		H_2(tmcomm_1,mr,mr_len);

		for(int i=0; i<PARAM_N*PARAM_D; i++){
			if(tmcomm_1[i]!=comm_1[i]) return false;
			
		}
		return true;

	}

	bool verifier2(poly_matrix mul_L_M_e,uint8_t* permute,poly_matrix permute_e,poly_matrix comm_1,poly_matrix comm_2){
		scalar* tmcomm_1 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );
		scalar* tmcomm_2 = malloc( PARAM_N*PARAM_D*sizeof(scalar) );

		int mr_len = ((PARAM_M-PARAM_D)*sizeof(int)) + ( PARAM_N*PARAM_D*sizeof(scalar) );
		uint8_t mr[mr_len];
		memcpy(mr,permute,((PARAM_M-PARAM_D)*sizeof(int)));
		memcpy(&mr[((PARAM_M-PARAM_D)*sizeof(int))],mul_L_M_e,PARAM_N*PARAM_D*sizeof(scalar) );

		H_2(tmcomm_1,mr,mr_len);

		H_2(tmcomm_2 , (uint8_t*)permute_e , PARAM_N*(PARAM_M-PARAM_D) *sizeof(scalar) ) ;

		for(int i=0; i<PARAM_N*PARAM_D; i++){
			if(tmcomm_2[i]!=comm_2[i]) return false;
			if(tmcomm_1[i]!=comm_1[i]) return false;
		}

		return true;

	}
	
bool Verify2( poly_matrix sig_1 , poly_matrix sig_2, poly_matrix sig_3, uint8_t *r, poly_matrix A, poly_matrix h , uint8_t * sd, uint8_t *m, int m_len , poly_matrix v , poly_matrix L ,poly_matrix M)
	{
	
	//double_scalar norm_val_sig_2 = norm_squared_mod( sig_2, PARAM_M );  

	//if( norm_val_sig_2 > PARAM_Q*PARAM_Q*(PARAM_N*PARAM_M) ){
	//		return false;
	//}
	
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	begin_timing = cpucycles_start();
	scalar* temp = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;

	mul_crt_poly_matrix(temp, A , sig_2 , PARAM_D,PARAM_M,1);

	bool flag = true ;

	for(int i=0 ; i< PARAM_D*PARAM_N  ;i++){
		if( temp[i] != sig_1[i]  ){
			//	flag = false;
			  break;	
		}
	}
	
	if(!flag)
		return false;

	scalar* e = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;

	SampleR_matrix_centered((signed_poly_matrix) e, (PARAM_M-PARAM_D), 1,PARAM_SIGMA);

	for(int i=0;i<PARAM_N*(PARAM_M-PARAM_D);i++){
		e[i] += PARAM_Q ;
		e[i] %= PARAM_Q ;
		}

	int rand_per[PARAM_M-PARAM_D] ;

	for(int i=PARAM_M-PARAM_D-1 ; i>=0 ; i--  ){
		rand_per[i] = i ;
	}

	for(int i=PARAM_M-PARAM_D-1 ; i>=0 ; i--  ){
		int q = rand()%( PARAM_M-PARAM_D );
		
		int temp = rand_per[q] ;
		rand_per[q] = rand_per[i] ;
		rand_per[i] = temp ;
	}
	
	//scalar* L = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	
	//salt(sd) ;

	//hash_com(L ,sd) ;

	//scalar* M = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));

	//salt(r) ;

	//Hash_com_message(M,m,m_len,r);

	scalar* add_L_M = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));

	for(int i= 0 ; i<PARAM_N * PARAM_D * (PARAM_M - PARAM_D)  ; i++){
		add_L_M[i] = L[i] ;
	}

	add_to_poly_matrix( add_L_M  , M, PARAM_D , (PARAM_M - PARAM_D)  );
	
	scalar* mul_L_M_e = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;
	
	mul_crt_poly_matrix(mul_L_M_e, add_L_M , e, PARAM_D,PARAM_M-PARAM_D,1 );

	scalar* comm_1 = malloc(PARAM_N*PARAM_D*sizeof(scalar));
	// concat of permute and mul_L_M_e
	int mr_len = ((PARAM_M-PARAM_D)*sizeof(int)) + ( PARAM_N*PARAM_D*sizeof(scalar) );
	uint8_t mr[mr_len];
	memcpy(mr,rand_per,((PARAM_M-PARAM_D)*sizeof(int)));
	memcpy(&mr[((PARAM_M-PARAM_D)*sizeof(int))],mul_L_M_e,PARAM_N*PARAM_D*sizeof(scalar) );

	H_2(comm_1,mr,mr_len);
	
	/////////////////////////////////////////////////////////////

	scalar* permute_e = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;

	for(int i= 0; i<  PARAM_M-PARAM_D ; i+= PARAM_N  ){
		for(int k=0 ; k<PARAM_N ; k++)
		permute_e[ i*PARAM_N+k ] = e[rand_per[i]*PARAM_N+k] ;
	}

	scalar* comm_2 = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;

	H_2(comm_2 , (uint8_t*)permute_e , PARAM_N*(PARAM_M-PARAM_D) *sizeof(scalar) ) ;

	scalar* add_v_e = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;
	
	for(int i= 0 ; i< PARAM_N*(PARAM_M-PARAM_D)  ;i++){
		add_v_e[i] = v[i] ;
	}
	
	add_to_poly_matrix( add_v_e, e , (PARAM_M - PARAM_D) , 1  );
	
	scalar* permute_v = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;
	
	for(int i= 0; i<  PARAM_M-PARAM_D ; i+= PARAM_N  ){
		for(int k=0 ; k<PARAM_N ; k++)
		permute_v[ i*PARAM_N +k] = v[rand_per[i]*PARAM_N+k] ;
	}

	
	scalar* permute_v_e = malloc(  PARAM_N*(PARAM_M-PARAM_D)*sizeof(scalar)) ;

	for(int i= 0; i<  PARAM_M-PARAM_D ; i+= PARAM_N  ){
		for(int k=0 ; k<PARAM_N ; k++)
			permute_v_e[ i*PARAM_N+k ] = add_v_e[rand_per[i]*PARAM_N+k] ;
	}

	scalar* comm_3 = malloc( PARAM_N*PARAM_D*sizeof(scalar) ) ;

	H_2(comm_3 , (uint8_t*)permute_v_e , PARAM_N*(PARAM_M-PARAM_D)  ) ;
	
	int ch = rand()%3;

	
	bool res = true;
	if(ch==0){
		// printf("Ch=0 testing\n") ;
		if(!verifier0(permute_e,permute_v,comm_2,comm_3)) res = false;
		//printf("Ch=0 failed\n") ;
		
	}
	else if(ch==1){
		//printf("Ch=1 testing\n") ;
		if(!verifier1((uint8_t*)rand_per,add_v_e,comm_1,add_L_M,h,sig_3)) res = false;
		//printf("Ch=1 failed\n") ;
		
	}
	else{
		//printf("Ch=2 testing\n") ;
		if(!verifier2(mul_L_M_e,(uint8_t*)rand_per,permute_e,comm_1,comm_2)) res = false;
		//printf("Ch=2 failed\n") ;
		
	}
	
	end_timing = cpucycles_stop();

	unsigned long long verify_cycles = end_timing-begin_timing;
	printf("Verify_cycles: ");
	printf("%llu",verify_cycles);
	printf("\n");
	
	if(res) return true;
	return false;
	
	}

	
	

/*
	Checks is the signature nu is valid for the message m, given the verification key A
*/
bool Verify(poly_matrix nu, uint8_t *r, poly_matrix A, uint8_t *m, int m_len)
	{
	// Construct a target u using the message and the salt
	scalar u_coeffs[PARAM_D * PARAM_N];
	poly_matrix u = u_coeffs;
	
	H(u, m, m_len, r);
	
	
	// Verify that A * nu = u mod q
	scalar prod_coeffs [ PARAM_N * PARAM_D ];
	poly_matrix prod = prod_coeffs;
	
	multiply_by_A(prod, A, (poly_matrix) nu);
	
	
	
	// Verify that nu has a small norm
	divide_by_2pow32(nu, PARAM_N * PARAM_M - 1);
	matrix_invntt(nu, PARAM_M, 1);
	freeze_poly(nu, PARAM_N * PARAM_M - 1);
	
	double_scalar norm_nu_squared = norm_squared((poly_matrix) nu, PARAM_M);
	double bound_squared = PARAM_T * PARAM_T * PARAM_ZETA * PARAM_ZETA * PARAM_N * PARAM_M;
	
	
	matrix_ntt(nu, PARAM_M, 1);
	
	if(!equals_poly(prod, u, PARAM_N * PARAM_D - 1) || (norm_nu_squared >= bound_squared))
		{
		return false;
		}
	
	return true;
	}
