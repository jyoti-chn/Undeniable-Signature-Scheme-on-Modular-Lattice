#include <stdbool.h>

#include "common.h"

void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp);

void Sign(poly_matrix nu, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len);

bool Verify(poly_matrix nu, uint8_t *r, poly_matrix A, uint8_t *m, int m_len);


void KeyGen2(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp ,poly_matrix h , uint8_t * sd, poly_matrix v ,poly_matrix L   );

void Sign2(poly_matrix sig_2, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len, poly_matrix h , poly_matrix v , poly_matrix sig_1, poly_matrix sig_3, poly_matrix M);
bool Verify2( poly_matrix sig_1 , poly_matrix sig_2, poly_matrix sig_3, uint8_t *r, poly_matrix A, poly_matrix h , uint8_t * sd, uint8_t *m, int m_len , poly_matrix v , poly_matrix L ,poly_matrix M);
bool verifier0(poly_matrix permute_e,poly_matrix permute_v,poly_matrix comm_2,poly_matrix comm_3);
bool verifier1(uint8_t* permute, poly_matrix add_v_e,poly_matrix comm_1, poly_matrix add_L_M, poly_matrix h, poly_matrix sig_3);
bool verifier2(poly_matrix mul_L_M_e,uint8_t* permute,poly_matrix permute_e,poly_matrix comm_1,poly_matrix comm_2);
