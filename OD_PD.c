#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <time.h>
#include <complex.h>
#include "taylor.h"

#define JJ _NUMBER_OF_JET_VARS_

#define uround 1e-16

#define OD_zero_TOL uround
#define EV_TOL OD_zero_TOL

#define SEC_IND 1

#define FORCE_JAC_PROJ 0

#define store_ALL  1
#define store_EVAL ( 1 || store_ALL)
#define store_EVEC ( 1 || store_ALL)
#define store_R    ( 1 || store_ALL)
#define store_Q    ( 1 || store_ALL)

#define verbose         0
#define verbose_none    0
#define verbose_A     ( 0 || verbose ) && !verbose_none
#define verbose_R     ( 0 || verbose ) && !verbose_none
#define verbose_Q     ( 0 || verbose ) && !verbose_none
#define verbose_QR    ( 0 || verbose ) && !verbose_none
#define verbose_perm  ( 1 || verbose ) && !verbose_none
#define verbose_wEVEC ( 0 || verbose ) && !verbose_none
#define verbose_EVEC  ( 0 || verbose ) && !verbose_none

MY_FLOAT s=77.27, w=0.161, q=8.375e-6, fff=1.0;

// 1e+5
// EIGENVALUES:

// lambda_1 = 1e-0.00000000000073989 
// lambda_2 = 1e-66.02752086534367493 
// lambda_3 = 1e-3244180.17328815767541528 


// EIGENVECTORS:

// w_2: 5.369091681891230e-01 8.436400566085421e-01 0.000000000000000e+00 , ||w_2||=1e66.027521, (3 iterations)
// EV_2: -8.281628667595103e-02 8.281753723517560e-02 -9.931176759012535e-01 
// Execution time: 6.057112 seconds


// 1e+2
// EIGENVALUES:

// lambda_1 = 1e0.0000000000 
// lambda_2 = 1e-66.0275208658 
// lambda_3 = 1e-3244180.1732873865 


// EIGENVECTORS:

// w_2: 5.369091681890859e-01 8.436400566085659e-01 0.000000000000000e+00 , ||w_2||=1e66.027521, (3 iterations)
// EV_2: -8.281628667593544e-02 8.281753723517477e-02 -9.931176759012550e-01 
// Execution time: 33.491480 seconds

double sign(double x){

    if (x<0.0) return -1.0;
    else return 1.0;

}

// void Householder(double *v, double *A, int rows, int k){

//     double alpha=0, norm=0;

//     for (i=s; i<S_N; i++) alpha+=A[m][i*JJ+t]*A[m][i*JJ+t];
//     alpha = sign(A[m][s*JJ+t])*sqrt(alpha);

//     for (i=0; i<S_N-s; i++) {
//         v[i] = A[m][(i+s)*JJ+t];
//         if (i==0) v[i]-=alpha;
//         norm+=v[i]*v[i];
//     }
//     norm=sqrt(norm);

//     for (i=0; i<S_N-s; i++) v[i]/=norm; 

// }

void PHTD(double **A, double **Q, int od){

    if (verbose_R || verbose_Q) printf("\nPHTD:\n\n"); 

    double alpha, norm, v[JJ], v_Q[JJ], M[JJ*JJ], R[JJ*JJ];

    int i,j;
    int m, s=0, t=0, S_0, S_N;

    int HH_ind=0;
    

    double temp_1, temp_2;
    //Applying the permutation matrix P·A_od and A_1·P^T (also apply it to the corresponding Q)

    if (verbose_perm){
        // printf("Permutation of cols of first matrix and rows of last matrix:\n");
        // printf("P:\n");
        // for (i=0; i<JJ; i++){
        //     for (j=0; j<JJ; j++){
        //         printf("%le ", Q[0][j+i*JJ]);
        //     }
        //     printf("\n");
        // }
        printf("A_1:\n");
        for (i=0; i<JJ; i++){
            for (j=0; j<JJ; j++){
                printf("%le ", A[0][j+i*JJ]);
            }
            printf("\n");
        }
        printf("A_%d:\n", od);
        for (i=0; i<JJ; i++){
            for (j=0; j<JJ; j++){
                printf("%le ", A[od-1][j+i*JJ]);
            }
            printf("\n");
        }
        printf("\n");
    }
    

    /*//Cyclic permutation 
    double temp_A[JJ*JJ], temp_Q[JJ*JJ];    
    for (i=0; i<JJ; i++) for (j=0; j<JJ; j++){ temp_A[j+i*JJ] = A[0][j+i*JJ]; temp_Q[j+i*JJ] = Q[0][j+i*JJ];}
    for (i=0; i<JJ; ++i){
        for (j=0; j<JJ; ++j){
            A[0][((JJ-1-SEC_IND)+j)%JJ+i*JJ] = temp_A[j+i*JJ];
            Q[0][((JJ-1-SEC_IND)+j)%JJ+i*JJ] = temp_Q[j+i*JJ];
        }
    }

    for (i=0; i<JJ; i++) for (j=0; j<JJ; j++){ temp_A[j+i*JJ] = A[od-1][j+i*JJ];}
    for (i=0; i<JJ; ++i){
        for (j=0; j<JJ; ++j){
            A[od-1][j+(((JJ-1-SEC_IND)+i)%JJ)*JJ] = temp_A[j+i*JJ];
        }
    }*/

    
    //Send SEC_IND to last 
    double temp_A, temp_Q;    
    for (i=0; i<JJ; i++){
        temp_A=A[0][SEC_IND+i*JJ];
        temp_Q=Q[0][SEC_IND+i*JJ];
        for (j=SEC_IND; j<JJ-1; j++){
            A[0][j+i*JJ] = A[0][(j+1)+i*JJ];
            // Q[0][j+i*JJ] = Q[0][(j+1)+i*JJ];
        }
        A[0][(JJ-1)+i*JJ] = temp_A;
        Q[0][(JJ-1)+i*JJ] = temp_Q;
    }

    for (j=0; j<JJ; j++){
        temp_A=A[od-1][j+SEC_IND*JJ];
        for (i=SEC_IND; i<JJ-1; i++){
            A[od-1][j+i*JJ] = A[od-1][j+(i+1)*JJ];
        }
        A[od-1][j+(JJ-1)*JJ] = temp_A;
    }

    if (verbose_perm){
        // printf("Permutation of cols of first matrix and rows of last matrix:\n");
        // printf("P:\n");
        // for (i=0; i<JJ; i++){
        //     for (j=0; j<JJ; j++){
        //         printf("%le ", Q[0][j+i*JJ]);
        //     }
        //     printf("\n");
        // }
        printf("A_1:\n");
        for (i=0; i<JJ; i++){
            for (j=0; j<JJ-1; j++){
                printf("%.17le ", A[0][j+i*JJ]);
            }
            printf("\n");
        }
        printf("A_%d:\n", od);
        for (i=0; i<JJ-1; i++){
            for (j=0; j<JJ; j++){
                printf("%.17le ", A[od-1][j+i*JJ]);
            }
            printf("\n");
        }
        printf("\n");
    }

    m=(HH_ind+1)%od;
    do{

        S_0=JJ; S_N=JJ;
        if (m==od-1)   S_N=JJ-1;
        else if (m==0) S_0=JJ-1;

        // printf("\n1. m=%d, s=%d, t=%d\n", m,s,t);

        alpha=0; norm=0;

        for (i=s; i<S_N; i++) alpha+=A[m][i*JJ+t]*A[m][i*JJ+t];
        alpha = sign(A[m][s*JJ+t])*sqrt(alpha);

        for (i=0; i<S_N-s; i++) {
            v[i] = A[m][(i+s)*JJ+t];
            if (i==0) v[i]-=alpha;
            norm+=v[i]*v[i];
        }
        norm=sqrt(norm);

        for (i=0; i<S_N-s; i++) v[i]/=norm; 

        // printf("alpha=%le, v=", alpha);
        // for (i=0; i<S_N-s; i++) printf("%.17le ", v[i]);
        // printf("\n");

        /*//Applying explicitly Q_{m+1} to Q (this is not backwards-stable ?? )
        //{
        //BUILDING Q_{m+1}^T from v
        for (i=0; i<S_N; i++){
            if (i<s) v_Q[i]=0.0;
            else v_Q[i]=v[i-s];
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, S_N, S_N, 1, -2.0, v_Q, 1, v_Q, 1, 0.0, M, JJ);
        for (i=0; i<S_N; i++){ M[i+i*JJ]+=1.0; }
        
        if (m+1==od){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, S_N, S_N, S_N, 1.0, Q[0], JJ, M, JJ, 0.0, R, JJ);
            for (i=0; i<S_N; i++){ for (j=0; j<S_N; j++){ Q[0][j+i*JJ]=R[j+i*JJ]; }}

            if (verbose_Q){
                printf("Q_%d\n", 1);
                for (i=0; i<S_N; i++){ for (j=0; j<S_N; j++){ printf("%le ", Q[0][j+i*JJ]); } printf("\n");}
            }
        }
        else{
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, JJ, JJ, JJ, 1.0, Q[m+1], JJ, M, JJ, 0.0, R, JJ);
            for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ Q[m+1][j+i*JJ]=R[j+i*JJ]; }}

            if (verbose_Q){
                printf("Q_%d\n", m+2);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", Q[m+1][j+i*JJ]); } printf("\n");}
            }
        }
        //}*/


        //Applying implicitly Q_{m+1}^T A_m and A_{m+1} Q_{m+1} (without using the matrix M. This is backwards-stable)
        for (i=0; i<S_N-s; i++){
            for (j=0; j<S_0; j++){
                M[j+i*JJ]=A[m][j+(i+s)*JJ];
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, S_0, S_N-s, 1.0, v, S_N-s, M, JJ, 0.0, R, JJ);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, S_N-s, S_0, 1, 2.0, v, 1, R, JJ, 0.0, M, JJ);


        for (i=0; i<S_N; i++){
            for (j=0; j<S_0; j++){
                if (i==s && j==t) A[m][j+i*JJ] = alpha;
                else if (i>s && j==t) A[m][j+i*JJ] = 0.0;
                else if (i>=s) A[m][j+i*JJ] -= M[j+(i-s)*JJ];
            }
        }

        if (verbose_R){
            printf("R_%d\n", m+1);
            for (i=0; i<S_N; i++){ for (j=0; j<S_0; j++){ printf("%.17le ", A[m][j+i*JJ]); } printf("\n"); }
        }

        if (m==HH_ind) t++;

        m++;
        if (m==od) m=0;

        // printf("\n2. m=%d, s=%d, t=%d\n", m,s,t);

        S_0=JJ; S_N=JJ;
        if (m==od-1)   S_N=JJ-1;
        else if (m==0) S_0=JJ-1;

        for (i=0; i<S_N; i++){
            for (j=0; j<S_0-s; j++){
                M[j+i*JJ]=A[m][(j+s)+i*JJ];
            }
        }

        /*cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, JJ, JJ, JJ, 1.0, A[m], JJ, M, JJ, 0.0, R, JJ);*/

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, S_N, 1, S_0-s, 1.0, M, JJ, v, 1, 0.0, R, 1);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, S_N, S_0-s, 1, 2.0, R, 1, v, S_0-s, 0.0, M, JJ);

        for (i=0; i<S_N; i++){
            for (j=0; j<S_0-s; j++){
                A[m][(j+s)+i*JJ] -= M[j+i*JJ];
            }
        }

        if (verbose_R){
            printf("R_%d\n", m+1);
            for (i=0; i<S_N; i++){ for (j=0; j<S_0; j++){ printf("%.17le ", A[m][j+i*JJ]); } printf("\n"); }
        }

        if (m==HH_ind) s++;

        if (m==od-1 && s==JJ-2) break;

    } while (s<JJ-1);

}

void Givens_coef(double *c, double *s, double *r, double a, double b){

    double abs_a=fabs(a), abs_b=fabs(b);
    double temp_1, temp_2;

    if (abs_b<OD_zero_TOL){
        if (abs_a<OD_zero_TOL) *c=1.0;
        else *c=a/abs_a;
        *s = 0;
        *r = abs_a;
    } else if (abs_b<OD_zero_TOL){
        *c = 0;
        *s = -b/abs_b;
        *r = abs_b;
    } else if (abs_a > abs_b){
        temp_1 = b/a;
        temp_2 = a/abs_a * sqrt(1+temp_1*temp_1);
        *c = 1/temp_2;
        *s = -temp_1/temp_2;
        *r = a * temp_2;
    } else {
        temp_1 = a/b;
        temp_2 = b/abs_b * sqrt(1+temp_1*temp_1);
        *s = -1/temp_2;
        *c = temp_1/temp_2;
        *r = b*temp_2;
    }

}

int search_PDEF(int *ind, int *pos, double **A, int OD){

    int od, jj;

    for (jj=0; jj<JJ-1; jj++){
        if (fabs(A[0][jj+(jj+1)*JJ])<OD_zero_TOL){
            *ind=0;
            *pos=jj;
            return 1;
        } else if (fabs(A[0][jj+(jj+1)*JJ])<uround) printf("WARNING. Possible deflation for ind %d in pos %d.\nWith OD_zero_TOL=%.15le, the value %le was rejected.\n", 0, jj, OD_zero_TOL, fabs(A[0][jj+(jj+1)*JJ]));
    }

    for (od=1; od<OD-1; od++){
        for (jj=0; jj<JJ; jj++){
            if (fabs(A[od][jj+jj*JJ])<OD_zero_TOL){
                *ind=od;
                *pos=jj;
                return 1;
            } else if (fabs(A[od][jj+jj*JJ])<uround) printf("WARNING. Possible deflation for ind %d in pos %d.\nWith OD_zero_TOL=%.15le, the value %le was rejected.\n", od, jj, OD_zero_TOL, fabs(A[od][jj+jj*JJ]));
        }
    }
    for (jj=0; jj<JJ-1; jj++){
        if (fabs(A[od][jj+jj*JJ])<OD_zero_TOL){
            *ind=od;
            *pos=jj;
            return 1;
        } else if (fabs(A[OD-1][jj+jj*JJ])<uround) printf("WARNING. Possible deflation for ind %d in pos %d.\nWith OD_zero_TOL=%.15le, the value %le was rejected.\n", od, jj, OD_zero_TOL, fabs(A[OD-1][jj+jj*JJ]));
    }

    return 0;

}

int PDEF(double **A, double **Q, int od){

    if (verbose_R || verbose_Q) printf("\n\nPDEF:\n\n"); 

    int i,j, m, m_next;
    int S_0, S_N=JJ-1;

    double G_c, G_s, G_r, G_a, G_b;
    double G_temp_A_1, G_temp_A_2, G_temp_Q_1, G_temp_Q_2;

    for (m=0; m<od; m++){

        if (m==od-1) S_N=JJ-2;

        for (S_0=0; S_0<S_N; S_0++){

            m_next=(m+1)%od;

            Givens_coef (&G_c, &G_s, &G_r, A[m][S_0+S_0*JJ], A[m][S_0+(S_0+1)*JJ]);

            A[m][S_0+S_0*JJ]  = G_r;
            A[m][S_0+(S_0+1)*JJ]  = 0.0;
            for (j=0; j<JJ; j++){
                if (j!=S_0){
                    G_temp_A_1 = A[m][j+S_0*JJ];
                    G_temp_A_2 = A[m][j+(S_0+1)*JJ];
                    A[m][j+S_0*JJ]  = G_c*G_temp_A_1-G_s*G_temp_A_2;
                    A[m][j+(S_0+1)*JJ]  = G_s*G_temp_A_1+G_c*G_temp_A_2;
                }
            }

            for (i=0; i<JJ; i++){  
                G_temp_A_1 = A[m_next][S_0+i*JJ];
                G_temp_A_2 = A[m_next][(S_0+1)+i*JJ];
                A[m_next][S_0+i*JJ] = G_c*G_temp_A_1-G_s*G_temp_A_2;
                A[m_next][(S_0+1)+i*JJ]  = G_s*G_temp_A_1+G_c*G_temp_A_2;
            }

            //Only for storing Q
            //{
            for (i=0; i<JJ; i++){  
                G_temp_Q_1 = Q[m_next][S_0+i*JJ];
                G_temp_Q_2 = Q[m_next][(S_0+1)+i*JJ];
                Q[m_next][S_0+i*JJ] = G_c*G_temp_Q_1-G_s*G_temp_Q_2;
                Q[m_next][(S_0+1)+i*JJ]  = G_s*G_temp_Q_1+G_c*G_temp_Q_2;
            }
            //}

            if (verbose_R){
                printf("R_%d\n", m+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", A[m][j+i*JJ]); } printf("\n"); }
            }

            if (verbose_Q){
                printf("Q_%d\n", m_next+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", Q[m_next][j+i*JJ]); } printf("\n");}
            }
            
            if (verbose_R){
                printf("R_%d\n", m_next+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", A[m_next][j+i*JJ]); } printf("\n"); }
            }

        }

    }

}

int PSD(double **A, double **Q, int od, int S_0, int S_N){

    if (verbose_R || verbose_Q) printf("\n\nPSD:\n\n"); 

    int i, j, m=0;

    double shift, G_c, G_s, G_r, G_a, G_b;

    double stop_condition[S_N-S_0-1];
    int stop=1;

    while (stop){

        // Use for shift = "last diagonal element"
        /*shift = A[0][(S_N-1)+(S_N-1)*JJ];
        G_a = shift;

        for (i=1; i<od; i++){
            shift*=A[i][(S_N-1)+(S_N-1)*JJ];
            G_a*=(A[i][(S_N-1)+(S_N-1)*JJ]/A[i][S_0+S_0*JJ]);
        } */

        //Use for shift = 0
        shift = 0.0;
        G_a = 0.0;

        /*//Use for shift = 1
        if (S_0+S_N == JJ){
            shift = 1.0;
            G_a = shift;

            for (i=1; i<od; i++){
                G_a/=A[i][S_0+S_0*JJ];
            } 
        } else {
            //Only use shift = 1 for the full matrix
            //Then use this for shift = 0
            shift = 0.0;
            G_a = 0.0;

            // //Then use this for shift = "first diagonal element"
            // shift = A[0][S_0+S_0*JJ];
            // G_a = A[0][S_0+S_0*JJ];

            // for (i=1; i<od; i++){
            //     shift*=A[i][S_0+S_0*JJ];
            //     // G_a*=(A[i][(S_N-2)+(S_N-2)*JJ]/A[i][S_0+S_0*JJ]);
            // } 
        }*/
        
        //The rest is common for all shifts

        G_a = A[0][S_0+S_0*JJ] - G_a;
        G_b = A[0][S_0+(S_0+1)*JJ];

        // printf("Shift: %le, G_a: %le, G_b:% le\n", shift, G_a, G_b);

        Givens_coef (&G_c, &G_s, &G_r,G_a, G_b);

        // printf("--> G_c: %le, G_s: %le, G_r:% le\n", G_c, G_s, G_r);

        double G_temp_A_1, G_temp_A_2, G_temp_Q_1, G_temp_Q_2;

        do{

            // printf("\nm=%d\n", m);

            A[m][S_0+S_0*JJ]  = G_r;
            A[m][S_0+(S_0+1)*JJ]  = 0.0;
            for (j=0; j<JJ; j++){
                if (j!=S_0){
                    G_temp_A_1 = A[m][j+S_0*JJ];
                    G_temp_A_2 = A[m][j+(S_0+1)*JJ];
                    A[m][j+S_0*JJ]  = G_c*G_temp_A_1-G_s*G_temp_A_2;
                    A[m][j+(S_0+1)*JJ]  = G_s*G_temp_A_1+G_c*G_temp_A_2;
                }
            }

            if (verbose_R){
                printf("R_%d\n", m+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", A[m][j+i*JJ]); } printf("\n"); }
            }

            m++;
            if (m==od) m=0;

            for (i=0; i<JJ; i++){  
                G_temp_A_1 = A[m][S_0+i*JJ];
                G_temp_A_2 = A[m][(S_0+1)+i*JJ];
                A[m][S_0+i*JJ] = G_c*G_temp_A_1-G_s*G_temp_A_2;
                A[m][(S_0+1)+i*JJ]  = G_s*G_temp_A_1+G_c*G_temp_A_2;
            }

            //Only for storing Q
            //{
            for (i=0; i<JJ; i++){  
                G_temp_Q_1 = Q[m][S_0+i*JJ];
                G_temp_Q_2 = Q[m][(S_0+1)+i*JJ];
                Q[m][S_0+i*JJ] = G_c*G_temp_Q_1-G_s*G_temp_Q_2;
                Q[m][(S_0+1)+i*JJ]  = G_s*G_temp_Q_1+G_c*G_temp_Q_2;
            }
            //}

            if (verbose_Q){
                printf("Q_%d\n", m+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", Q[m][j+i*JJ]); } printf("\n");}
            }
            
            if (verbose_R){
                printf("R_%d\n", m+1);
                for (i=0; i<JJ; i++){ for (j=0; j<JJ; j++){ printf("%le ", A[m][j+i*JJ]); } printf("\n"); }
            }

            // printf("G_a: %le, G_b:% le\n", A[m][S_0+S_0*JJ], A[m][S_0+(S_0+1)*JJ]);

            Givens_coef (&G_c, &G_s, &G_r, A[m][S_0+S_0*JJ], A[m][S_0+(S_0+1)*JJ]);

            // printf("--> G_c: %le, G_s: %le, G_r:% le\n", G_c, G_s, G_r);

        }while(m!=0);

        for (i=S_0; i<S_N-1; i++){
            stop_condition[i-S_0] = fabs(A[0][i+(i+1)*JJ]);
        }

        stop = (stop_condition[0] > OD_zero_TOL);
        for (i=1; i<S_N-S_0-1; i++){
            stop = stop && (stop_condition[i] > OD_zero_TOL); 
        }

        for (i=0; i<S_N-S_0-1; i++){
            stop = (stop_condition[i] > OD_zero_TOL); 
            if (stop == 0) return i; 
        }

    }


}


int main(int argc, char *argv[]){

    int jj, cc, i, od;
    int OD;

    char OD_num_name[] = "ORBIT DECOMP/num_OD.txt", OD_A_name[] = "ORBIT DECOMP/A.bin", OD_EVAL_name[] = "ORBIT DECOMP/EVAL.bin", OD_EVEC_name[] = "ORBIT DECOMP/EVEC.bin", OD_R_name[] = "ORBIT DECOMP/R.bin", OD_Q_name[] = "ORBIT DECOMP/Q.bin";
    FILE *OD_num_txt, *OD_A_bin, *OD_EVAL_bin, *OD_EVEC_bin, *OD_R_bin, *OD_Q_bin;

    OD_num_txt = fopen(OD_num_name, "r");
    if (OD_num_txt==NULL) exit(0);
    fscanf(OD_num_txt, "%d", &OD);

    OD_A_bin = fopen(OD_A_name, "rb");
    if (OD_A_bin==NULL) exit(0);

    double **OD_JAC = malloc(OD*sizeof(double *));
    for (od=0; od<OD; od++){
        OD_JAC[od] = malloc(JJ*JJ*sizeof(double));
        fread(OD_JAC[od], sizeof(double), JJ*JJ, OD_A_bin);
    }

    if (FORCE_JAC_PROJ){
        for (jj=0; jj<JJ; jj++){
            OD_JAC[0][SEC_IND+jj*JJ] = 0.0;
            OD_JAC[OD-1][jj+SEC_IND*JJ] = 0.0;
        }
    }

    if (verbose_A){  
        for (od=0; od<OD; od++){
            printf("A_%d:\n", od+1);
            for (jj=0; jj<JJ; jj++){
                for (cc=0; cc<JJ; cc++){
                    printf("%.17le ", OD_JAC[od][cc+jj*JJ]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    clock_t CLOCK_begin = clock();

    if (store_EVAL) OD_EVAL_bin = fopen(OD_EVAL_name, "w");
    if (store_EVEC) OD_EVEC_bin = fopen(OD_EVEC_name, "w");
    if (store_R)    OD_R_bin = fopen(OD_R_name, "w");
    if (store_Q)    OD_Q_bin= fopen(OD_Q_name, "w");

    double **OD_Q;

    //Perform the periodic Schur decomposition (we assume that all eigenvalues are real, so all matrices are diagonal)
    OD_Q = malloc(OD*sizeof(double *));
    for (od=0; od<OD; od++) {
        OD_Q[od] = malloc(JJ*JJ*sizeof(double));
        for (jj=0; jj<JJ; jj++){
            for (cc=0; cc<JJ; cc++){
                if (jj==cc) OD_Q[od][cc+jj*JJ] = 1.0;
                else  OD_Q[od][cc+jj*JJ] = 0.0;
            }
        }
    }
    //Make the Q_0 a (JJ-1)x(JJ-1) matrix
    OD_Q[0][SEC_IND+SEC_IND*JJ]=0.0;

    PHTD(OD_JAC, OD_Q, OD);

    int found_PDEF, PDEF_ind, PDEF_pos;
    found_PDEF = search_PDEF(&PDEF_ind, &PDEF_pos, OD_JAC, OD);
    printf("found PDEF? ");

    if (found_PDEF){
        printf("Yes --> ind %d, pos %d\n", PDEF_ind, PDEF_pos);
        if (PDEF_ind==0){
            if (PDEF_pos==0) PSD(OD_JAC, OD_Q, OD, 1, JJ);
            if (PDEF_pos==1) PSD(OD_JAC, OD_Q, OD, 0, JJ-1);
        } else {
            if (PDEF_pos==JJ-1) PDEF(OD_JAC, OD_Q, PDEF_ind+1);
            else exit(0);
        }
    }else{
        printf("No\n");
        int PSD_partition;
        PSD_partition = PSD(OD_JAC, OD_Q, OD, 0, JJ);
        if (PSD_partition == 0) PSD(OD_JAC, OD_Q, OD, 1, JJ);
        else if (PSD_partition == 1) PSD(OD_JAC, OD_Q, OD, 0, JJ-1);
    } 


    // PDEF(OD_JAC, OD_Q, OD);

    // PSD(OD_JAC, OD_Q, OD, 0, JJ-1);

    

    if (verbose_QR){
        printf("\n\n");
        for (od=0; od<OD; od++){
            printf("Q_%d\n", od+1);
            for (jj=0; jj<JJ; jj++){
                for (cc=0; cc<JJ; cc++){
                    printf("%.17le, ", OD_Q[od][cc+jj*JJ]);
                }
                printf("\n");
            }
        }
        printf("\n");
        for (od=0; od<OD; od++){
            printf("R_%d\n", od+1);
            for (jj=0; jj<JJ; jj++){
                for (cc=0; cc<JJ; cc++){
                    printf("%.17le, ", OD_JAC[od][cc+jj*JJ]);
                }
                printf("\n");
            }
        }
    }

    //Compute the eigenvalues (in fact, as a power of 10, adding the log10's of the factors)
    printf("\nEIGENVALUES:\n\n");

    double WR[JJ], W_sign[JJ], W_temp;
    int W_zero_factors;
    for (jj=0; jj<JJ; jj++) W_sign[jj]=1.0;

    for (jj=0; jj<JJ; jj++){
        WR[jj]=0.0;
        W_zero_factors=0;
        for (od=0; od<OD; od++){
            W_temp = OD_JAC[od][jj+jj*JJ];
            if (fabs(W_temp)>1e-200)  WR[jj]+=log10(fabs(W_temp));
            else W_zero_factors++; 

            if (W_temp<0) W_sign[jj]*=-1;
        } 

        if (W_zero_factors==0){
            if (W_sign[jj]>0) printf("lambda_%d = 1e%.17lf \n", jj+1, WR[jj]);
            else printf("lambda_%d = -1e%.17lf \n", jj+1, WR[jj]);
        } else {
            printf("lambda_%d = 0", jj+1);
            if (jj==JJ-1) {
                printf("\t(lambda_%d has %d zero factors. Without them, it would be ", jj+1, W_zero_factors);
                if (W_sign[jj]>0) printf("lambda_%d=-1e%.17lf )\n", jj+1, WR[jj]);
                else printf("lambda_%d=1e%.17lf )\n", jj+1, WR[jj]);
            }
            else printf("WARNING, lambda_%d has %d zero factors (discarded)\n", jj+1, W_zero_factors);
        }

    }




    //Compute the eigenvectors
    printf("\n\nEIGENVECTORS:\n\n"); 

    double **EV = malloc(OD*sizeof(double *));
    for (od=0; od<OD; od++){ EV[od] = malloc(JJ*sizeof(double)); }

    double /*EV_last[JJ],*/ EV_temp[JJ], EV_temp_norm, EV_exp_norm;
    int EV_ind=1, it_EV=0, EV_stop_next_it=0;

    // for (jj=0; jj<JJ; jj++) {
    //     if (jj<=EV_ind) EV[0][jj]=1.0;
    //     else EV[0][jj]=0.0;
    // }

    // if (verbose_wEVEC || 1){
    //     printf("w_1^(0)= ");
    //     for (jj=0; jj<JJ; jj++) printf("%le ", EV[0][jj]);
    //     printf("\n");
    // }

    // lapack_int LAPACK_info, *LAPACK_ipiv;
    // LAPACK_ipiv = (lapack_int *)malloc(JJ*sizeof(lapack_int)) ;

    // double temp_JAC[JJ*JJ];
    // do{
    //     // printf("it_EV=%d\n", it_EV);
    //     EV_exp_norm=0.0;

    //     od=0;
    //     do{

    //         for (jj=0; jj<JJ; jj++) EV_temp[jj]=EV[od][jj];
    //         for (jj=0; jj<JJ; jj++) for (cc=0; cc<JJ; cc++) temp_JAC[cc+jj*JJ] = OD_JAC[od][cc+jj*JJ];
    //         if (od==OD-1){
    //             LAPACK_info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', JJ-1, JJ, 1, temp_JAC, JJ, EV_temp, 1);
    //             if (LAPACK_info != 0) { printf("LAPACKE_dtrtrs error %d. (At %d-eigenvector computation).\n", LAPACK_info, od+1); exit(0); } 
    //         } else if (od==0){
    //             LAPACK_info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', JJ-1, JJ-1, 1, temp_JAC, JJ, EV_temp, 1);
    //             if (LAPACK_info != 0) { printf("LAPACKE_dtrtrs error %d. (At %d-eigenvector computation).\n", LAPACK_info, od+1); exit(0); }
    //             EV_temp[JJ-1]=0.0; //Remove artifact from LAPACK_dgels output
    //         } else {
    //             LAPACK_info = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 'U', 'N', 'N', JJ, 1, temp_JAC, JJ, EV_temp, 1);
    //             if (LAPACK_info != 0) { printf("LAPACKE_dtrtrs error %d. (At %d-eigenvector computation).\n", LAPACK_info, od+1); exit(0); } 
    //         }

    //         // for (jj=EV_ind+1; jj<JJ; jj++){
    //         //     if (fabs(EV_temp[jj])>1e-200){
    //         //         printf("Artifact removed w_%d^(%d)=", (od+2)%OD, it_EV);
    //         //         for (cc=0; cc<JJ; cc++) printf("%le ", EV_temp[cc]);
    //         //         printf("\n");

    //         //         EV_temp[jj]=0.0;

    //         //         printf("\t     --> w_%d^(%d)=", (od+2)%OD, it_EV);
    //         //         for (cc=0; cc<JJ; cc++) printf("%le ", EV_temp[cc]);
    //         //         printf("\n");
    //         //     } 
                
    //         // }

    //         od++;
    //         if (od==OD) od=0;

    //         EV_temp_norm = 0.0;
    //         for (jj=0; jj<JJ; jj++) EV_temp_norm+=EV_temp[jj]*EV_temp[jj];
    //         EV_temp_norm = sqrt(EV_temp_norm);

    //         for (jj=0; jj<JJ; jj++) EV[od][jj]=EV_temp[jj]/EV_temp_norm;

    //         EV_exp_norm += log10(EV_temp_norm);

    //         if (verbose_wEVEC || (od-1+OD)%OD==0 || (od-1+OD)%OD==OD-1 || (od-1+OD)%OD==OD-2){
    //             printf("w_%d^(%d)= ", od+1, it_EV);
    //             for (jj=0; jj<JJ; jj++) printf("%le, ", EV[od][jj]);
    //             printf("\n");
    //         }


    //     }while (od!=0);

    //     if (it_EV>0){
    //         EV_error = 0.0;
    //         for (jj=0; jj<JJ; jj++){ 
    //             EV_error+=fabs(sign(EV[0][0])*EV[0][jj]-sign(EV_last[0])*EV_last[jj]); //To avoid errors, we must consider the case of vectors w/ opposite sign 
    //         }
                        
    //         if (EV_error<EV_TOL) EV_stop_next_it++;
    //         // printf("error E_%d=%le, it_EV=%d, EV_stop_next_it=%d\n", EV_ind+1, EV_error, it_EV, EV_stop_next_it);
    //     }
    //     for (jj=0; jj<JJ; jj++){ EV_last[jj]=EV[0][jj]; }

    //     if (EV_stop_next_it>1) EV_stop_next_it++;
    //     it_EV++;

    // }while(EV_stop_next_it<2);

    // free(LAPACK_ipiv);






    double **EV_last = malloc(OD*sizeof(double *));
    for (od=0; od<OD; od++){
        EV_last[od] = malloc(JJ*sizeof(double));
        for (jj=0; jj<JJ; jj++){
            /*if (jj<=EV_ind)*/ EV[od][jj] = 1.0;
            //else EV[od][jj] = 0.0;

            EV_last[od][jj] = EV[od][jj];
        }
    } 
    EV[0][JJ-1]=0.0; EV_last[0][JJ-1]=0.0;

    double EV_error, EV_total_difference_error, EV_total_relative_diff_error, EV_total_residual_error;
    int EV_difference_error_count, EV_relative_diff_error_count, EV_residual_error_count;
    double EV_temp_double, EV_prod_R_w[JJ], EV_sign;

    do{
        EV_difference_error_count=0; EV_relative_diff_error_count=0; EV_residual_error_count=0;
        EV_total_difference_error=0.0; EV_total_relative_diff_error=0.0; EV_total_residual_error=0.0;

        od=0;
        for (i=0; i<OD; i++){

            int last_od = ((od-1)+OD)%OD;
            int next_od = (od+1)%OD;

            /*jj=JJ-1;
            EV[od][jj]=(OD_JAC[last_od][jj+jj*JJ]*EV[last_od][jj])/OD_JAC[last_od][EV_ind+EV_ind*JJ];
            for (jj=0; jj<JJ-1; jj++){
                EV_prod_R_w[jj]=0.0;
                for (cc=jj+1; cc<JJ; cc++) EV_prod_R_w[jj]+=OD_JAC[od][cc+jj*JJ]*EV[od][cc];
                EV[od][jj]=(OD_JAC[od][EV_ind+EV_ind*JJ]*EV[(od+1)%OD][jj]-EV_prod_R_w[jj])/OD_JAC[od][jj+jj*JJ];
            }*/


            // for (jj=0; jj<JJ-1; jj++){
            //     EV_prod_R_w[jj]=0.0;
            //     EV_prod_R_w[jj]+=OD_JAC[last_od][jj+jj*JJ]*EV_last[last_od][jj];
            //     for (cc=jj+1; cc<JJ-1; cc++) EV_prod_R_w[jj]+=OD_JAC[last_od][cc+jj*JJ]*EV_last[last_od][cc];
            //     EV[od][jj]=EV_prod_R_w[jj]/OD_JAC[last_od][EV_ind+EV_ind*JJ];
            // }

            // if (od==0){
            //     for (jj=JJ-2; jj>=0; jj--){
            //         EV_prod_R_w[jj]=0.0;
            //         for (cc=jj+1; cc<JJ-1; cc++) EV_prod_R_w[jj]+=OD_JAC[od][cc+jj*JJ]*EV[od][cc];
            //         EV[od][jj]=(OD_JAC[od][EV_ind+EV_ind*JJ]*EV_last[next_od][jj]-EV_prod_R_w[jj])/OD_JAC[od][jj+jj*JJ];
            //     }
            // } else {
            //     for (jj=JJ-1; jj>=0; jj--){
            //         EV_prod_R_w[jj]=0.0;
            //         for (cc=jj+1; cc<JJ; cc++) EV_prod_R_w[jj]+=OD_JAC[od][cc+jj*JJ]*EV[od][cc];
            //         EV[od][jj]=(OD_JAC[od][EV_ind+EV_ind*JJ]*EV_last[next_od][jj]-EV_prod_R_w[jj])/OD_JAC[od][jj+jj*JJ];
            //     }
            // }

            if (od==0){
                for (jj=JJ-2; jj>=0; jj--){
                    EV_prod_R_w[jj]=0.0;
                    for (cc=jj; cc<JJ-1; cc++) EV_prod_R_w[jj]+=OD_JAC[last_od][cc+jj*JJ]*EV_last[last_od][cc];
                    EV[od][jj]=(EV[od][jj]-EV_prod_R_w[jj])/OD_JAC[last_od][EV_ind+EV_ind*JJ];
                }
            } else {
                for (jj=JJ-1; jj>=0; jj--){
                    EV_prod_R_w[jj]=0.0;
                    for (cc=jj; cc<JJ; cc++) EV_prod_R_w[jj]+=OD_JAC[last_od][cc+jj*JJ]*EV_last[last_od][cc];
                    EV[od][jj]=(EV[od][jj]-EV_prod_R_w[jj])/OD_JAC[last_od][EV_ind+EV_ind*JJ];
                }
            }

            // EV_temp_double=0.0;
            // for (jj=0; jj<JJ; jj++){
            //     EV_temp_double+=EV[od][jj]*EV[od][jj];
            // }
            // EV_temp_double=sqrt(EV_temp_double);
            // for (jj=0; jj<JJ; jj++) EV[od][jj]/=EV_temp_double;

            if (verbose_wEVEC /*|| od==0 || od==1 || od==OD-1*/){
                printf("w_%d^(%d)= ", od+1, it_EV);
                for (jj=0; jj<JJ; jj++) printf("%.17le ", EV[od][jj]);
                printf("\n");
            }


            if (it_EV>0){
                // cblas_dgemv(CblasRowMajor, CblasNoTrans, JJ, JJ, 1.0, OD_JAC[od], JJ, EV[od], 1, 0, EV_temp, 1);
                // EV_error=0;
                // for (jj=0; jj<JJ; jj++){
                //     EV_temp_double = EV_temp[jj]-OD_JAC[od][EV_ind+EV_ind*JJ]*EV_last[(od+1)%OD][jj];
                //     EV_error+=EV_temp_double*EV_temp_double;
                // }
                // EV_error=EV_temp_double*EV_temp_double;
                // EV_diff_total_error+=EV_error;
                // EV_error = sqrt(EV_error);


                for (jj=0; jj<JJ; jj++){
                    // EV_prod_R_w[jj]=0.0;
                    // for (cc=jj; cc<JJ; cc++) EV_prod_R_w[jj]+=OD_JAC[od][cc+jj*JJ]*EV[od][cc];
                    EV_error = OD_JAC[last_od][EV_ind+EV_ind*JJ]*EV[od][jj]-EV_prod_R_w[jj];
                    EV_total_residual_error+=EV_error*EV_error;
                    if (fabs(EV_error)>EV_TOL) EV_residual_error_count++;
                }

                // EV_temp[jj] = OD_JAC[od][jj+jj*JJ]*EV[od][jj];
                // EV_error = fabs(EV_temp[jj]-OD_JAC[od][EV_ind+EV_ind*JJ]*EV_last[(od+1)%OD][jj]);
                // EV_total_residual_error+=EV_error*EV_error;
                // if (EV_error>EV_TOL) EV_residual_error_count++; //Maybe should use square of TOL instead of computing sqrt of ERROR.
            

                // EV_error=0.0;
                // for (jj=0; jj<JJ; jj++) {
                //     EV_temp_double = (EV[od][jj]-EV_last[od][jj]);
                //     EV_error+=EV_temp_double*EV_temp_double;
                // }
                // EV_total_error+=EV_error;
                // EV_error = sqrt(EV_error);

                EV_temp_double = 0.0;
                for (jj=0; jj<JJ; jj++){
                    EV_error = fabs(EV[od][jj]-EV_last[od][jj]);
                    if (EV_error>EV_TOL) EV_difference_error_count++; //Maybe should use square of TOL instead of computing sqrt of ERROR.
                    EV_total_difference_error+=EV_error*EV_error;

                    if (fabs(EV[od][jj])>1e-200 && EV_error/fabs(EV[od][jj])>EV_temp_double) EV_temp_double=EV_error/fabs(EV[od][jj]);
                }
                if (EV_temp_double>EV_total_relative_diff_error) EV_total_relative_diff_error=EV_temp_double;
            }

            od=(od-1+OD)%OD;
        }

        for (od=0; od<OD; od++) for (jj=0; jj<JJ; jj++) EV_last[od][jj] = EV[od][jj];

        it_EV++;

        EV_total_residual_error=sqrt(EV_total_residual_error/(double)OD);
        EV_total_difference_error=sqrt(EV_total_difference_error/(double)OD);

        // printf("RESID: num=%d, acc=%le\t DIFF: num=%d, acc=%le\t REL_DIFF=%le\t\n", EV_residual_error_count, EV_total_residual_error, EV_difference_error_count, EV_total_difference_error, EV_total_relative_diff_error);

    }while (it_EV<2 || EV_residual_error_count>0 && EV_total_residual_error>EV_TOL && EV_difference_error_count>0 && EV_total_difference_error>EV_TOL);


    //-2.58601525075823247e+02



    printf("\n\nw_%d^(%d)= ", OD, it_EV);
    for (jj=0; jj<JJ; jj++) printf("%.17le ", EV[OD-1][jj]);
    printf("\t (lambda_%d=%.17le)\n\n", OD, OD_JAC[OD-1][EV_ind+EV_ind*JJ]);

    printf("w_%d: ", EV_ind+1);
    for (jj=0; jj<JJ; jj++) printf("%.15le ", EV[0][jj]);
    printf(", ||w_%d||=1e%lf, (%d iterations)\n", EV_ind+1, EV_exp_norm, it_EV);

    
    //Compute the product v=Qw to recover the eigenvector
    cblas_dgemv(CblasRowMajor, CblasNoTrans, JJ, JJ, 1.0, OD_Q[0], JJ, EV[0], 1, 0, EV_temp, 1);
    fwrite(EV_temp, sizeof(double), JJ, OD_EVEC_bin);

    printf("EV_%d: ", EV_ind+1);
    for (jj=0; jj<JJ; jj++) printf("%.17le, ", EV_temp[jj]);
    printf("\n");
    
    free(EV[0]);
    for (od=1; od<OD; od++){
        if (store_EVEC || verbose_EVEC) cblas_dgemv(CblasRowMajor, CblasNoTrans, JJ, JJ, 1.0, OD_Q[od], JJ, EV[od], 1, 0, EV_temp, 1);
        if (store_EVEC) fwrite(EV_temp, sizeof(double), JJ, OD_EVEC_bin);

        if (verbose_EVEC){
            printf("EV_%d_(%d): ", EV_ind+1, od+1);
            for (jj=0; jj<JJ; jj++) printf("%.17le, ", EV_temp[jj]);
            printf("\n");
        }

        free(EV[od]);
    }
    free(EV);

    for (od=0; od<OD; od++){
        if (store_R) fwrite(OD_JAC[od], sizeof(double), JJ*JJ, OD_R_bin);
        if (store_Q) fwrite(OD_Q[od], sizeof(double), JJ*JJ, OD_Q_bin);

        if (store_EVAL){
            for (jj=0; jj<JJ; jj++) EV_temp[jj] = OD_JAC[od][jj+jj*JJ];
            fwrite(EV_temp, sizeof(double), JJ, OD_EVAL_bin);
        } 

        free(OD_JAC[od]);
        free(OD_Q[od]);
    }
    free(OD_JAC);
    free(OD_Q);

    clock_t CLOCK_end = clock();
    printf("\nExecution time: %lf seconds\n\n\n", (double)(CLOCK_end-CLOCK_begin)/CLOCKS_PER_SEC);

    return 0;
}