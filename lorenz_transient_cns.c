#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


#include<limits.h>
#include<float.h>

#include "gmp.h"
#include "mpfr.h"

#define H 2
#define R 20
#define ENDTIME 200
#define TAYLOR 20
#define PREC 200
#define SAMPLENUMBER 200
#define OUTPUT 0

double CNS_lorenz_transient(int taylor, double endTime, mpfr_t xIni, mpfr_t yIni, mpfr_t zIni, int tag, int bool_output);
void sample_1d(mpfr_t xSample[], mpfr_t xStart, mpfr_t xEnd, int sampleNumber);


mpfr_t x[TAYLOR+1], y[TAYLOR+1], z[TAYLOR+1];
mpfr_t X, Y, Z, r, error_r, error_set, sigma, b, h_mpfr, X_stable, Y_stable, Z_stable, X_check, Y_check, Z_check;
mpfr_t save_mul1, save_mul2, save_mul3, save_mul4, save_mul5, save_mul6, save_mul7, save_mul8, save_minus, saver1, saver2;

mpfr_t X_start, X_end, Y_start, Y_end, Z_start, Z_end;
mpfr_t X_sample[SAMPLENUMBER], Y_sample[SAMPLENUMBER], Z_sample[SAMPLENUMBER];

int main(int argc, char *argv[])
{
    int i, max_i, n, check;
    int N_taylor_change;
    FILE *lf;
    time_t t;
    int seed;
    double T_escape[SAMPLENUMBER];
/*-------------------------------------------------------------------------------------*/
    mpfr_inits2(PREC, X, Y, Z, r, error_r, error_set, sigma, b, h_mpfr, X_stable, Y_stable, Z_stable, X_check, Y_check, Z_check, save_mul1, save_mul2, save_mul3, save_mul4, save_mul5, save_mul6, save_mul7, save_mul8, save_minus, saver1, saver2, X_start, X_end, Y_start, Y_end, Z_start, Z_end, NULL);
    for (i=0; i<=TAYLOR; i++)
    {
        mpfr_init2(x[i], PREC);
        mpfr_init2(y[i], PREC);
        mpfr_init2(z[i], PREC);
    }
    for (i=0; i<SAMPLENUMBER; i++)
    {
        mpfr_init2(X_sample[i], PREC);
        mpfr_init2(Y_sample[i], PREC);
        mpfr_init2(Z_sample[i], PREC);
    }
    mpfr_set_str(sigma, "10.0", 10, GMP_RNDN);
    mpfr_set_str(b, "8.0", 10, GMP_RNDN);
    mpfr_div_si(b, b, 3, GMP_RNDN);
    mpfr_set_d(r, R, GMP_RNDN);
    mpfr_sub_si(Z_stable, r, 1, GMP_RNDN);
    mpfr_mul(X_stable, Z_stable, b, GMP_RNDN);
    mpfr_sqrt(X_stable, X_stable, GMP_RNDN);
    mpfr_set(Y_stable, X_stable, GMP_RNDN);
/*-------------------------------------------------------------------------------------*/
    //Input of X, Y, Z, start and end points
    mpfr_set_str(X_start, "-10.0", 10, GMP_RNDN);
    mpfr_set_str(X_end, "10.0", 10, GMP_RNDN);
    mpfr_set_str(Y_start, "-10.0", 10, GMP_RNDN);
    mpfr_set_str(Y_end, "10.0", 10, GMP_RNDN);
    mpfr_set_str(Z_start, "-10.0", 10, GMP_RNDN);
    mpfr_set_str(Z_end, "10.0", 10, GMP_RNDN);
/*-------------------------------------------------------------------------------------*/
    char filename_lf[30] = { };
    sprintf(filename_lf, "lifetime.txt");
    lf = fopen(filename_lf, "w");
    sample_1d(X_sample, X_start, X_end, SAMPLENUMBER);
    sample_1d(Y_sample, Y_start, Y_end, SAMPLENUMBER);
    sample_1d(Z_sample, Z_start, Z_end, SAMPLENUMBER);
    for (i = 0; i < SAMPLENUMBER; i++) {
        mpfr_set(X, X_sample[i], GMP_RNDN);
        mpfr_set(Y, Y_sample[i], GMP_RNDN);
        mpfr_set(Z, Z_sample[i], GMP_RNDN);
        mpfr_fprintf(lf, "%1.16Re, ", X);
        mpfr_fprintf(lf, "%1.16Re, ", Y);
        mpfr_fprintf(lf, "%1.16Re, ", Z);
        T_escape[i] = CNS_lorenz_transient(TAYLOR, ENDTIME, X, Y, Z, i, OUTPUT);
        fprintf(lf, "%f\n", T_escape[i]);
    }
    mpfr_clears(X, Y, Z, r, error_r, X_start, X_end, error_set, sigma, b, h_mpfr, X_stable, Y_stable, Z_stable, X_check, Y_check, Z_check, save_mul1, save_mul2, save_mul3, save_mul4, save_mul5, save_mul6, save_mul7, save_mul8, save_minus, saver1, saver2, NULL);
    for (i=0; i<=TAYLOR; i++)
    {
        mpfr_clear(x[i]);
        mpfr_clear(y[i]);
        mpfr_clear(z[i]);
    }
    fclose(lf);
    return 0;
}

void sample_1d(mpfr_t xSample[], mpfr_t xStart, mpfr_t xEnd, int sampleNumber)
{
    int i;
    for (i = 0; i < sampleNumber; i++) {
        mpfr_set_zero(xSample[i], 1);
        mpfr_set(X, xStart, GMP_RNDN);
        mpfr_sub(save_minus, xEnd, xStart, GMP_RNDN);
        mpfr_div_si(save_minus, save_minus, sampleNumber, GMP_RNDN);
        mpfr_mul_si(save_minus, save_minus, i, GMP_RNDN);
        mpfr_add(xSample[i], xStart, save_minus, GMP_RNDN);
    }
}

double CNS_lorenz_transient(int taylor, double endTime, mpfr_t xIni, mpfr_t yIni, mpfr_t zIni, int tag, int bool_output)
{
   
    int i, j, k, step, endStep, check;
    int out_check;
    double h;
    FILE *fp;
    char filename[30] = { };
    sprintf(filename, "trajectory_%d_%d_%d_%d", taylor, PREC, H, tag);
    fp=fopen(filename,"w");
    step = 0;
    check = 1;
    out_check = pow(10, H);
    h = 1.0/out_check;
    
    mpfr_set_d(h_mpfr, h, GMP_RNDN);
    endStep = floor(endTime / h);
    if (bool_output == 1) {
        mpfr_fprintf(fp, "%f, ", +0.0);
        mpfr_fprintf(fp, "%1.16Re, ", xIni);
        mpfr_fprintf(fp, "%1.16Re, ", yIni);
        mpfr_fprintf(fp, "%1.16Re\n", zIni);
    }
    for (step=1; step < endStep; step++)
    {
        for (i=0; i<=taylor; i++)
        {
            mpfr_set_d (x[i], 0.0, GMP_RNDN);
            mpfr_set_d (y[i], 0.0, GMP_RNDN);
            mpfr_set_d (z[i], 0.0, GMP_RNDN);
        }
        mpfr_set(x[0], xIni, GMP_RNDN);
        mpfr_set(y[0], yIni, GMP_RNDN);
        mpfr_set(z[0], zIni, GMP_RNDN);
        mpfr_set_d(xIni, 0.0, GMP_RNDN);
        mpfr_set_d(yIni, 0.0, GMP_RNDN);
        mpfr_set_d(zIni, 0.0, GMP_RNDN);
        
        mpfr_sub(x[1], y[0], x[0], GMP_RNDN);
        mpfr_mul(x[1], x[1], sigma, GMP_RNDN);
        
        mpfr_mul(save_mul1, x[0], z[0], GMP_RNDN);
        mpfr_mul(y[1], x[0], r, GMP_RNDN);
        mpfr_sub(y[1], y[1], save_mul1, GMP_RNDN);
        mpfr_sub(y[1], y[1], y[0], GMP_RNDN);
        
        mpfr_mul(save_mul2, z[0], b, GMP_RNDN);
        mpfr_mul(z[1], x[0], y[0], GMP_RNDN);
        mpfr_sub(z[1], z[1], save_mul2, GMP_RNDN);
        
        for (k = 1; k < taylor; k++)
        {
            mpfr_sub(x[k+1], y[k], x[k], GMP_RNDN);
            mpfr_mul(x[k+1], x[k+1], sigma, GMP_RNDN);
            mpfr_div_ui(x[k+1],x[k+1],(k+1),GMP_RNDN);
            
            mpfr_mul(y[k+1], x[k], r, GMP_RNDN);
            mpfr_sub(y[k+1], y[k+1], y[k], GMP_RNDN);
            mpfr_mul(save_minus, z[k], b, GMP_RNDN);
            
            mpfr_set_d(saver1, 0.0, GMP_RNDN);
            mpfr_set_d(saver2, 0.0, GMP_RNDN);
            for (i = 0; i <= k; i++){
                mpfr_mul(save_mul4, x[k-i], z[i], GMP_RNDN);
                mpfr_mul(save_mul5, x[k-i], y[i], GMP_RNDN);
                mpfr_add(saver1, saver1, save_mul4, GMP_RNDN);
                mpfr_add(saver2, saver2, save_mul5, GMP_RNDN);
            }
            mpfr_sub(y[k+1], y[k+1], saver1, GMP_RNDN);
            mpfr_sub(z[k+1], saver2, save_minus, GMP_RNDN);
            mpfr_div_ui(y[k+1],y[k+1],(k+1),GMP_RNDN);
            mpfr_div_ui(z[k+1],z[k+1],(k+1),GMP_RNDN);
        }
        for (i = 1; i <= taylor; i++)
        {
            mpfr_set(save_mul6, x[i], GMP_RNDN);
            mpfr_set(save_mul7, y[i], GMP_RNDN);
            mpfr_set(save_mul8, z[i], GMP_RNDN);
            for (j=1; j<=i; j++) {
                mpfr_mul(save_mul6, save_mul6, h_mpfr, GMP_RNDN);
                mpfr_mul(save_mul7, save_mul7, h_mpfr, GMP_RNDN);
                mpfr_mul(save_mul8, save_mul8, h_mpfr, GMP_RNDN);
            }
            mpfr_add(xIni, xIni, save_mul6, GMP_RNDN);
            mpfr_add(yIni, yIni, save_mul7, GMP_RNDN);
            mpfr_add(zIni, zIni, save_mul8, GMP_RNDN);
        }
        mpfr_add(xIni, xIni, x[0], GMP_RNDN);
        mpfr_add(yIni, yIni, y[0], GMP_RNDN);
        mpfr_add(zIni, zIni, z[0], GMP_RNDN);
        mpfr_abs(X_check, xIni, GMP_RNDN);
        mpfr_abs(Y_check, yIni, GMP_RNDN);
        mpfr_set(Z_check, zIni, GMP_RNDN);
        mpfr_sub(X_check, X_stable, X_check, GMP_RNDN);
        mpfr_sub(Y_check, Y_stable, Y_check, GMP_RNDN);
        mpfr_sub(Z_check, Z_stable, Z_check, GMP_RNDN);
        mpfr_mul(X_check, X_check, X_check, GMP_RNDN);
        mpfr_mul(Y_check, Y_check, Y_check, GMP_RNDN);
        mpfr_mul(Z_check, Z_check, Z_check, GMP_RNDN);
        mpfr_add(X_check, X_check, Y_check, GMP_RNDN);
        mpfr_add(X_check, X_check, Z_check, GMP_RNDN);
        mpfr_sqr(X_check, X_check, GMP_RNDN);
        check = mpfr_cmp_si(X_check, 1);
        if (bool_output == 1) {
            if (step % (out_check/100) == 0) {
                mpfr_fprintf(fp, "%f, ",step*h);
                mpfr_fprintf(fp, "%1.16Re, ", xIni);
                mpfr_fprintf(fp, "%1.16Re, ", yIni);
                mpfr_fprintf(fp, "%1.16Re\n", zIni);
            }
            
        }
        if (check != 1) break;
    }
    if (bool_output == 2) {
        mpfr_fprintf(fp, "%f, ",step*h);
        mpfr_fprintf(fp, "%1.16Re, ", xIni);
        mpfr_fprintf(fp, "%1.16Re, ", yIni);
        mpfr_fprintf(fp, "%1.16Re\n", zIni);
    }
    fclose(fp);
    return step * h;
}
