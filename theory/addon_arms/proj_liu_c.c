/* proj_liu_c.c -- compiled copy of proj_A_liu() in proj_liu.R (same definition, same zero rule).
 *
 *   A_ij = sum_l A_ijl,  A_ijl = (pi - angle(X_i - X_l, X_j - X_l)) / (2 pi),
 *   A_ijl = 1/2 when exactly one of X_i - X_l, X_j - X_l is the zero vector, 1 when both are.
 *
 * X is the n x p matrix (column-major) AFTER standardisation (done in R). Only the upper triangle
 * is computed and mirrored. Checked against proj_A_liu() in check_proj_liu_c.R.
 * Build:  R CMD SHLIB proj_liu_c.c
 */
#include <math.h>
#include <stdlib.h>
#include <R.h>

void proj_A_liu_c(double *X, int *pn, int *pp, double *peps, double *A)
{
    int n = *pn, p = *pp, i, j, k, l;
    double eps = *peps, c, s, nm;
    const double PIV = 3.14159265358979323846;
    double *U = (double *) R_alloc((size_t) n * p, sizeof(double));
    int *z = (int *) R_alloc((size_t) n, sizeof(int));
    for (i = 0; i < n * n; i++) A[i] = 0.0;
    for (l = 0; l < n; l++) {
        for (i = 0; i < n; i++) {                 /* unit difference vectors, row-major in U */
            s = 0.0;
            for (k = 0; k < p; k++) { c = X[i + (size_t) k * n] - X[l + (size_t) k * n]; U[(size_t) i * p + k] = c; s += c * c; }
            nm = sqrt(s);
            z[i] = (nm <= eps);
            if (!z[i]) for (k = 0; k < p; k++) U[(size_t) i * p + k] /= nm;
        }
        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                double K;
                if (z[i] && z[j]) K = 1.0;
                else if (i == j) K = 0.5;                 /* angle 0 exactly */
                else if (z[i] || z[j]) K = 0.5;
                else {
                    c = 0.0;
                    for (k = 0; k < p; k++) c += U[(size_t) i * p + k] * U[(size_t) j * p + k];
                    if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
                    K = (PIV - acos(c)) / (2.0 * PIV);
                }
                A[i + (size_t) j * n] += K;
            }
        }
    }
    for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) A[j + (size_t) i * n] = A[i + (size_t) j * n];
}
