#include <cmath>
#include "geq.h"

void intpol(const double *, double, double, double *);

void condit(double *v, double rx, double zx, int itp, double **b, int i, int j) {

    double var[3];
    double x = (rx - Rmin)/dr;
    double y = (zx - Zmin)/dz;

    intpol(v, x, y, var);
    b[i][j] = var[itp];

}

double p(double, double, double);
double dp(double, double, double);

void intpol(const double *v, double x, double y, double * w) {

    double af[3], adf[3];
    double a1, a2;

    double dx = x - floor(x + 0.5);
    double dy = y - floor(y + 0.5);
    int m = (int) floor(x - 0.5);
    int n = (int) floor(y + 0.5);
    int mni = m + n * Mr;

    if ( ! (m < 0 || m > Mr - 2 || n < 0 || n > Nz -1)) {
        if (n == 0) {
            for (int i = 0; i < 3; i++) {
                int ki = mni;
                a1 = 2. * (v[ki + Mr] - v[ki]);
                a2 = 0.;
                af[i] = p(dy, a1, a2) + v[ki];
                adf[i] = dp(dy, a1, a2);
            }
        }
        else {
            for (int i = 0; i < 3; i++) {
                int ki = mni + i;
                a1 = v[ki + Mr] + v[ki - Mr] - 2.0 * v[ki];
                a2 = v[ki + Mr] - v[ki - Mr];
                af[i] = p(dy, a1, a2) + v[ki];
                adf[i] = dp(dy, a1, a2);
            }
        }

        a1 = af[0] + af[2] - 2. * af[1];
        a2 = af[2] - af[0];
        w[0] = p(dx, a1, a2) + af[1];
        w[2] = dp(dx, a1, a2);
        a1 = adf[0] + adf[2] - 2. * adf[1];
        a2 = adf[2] - adf[0];
        w[1] = p(dx, a1, a2) + adf[1];

    }

}

double p(double x, double a1, double a2) {
    return (a1 * x + a2) * x / 2. + 0.125 * a1;
}

double dp(double x, double a1, double a2) {
    return a1 * x + a2 / 2.;
}