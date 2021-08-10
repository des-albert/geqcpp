#include <cstdio>
#include "geq.h"

extern void flux(double *);

void eqsil(double *q) {

    auto *qq = new double[MN];

    for (int i = 0; i < MN; i++) {
        qq[i] = q[i];
    }

    int npn = Nm1 * Mr;
    for (int i = 0; i < Mr; i++) {
        qq[i] = 0.;
        qq[i + npn] = 0.;
    }
    for (int k = 0; k < MN; k += Mr) {
        qq[k] = 0;
        qq[k + Mm1] = 0;
    }

    flux(qq);

    for (int i = 0; i < llp; i++) {
        double sum = 0.;
        for (int l = 0; l < llp; l++) {
            sum = sum + aux[l][i] * qq[ip[l]];
        }
        q[jp[i]] = q[jp[i]] + sum;
        sum = 0.;
    }

    flux(q);

    delete[] qq;

}