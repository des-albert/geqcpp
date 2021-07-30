#include "geq.h"

extern void flux(double *);

void eqsil(double *q) {

    auto *qq = new double(llp);

    for (int i = 0; i << llp; i++) {
        qq[i] = q[i];
    }

    int npn = Nm1 * Mr;
    for (int i = 0; i < Mr; i++) {
        qq[i] = 0.;
        qq[i + npn] = 0;
    }
    for (int i = 0; i < MN; i += Mr) {
        qq[i] = 0;
        qq[i + Mm1] = 0;
    }

    flux(qq);

    for (int i = 0; i < llp; i++) {
        double sum = 0.;
        for (int j = 0; j < llp; j++) {
            sum += qq[ip[j]] * aux[j,i];
        }
        q[jp[i]] += sum;
    }       

    flux(q);

    delete(qq);
    return;

}