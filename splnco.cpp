#include <algorithm>
#include "geq.h"

using namespace std;

void splnco(double *psi) {

    double a[8];

    double a1 = -4.0;

    for (int i = 0; i < 4; i++) {
        a[i] = 1. / a1;
        a1 = a1 - 2.0;
        a[7 - i] = 1. / a1;
        a1 = a1 * a1;
    }
    /*

       Control for inner loop on columns

     */

    int jh1 = 1;
    int lend = Mr * Nm1 + 1;
    int ls = Mr;
    int jh2 = min(16, Mm1 / 2);
    int il = Mm1 - 1;

    L20:
    int k = 0;
    int jh = jh1;
    int mode = 2;

    L30:
    int j = 2 * jh;

    for (int l = 1; l <= lend; l += ls) {
        int init = l + jh * mode;
        int iend = l + il;
        for (int i = init - 1; i < iend; i += j) {
            psi[i] = psi[i] + a[k] * (psi[i + jh] + psi[i - jh] - 2. * psi[i]);
        }
    }
    k += 1;

    if (mode == 2) {
        jh = 2 * jh;
        if (jh < jh2) {
            goto L30;
        }
        mode = 1;
        if (k == 4) {
            jh = jh / 2;
            if (jh >= jh1)
                goto L30;
            if (jh1 == Mr) {
                goto L50;
            }
            jh1 = Mr;
            lend = Mr;
            ls = 1;
            jh2 = min(16 * Mr, Mr * Nm1 / 2);
            il = Mr * (Nm1 - 1);
            goto L20;
        }
        k = 7 - k;
        goto L30;
    } else {
        jh = jh / 2;
        if (jh >= jh1)
            goto L30;
        if (jh1 == Mr) {
            goto L50;
        }

        /*
                Control for inner loop on rows
        */

        jh1 = Mr;
        lend = Mr;
        ls = 1;
        jh2 = min(16 * Mr, Mr * Nm1 / 2);
        il = Mr * (Nm1 - 1);
        goto L20;

    }
    L50:
    return;
}

