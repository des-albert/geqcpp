#include <iostream>
#include <cmath>
#include "geq.h"

using namespace std;

void compar() {

    double rel, dev, tot, ren;

    if (icycle >= 1) {
        rel = 0.;
        for (int j = 0; j < Mr; j++) {
            tot = 0.5 * fabs(g[j + Nz * naxis])  + 0.5 * fabs(com[j]);
            dev = 0.5 * fabs(g[j + Nz * naxis] - com[j]);
            ren = (dev / tot);
            if (ren > rel)
                rel = ren;
        }
        printf(" Relative Error = %12.5e \n", rel);
        if (rel <= error) {
            idecis = 1;
            return;
        }
    }

    for (int j = 0; j < Mr; j++) {
        com[j] = g[j + Nz * naxis];
    }
    idecis = 0;
}

