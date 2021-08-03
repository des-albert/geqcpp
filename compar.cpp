#include <algorithm>
#include <iostream>
#include "geq.h"

using namespace std;

void compar() {

    double rel, dev, tot, ren;

    if (icycle > 1) {
        rel = 0.;
        for (int j = 0; j < Mr; j++) {
            tot = abs(g[j * Nz + naxis]) /2. + abs(com[j]) /2.;
            dev = (abs(g[j * Nz + naxis]) - com[j]) / 2.;
            ren = (dev/tot);
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
        com[j] = g[j * Nz + naxis];
    }
    idecis = 0;
}

