#include <cmath>
#include <iostream>
#include "geq.h"

using namespace std;

void startt() {

    double cp = 3. * (totcurr) / (20. * dr * dz) / fabs((ndes - naxis));

    cout << "Start Current : " << totcurr << endl;

    for (int i = 0; i < Mr; i++) {
        for (int j = 0; j < Nz; j++) {
            g[i * Nz + j] = 0.;
        }
    }

    int nmin = 2 * naxis - ndes;
    int njmax = ndes;
    if (njmax <= nmin) {
        njmax = nmin;
        nmin = ndes;
    }

    for (int j = nmin; j < njmax; j++) {
        for (int i = 0; i < 5; i++) {
            g[(jaxis - 2 + i) + j * Mr] =
                    cp * (1. - pow((double) (j - naxis) / (double) (ndes - naxis), 2.0)) * dz * dz;
        }
    }

}

