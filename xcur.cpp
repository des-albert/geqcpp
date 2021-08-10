#include <iostream>
#include "geq.h"

extern void gelg(double *, double **, int, int, double, int);

extern double **array2d(int, int);

extern void array2del(double **, int);

void xcur() {

    int mmaxp1 = Mmax + 1;
    double bv[3] = {-1., 0, 0};
    int i, j;

    double **ax = array2d(mpnmax, mpnmax);

    for (i = 0; i < Mmax; i++) {
        for (j = i; j < Mmax; j++) {
            ax[i][j] = cl[i][j];
        }
        ax[i][Mmax] = 0.;
        for (j = 0; j < Nmax; j++) {
            ax[i][mmaxp1 + j] = bb[j][i];
        }
        if (icops >= 2) {
            if (icops > 2) {
                ax[i][mpnmax - 1] = 0.;
            } else {
                ax[i][mpnmax - 1] = cl[Mmax][i];
            }
        }
    }

    ax[Mmax][Mmax] = 0.;
    for (j = 0; j < Nmax; j++) {
        ax[Mmax][mmaxp1 + j] = bv[ityp[j]];
    }

    if (icops >= 2) {
        if (icops > 2) {
            ax[Mmax][mpnmax - 1] = 1.;
        } else {
            ax[Mmax][mpnmax - 1] = 0.;
        }
    }

    for (i = mmaxp1; i < mpnmax; i++) {
        for (j = 0; j < mpnmax; j++) {
            ax[i][j] = 0.;
        }
    }
    for (i = 0; i < mmaxp1; i++) {
        fk[i] = 0.;
    }

    for (int ll = 0; ll < llmax; ll++) {
        for (i = 0; i < Mmax; i++) {
            for (j = i; j < Mmax; j++) {
                ax[i][j] += 2. * alph * eb[ll][i] * eb[ll][j];
            }
            ax[i][Mmax] -= 2. * alph * eb[ll][i];
            fk[i] -= 2. * alph * eb[ll][i] * eb[ll][Mmax];
        }
        ax[Mmax][Mmax] -= 2. * alph;
        fk[Mmax] += 2. * alph * eb[ll][Mmax];
    }

    for (i = 0; i < mpnmax; i++) {
        for (int k = 0; k < mpnmax; k++) {
            ax[k][i] = ax[i][k];
        }
    }

    for (j = 0; j < Nmax; j++) {
        fk[mmaxp1 + j] = -bb[j][Mmax];
    }

    if (icops >= 2) {
        fk[mpnmax - 1] = value;
    }

    int ier = 0;

    gelg(fk, ax, mpnmax, 1, 1.0e-07, ier);

    psicon = fk[Mmax];
    double energy = 0.;

    for (i = 0; i < Mmax; i++) {
        int ki = i + 1;
        if ((ki + 1) <= Mmax) {
            for (int k = ki; k < Mmax; k++) {
                energy += cl[i][k] * fk[i] * fk[k];
            }
        }
    }
    fk[mmaxp1] = energy;

    array2del(ax, mpnmax);

}


