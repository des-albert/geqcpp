#include "geq.h"

extern void gauss(double *, double **, int);
extern double **array2d(int, int);
extern void array2del(double **, int);

void xcur() {

    int mmaxp1 = Mmax + 1;

    double **ax = array2d(mpnmax, mpnmax);
    double bv[3] = {-1., 0, 0};

    for (int i = 1; i < Mmax; i ++) {
        for(int j = i; j < Mmax; j++) {
            ax[i][j] = cl[i][j];
        }
        ax[i][Mmax] = 0.;
        for (int j = 1; j < Nmax; j++) {
            ax[i][mmaxp1 + j] = bb[j][i];
        }
        if (icops >= 2) {
            if (icops > 2) {
                ax[i][Mmax] = 0.;
            }
            else {
                ax[i][Mmax] = cl[Mmax][i];
            }
        }
    }

    ax[Mmax][Mmax] = 0.;
    for (int j = 0; j < Nmax; j++) {
        ax[Mmax][mmaxp1 + j] = bv[ityp[j]];
    }

    if (icops >= 2) {
        if (icops > 2) {
            ax[Mmax][mpnmax - 1] = 1.;
        }
        else {
            ax[Mmax][mpnmax - 1] = 0.;
        }
    }

    for (int i = mmaxp1; i < mpnmax; i++) {
        for (int j = 0; j < mpnmax; j++) {
            ax[i][j] = 0.;
        }
    }
    for (int i = 0; i < mmaxp1; i++) {
        fk[i] = 0.;
    }

    for (int ll = 0; ll < Mmax; ll++) {
        for (int i = 0; i < Mmax; i++) {
            for (int j = 0; j < Mmax; j++) {
                ax[i][j] += +2. * alph * eb[ll][j];
            }
            ax[i][Mmax] -= 2. * alph * eb[ll][i];
            fk[i] -= 2. * alph * eb[ll][Mmax];
        }
        ax[Mmax][Mmax] -= 2. * alph;
        fk[Mmax] += 2. * alph * eb[ll][Mmax];
    }

    for (int i = 0; i < mpnmax; i++) {
        for (int k = 0; k < mpnmax; k++) {
            ax[k][i] = ax[i][k];
        }
    }

    for (int j = 0; j < Nmax; j++) {
        fk[mmaxp1 + j] = - bb[j][Mmax];
    }

    if (icops >= 2) {
        fk[mpnmax - 1] = value;
    }

    gauss(fk, ax, mpnmax);

    psicon = fk[Mmax];
    double energy = 0.;

    for (int i = 0; i < Mmax; i++) {
        int ki = i + 1;
        if ((ki + 1) <= Mmax) {
            for (int k = ki; k < Mmax; k++ ) {
                energy += cl[i][k] * fk[k];
            }
        }
    }
    fk[mmaxp1] = energy;

    array2del(ax, mpnmax);

}


