#include <cmath>
#include <iostream>
#include <fstream>

#define EXTERN

#include "geq.h"

using namespace std;

extern void bndmat();
extern void eqsil(double *q);
extern double gfl(double rv, double rst, double zv, double del);

int main() {

    double Offset, El, Rxpsn, Zxpsn, tri, elxp, trixp, ang, al1, al2, rc1, rc2, anga, angb, ang1, ang2;
    double rac, zac, exc, Ra[Ng][Mc], Za[Ng][Mc], Ex[Ng][Mc], Rl[Ng][Mc];
    double Rc[Nmax], Zc[Nmax], Rcc[llmax], Zcc[llmax], *expsi, *fool;
    int ic[Mc] = {};

    cout << "Garching Tokamak Equlibrium" << endl;

    int Nexp = 6;
    Mr = 1 + (1 << Nexp);
    Nz = Mr;
    MN = Mr * Nz;
    Mm1 = Mr - 1;
    Nm1 = Nz - 1;
    llp = 2 * (Mr + Nz) - 8;
    pi = 3.141592653589793238462643383279502884197;

    ifstream fin;

    fin.open("solver.dat");

    fin >> Rmin >> Rmax >> Zmin >> Zmax >> error >> meshfg >> mprfg;
    fin >> Rmpl >> Offset >> Apl >> El >> tri >> Rxpsn >> Zxpsn;


    if (meshfg > 0) {
        Rmin = Rmpl - 32.0 * Apl / 20.0;
        Rmax = Rmpl + 32.0 * Apl / 20.0;
        Zmin = Offset - 32.0 * (Offset + abs(Zxpsn)) / 20.0;
        Zmax = Offset + 32.0 * (Apl * El) / 20.0;
    }

    R = new double[Mr];
    Z = new double[Nz];
    ip = new int[llp];
    jp = new int[llp];
    ityp = new int[6];
    aux = new double[llp * llp];


    dr = (Rmax - Rmin) / (double) Mm1;
    dz = (Zmax - Zmin) / (double) Nm1;

    for (int i = 0; i < Mr; i++) {
        R[i] = Rmin + dr * (double) i;
    }
    for (int j = 0; j < Nz; j++) {
        Z[j] = Zmin + dz * (double) j;
    }

    alpha = (Rmax - Rmin) / (Rmax + Rmin);
    sh = dz / dr;
    ss = sh * sh;

    cout << "Rmin  = " << Rmin << " Rmax = " << Rmax << " dr = " << dr << endl;
    cout << "Zmin  = " << Zmin << " zmax = " << Zmax << " dz = " << dz << endl;
    cout << "alpha = " << alpha << "   sh = " << sh << " ss = " << ss << endl;

    bndmat();

/*
    Read in Poloidal Field Coil Data
*/

    int k = 0;

    while (true) {
        ic[k] = -1;
        while (true) {
            fin >> rac >> zac >> exc;

            if (rac <= 0.0) goto L10;
            ic[k] = ic[k] + 1;
            Ra[ic[k]][k] = rac;
            Za[ic[k]][k] = zac;
            Ex[ic[k]][k] = exc;
            Rl[ic[k]][k] = 1.0e-20 * rac;

            k += 1;
        }

        L10:
        if (ic[k] < 0) goto L20;
    }

    L20:
    int Mmax = k;

    if (mprfg != 0) {
        cout << "Conductor groups available for optimization" << endl;
        for (int j = 0; j < Mmax; k++) {
            cout << " Group : " << k + 1 << endl;
            for (int i = 0; i <= ic[j]; i++) {
                printf(" %7.3f   %7.3f   %7.3f \n", Ra[i][j], Za[i][j], Ex[i][j]);
            }
        }
    }

    /*
        Conditions to be satisfied by resulting equilibrium
    */

    elxp = abs(Offset - Zxpsn) / Apl;
    trixp = (Rmpl - Rxpsn) / Apl;

    for (int j = 0; j < 3; j++) {
        ityp[j] = j;
        Rc[j] = Rxpsn;
        Zc[j] = Zxpsn;
    }

    ityp[3] = 1;
    Rc[3] = Rmpl - Apl;
    Zc[3] = Offset;
    ityp[4] = 1;
    Rc[4] = Rmpl + Apl;
    Zc[4] = Offset;
    ityp[5] = 1;
    Rc[5] = Rmpl - Apl * tri;
    Zc[5] = Offset + Apl * El;

    for (int j = 0; j < 4; j++) {
        ang = (pi * double(j + 1)) / 10.0;
        Rcc[j] = Rmpl + Apl * cos(ang + tri * sin(ang));
        Zcc[j] = Offset + El * Apl * sin(ang);
        Rcc[j + 4] = Rmpl + Apl * cos(ang + pi / 2.0 + tri * sin(ang + pi / 2.0));
        Zcc[j + 4] = Offset + El * Apl * sin(ang + pi / 2.0);
    }

    for (int j = 0; j < 4; j++) {
        al1 = Apl * (((1.0 + trixp) * (1.0 + trixp)) + elxp * elxp) / (2.0 * (1.0 + trixp));
        al2 = Apl * (((1.0 - trixp) * (1.0 - trixp)) + elxp * elxp) / (2.0 * (1.0 - trixp));
        anga = atan(2.0 * elxp * (1.0 + trixp) / (elxp * elxp - (1.0 + trixp) * (1.0 + trixp)));
        angb = atan(2.0 * elxp * (1.0 - trixp) / (elxp * elxp - (1.0 - trixp) * (1.0 - trixp)));
        rc1 = Rmpl + Apl - al1;
        rc2 = Rmpl - Apl + al2;
        ang1 = anga * (j + 1) / 5.0;
        ang2 = angb * (j + 1) / 5.0;
        Rcc[j + 8] = rc1 + al1 * cos(ang1);
        Zcc[j + 8] = Offset - al1 * sin(ang1);
        Rcc[j + 12] = rc2 - al2 * cos(ang2);
        Zcc[j + 12] = Offset - al2 * sin(ang2);
    }

    if (mprfg >= 0) {

        cout << "  Single null case: boundary points " << endl;

        for (int j = 0; j < llmax; j++) {
            printf(" %7.3f   %7.3f\n", Rcc[j], Zcc[j]);
        }
    }

    expsi = new double(MN);
    fool = new double(MN);

    int icl, nof, jn;
    for (int kk = 0; kk < Mmax; kk++) {
        icl = ic[kk];
        for (int i = 0; i < Nz; i++) {
            nof = i * Mr;
            for (int j = 0; j < Mr; j++) {
                jn = nof + j;
                expsi[jn] = 0.;
            }
        }

        for (int i = 0; i < icl; i++) {
            if ( ! (((Zmax - Za[i][kk]) * (Zmin - Za[i][kk]) <= 0.) && ((Rmax - Ra[i][kk]) * (Rmin - Ra[i][kk]) <= 0.)) ) {
                for (int l = 0; l < Nz; l += Nm1) {
                    nof = l * Mr;
                    for (int j = 0; j < Mr; j++) {
                        jn = nof + j;
                        expsi[jn] += Ex[i][kk] * gfl(R[j], Ra[i][kk], Z[l] - Za[i][kk], 0.);
                    }
                }

                for (int l = 0; l < Nz; l++) {
                    nof = l * Mr;
                    for (int j = 0; j < Mr; j+= Mm1) {
                        jn = nof + j;
                        expsi[jn] += Ex[i][kk] * gfl(R[j], Ra[i][kk], Z[l] - Za[i][kk], 0.);
                    }
                }
            }
        }

        eqsil(expsi);

    }


    fin.close();

}
