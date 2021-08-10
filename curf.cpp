#include <cstdio>
#include <cmath>
#include "geq.h"

extern void topol();

double pp(double, double);
double ff(double, double);
double ppin(double, double);
double ffin(double, double);

int lcod = 0;

void curf() {
    int j, k, n, ipoi;
    double dint[10], s[10], beh[10];
    double rax, zax, xalp, xalp1, xalp2, xalpr2, xalpz2, rzy, gn, xind, cjm, c0pr, c0btr;

    if (mprfg != 0)
        printf(" psi saddle point                = %12.5f \n", psicon);

    topol();

    for (j = 0; j < 10; j++) {
        dint[j] = 0;
    }

    fmaxa = 0.;
    for (j = 0; j < Mr; j++) {
        for (n = 0; n < Nz; n++) {
            if (g[j + n * Nz] > fmaxa) {
                fmaxa = g[j + n * Nz];
                jaxis = j;
                naxis = n;
                rax = R[jaxis];
                zax = Z[naxis];
            }
        }
    }

    fab = fmaxa + psicon;
    if (mprfg != 0) {
        printf (" Magnetic Axis radius = %12.5f  height =  %12.5f psi =  %12.5f \n", rax, zax, fab);
    }
    ipoi = 0;
    idol = 0;
    xalp = 0.;

    for (j = 0; j < Mr; j++) {
        rzy = R[j];
        xalp1 = 0;
        for (k = 0; k < 10; k ++) {
            beh[k] = 0.;
        }
        for (n = 0; n < Nz; n++) {
            xalp2 = 0.;
            if ( g[j + n * Nz] > 0. ) {
                if ( ((j > 0) && (j < (Mr - 1) )) && ( (n > 0) && (n < (Nz - 1) ) ) ) {
                    xalpr2 = pow( (g[j + 1 + n * Nz] - g[j - 1 + n * Nz]) / (2. * dr), 2.);
                    xalpz2 = pow( (g[j + (n + 1) * Nz] - g[j + (n - 1) * Nz]) / (2. * dz), 2.);

                    if (g[j - 1 + n * Nz] <= 0. ) {
                        xalpr2 = pow((g[j + 1 + n * Nz] - g[j + n * Nz])/dr, 2.0);
                    }
                    if (g[j + 1 + n * Nz] <= 0. ) {
                        xalpr2 = pow((g[j + n * Nz] - g[j - 1 + n * Nz])/dr, 2.0);
                    }
                    if (g[j + (n - 1) * Nz] <= 0. ) {
                        xalpz2 = pow((g[j + (n + 1) * Nz] - g[j + n * Nz])/dz, 2.0);
                    }
                    if (g[j + (n + 1) * Nz] <= 0. ) {
                        xalpz2 = pow((g[j + n * Nz] - g[j + (n - 1) * Nz])/dz, 2.0);
                    }

                    xalp2 = xalp2 + (xalpr2 + xalpz2)/rzy;
                    ipoi += 1;
                }
                else {
                    idol += 1;
                }

                gn = g[j + n * Nz] / fmaxa;
                s[0] = pp(gn, alfac) * rzy;
                s[1] = ff(gn, alfac) / rzy;
                s[2] = 1. / pow(rzy, 2.);
                s[3] = ppin(gn, alfac) * fmaxa;
                s[4] = ffin(gn, alfac) *fmaxa / pow (rzy, 2.);
                s[5] = 1.;
                s[6] = 2. * pi * rzy * ppin(gn, alfac) * fmaxa;
                s[7] = 2. * pi * rzy;
                s[8] = 0.;
                s[9] = 0.;
                for (k = 0; k < 10; k++) {
                    beh[k] += s[k];
                }
            }
        }

        xalp += xalp1;
        for (k = 0; k < 10; k++) {
            dint[k] += beh[k];
        }
    }

    for (k = 0; k < 10; k++) {
        dint[k] = dint[k] * dz * dr;
    }
    xalp = 2. * pi *xalp * dr * dz;
    xind = 2. * xalp / (totcurr * totcurr * Rmpl);
    cjm = 0;
    c0pr = betapol * pow(totcurr, 2.)/(8. * pi *dint[3]);
    c0btr = (totcurr - c0pr*dint[0]) / dint[1];

    if (mprfg != 0) {
        printf(" c0pr = %12.5f  c0btr = %12.5f \n", c0pr, c0btr);
        printf(" Plasma Area = %12.5f  Plasma Volume = %12.5f \n\n",dint[5], dint[7]);
    }

    for (j = 0; j < Mr; j++ ) {
        gn = g[j + naxis * Nz] / fmaxa;
        pr[j] = c0pr * ppin(gn, alfac) * fmaxa;
        bt2[j] = c0btr * ffin(gn, alfac) * fmaxa / pow(R[j], 2.);
    }
    for (j = 0; j < Mr; j++ ) {
        for (n = 0; n < Nz; n++) {
            gn = g[j + n * Nz] / fmaxa;
            g[j + n * Nz] = (c0pr * pp(gn, alfac) * pow(R[j], 2.) + c0btr * ff(gn, alfac)) * pow(dz, 2.);
        }
    }

    pint = c0pr * dint[3];
    pintvo = c0pr * dint[6];
    bt2int = c0btr * dint[4];
    betap = 2. * pint / (pow(totcurr, 2.) * dint[5]);
    betapvo = 2. * pintvo / (pow(totcurr, 2.) * dint[7]);
    betat = dint[2];

    for (j = 0; j < Mr; j++) {
        cjt[j] = g[j + naxis * Nz] / R[j];
        if ( fabs(cjt[j]) > cjm) {
            cjm = fabs(cjt[j]);
        }
    }
    if (cjm < pow(dz, 2.)) {
        cjm = 1.;
    }
    for (j = 0; j < Mr; j++) {
        cjt[j] = cjt[j]/cjm;
    }


}

double pp(double x, double alfa) {
        return pow(x, alfa) * (2. - pow(x, alfa));
}

double ff(double x, double alfa) {
    return pow(x, alfa) * (2. - pow(x, alfa));
}

double ppin(double x, double alfa) {
    return (2. / (alfa + 1.)) * pow(x, alfa + 1.) - (1. / (2. * alfa + 1.)) * pow(x, 2. * alfa + 1.);
}

double ffin(double x, double alfa) {
    return (2. / (alfa + 1.)) * pow(x, alfa + 1.) - (1. / (2. * alfa + 1.)) * pow(x, 2. * alfa + 1.);
}