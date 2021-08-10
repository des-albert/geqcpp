#include "geq.h"

double f1(double, double, double);

void saddle() {

    double hk, zfrr, zfzz, zfrz;
    int i, k, l, ko, io;
    bool sg, gg;

    irsp = 0;
    for (i = 1; i < Mm1; i++) {
        for (k = 1; k < Nm1; k++) {
            hk = g[k + i * Mr];
            if (!(irsp != 0 && psicon >= hk)) {
                zfrr = g[k + 1 + i * Mr] - 2. * hk + g[k - 1 + i * Mr];
                zfzz = g[k + (i + 1) * Mr] - 2. * hk + g[k + (i - 1) * Mr];
                zfrz = g[k + 1 + (i + 1) * Mr] + g[k - 1 + (i - 1) * Mr] - g[k + 1 + (i - 1) * Mr] -
                       g[k - 1 + (i + 1) * Mr];
                if (16. * zfrr * zfzz - zfrz * zfrz < 0.) {
                    sg = false;
                    gg = false;
                    ko = 1;
                    io = -2;
                    for (l = 1; l <= 4; l++) {
                        if (l < 4)
                            io += 1;
                        if (l == 4)
                            ko = 0;
                        if (f1(g[k + ko + (i + io) * Mr], g[k - ko + (i - io) * Mr], hk) > 0.) {
                            if (g[k + ko + (i + io) * Mr] != hk) {
                                if (g[k + ko + (i + io) * Mr] > hk) {
                                    sg = true;
                                } else {
                                    gg = true;
                                }
                                if (sg && gg)
                                    goto L10;
                            }
                        }
                    }
                    goto L20;
                    L10:
                    irsp = k;
                    izsp = i;
                    psicon = hk;
                }
            }
            L20:
            continue;
        }
    }

}

double f1(double a, double b, double hk) {
    return (a - hk) * (b - hk);

}