#include <cstdio>
#include "geq.h"

extern int lcod;

void topol() {

    int lsym, nmx[2], j, l, n, jma, jmax, jnu, ji, jmin, j1, j2, jmi;
    int ndi[2] = {1, -1};

    lsym = 0;
    nmx[1] = 0;
    if (naxis != 1)
        lsym = 1;

    for (j = 0; j < Mr; j++) {
        for (n = 0; n < Nz; n++) {
            g[j + n * Nz] -= psicon;
            if (g[j + n * Nz] <= 0.) {
                g[j + n * Nz] = 0.;
            }
        }
    }

    n = naxis;
    jma = jaxis;

    if (g[jma + n * Nz] == 0.) {
        lcod = 5;
        printf("Plasma disappeared of region, probably squeezed of f at center\n");
    } else {
        for (l = 0; l <= lsym; l++) {
            jma = jaxis;
            n = naxis;
            for (j = jma; j < Mr; j++) {
                if (g[j + n * Nz] == 0.) {
                    goto L10;
                }
            }
            goto L100;
            L10:
            jmax = j - 1;
            jnu = Mr - jmax - 1;

            for (ji = jnu; ji < Mr; ji++) {
                j = Mr - ji - 1;
                if (g[j + n * Nz] == 0.) {
                    goto L20;
                }
            }
            goto L110;

            L20:
            jmin = j + 1;

            while (true) {
                j1 = jmin - 1;
                j2 = jmax + 1;
                for (j = 0; j <= j1; j++) {
                    g[j + n * Nz] = 0.;
                }

                for (j = j2; j < Mr; j++) {
                    g[j + n * Nz] = 0.;
                }

                if (jmax <= jmin)
                    goto L90;

                n += ndi[l];
                if ( ! ( (n < (Nz - 1)) && (n > 0) ) ) {
                    goto L140;
                }
                if (g[jmax + n * Nz] == 0.) {
                    jnu = Mr - jmax - 1;
                    for (ji = jnu; ji < Mr; ji++) {
                        j = Mr - ji - 1;
                        if (g[j + n * Nz] != 0.)
                            goto L30;
                    }
                    goto L70;

                    L30:
                    jmax = j;
                } else {
                    jma = jmax;
                    for (j = jma; j < Mr; j++) {
                        if (g[j + n * Nz] == 0.) {
                            goto L40;
                        }
                    }
                    goto L120;

                    L40:
                    jmax = j - 1;
                }

                jmi = jmin;
                if (g[jmin + n * Nz] == 0.) {
                    for (j = jmi; j < Mr; j++) {
                        if (g[j + n * Nz] != 0.) {
                            goto L50;
                        }
                    }
                    goto L80;

                    L50:
                    jmin = j;
                } else {
                    jnu = Mr - jmin - 1;
                    for (ji = jnu; j < Mr; ji++) {
                        j = Mr - ji - 1;
                        if (g[j + n * Nz] == 0.) {
                            goto L60;
                        }
                    }
                    goto L130;
                    L60:
                    jmin = j + 1;
                }
            }

            L70:
            n -= ndi[l];
            goto L90;

            L80:
            n -= 1;

            L90:
            nmx[l] = n + ndi[l];
        }

        for (n = 0; n < Nz; n++) {
            if ( (nmx[0] - n) * (nmx[1] - n) >= 0.) {
                for (j = 0; j < Mr; j++) {
                    g[j + n * Nz] = 0.;
                }
            }
        }

        return;

        L100:
        printf("Plasma runs out of grid on outside at axis\n");
        lcod = 1;
        goto L150;

        L110:
        printf("Plasma runs out of grid on inside at axis\n");
        lcod = 6;
        goto L150;

        L120:
        printf("Plasma runs out of grid on outside for n = %3i \n", n);
        lcod = 3;
        goto L150;

        L130:
        printf("Plasma runs out of grid on inside for n = %3i \n", n);
        lcod = 4;
        goto L150;
        L140:
        printf("Plasma runs out on top immersing probably top conductor\n");
        lcod = 2;
    }

    L150:
    printf(" Case abandoned");
}

