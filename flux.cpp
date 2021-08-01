#include "geq.h"
#include <cstdio>
#include <cmath>

void flux(double *psi) {

    auto *tc = new double[llp];
    auto *a = new double[llp];

    int nn = Nm1;
    int iu = Mm1 - 1;
    int itc = 2 * Nm1;
    tc[itc + nn/2 - 1] = 0.;
    int lo = nn/2;

L10: 
    int l = lo/2;
    tc[itc + l - 1] = sqrt(2.0 + tc[itc + lo - 1]);
    lo = l;
L20:
    tc[itc + nn - l - 1] = -tc[itc + l - 1];
    l = l + 2*lo;    

    int sw =  (2*l/nn) * (2*lo - 3); 
    
    if (sw > 0) {
        goto L10;
    }
    else if (sw == 0) {
        tc[itc + l - 1] = (tc[itc + l + lo - 1] + tc[itc + l - lo - 1])/tc[itc + lo - 1];
        goto L20;
    }
    else if (sw < 0) {
        for (int i = 0; i < iu; i++) {
            double d = alpha/((double)Mm1 + (double)(2*(i + 1) - Mm1) * alpha);
            tc[i] = ss/(1. - d);
            tc[i + Mm1] = ss/(1. + d);
        }
    }

    lo = nn/2;
    int j1 = 1 + Mr * (Nm1 / nn);
    int ju = Nm1 * Mr;

    for (int i = j1 - 1; i < ju; i += Mr) {
        psi[i + 1] = psi[i + 1] + tc[0] * psi[i];
        psi[iu + i] =  psi[iu + i] + tc[iu + Mm1 - 1] * psi[i + Mm1];
    }

    a[Mm1 - 1] = 0.;
    a[2*Mm1 - 1] = 0.;
    int mode = 2;
    int is = -1;

L30:
    int li = 2*lo;  

    int iphase = 2*mode - li/nn;
    int jd = Mr * nn/li;
    int jh = jd/2;
    int jt = jd + jh;
    int ji = 2*jd;
    int jo = 1 + jd*mode*(1 - is)/2;

    for (int j = jo; j <= ju; j+= ji) {
        j1 = j + 1;
        int jdm = jd * is;
        int jhm = jh * is;
        int jtm = jt * is;
        int jiu = j + iu;

        switch (iphase) {
            case 1:
                for (int i = j1 - 1; i < jiu; i++) {
                    a[i - j] = psi[i] + psi[i + jd] + psi[i + jdm];
                    psi[i] = 0.;
                }
                break;
            case 2:
                for (int i = j1 - 1; i < jiu; i++) {
                    a[i - j] = psi[i] + psi[i + jd] + psi [i + jdm];
                    psi[i] = (psi[i] - psi[i + jh] - psi[i + jhm]) / 2.;
                } 
                break;
            case 3:
                for (int i = j1 - 1; i < jiu; i++) {
                    a[i - j] = 2.* psi[i];
                    psi[i] = psi[i + jd] + psi[i + jdm];
                }
                break;
            case 4:
                for (int i = j1 - 1; i < jiu; i++) {
                    double d = psi[i] - psi[i + jt] - psi[i + jtm];
                    psi[i] = psi[i] - psi[i + jh] - psi[i + jhm] + psi[i + jd] + psi[i + jdm];
                    a[i - j] = d + psi[i];  
                }
                break;
            default:
                break;
        }
        

        for (l = lo; l <= nn; l += li) {
            double d = 2. - tc[itc + l - 1];
            for (int i = 0; i < iu; i++) {
                int k = Mm1 - i - 2;
                double b = 1. / (d + tc[k] + tc[k + Mm1] * (1. - a[k + Mr]));
                a[k + Mm1] = tc[k] * b;
                a[k] = (a[k] + tc[k + Mm1] * a[k + 1] )* b;
            }
            for (int i = 1; i < iu; i++) {
                a[i] = a[i + Mm1] * a[i - 1] + a[i];
            }
        }

        for (int i = j1 - 1; i < jiu; i++) {
            psi[i] = psi[i] + a[i - j];
        }
        is = -1;

    }  

    if (iphase == 1) {
        goto L40;
    }
    else if (iphase == 2) {
        lo = 2 * lo;
        if (lo < nn)
            goto L30;
    }
    else if (iphase == 3 || iphase == 4) {
        lo = lo/2;
        if (lo == 1) 
            mode = 1;
        goto L30;    
    }    

L40:
    delete[] tc;
    delete[] a;
    
}