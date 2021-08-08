#include <cmath>
#include "geq.h"

using namespace std;

double gfl(double, double, double, double);
double k(double, double);
double e(double, double);

void bndmat() {

    auto *rt = new double[llp];
    auto *zt = new double[llp];

    double ar = 2.0 / (double) Mm1;
    double r0 = 1.0 / alpha;
    double az = sh * ar;
/*
    Index vector interior neighbours to boundary, bottom-top
*/

    int j1 = Mr;
    int j2 = Mr * (Nz - 2);
    int kp = 0;

    for (int i = 2; i <= Mm1; i++) {
        *(ip + kp) = i + j1 - 1;
        *(ip + kp + 1) = i + j2 - 1;
        kp += 2;
    }

/*
    Left - Right
*/

    for (int j = Mr; j <= j2; j += Mr) {
        ip[kp] = j + 1;
        ip[kp + 1] = j + Nm1 - 1;
        kp += 2;
    }

/*
     r - coordinates of boundary points, index vector
*/

    double ra = r0 - 1.;
    double rb = r0 + 1.;
    double za = 0.;
    double zb = az * (double) Nm1;
    int m2 = Mr - 2;
    j2 = Mr * Nm1;

/*
    Bottom - Top
*/
    int nl = 0;
    for (int i = 1; i <= m2; i++) {
        rt[nl] = ra + ar * (double) i;
        rt[nl + 1] = rt[nl];
        zt[nl] = za;
        zt[nl + 1] = zb;
        jp[nl] = i;
        jp[nl + 1] = i + j2;
        nl += 2;
    }

/*
    Left - Right
*/

    int n2 = Nz - 2;
    for (int j = 1; j <= n2; j++) {
        rt[nl] = ra;
        rt[nl + 1] = rb;
        zt[nl] = az * (double) j;
        zt[nl + 1] = zt[nl];
        jp[nl] = j * Mr;
        jp[nl + 1] = (j + 1) * Mr - 1;
        nl += 2;
    }

/*
    Matrix elements

*/
    double arh = 0.5 * ar ;

    for (int i = 0; i < kp; i++) {

/*
    Bottom - Top
*/

        nl = 0;
        for (int j = 0; j < m2; j++) {
            aux[nl][i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], ar) / (sh * rt[nl]);
            aux[nl + 1][i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], ar)/ (sh * rt[nl + 1]);
            nl += 2;
        }
/*
    Left - Right
*/
        for (int j = 0; j < n2; j++) {
            aux[nl][i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], az) * sh / (rt[nl] + arh);
            aux[nl + 1][i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], az)* sh /(rt[nl + 1] - arh);
            nl += 2;
        }

    }

    delete[] rt;
    delete[] zt;

}

double gfl(double rv, double rst, double zv, double del)
{

    double xdl;
    double ak = 4.0 * rv * rst /((rv + rst)*(rv + rst) + zv*zv);
    double x = ((rv - rst)*(rv - rst) + zv*zv)/((rv + rst)*(rv + rst) + zv*zv);
    if (x == 0.) {
        xdl = 2.0 *(log(del/(4. * rv)) - 1.);
    }
    else {
        xdl = log(x);
    }
    return sqrt(rv * rst/ak) * ((1. - 0.5 * ak)*k(x, xdl) - e(x, xdl)) / pi;

}

double p1(double x) {
    return (((.01736506451*x + .04757383546)*x + .06260601220)*x + .44325141463)*x + 1.0;
}
double p2(double x) {
    return (((.00526449639*x + .04069697526)*x + .09200180037)*x + .24998368310)*x;
}
double p3(double x) {
    return (((.01451196212*x + .03742563713)*x + .03590092383)*x + .09666344259)*x + 1.38629436112;
}
double p4(double x) {
    return (((.00441787012*x + .03328355346)*x + .06880248576)*x + .12498593597)*x + .5;
}
double e(double x, double xdl) {
    return p1(x) - p2(x) * xdl;
}
double k(double x, double xdl) {
    return p3(x) - p4(x) * xdl;
}

