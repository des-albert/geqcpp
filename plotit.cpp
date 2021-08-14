#include "geq.h"
#include "/usr/local/dislin/include/discpp.h"


void plotit() {

    Dislin d;

    float plev;
    float rf[MN], rr[Mr], rz[Nz], cjr[Mr], p[Mr];

    for (int i = 0; i < Mr; i++) {
        rr[i] = (float) R[i];
        cjr[i] = (float) cjt[i];
        p[i] = (float) pr[i];
    }
    for (int i = 0; i < Nz; i++) {
        rz[i] = (float) Z[i];
    }
    for (int i = 0; i < Mr; i++) {
        for (int j = 0; j < Nz; j++) {
            rf[i * Mr + j] = (float) g[i + Nz * j];
        }
    }

    d.setpag("ps4l");
    d.disini();
    d.complx();
    d.pagera();
    d.pagfll(255);
    d.color("blue");
    d.axslen(1120,1600);
    d.graf(0., 14., 0., 2., -10., 10., -10., 2.);
    d.xaxgit();
    d.color("green");
    d.rlrec(1.3135,6.0605,.749,1.979);
    d.rlrec(1.3135,4.0325,.749,1.979);
    d.rlrec(1.3135,2.0035,.749,1.979);
    d.rlrec(1.3135,-0.0245,.749,1.979);
    d.rlrec(1.3135,-2.0535,.749,1.979);
    d.rlrec(1.3135,-4.0815,.749,1.979);
    d.rlrec(3.47,8.045,.968,.976);
    d.rlrec(7.9845,6.8275,.649,.595);
    d.rlrec(11.581,3.8275,.708,1.125);
    d.rlrec(11.5805,-1.6805,.649,1.125);
    d.rlrec(7.985,-6.2575,.82,.945);
    d.rlrec(3.4705,-7.069,1.633,.976);
    d.color("blue");

    for (int i = 2; i < 16; i++) {
        plev = (float)(psicon + (i-1)*(fab - psicon)/15.);
        d.contur(rr, Mr, rz, Nz, rf, plev);
    }
    d.color("red");
    d.contur(rr, Mr, rz, Nz, rf, (float) psicon);
    d.endgrf();
    d.disfin();

    /*
        Current and Pressure profile
    */

    d.disini();
    d.complx();
    d.pagera();
    d.pagfll(255);
    d.color("blue");
    d.polcrv("linear");
    d.titlin("Current Density & Pressure", 1);
    d.graf(3., 10., 3., 1., 0., 2., 0., 1.);
    d.title();
    d.name("Major Radius (m)","X");
    d.color("red");
    d.curve(rr, cjr, Mr);
    d.color("green");
    d.curve(rr, p, Mr);
    d.disfin();
}

