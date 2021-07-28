#include <cmath>
#include <iostream>
#include <fstream>

#define EXTERN

#include "geq.h"

using namespace std;

extern void bndmat();

int main()
{   
    Mc = 15;
    Ng = 1;
    Nmax = 6;
    
    double Offset, El, Rxpsn, Zxpsn, tri, elxp, trixp;
    double rac, zac, exc, Ra[Ng][Mc], Za[Ng][Mc], Ex[Ng][Mc], Rl[Ng][Mc];
    double Rc[Nmax], Zc[Nmax];
    int ic[Mc];

    cout << "Garching Tokamak Equlibrium" << endl;

    int Nexp = 6;
    Mr = 1 + (1<<Nexp);
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
    
    
    if (meshfg > 0 ) {
        Rmin = Rmpl - 32.0 *Apl / 20.0;
        Rmax = Rmpl + 32.0 *Apl / 20.0;
        Zmin = Offset - 32.0 * (Offset + abs(Zxpsn)) / 20.0;
        Zmax = Offset + 32.0 * (Apl*El) / 20.0;
    }

    R = new double[Mr];
    Z = new double[Nz];
    ip = new int[llp];
    jp = new int[llp];
    ityp = new int[6];
    aux = new double[llp * llp];


    dr = (Rmax - Rmin) / (double) Mm1;
    dz = (Zmax - Zmin) / (double) Nm1 ;

    for (int i = 0; i < Mr; i++) 
    {
       R[i] = Rmin + dr * (double)i; 
    }
    for (int j = 0; j < Nz; j++) 
    {
       Z[j] = Zmin + dz * (double)j; 
    }  

    alpha = (Rmax - Rmin) / (Rmax + Rmin);
    sh = dz / dr;
    ss = sh * sh;

    cout << "Rmin  = " << Rmin <<  " Rmax = " << Rmax << " dr = " << dr << endl;
    cout << "Zmin  = " << Zmin <<  " zmax = " << Zmax << " dz = " << dz << endl;
    cout << "alpha = " << alpha << "   sh = " << sh   << " ss = " << ss << endl;

    bndmat();  

 
/*
    Read in Poloidal Field Coil Data

*/

    int k = 0;

    while (true) {
        ic[k] = -1;
        cout << " ic " << ic[k] << endl;
        while (true) {
            fin >> rac >> zac >> exc;

            if ( rac <= 0.0) goto L10;
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

        if(mprfg != 0) {
            cout << "Conductor groups available for optimization" << endl;
            for (int k = 0; k < Mmax; k++) {
                cout << " Group : " << k + 1 << endl;
                for (int i = 0; i <= ic[k]; i++)  {
                    printf (" %7.3f   %7.3f   %7.3f \n", Ra[i][k], Za[i][k], Ex[i][k]);
                }
            }
        }

    /*
        Conditions to be satisfied by resulting equilibrium
    */    

    elxp = abs(Offset - Zxpsn)/Apl;
    trixp = (Rmpl - Rxpsn)/Apl;

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
    Rc[5] = Rmpl - Apl*tri;
    Zc[5] = Offset + Apl*El;

    fin.close();

}