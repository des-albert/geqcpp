#ifndef EXTERN
#define EXTERN extern
#endif

#define Mc 15
#define Ng 1
#define Nmax 6
#define llmax 16

EXTERN int Mr, Nz, MN, Mm1, Nm1, llp, meshfg, mprfg, *ip, *jp, *ityp;
EXTERN int ndes, jaxis, naxis, icycle, idecis, Mmax, mpnmax, icops;
EXTERN double alpha, sh, pi, **aux, dr, dz, ss, *R, *Z, value;
EXTERN double Rmax, Rmin, Zmax, Zmin, error, Rmpl, Zmpl, Apl;
EXTERN double totcurr, betapol, alfac, alph, *com, psicon;
EXTERN double **bb, **eb, **cl, *fk, *g;

