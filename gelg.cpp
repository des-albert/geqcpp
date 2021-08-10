#include <cmath>

using namespace std;

void gelg(double *rr, double **a, int m, int n, double eps, int ier) {

    double piv, pivi, tb;
    int i, ii, j, k, l, ll, mm, nm, ist;

    ier = 0;
    piv = 0.;
    mm = m * m;
    nm = n * m;
    auto *ag = new double[mm];

    l = 0;
    for (j = 0; j < m; j++) {
        for (i = 0; i < m; i++) {
            ag[l++] = a[i][j];
        }
    }
    /*
     Search for greatest element in matrix
    */

    for (l = 0; l < mm; l++) {
        tb = fabs(ag[l]);
        if (tb > piv) {
            piv = tb;
            i = l + 1;
        }
    }
    double tol = eps * piv;

    /*
     a(i) is pivot element. piv contains the absolute value of a(i).

     Start elimination loop
    */

    int lst = 1;

    for (k = 1; k <= m; k++) {
        if (piv <= 0) goto Lend;
        if (ier == 0) {
            if (piv <= tol) {
                ier = k - 1;
            }
        }
        pivi = 1.0 / ag[i - 1];
        j = (i - 1) / m;
        i = i - j * m - k;
        j = j + 1 - k;

        /*
            i+k is row-index, j+k column-index of pivot element
            pivot row reduction and row interchange in right-hand side
        */

        for (l = k - 1; l < nm; l += m) {
            ll = l + i;
            tb = pivi * rr[ll];
            rr[ll] = rr[l];
            rr[l] = tb;
        }
        /*
          Is elimination terminated
        */

        if (k >= m) goto Lback;

        /*
          Column interchange in matrix a
        */

        int lend = lst + m - k;
        if (j > 0) {
            ii = j * m;
            for (l = lst - 1; l < lend; l++) {
                tb = ag[l];
                ll = l + ii;
                ag[l] = ag[ll];
                ag[ll] = tb;
            }
        }

        for (l = lst - 1; l < mm; l += m) {
            ll = l + i;
            tb = pivi * ag[ll];
            ag[ll] = ag[l];
            ag[l] = tb;
        }
        /*
        Save column interchange information
        */

        ag[lst - 1] = j;

        piv = 0.;
        lst = lst + 1;
        j = 0;
        for (ii = lst; ii <= lend; ii++) {
            pivi = -ag[ii - 1];
            ist = ii + m;
            j = j + 1;
            for (l = ist - 1; l < mm; l += m) {
                ll = l - j;
                ag[l] = ag[l] + pivi * ag[ll];
                tb = fabs(ag[l]);
                if (tb > piv) {
                    piv = tb;
                    i = l + 1;
                }
            }
            for (l = k - 1; l < nm; l += m) {
                ll = l + j;
                rr[ll] = rr[ll] + pivi * rr[l];
            }
        }
        lst = lst + m;
    }

    Lback:

    ist = mm + m;
    lst = m + 1;

    for (i = 2; i <= m; i++) {
        ii = lst - i;
        ist = ist - lst;
        l = ist - m;
        l = (int) floor(ag[l - 1] + 0.5);
        for (j = ii - 1; j < nm; j += m) {
            tb = rr[j];
            ll = j + 1;
            for (k = ist - 1; k < mm; k += m) {
                ll = ll + 1;
                tb = tb - ag[k] * rr[ll - 1];
            }
            k = j + l;
            rr[j] = rr[k];
            rr[k] = tb;
        }
    }
    delete[] ag;
    return;
    Lend:
    ier = -1;

}


