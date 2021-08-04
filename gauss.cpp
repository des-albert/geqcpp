#include "geq.h"
#include <algorithm>

using namespace std;

double **array2d(int, int);

void array2del(double **, int);

void gauss(double *b, double **a, int n) {

    double **ag = array2d(n, n + 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ag[i][j] = a[i][j];
        }
        ag[i][n] = b[i];
    }

    // Pivoting

    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            if (abs(ag[i][i]) < abs(ag[k][i])) {
                for (int j = 0; j <= n; j++) {
                    double temp = ag[i][j];
                    ag[i][j] = ag[k][j];
                    ag[k][j] = temp;
                }
            }
        }
    }

    //  Gauss elimination

    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            double t = ag[k][i] / ag[i][i];
            for (int j = 0; j <= n; j++) {
                ag[k][j] = ag[k][j] - t * ag[i][j];
            }
        }
    }

    // Back substitution

    for (int i = n - 1; i >= 0; i--) {
        b[i] = ag[i][n];
        for (int j = i + 1; j < n; j++) {
            if (j != i)
                b[i] = b[i] - ag[i][j] * b[j];
        }
        b[i] = b[i] / ag[i][i];
    }

    array2del(ag, n);
}

