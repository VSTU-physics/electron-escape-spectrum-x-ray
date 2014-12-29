#include <ctime>
#include <cstdlib>
#include "plots.h"
#include "spectrum.h"
#include "INIReader.h"


int main()
{
    #ifdef __WIN__
        system("chcp 65001");
    #endif // __WIN__
    srand (time(NULL));
    INIReader reader("config.ini");

    if (reader.ParseError() < 0) {
        printf("Can't load 'config.ini'\n");
        return 1;
    }

    int N = 1000;
    int M = 5000;

    int z;
    bool A, P, T;
    z = reader.GetInteger("analytical", "z", 0);
    if (z)
    {
        double l = reader.GetReal("analytical", "l", 0.001);
        double ecut = reader.GetReal("analytical", "ecut", 400);
        analytical(z, M, N, l, ecut);
        A = true;
    }

    z = reader.GetInteger("table", "z", 0);
    if (z)
    {
        double l = reader.GetReal("table", "l", 0.001);
        double ecut = reader.GetReal("table", "ecut", 400);
        table(z, M, N, l, ecut);
        T = true;
    }

    z = reader.GetInteger("approximation", "z", 0);
    if (z)
    {
        double l = reader.GetReal("approximation", "l", 0.001);
        double ecut = reader.GetReal("approximation", "ecut", 400);
        approximation(z, M, N, l, ecut);
        P = true;
    }

    z = reader.GetInteger("monte-carlo", "z", 0);
    if (z)
    {
        double l = reader.GetReal("monte-carlo", "l", 0.00001);
        double s = reader.GetReal("monte-carlo", "s", 0.0001);
        double ecut = reader.GetReal("monte-carlo", "ecut", 10);
        int nparticles = reader.GetInteger("monte-carlo", "particles", 400);
        int ntimes = reader.GetInteger("monte-carlo", "times", 400);
        int N = 100;
        monte_carlo(z, nparticles, ntimes, N, ecut, s, l);
        plot_mc();
    }

    plot_analytics(A, P, T);
    return 0;
}
