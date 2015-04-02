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
    INIReader reader("analytics.ini");

    if (reader.ParseError() < 0) {
        printf("Can't load 'analytics.ini'\n");
        return 1;
    }

    int N = 1000;
    int M = 5000;

    int z;
    bool A = false, P = false, T = false;
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

    plot_analytics(A, P, T);
    return 0;
}
