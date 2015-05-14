

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
    INIReader reader("mc.ini");

    if (reader.ParseError() < 0) {
        printf("Can't load 'mc.ini'\n");
        return 1;
    }

    auto z = reader.GetInteger("monte-carlo", "z", 0);
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
    return 0;
}
