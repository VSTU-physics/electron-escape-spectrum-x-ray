#include "plots.h"

unsigned wxt = 0;

void plot_analytics(bool A, bool P, bool T)
{
    if (!(A || P || T))
        return;

    FILE* fd = popen("gnuplot -p", "w");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "set title 'Зависимость l_tr(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:2 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:2 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:2 lw 2 with lines title 'T'\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "set title 'Зависимость dE/dS(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:3 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:3 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:3 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "set title 'Спектр n(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:4 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:4 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:4 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "set title 'Зависимость пробега R(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:5 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:5 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:5 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "set title 'Граничное условие bc(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:6 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:6 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:6 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    pclose(fd);
}

void plot_mc()
{
    FILE* fd = popen("gnuplot -p", "w");
    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Спектр n(E)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 1:2 lw 2 with lines\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость числа вышедших электронов от пробега n(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:4 lw 2 with lines\n");

    fprintf(fd, "set terminal wxt %d\n", wxt++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость средней энергии вышедших электронов от пробега E(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:5 lw 2 with lines\n");

    pclose(fd);
}
