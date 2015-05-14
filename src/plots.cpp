/*!
    \file
    В этом файле описаны функции для построения графиков.

    Графики строятся при помощи программы Gnuplot. В Windows используется Wxt,
    в остальных ОС - qt терминал.
*/

#include "plots.h"


// Gnuplot term: wxt for windows, qt for linux
#ifdef __WIN__
    char term[] = "wxt";
#else
    char term[] = "qt";
#endif // __WIN__

unsigned n = 0;

/*!
    \brief
    plot_analytics строит графики по данным, полученным каким либо аналитическим
    методом: чисто аналитическим, аппроксимацией или табличным. Данные для
    построения берутся из текстовых файлов, созданных на этапе рассчёта каждым
    из методов.

    \details
    \param[in] A Строить ли данные для чисто аналитического метода
    \param[in] P Строить ли данные для аппроксимационного метода
    \param[in] T Строить ли данные для табличного метода
*/
void plot_analytics(bool A, bool P, bool T)
{
    if (!(A || P || T))
        return;

    FILE* fd = popen("gnuplot -p", "w");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "set title 'Зависимость l_tr(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:2 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:2 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:2 lw 2 with lines title 'T'\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "set title 'Зависимость dE/dS(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:3 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:3 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:3 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "set title 'Спектр n(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:4 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:4 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:4 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "set title 'Зависимость пробега R(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:5 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:5 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:5 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "set title 'Граничное условие bc(E)'\n");
    fprintf(fd, "plot\\\n");
    if(A) fprintf(fd, "'data_a.dat' using 1:6 lw 2 with lines title 'A',\\\n");
    if(P) fprintf(fd, "'data_p.dat' using 1:6 lw 2 with lines title 'P',\\\n");
    if(T) fprintf(fd, "'data_t.dat' using 1:6 lw 2 with lines title 'T'\n");
    fprintf(fd, "\n");

    pclose(fd);
}

/*!
    \brief
    plot_mc строит графики по данным, полученным методом Монте-Карло. Данные
    берутся из соответствующего текстового файла.
*/
void plot_mc()
{
    FILE* fd = popen("gnuplot -p", "w");
    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Спектр n(E)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 1:2 lw 2 with lines\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость числа вышедших электронов от пробега n(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:4 lw 2 with lines\n");

    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость средней энергии вышедших электронов от пробега E(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:5 lw 2 with lines\n");

    pclose(fd);
}

/*!
    \brief
    plot_k  строит график интегральной функции выхода
    \f[
        K(t) = \frac{\sum_i P_i \int_0^t q_i(z) dz}
                    {\sum_i P_i \int_0^\infty q_i(z) dz},
    \f]
    где \f$ P_i \f$ -- вероятность рождения электрона \f$ i \f$-ой группы,
    \f$ q_i(z) \f$ -- плотность вероятности выхода этого электрона с глубины
    \f$ z \f$. Другими словами, это отношение числа электронов, вылетающих из
    образца с конечной толщиной к числу электронов, вылетающих из
    полубесконечного образца.
*/
void plot_k()
{
    FILE* fd = popen("gnuplot -p", "w");
    fprintf(fd, "set terminal %s %d\n", term, n++);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Интегральная функция выхода K(t)'\n");
    fprintf(fd, "plot 'Kt.dat' lw 2 with lines smooth bezier\n");
    pclose(fd);
}

