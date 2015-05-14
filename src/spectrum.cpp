/*!
    \file
    Содержит функции для поиска необходимых массивов для последующего
    построения графиков и обработки.
*/
#include <cstdio>
#include <cmath>
#ifdef __WIN__
    #include <windows.h>
#endif // __WIN__
#include "physics.h"
#include "calculations.h"
#include "parse.h"

/*!
    \brief Выводит данные в консоль, проверяет, что всё загружено
    правильно.

    \param[in] Z Номер элемента
    \param[in] s Структура с данными по веществу
*/
void check_data(int Z, subst_t s/*, approx_t ap*/)
{
    printf("Загружены данные для Z = %d:\n", Z);
    printf("Свойства материала:\n");
    printf("Зарядовое число: %d\n", s.Z);
    printf("Плотность: %f\n", s.rho);
    printf("Молярная масса: %f\n", s.M);
    printf("Работа выхода: %f\n", s.U0);
    printf("\n");
    printf("Характеристики электронов:\n");
    printf("Кол-во типов: %d\n", s.N);
    printf("Энергии: [");
    for (int i = 0; i < s.N; i++) {
        printf(" %f", s.E[i]);
    }
    printf(" ]\n");
    printf("Вероятности: [");
    for (int i = 0; i < s.N; i++) {
        printf(" %f", s.P[i]);
    }
    printf(" ]\n");
    //printf("Для метода аппроксимации:\n");
    //printf("Кол-во типов: %d\n", ap.N);
    //printf("Энергии: [");
    //for (int i = 0; i < ap.N; i++) {
        //printf(" %f", ap.E[i]);
    //}
    //printf(" ]\n");
    //printf("alpha: [");
    //for (int i = 0; i < ap.N; i++) {
        //printf(" %f", ap.alpha_i[i]);
    //}
    //printf(" ]\n");
    //printf("d1: [");
    //for (int i = 0; i < ap.N; i++) {
        //printf(" %f", ap.d1_i[i]);
    //}
    //printf(" ]\n");
    //printf("R0: [");
    //for (int i = 0; i < ap.N; i++) {
        //printf(" %f", ap.R0_i[i]);
    //}
    //printf(" ]\n");
}

/*!
    \brief Решает задачу методом решения дифференциального уравнения
    на заданной сетке.

    \param[in] s Информация о веществе
    \param[in] N Размер сетки по координате \f$ z \f$
    \param[in] z Массив, содержащий сетку по \f$ z \f$
    \param[in] M Размер сетки по энергии \f$ E \f$
    \param[in] E Сетка по \f$ E \f$
    \param[in] ltrs Массив \f$ \lambda_{tr}(E) \f$
    \param[in] epss Массив \f$ \bar{\varepsilon} \f$
    \param[in] I1s Массив интегралов из граничных условий \f$ I_1(E) \f$
    \param[in] I2s Массив интегралов из граничных условий \f$ I_2(E) \f$
    \param[out] fs Рассчитанный спектр \f$ f(E) \f$

    Решается дифференциальное уравнение следующего вида:
    \f[
        \frac{\partial f}{\partial E} =
        - \frac{\lambda_{tr}(E)}{3\bar{\varepsilon}(E)} \frac{\partial^2 f}{\partial z^2} -
        \sum\limits_{i = 1}^N P_i \delta(E - E_i)
    \f]
    При граничных условиях:
    \f{align*}{
        & \frac{I_1(E)}{2 - 3I_2(E)} f -
        \frac{\lambda_{tr}}{3} \frac{\partial f}{\partial z}
        \Bigg|_{z = 0} = 0 \\
        & f \Bigg|_{z = \infty} = 0
    \f}
    Так как численным методом без преобразований использовать бесконечность
    невозможно, схема решает уравнение от 0 до некоторого \f$ z_{max} \f$. Следует
    понимать, что выбор \f$ z_{max} \f$ один из самых трудных вопросов задачи.
    При большом \f$ z_{max} \f$ шаг сетки оказывается слишком большим, коэффициенты велики
    и метод даёт большую ошибку. При малом \f$ z_{max} \f$ мы не можем утверждать, что
    \f$ z_{max} \approx \infty \f$, то есть решается совершенно другая задача не эквивалентная
    первой. В этом состоит ограничение метода конечных разностей. \f$ z_{max} \f$ нужно выбирать таким,
    чтобы при его небольшом увеличении решение не меняло своего характера. В данной работе выбор \f$ z_{max} \f$
    не проводился. Было выбрано значение, при котором результаты похожи на правду.

    Вторая проблема, которая возникла в ходе работы, - это использование \f$ \delta \f$-функций.
    Сингулярности \f$ \delta \f$-функций приводят к тому, что, заменяя \f$ \delta \f$-функцию её
    приближённым значением, требования к сетке сильно возрастают, и с её уплотнением начинает
    расти решение  \f$ f \f$. Чтобы избежать таких эффектов, использовался стандартный метод решения
    задач с \f$ \delta \f$-функциями. Точное решение задачи без источников можно искать на участках
    \f$ (E_i, E_{i+1}) \f$. Источник в форме \f$ \delta \f$-функции удет приводить к тому, что при
    пересечении \f$ E_i \f$ функция испытывает разрыв. Величина разрыва:
    \f{gather*}{
        \lim\limits_{\eta \to 0} \Big[ f(E_i + \eta) - f(E_i - \eta) \Big] = \lim\limits_{\eta \to 0}
        \int\limits_{E_i - \eta}^{E_i + \eta} \Bigg[- \frac{\lambda_{tr}(E)}{3\bar{\varepsilon}(E)}
        \frac{\partial^2 f}{\partial z^2} - \sum\limits_{j = 1}^N P_j \delta(E - E_j) \Bigg] dE
        = \\ =
        - \frac{\lambda_{tr}(E_i)}{3\bar{\varepsilon}(E_i)}
        \frac{\partial^2 f(E_i, z)}{\partial z^2} \lim\limits_{\eta \to 0} 2\eta
        - P_i = - P_i
    \f}
    Массив энергий E[i] упорядочен по убыванию. Поэтому в точках разрыва
    \f[
        \lim\limits_{\eta \to 0} \Big[ f(E_i - \eta) - f(E_i + \eta) \Big] = P_i
    \f]
    \code
        for (int k = 0; k < s.N; k++)
        {
            if (fabs(E[i] - s.E[k]) < fabs(dE/2))
            {
                for (int j = 0; j < N; j++) u[j] += s.P[k];
            }
        }
    \endcode
*/
void solve(/*const char* fname,*/ subst_t s, int N, double* z,
           int M, double* E, double* ltrs, double* epss,
           double* I1s, double* I2s, double* fs)
{
    //FILE* fd = fopen(fname, "w");
    double* tmp;
    double* u = new double[N];
    double* up = new double[N];
    for (int i = 0; i < N; i++) u[i] = 0;

    for (int i = 1; i < M; i++)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        double dE = E[i] - E[i - 1];

        spe(u, up, z, N, dE,
            - 1./3 * ltrs[i] / epss[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        for (int k = 0; k < s.N; k++)
        {
            if (fabs(E[i] - s.E[k]) < fabs(dE/2))
            {
                for (int j = 0; j < N; j++) u[j] += s.P[k];
            }
        }
        fs[i] = 3 * I1s[i] / (2 - 3 * I2s[i]) * u[0]; // *
        //if (i % 100 == 0)
        //{
            //for (int j = 0; j < N; j += 10)
            //{
                //fprintf(fd, "%e %e %e\n", E[i], z[j], u[j]);
            //}
            //fprintf(fd, "\n");
        //}
    }
    //fclose(fd);
}

/*!
    \brief Проводит расчёт коэффициентов аналитическим методом, решает
    уравнение и готовит данные для последующей обработки.

    \param[in] Z Номер элемента
    \param[in] M Число узлов сетки по \f$ E \f$
    \param[in] N Число узлов сетки по \f$ z \f$
    \param[in] l Максимальная длина промежутка \f$ z_{max} \f$
    \param[in] Emin Минимальная энергия (ограничение требуется, так
    как коэффициенты ведут себя неадекватно при \f$ E \to 0 \f$)


*/
void analytical(int Z, int M, int N, double l, double Emin)
{
    //Грузим данные о веществе
    subst_t s;
    load_subst(Z, "data/subst.dat", s);
    double Emax = 1.1 * s.E[0];

    //Готовим сетки для расчётов
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    //Выделяем память под массивы
    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    //Ищем коэффициенты по аналитическим формулам
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        ltrs[i] = l_tr(s, E[i]);
        epss[i] = eps(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    //Решаем уравнение
    solve(s, N, z, M, E, ltrs, epss, I1s, I2s, fs);

    //Находим средний пробег в зависимости от энергии
    rs[M-1] = 0;
    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }

    //Выводим данные в файл
    FILE* fd = fopen("data_a.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);

    //Чистим память
    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}

/*!
    \brief Проводит расчёт коэффициентов с помощью данных из таблиц института Иоффе, решает
    уравнение и готовит данные для последующей обработки.

    \param[in] Z Номер элемента
    \param[in] M Число узлов сетки по \f$ E \f$
    \param[in] N Число узлов сетки по \f$ z \f$
    \param[in] l Максимальная длина промежутка \f$ z_{max} \f$
    \param[in] Emin Минимальная энергия

*/
void table(int Z, int M, int N, double l, double Emin)
{
    //Грузим данные о веществе
    subst_t s;
    load_subst(Z, "data/subst.dat", s);
    double Emax = 1.1 * s.E[0];

    //Готовим сетки для расчётов
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    //Выделяем память под массивы
    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    //Загружаем данные из таблиц института Иоффе
    load_ltr(ltrs, E, M, "data/", s);
    load_eave(epss, E, M, "data/", s);

    //Вычисляем оставшиеся коэффициенты по аналитическим формулам
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    //Ищем средний пробег в зависимости от энергии
    rs[M-1] = 0;
    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }

    //Решаем уравнение
    solve(s, N, z, M, E, ltrs, epss, I1s, I2s, fs);

    //Выводим данные в файл
    FILE* fd = fopen("data_t.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);

    //Чистим память
    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}

/*!
    \brief Проводит расчёт коэффициентов с помощью аппроксимирующих формул для данных
    из таблиц института Иоффе, решает
    уравнение и готовит данные для последующей обработки.

    \param[in] Z Номер элемента
    \param[in] M Число узлов сетки по \f$ E \f$
    \param[in] N Число узлов сетки по \f$ z \f$
    \param[in] l Максимальная длина промежутка \f$ z_{max} \f$
    \param[in] Emin Минимальная энергия

*/
void approximation(int Z, int M, int N, double l, double Emin)
{
    //Грузим данные о веществе и аппроксимации
    subst_t s;
    approx_t ap;
    load_subst(Z, "data/subst.dat", s);
    load_approx(Z, "data/approx.dat", ap);
    double Emax = 1.1 * s.E[0];

    //Готовим сетки
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    //Выделяем память под массивы
    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    //Ищем нужные коэффициенты по аппроксимирующим формулам
    le_approx(ltrs, epss, E, M, ap, s);

    //Ищем оставшиеся коэффициенты по аналитическим формулам
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    //Ищем средний пробег в зависимости от энергии
    rs[M-1] = 0;
    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }

    //Решаем уравнение
    solve(s, N, z, M, E, ltrs, epss, I1s, I2s, fs);

    //Выводим данные в файл
    FILE* fd = fopen("data_p.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);

    //Чистим память
    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}

/*!
    \brief Рассчитывает интегральную функцию выхода \f[ K(t) \f].

    \param[in] Z Номер элемента
    \param[in] M Число узлов сетки по \f$ E \f$
    \param[in] N Число узлов сетки по \f$ z \f$
    \param[in] l Максимальная длина промежутка \f$ z_{max} \f$
    \param[in] Emin Минимальная энергия

*/
void quit_function(int Z, int M, int N, double l, double Emin)
{
    //Грузим данные
    subst_t s;
    load_subst(Z, "data/subst.dat", s);
    double Emax = 1.1 * s.E[0];

    //Готовим сетки
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    //Выделяем память
    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    //Ищем коэффициенты аналитическим методом
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        ltrs[i] = l_tr(s, E[i]);
        epss[i] = eps(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    int L = 20;
    double* Kt = new double[L];
    for (int k = 1; k<=L; k++)
    {
        for (int i = 0; i<N; i++)
            z[i] = i * k * l / (N - 1) / L;
        solve(s, N, z, M, E, ltrs, epss, I1s, I2s, fs);
        Kt[k-1] = 0;
        for (int i = 0; i < M; ++i)
            Kt[k-1] += fs[i];
    }

    FILE* fd = fopen("Kt.dat", "w");
    fprintf(fd, "%f %f\n", .0, .0);
    for (int k = 0; k<L; k++)
        fprintf(fd, "%f %f\n", l*(k+1)/L, Kt[k] / Kt[L-1]);
    fclose(fd);

    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] bs;
    delete [] fs;
    delete [] E;
    delete [] Kt;
}

/*!
    \brief Реализация метода Монте-Карло. Функция рассчитывает число частиц с энергией в данной ячейке сетки по энергиям, покинувшим образец и
    число частиц с данным пробегом, покинувшим образец.

    \param[in] Z Номер элемента
    \param[in] nparticles Число частиц
    \param[in] ntimes Число столкновений
    \param[in] N Число узлов сетки по энергиям
    \param[in] Emin Минимальная энергия
    \param[in] Smax Максимальный пройденный путь до выхода из образца
    необходим для сетки по \f$ S \f$
    \param[in] lmax Образец представляет собой куб со стороной \f$ l_{max} \f$

    Функция использует следующие результаты теории вероятностей.
    Пусть известна плотность распределения \f$ f(x) \f$, и мы можем получить величину \f$ y \f$ равномерно распределённую
    в промежутке \f$ [0, l] \f$. Если определить функцию \f$ G \f$, такую что:
    \f[
        G : y \in [0, l] \rightarrow x \in [-\infty, \infty],
    \f]
    то этой функцией будет:
    \f[
        x = G(y) = F^{-1}(y),
    \f]
    где \f$F(x)\f$ - первообразная функции \f$f(x)\f$, обладающая следующими свойствами:
    \f[
        F(-\infty) = 0,
    \f]
    \f[
        F(\infty) = l.
    \f]
    Очевидно, что \f$ F(x) = lP(x) \f$, где \f$ P(x) \f$ - функция распределения величины \f$ x \f$
    \f[
        P(x) = \int\limits_{-\infty}^x f(\xi) d\xi.
    \f]

    Будем считать что электроны, появивишиеся в результате Оже-процессов, двигаются прямолинейно и равномерно
    до столкновения. В качестве длины пути такого электрона с энергией \f$ E \f$ будем использовать длину
    \f$ \lambda \f$:
    \f[
        \lambda^{-1} = \lambda_{el}^{-1} + \lambda_{in}^{-1}.
    \f]
    Также полагаем, что в результате упругих столкновений энергия электрона не меняется, но меняется направление
    движения, а в результате неупругих столкновений направление движения не меняется, но меняется энергия и эти
    процессы протекают независимо друг от друга. После упругого столкновения случайными величинами являются два угла:
    угол \f$ \theta' \f$ между новым направлением движения и старым направлением движения и угол \f$ \varphi' \f$
    между проекцией нового направления
    движения на плоскость перпендикулярную старому направлению движения и проекцией оси \f$ x \f$ на эту же плоскость.
    Сечение упругого рассеяния зависит только от \f$ \theta' \f$ и \f$ E \f$. Как следствие функция распределения
    угла \f$ \varphi' \f$ равномерна на участке от \f$ 0 \f$ до \f$ 2\pi \f$.
    Функция распределения по \f$ \theta' \f$:
    \f[
        P(\theta', E) =
        \frac{\int\limits_0^{\theta'} \sigma(\theta, E) d\theta}
        {\int\limits_0^{\pi} \sigma(\theta, E) d\theta}.
    \f]
    Выражаем её приближённо:
    \f[
        P(\theta', E) \approx \frac{2}{\pi} \arctan (\alpha(E) \theta').
    \f]
    И для угла:
    \f[
        \theta' = \frac{1}{\alpha(E)} \tan \left( \frac{\pi X}{2} \right),
    \f]
    здесь \f$ X \f$ - равномерно-распределённая величина на промежутке \f$ [0, 1] \f$.
    Функция распределения по энергиям неупругого рассеяния \f$ Q \f$:
    \f[
        P(Q, E) = \frac{\int\limits_0^{Q} \frac{d\lambda^{-1}_{in}}{d Q} dQ}
        {\int\limits_0^{E/2} \frac{d\lambda^{-1}_{in}}{d Q} dQ}
        \approx \frac{2}{\pi} \arctan \left(\beta(E) \frac{Q}{E}\right).
    \f]
    и
    \f[
        Q = \frac{E}{\beta(E)} \tan \left( \frac{\pi X}{2} \right).
    \f]
    Параметры \f$ \alpha(E) \f$ и \f$ \beta(E) \f$ находятся по данным из таблиц Иоффе.
    Требовалось также учесть отражение на границе. Для этого использовалась функция reflection.
    Существует ненулевая работа выхода из образца \f$ U_0 \f$. Если энергия частицы \f$ E \cos^2 \theta < U_0 \f$,
    то частица образец не покинет, в противном случае, можно рассматривать отражение частицы, как отражение
    от барьера высотой \f$ U_0 \f$. В этом случае действуют вероятностные законы. Которые и используются в функции
    reflection.
    Осталось только привести формулы, по которым можно рассчитать новые углы \f$ \theta_n \f$ и \f$ \varphi_n \f$, по
    известным старым углам \f$ \theta \f$ и \f$ \varphi \f$ и новым углам относительно старых \f$ \theta' \f$ и \f$ \varphi' \f$:
    \f[
        \tan \varphi_n =
        \frac{
        \sin (\theta) \sin (\varphi) \cos (\theta') -
        \cos (\varphi) \sin (\theta') \cos (\varphi') +
        \cos(\theta) \sin(\varphi) \sin(\theta') \sin(\varphi')
        }
        {
        \sin(\theta) \cos(\varphi) \cos(\theta') +
        \sin(\varphi) \sin(\theta') \cos(\varphi') +
        \cos(\theta) \cos(\varphi) \sin(\theta')  \sin(\varphi')
        }
    \f]
    \f[
        \cos \theta_n = \cos(\theta') \cos(\theta) - \sin(\theta') \sin(\theta) \sin(\varphi')
    \f]
    В чём возможные ошибки метода?
    1. Неправильно выбрана аппроксимирующая функция для вероятностей (углы \f$ \theta \f$ лежат в промежутке
    от 0 до \f$ \pi \f$, а при данной аппроксимации от 0 до \f$ \infty \f$, анлогично и для энергий).
    2. Не учитывается форма распределений по пробегам (а логарифм (log(rand)), вообще говоря, использовать нельзя).
    3. Соотношения для отражений не проверялись, а взяты как есть из статьи.

*/
void monte_carlo(int Z, int nparticles, int ntimes, int N, double Emin, double Smax, double lmax)
{
    double *E, *beta, *alpha, beta_t, alpha_t, theta_s, phi_s, theta_t, phi_t;
    E = new double[N];
    beta = new double[N];
    alpha = new double[N];

    subst_t s;
    load_subst(Z, "data/subst.dat", s);
    double Emax = 1.1 * s.E[0];
    double dE = (Emax - Emin) / (N - 1);
    for (int i = 0; i < N; i++)
        E[i] = Emin + dE * i;

    double inv_lambda_el, *inv_lambda_in;
    inv_lambda_in = new double[N];

    load_mc_elastic(alpha, E, N, inv_lambda_el, "data/", s);
    load_mc_inelastic(beta, E, N, inv_lambda_in, "data/", s);

    double lambda;

    particle_t p;

    double x, y, z, dS;
    int k;
    dS = Smax / N;

    double *n_E, *n_S, *E_S;
    double E_exit;
    n_E = new double[N]; for (int i = 0; i<N; i++) n_E[i] = 0;
    n_S = new double[N]; for (int i = 0; i<N; i++) n_S[i] = 0;
    E_S = new double[N]; for (int i = 0; i<N; i++) E_S[i] = 0;
    int general_sum = 0;

    for (int i = 0; i<nparticles; i++)
    {
        bool stop = true;
        double S = 0;
        p.x = - random(lmax/2) + lmax;
        p.y = - random(lmax/2) + lmax;
        p.z = random(lmax);
        p.theta = random(M_PI);
        p.phi = random(2*M_PI);
        double rand = random(1);
        p.E = auger_source(s, rand);
        for (int j = 0; j<ntimes && stop; j++)
        {
            k = (int)((p.E - Emin) / dE);
            lambda = 1 / (inv_lambda_el + inv_lambda_in[k]);

            rand = random(1);
            x = p.x - lambda * log(rand) * sin(p.theta) * cos(p.phi);
            y = p.y - lambda * log(rand) * sin(p.theta) * sin(p.phi);
            z = p.z - lambda * log(rand) * cos(p.theta);
            S -= lambda * log(rand);
            if (z<0)
            {
                if (reflection(p.E, p.theta, s.U0))
                {
                    p.x = p.x + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * sin(p.theta) * cos(p.phi);
                    p.y = p.y + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * sin(p.theta) * sin(p.phi);
                    p.z = p.z + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * cos(p.theta);
                    p.theta = M_PI - p.theta;
                    p.phi = M_PI + p.phi;
                } else {
                    E_exit = p.E - s.U0;
                    //theta_exit = sqrt(1 - pow(sin(p.theta), 2) / (1 - s.U0 / p.E));
                    k = (int)((E_exit - Emin) / dE);
                    if (k < N) n_E[k] += 1;
                    k = (int) (S/dS);
                    //printf("%d %d %e %e\n", general_sum, k, p.E, E_S[k]);
                    if (k < N)
                    {
                        n_S[k] += 1;
                        E_S[k] += E_exit;
                    }
                    stop = 0;
                    general_sum ++;
                    //
                }
            } else {
                p.x = x;
                p.y = y;
                p.z = z;
            }

            k = (int)((p.E - Emin) / dE);
            if (random(1) < inv_lambda_in[k] * lambda)
            {
                beta_t = linterp(p.E, E[k], beta[k], E[k+1], beta[k+1]);
                p.E = p.E*(1 - 1 / beta_t * tan(M_PI * random(1) / 2));
                if (p.E < Emin)
                {
                    stop = 0;
                }
            } else {
                alpha_t = linterp(p.E, E[k], alpha[k], E[k+1], alpha[k+1]);
                theta_s = 1 / alpha_t * tan(M_PI * random(1) / 2);
                phi_s = 2*M_PI*random(1);
                phi_t = atan2(sin(p.theta)*sin(p.phi)*cos(theta_s) -
                              cos(p.phi)*sin(theta_s)*cos(phi_s) +
                              cos(p.theta)*sin(p.phi)*sin(theta_s)*sin(phi_s),
                              sin(p.theta)*cos(p.phi)*cos(theta_s) +
                              sin(p.phi)*sin(theta_s)*cos(phi_s) +
                              cos(p.theta)*cos(p.phi)*sin(theta_s)*sin(phi_s));
                theta_t = acos(cos(theta_s)*cos(p.theta) - sin(theta_s)*sin(p.theta)*sin(phi_s));
                p.theta = theta_t;
                p.phi = phi_t;
            }
        }
    }

    FILE *fd;
    fd = fopen("results_mc.dat", "w");
    for (int i = 0; i<N; i++)
    {
        fprintf(fd, "%e %e %e %e %e\n",
                E[i],
                n_E[i],
                (i+1)*dS,
                n_S[i],
                ((int) n_S[i] == 0) ? 0 : E_S[i] / n_S[i]);
    }
    fclose(fd);

    delete [] E;
    delete [] beta;
    delete [] alpha;
    delete [] n_E;
    delete [] n_S;
    delete [] E_S;
}
