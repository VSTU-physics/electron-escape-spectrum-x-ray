/*!
    \mainpage

    \section problem Постановка задачи

    Требуется решить уравнение:
    \f{equation}{
    \begin{gathered}
        \vec{\Omega}\cdot\nabla F(\vec{r}, \vec{\Omega}, E) =
        \frac{\partial}{\partial E}[\bar{\varepsilon}(E)F(\vec{r}, \vec{\Omega}, E)] + \\ +
        \int\limits_{4\pi} \lambda^{-1}_{el}(E, \vec{\Omega}', \vec{\Omega}) \big[ F(\vec{r}, \vec{\Omega}', E) - F(\vec{r}, \vec{\Omega}, E) \big] d\Omega' +
        Q(\vec{r}, \vec{\Omega}, E),
    \end{gathered}
    \f}
    где
    \f$\vec{\Omega}\f$ -- единичный вектор;
    \f$F(\vec{r}, \vec{\Omega}, E)\f$ --  плотность потока частиц в
            окрестности вектора \f$\vec{r}\f$, движущуюся в телесном угле \f$d\Omega\f$ в
            направлении \f$\vec{\Omega}\f$, с энергией в диапазоне \f$[E, E+dE]\f$;
    \f$\bar{\varepsilon}(E)\f$ -- средние потери на единице пути:
        \f{equation}{\label{eq: sharp}
            \bar{\varepsilon}(E) = \int\limits_{0}^{E/2} \varepsilon \frac{\partial \lambda_{in}^{-1}(E,
            \varepsilon)}{\partial \varepsilon} d\varepsilon;
        \f}
    \f$\lambda_{in}^{-1}(E, \varepsilon)\f$ -- вероятность потерять на единице длины энергию \f$\varepsilon\f$ (вероятность неупругого столкновения);
    \f$\lambda^{-1}_{el}(E, \vec{\Omega}', \vec{\Omega})\f$ -- вероятность изменения направления движения с \f$\vec{\Omega}'\f$ на \f$\vec{\Omega}\f$ на единице пути (вероятность упругого столкновения);
    \f$Q(\vec{r}, \vec{\Omega}, E)\f$ -- функция источника, для Оже-процессов имеет вид:
        \f{equation}{
            Q(\vec{r}, \vec{\Omega}, E) = \frac{1}{4\pi} \sum_{i} P_i\delta(E - E_i);
        \f}
    \f$P_i\f$ -- вероятность рождения Оже-электрона с энергией \f$E_i\f$.

    \f$P_1\f$ приближение:
    \f{gather}{\label{eq: P1}
        F(\vec{r}, \vec{\Omega}, E) = \frac{1}{4\pi} [ F_0(\vec{r}, E) + 3 \vec{\Omega}\cdot\vec{F_1}(\vec{r}, E) ].
    \f}

    Здесь:
    \f{align}{
        & F_0(\vec{r}, E) = \int\limits_{4\pi} F(\vec{r}, \vec{\Omega}, E) d\Omega, \\
        & \vec{F}_1(\vec{r}, E) = \int\limits_{4\pi} \vec{\Omega} F(\vec{r}, \vec{\Omega}, E) d\Omega.
    \f}

    Подставляя это выражение в основное уравнение, получим:
    \f{equation}{\label{eq:P1is}
    \begin{gathered}
        \vec{\Omega}\cdot\nabla [ F_0(\vec{r}, E) + 3 \vec{\Omega}\cdot\vec{F_1}(\vec{r}, E) ] =
        \frac{\partial}{\partial E}[\bar{\varepsilon}(E)[ F_0(\vec{r}, E) + 3 \vec{\Omega}\cdot\vec{F_1}(\vec{r}, E) ]] + \\ +
        \int\limits_{4\pi} \lambda^{-1}_{el}(E, \vec{\Omega}', \vec{\Omega}) 3 \big[ \vec{\Omega}' - \vec{\Omega}\big] \cdot \vec{F}_1 (\vec{r}, E) d\Omega' +
        \sum_{i} P_i\delta(E - E_i).
    \end{gathered}
    \f}

    Для упругого столкновения \f$\lambda^{-1}_{el}(E, \vec{\Omega}', \vec{\Omega})\f$ зависит только от угла между \f$\vec{\Omega}'\f$, \f$\vec{\Omega}\f$ -- \f$\theta\f$. Рассмотрим детально интеграл:
    \f{equation}{
    \begin{gathered}
        \int\limits_{4\pi} \lambda^{-1}_{el}(E, \vec{\Omega}', \vec{\Omega}) 3 \big[ \vec{\Omega}' - \vec{\Omega}\big] d\Omega' =
        \int\limits_{4\pi} \lambda^{-1}_{el}(E, \theta) 3 \big[ \vec{\Omega}'_{\perp} + \vec{\Omega} \cos \theta - \vec{\Omega}\big] \sin d\theta d\phi = \\
        = [\text{ интеграл от } \vec{\Omega}'_{\perp} \text{ по } \phi ~ = 0] =
        6\pi\int\limits_{0}^\pi \lambda^{-1}_{el}(E, \theta) \big[ \cos \theta - 1\big] \sin \theta d\theta \vec{\Omega} = \\ =
        - 3 \lambda^{-1}_{tr} \vec{\Omega},
    \end{gathered}
    \f}
    где обозначено:
    \f{equation}{
        \lambda^{-1}_{tr}(E) = 2\pi\int\limits_{0}^\pi \lambda^{-1}_{el}(E, \theta) \big[ 1 -\cos \theta \big] \sin \theta d\theta.
    \f}
    \f$\lambda_{tr}\f$ носит название транспортной длины.

    Далее рассмотрим выражение
    \f{equation}{
    \begin{gathered}
        \vec{\Omega}\cdot\nabla [\vec{\Omega}\cdot\vec{F_1}] =
        \vec{\Omega}\cdot[
        \vec{\Omega} \times \mathrm{rot\,} \vec{F_1} +
        (\vec{\Omega}\cdot\nabla) \vec{F_1} + 0
        ] = \vec{\Omega}(\vec{\Omega}\cdot\nabla)\cdot \vec{F_1} = \\
        = (\vec{\Omega}\times(\vec{\Omega}\times\nabla) + \nabla (\vec{\Omega}\cdot\vec{\Omega}))\cdot\vec{F_1} =
        (\vec{\Omega}\times(\vec{\Omega}\times\nabla) + \nabla)\cdot\vec{F_1}
    \end{gathered}
    \f}

    \b По неизвестной причине мы выкидываем \f$\vec{\Omega}\times(\vec{\Omega}\times\nabla)\f$. \b

    Направление \f$\vec{\Omega}\f$ произвольно, поэтому из линейной независимости в \f$\eqref{eq:P1is}\f$ следует:
    \f{align}{
        & \vec{F}_1 = -\frac{\lambda_{tr}}{3}\nabla F_0 + \lambda_{tr} \frac{\partial \bar{\varepsilon} \vec{F_1}}{\partial E}, \label{eq: P1-1}\\
        & \nabla\cdot\vec{F}_1 = \frac{\partial \bar{\varepsilon}F_0}{\partial E} + \sum_{i} P_i\delta(E - E_i). \label{eq: P1-2}
    \f}

    В диффузионном приближении пренебрегают последним слагаемым в \f$\eqref{eq: P1-1}\f$, тогда из \f$[\ref{eq: P1-1}-\ref{eq: P1-2}]\f$ следует:
    \f{equation}{
        \frac{\partial(\bar{\varepsilon}F_0)}{\partial E} + \sum_{i} P_i\delta(E - E_i) = -\frac{\lambda_{tr}}{3}\Delta F_0
    \f}


    \section analitics Аналитический метод

    В этом случае \f$\lambda_{tr}^{-1}\f$ считается аналитически:
    \f{equation}{
        \lambda_{tr}^{-1} =
        2 \pi \frac{\rho N_A}{M} \int\limits_{0}^{\pi}
        \frac{d\sigma_{el}}{d\Omega} (1 - \cos \theta) \sin \theta d \theta.
    \f}

    Модификация (кажется Бете) формулы для дифференциального сечения рассеяния (Резерфорда):
    \f{equation}{
        \frac{d\sigma_{el}}{d\Omega} =
        \frac{(T+1)^2 r_e^2}{T^2(T+2)^2} \frac{Z(Z+1)}{(1 - \cos \theta + 2 \eta)^2},
    \f}
    где
    \f$T\f$ - кинетическая энергия электрона, отнесённая к энергии покоя:
        \f{equation}{
            T = \frac{E}{mc^2};
        \f}
    \f$r_e\f$ -- классический радиус электрона:
        \f{equation}{
            r_e = \frac{e^2}{mc^2} (\text{В СГС}) = \frac{1}{4 \pi \varepsilon_0}\frac{e^2}{mc^2} (\text{В СИ});
        \f}
    \f$\eta\f$ -- поправка:
        \f{equation}{
            \eta = \frac{1}{4}
            \Big[
                \frac{\alpha Z^{1/3}}{0{,}885}
            \Big]^2
            \frac{1}{T(T+2)} \left(1{,}13 + 3{,}76 \alpha^2 Z^2 \frac{(T+1)^2}{T(T+2)} \right);
        \f}
    \f$\alpha\f$ -- постоянная тонкой структуры:
        \f{equation}{
            \alpha = \frac{e^2}{\hbar c} (\text{В СГС}) = \frac{1}{4 \pi \varepsilon_0} \frac{e^2}{\hbar c} (\text{В СИ});
        \f}
    \f$\rho\f$, \f$N_A\f$, \f$M\f$ -- плотность вещества, число Авогадро, молярная масса.

    После интегрирования:
    \f{equation}{
        \begin{gathered}
            \lambda_{tr}^{-1} =
            2 \pi \frac{\rho N_AZ(Z+1)}{M}
            \frac{(T+1)^2 r_e^2}{T^2(T+2)^2}
            \int\limits_{0}^{\pi}
            \frac{1}{(1 - \cos \theta + 2 \eta)^2} (1 - \cos \theta) \sin \theta d \theta
            = \\ =
            \left[
                \int\limits_{0}^{2} \frac{1}{(x + 2 \eta)^2} x\, dx =
                - \int\limits_{0}^{2}  x\, d\frac{1}{x + 2 \eta} =
                -\frac{2}{2 + 2 \eta} + \int\limits_{0}^{2} \frac{1}{x + 2 \eta} dx =
                -\frac{1}{1 + \eta} + \ln \frac{2 +  2 \eta}{2 \eta}
            \right]
            = \\ =
            2 \pi \frac{\rho N_AZ(Z+1)}{M}
            \frac{(T+1)^2 r_e^2}{T^2(T+2)^2}
            \left(
                \ln \left( \frac{1}{\eta} + 1 \right) - \frac{1}{\eta + 1}
            \right).
        \end{gathered}
    \f}

    \f$\bar{\varepsilon}\f$ находится по формуле Бете:
    \f{equation}{
        \bar{\varepsilon} =
        2 \pi \frac{\rho N_AZ}{M}
        \frac{(T+1)^2 r_e^2}{T(T+2)}
        \left(
            2 \ln
            \frac{1{,}166\,T}{J/mc^2}
        \right),
    \f}
    где
    \f{equation}{
        J [\text{эВ}] =
        \begin{cases}
            13{,}6 Z & Z<10, \\
            (9{,}76 + 58{,}8 Z^{-1{,}19}) Z & Z\geq 10.
        \end{cases}
    \f}

    \section table Табличный метод

    Формулы \f$\eqref{eq: sharp}\f$ и \f$\eqref{eq: ltr}\f$, данные из таблиц Иоффе.
*/

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
    bool A = false, P = false, T = false, K = false;
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

    z = reader.GetInteger("quit", "z", 0);
    if (z)
    {
        double l = reader.GetReal("quit", "l", 0.001);
        double ecut = reader.GetReal("quit", "ecut", 400);
        quit_function(z, M, N, l, ecut);
        K = true;
    }

    plot_analytics(A, P, T);

    if (K)
        plot_k();

    return 0;
}
