/*!
    \file
    Файл содержит функции для обработки файлов данных.
*/
#include "parse.h"

/*!
    \brief Загружает параметры вещества из файла вида:
    \code
    #Z,label,M,rho,U0,N,energies,probabilities
    14 Si 28.090 2.3210 18.470 3 [7.90 6.20 1.60] [0.08 0.92 0.88]
    32 Ge 72.600 5.3500 12.500 3 [7.90 6.83 1.16] [0.12 0.88 0.86]
    \endcode

    \details
    \param[in] Z Номер элемента в периодической системе.
    \param[in] ch Имя файла из директории запуска.
    \param[out] subs Структура в которую записываем параметры вещества.

    Первая строка файла пропускается
    \code
        fscanf(fd, "%s\n", chrm);
    \endcode
    В первой строке не должно быть пробелов! Остальные строки последовательно считываются в структуру subs,
    пока не будет найден нужный элемент. Если элемент не найден на стандартный поток ввода-вывода (консоль)
    подаётся сигнал, что элемент не найден.
*/

void load_subst(int Z, const char* ch, subst_t &subs){
	char chrm[50];
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", ch);
        return;
    };
	fscanf(fd, "%s\n", chrm);

    subs.E = new double[1];
    subs.P = new double[1];
    do
	{
        fscanf(fd, "%d %s %lf %lf %lf %d\n", &subs.Z, subs.atom, &subs.M, &subs.rho, &subs.U0, &subs.N);

        delete[] subs.E;
        delete[] subs.P;

        subs.E = new double[subs.N];
        subs.P = new double[subs.N];
		fscanf(fd, " [");
		for (int i = 0; i<subs.N; i++)
		{
			fscanf(fd, "%lf", &subs.E[i]);
            subs.E[i] *= 1000;
		}
		fscanf(fd, "]");
		fscanf(fd, " [");
		for (int i = 0; i<subs.N; i++)
		{
			fscanf(fd, "%lf", &subs.P[i]);
		}
		fscanf(fd, "]\n");
		if ((feof(fd)!=0)&&(subs.Z!=Z))
		{
			printf("Element is absent at list\n");
			return;
		}
	} while (subs.Z!=Z);
	fclose(fd);
};

/*!
    \brief Находит в нужных точках \f$ E_i \f$ обратную транспортную длину для упругих
    столкновений, интерполируя табличные данные из файла вида:
    \code
    33 theta-points
    49 E-points
    .1 .5 1.0 2.0 3.0 ... 180.0
    E[eV],Differential_cross_section[cm**2/str],Total_cross_section[cm**2]
    .3000E+05 .5231E-14 .4612E-14 .3225E-14 ... |.2074E-16
    .2900E+05 .5210E-14 .4612E-14 .3255E-14 ... |.2130E-16
    ...
    \endcode

    В первой и второй строке должны присутствовать слова theta-points и E-points и
    после них не должно быть пробелов. В четвёртой строке также не должно быть пробелов.
    33 в данном случае означает количество углов \f$ \theta_i \f$ {.1, .5, 1.0, 2.0, 3.0, ..., 180.0}.
    49 означает размер массива энергий \f$ E_i \f$ {.3000E+05, .2900E+05, ...}. Полное сечение (|.2074E-16,
    |.2130E-16 и т. д.) не используется.
    \f[
        \lambda^{-1}_{tr}
        =
        \frac{2 \pi N_a \rho}{M} \int\limits_0^\pi \sigma(\theta) (1 - \cos \theta) \sin \theta \, d\theta
    \f]

    \param[out] ltr Искомый массив транспортной длины \f$ \lambda_{tr} (E_i) \f$ .
    \param[in] E Энергии, в которых мы ищем \f$ \lambda_{tr} \f$
    \param[in] N Размер массива энергий.
    \param[in] ch Имя файла (из директории запуска), из которого берём данные для обработки.
    \param[in] s Структура, содержащая параметры вещества.
*/
void load_ltr(double *ltr, double *E, int N, const char* ch, subst_t s)
{
	char chrm[500], filename[50];
	sprintf(filename, "%s%d_el.dat", ch, s.Z);
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", filename);
        return;
    };
	int theta_l, e_l;
	fscanf(fd, "%d %s\n", &theta_l, chrm);
	fscanf(fd, "%d %s\n", &e_l, chrm);
	double theta[theta_l], dltr[theta_l], E_points[e_l], ltr_points[e_l];
	for (int j = 0; j<theta_l; j++)
	{
		fscanf(fd, "%lE", &theta[j]);
		theta[j]*=M_PI/180;
	}
    fscanf(fd, "%s\n", chrm);
    for (int i = 0; i<e_l; i++)
    {
        fscanf(fd, "%lE", &E_points[e_l - 1 - i]);
        for (int j = 0; j<theta_l; j++)
        {
            fscanf(fd, "%lE", &dltr[j]);
            dltr[j]*=s.rho*Na/s.M*(1 - cos(theta[j]))*sin(theta[j]);
        }
        fscanf(fd, "%s\n", chrm);
        ltr_points[e_l - 1 - i] = 1./(2*M_PI*int_cubic_spline(theta[0], theta[theta_l-1], theta, dltr, theta_l));
    }
    eval_cubic_spline(E, ltr, N, E_points, ltr_points, e_l);
    fclose(fd);
};

/*!
    Находит массив значений потерь энергии на единице длины \f$ \bar{\varepsilon}(E_i) \f$, интерполируя их из файла.

    Файл имеет вид:
    \code
    E(eV)
    Q(eV)
    dlambda(Q)**(-1)/dQ(cm*eV)**(-1)
    .5000E+01
    .1000E-01 .1030E-01 ... .4630E+01
    .1546E+04 .1693E+04 ... .0000E+00
    .1000E+02
    .1000E-01 .1030E-01 ... .9630E+01
    .8384E+03 .9191E+03 ... .0000E+00
    ...
    \endcode
    Первый ряд представляет собой энергию. Второй ряд - потери энергии. Третий
    ряд плотность вероятности потери энергии на единице длины.
    \f[
        \bar{\varepsilon} = \int\limits_0^{E/2} Q \frac{d \lambda^{-1}(Q)}{d Q}\, dQ
    \f]    \param[out] eave Искомый массив \f$ \bar{\varepsilon}(E_i) \f$ .
    \param[in] E Энергии, в которых мы ищем \f$ \bar{\varepsilon} \f$
    \param[in] N Размер массива энергий.
    \param[in] ch Имя файла (из директории запуска), из которого берём данные для обработки.
    \param[in] s Структура, содержащая параметры вещества.
*/
void load_eave(double *eave, double *E, int N, const char* ch, subst_t s)
{
	char chr, chrm[500], filename[50];
	sprintf(filename, "%s%d_in.dat", ch, s.Z);
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", filename);
        return;
    };
	fscanf(fd, "%s\n", chrm);
	fscanf(fd, "%s\n", chrm);
	fscanf(fd, "%s\n", chrm);
	double Q_points[500], dW_points[500], E_points[100], e_sharp[100];
	int i = 0;
	do
	{
		fscanf(fd, "%lE", &E_points[i]);
		int j = 0;
		do
		{
			fscanf(fd, "%lE", &Q_points[j]);
			fscanf(fd, "%c", &chr);
			j++;
		} while (chr!='\n');
		j = 0;
		do
		{
			fscanf(fd, "%lE", &dW_points[j]);
			fscanf(fd, "%c", &chr);
			j++;
		} while (chr!='\n');
		e_sharp[i] = 0.;

		for (int l = 1; l<j-1; l++)
		{
			e_sharp[i] += 0.5*(Q_points[l] - Q_points[l-1])*(Q_points[l]*dW_points[l] + Q_points[l-1]*dW_points[l-1]);
		}
		i++;
	} while (feof(fd)==0);
	fclose(fd);
	eval_cubic_spline(E, eave, N, E_points, e_sharp, i-1);
};

/*!
    \brief Читаем данные из файла для метода аппроксимации.

    Файл имеет вид:
    \code
        #Z,N,E[],alphai[],d1i[],R0i[nm]
        32 3 [7.90 6.83 1.16] [1.366 1.290 1.2541] [5.79 5.96 17.3] [726.1 606.1 37] Ge
        14 3 [7.90 6.20 1.60] [1.652 1.6257 1.4236] [2.1972 2.1516 1.5023] [878.8 587.93 71.16] Si
    \endcode
    Первая строка не должна содержать пробелов.
    Вторая строка содержит последовательно данные по Оже электронам: номер элемента в периодической системе,
    число Оже-электронов, энергии Оже-электронов, коэффициенты для аппроксимации \f$ \bar{\varepsilon} \f$.
    \param[in] Z Номер элемента в периодической системе
    \param[in] ch Имя файла
    \param[out] app Структура, в которую считываем данные
*/
void load_approx(int Z, const char* ch, approx_t &app)
{
	char chrm[50];
	int z;
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", ch);
        return;
    };
	fscanf(fd, "%s\n", chrm);
	app.E = new double[1];
	app.alpha_i = new double[1];
    app.R0_i = new double[1];
    app.d1_i = new double[1];
    do
	{
        fscanf(fd, "%d ", &z);
		fscanf(fd, "%d", &app.N);
		delete [] app.E;
        delete [] app.alpha_i;
        delete [] app.R0_i;
        delete [] app.d1_i;
		app.E = new double[1];
        app.alpha_i = new double[1];
        app.R0_i = new double[1];
        app.d1_i = new double[1];
		fscanf(fd, " [");
		for (int i = 0; i<app.N; i++)
		{
			fscanf(fd, "%lf", &app.E[i]);
            app.E[i] *= 1000;
		}
		fscanf(fd, "]");
		fscanf(fd, " [");
		for (int i = 0; i<app.N; i++)
		{
			fscanf(fd, "%lf", &app.alpha_i[i]);
		}
		fscanf(fd, "]");
        fscanf(fd, " [");
		for (int i = 0; i<app.N; i++)
		{
			fscanf(fd, "%lf", &app.d1_i[i]);
		}
		fscanf(fd, "]");
        fscanf(fd, " [");
		for (int i = 0; i<app.N; i++)
		{
			fscanf(fd, "%lf", &app.R0_i[i]);
            app.R0_i[i] *= 1e-7;
		}
		fscanf(fd, "]");
		fscanf(fd, "%s\n", chrm);
		if ((feof(fd)!=0)&&(z!=Z))
		{
			printf("Element is absent at list\n");
			return;
		}
	} while (z!=Z);
	fclose(fd);
};

/*!
    \brief Аппроксимация \f$ \lambda_{tr} \f$ и \f$ \bar{\varepsilon} \f$ на
    промежутках \f$ (E_{i+1}, E_i]\f$.

    Аппроксимация на промежутке \f$ (E_{i+1}, E_i] \f$ проводится по законам:
    \f[
        \bar{\varepsilon}(E') = \frac{E_i}{\alpha_i R_{0i}} \left( \frac{E'}{E_i} \right)^{1 - \alpha_i}
    \f]
    \f[
        \lambda_{tr}(E') = \frac{R_{0i}}{d_{1i}} \left( \frac{E'}{E_i} \right)^{\alpha_i}
    \f]
    \param[out] ltr Массив \f$ \lambda_{tr}(E'_i) \f$
    \param[out] eps Массив \f$ \bar{\varepsilon}(E'_i) \f$
    \param[in] E Массив \f$ E'_i \f$
    \param[in] N Размер массивов \f$ E'_i \f$, \f$ \lambda_{tr}(E'_i) \f$, \f$ \bar{\varepsilon}(E'_i) \f$
    \param[in] ap Содержит параметры аппроксимации
    \param[in] s Содержит параметры вещества
*/
void le_approx(double *ltr, double *eps, double *E, int N, approx_t ap, subst_t s)
{
    int k;
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<ap.N; j++)
        {
            if ((E[i] >= ap.E[j+1])&&(E[i] < ap.E[j]))
            {
                k = j;
            }
        }
        if (E[i] >= ap.E[0])
        {
            k = 0;
        }
        if (E[i] < ap.E[ap.N - 1])
        {
            k = ap.N - 1;
        }
        eps[i] = ap.E[k] / ap.alpha_i[k] / ap.R0_i[k] * pow(E[i] / ap.E[k], 1 - ap.alpha_i[k]);
        ltr[i] = ap.R0_i[k] / ap.d1_i[k] * pow(E[i] / ap.E[k], ap.alpha_i[k]);
    }
}

/*!
    \brief Функция используется в методе Монте-Карло, ищет параметры \f$ \alpha(E_i) \f$ и \f$ \lambda^{-1}_{el} \f$.

    Так как массивы данных велики, хранить их в памяти, выполнять многомерную интерполяцию
    неудобно, поэтому использовалась аппроксимация вероятности упругого рассеяния
    в углы \f$ (0, \theta) \f$ выражением:
    \f[
        P = 1/\Sigma\int\limits_0^\theta \sigma(\theta) d\theta \approx 2/\pi \arctan (\alpha \theta)
    \f]
    В этом случае задача свелась к поиску \f$ \alpha(E_i) \f$. Построив график по опытным данным,
    можно выбрать точку, для которой погрешности минимальны в широком диапазоне углов. Методом научного
    тыка было установлено, что эта точка приходится на \f$ P \approx 0.8 \f$. Таким образом, в функции
    проводится поиск полного сечения \f$ \Sigma(E_i) \f$ и частичных сечений \f$ \Sigma(E_i, \theta_k) \f$, по полному сечению
    определяется
    \f[ \lambda^{-1}_{el} = \frac{2 \pi \rho N_a \Sigma}{M} \f]
    Затем для данного сечения проводится поиск \f$ \alpha(E_i) \f$.

    \param[out] alpha Массив значений \f$ \alpha(E_i) \f$
    \param[in] E Массив \f$ E_i \f$
    \param[in] N Размер массивов
    \param[out] inv_lambda_el \f$ \lambda^{-1}_{el} \f$
    \param[in] ch Имя файла
    \param[in] s Информация о веществе
*/
void load_mc_elastic(double *alpha, double *E, int N, double &inv_lambda_el, const char* ch, subst_t s)
{
    char chrm[500], filename[50];
	sprintf(filename, "%s%d_el.dat", ch, s.Z);
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", filename);
        return;
    };

    int theta_l, e_l;
	fscanf(fd, "%d %s\n", &theta_l, chrm);
	fscanf(fd, "%d %s\n", &e_l, chrm);

    double *theta, *E_points, *alpha_points, *sigma, dsigma, horror;
    theta = new double[theta_l];
    sigma = new double[theta_l];
    E_points = new double[e_l];
    alpha_points = new double[e_l];

    /*
    Метод конечно ошибочный, но зато не нужно массивы хранить и интерполяцией заниматься
    \alpha - массив коэффициентов для приближённого выражения
    1/\Sigma\int\limits_0^\theta \sigma(\theta) d\theta \approx 2/\pi \arctan (\alpha \theta)
    */

	for (int j = 0; j<theta_l; j++)
	{
		fscanf(fd, "%lE", &theta[j]);
		theta[j]*=M_PI/180;
	}
	fscanf(fd, "%s\n", chrm);
	for (int i = 0; i<e_l; i++)
	{
		fscanf(fd, "%lE", &E_points[e_l - 1 - i]);
        sigma[0] = 0;
		for (int j = 1; j<theta_l; j++)
		{
			fscanf(fd, "%lE", &dsigma);
            sigma[j] = sigma[j-1] + dsigma * sin(theta[j]) * (theta[j] - theta[j-1]);
		}
        inv_lambda_el = 2*M_PI*s.rho/s.M*Na*sigma[theta_l - 1];
		fscanf(fd, " |%lE\n", &horror);
        for (int j = 0; j<theta_l; j++)
		{
			sigma[j] /= sigma[theta_l - 1];
            if (sigma[j] > 0.8)
            {
                alpha_points[e_l - 1 - i] = 1 / theta[j] * tan(M_PI * sigma[j] / 2);
                break;
            }
		}
	}
    eval_cubic_spline(E, alpha, N, E_points, alpha_points, e_l);
	fclose(fd);
};

/*!
    \brief Функция используется в методе Монте-Карло, ищет параметры \f$ \beta(E_i) \f$ и \f$ \lambda^{-1}_{in}(E_i) \f$.

    Аппроксимация вероятности неупругого рассеяния
    в углы \f$ (0, \theta) \f$ выражением:
    \f[
        P(E_i, Q) = \frac{\int\limits_0^{Q} \frac{d\lambda^{-1}_{in}}{d Q} dQ}
        {\int\limits_0^{E/2} \frac{d\lambda^{-1}_{in}}{d Q} dQ}
        \approx 2/\pi \arctan \left(\beta \frac{Q}{E_i}\right)
    \f]
    Задача сводится к поиску \f$ \beta(Q_j) \f$. Построив график по опытным данным,
    можно выбрать точку, для которой погрешности минимальны в широком диапазоне
    углов (\f$ P \approx 0.8 \f$).

    \param[out] beta Массив значений \f$ \beta(Q_i) \f$
    \param[in] E Массив \f$ E_i \f$
    \param[in] N Размер массивов
    \param[out] inv_lambda_in Массив \f$ \lambda^{-1}_{in}(E_i) \f$
    \param[in] ch Имя файла
    \param[in] s Информация о веществе
*/
void load_mc_inelastic(double *beta, double *E, int N, double *inv_lambda_in, const char* ch, subst_t s)
{
	char chr, chrm[500], filename[50];
	sprintf(filename, "%s%d_in.dat", ch, s.Z);
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", filename);
        return;
    };
    /*
    Метод конечно ошибочный, но зато не нужно массивы хранить и интерполяцией заниматься
    \beta - массив коэффициентов для приближённого выражения
    \approx 2/\pi \arctan (\alpha Q/E)
    */
	fscanf(fd, "%s\n", chrm);
	fscanf(fd, "%s\n", chrm);
	fscanf(fd, "%s\n", chrm);
	double Q_points[500], dW_points[500], E_points[100], e_sharp[500], beta_points[100], inv_l_points[100];
	int i = 0, imax;
	do
	{
		fscanf(fd, "%lE", &E_points[i]);
		int j = 0;
		do
		{
			fscanf(fd, "%lE", &Q_points[j]);
			fscanf(fd, "%c", &chr);
			j++;
		} while (chr!='\n');
		j = 0;
		do
		{
			fscanf(fd, "%lE", &dW_points[j]);
			fscanf(fd, "%c", &chr);
			j++;
		} while (chr!='\n');
		int jmax = j - 1;
		e_sharp[0] = 0.;

		for (int l = 1; l<jmax; l++)
		{
			e_sharp[l] = e_sharp[l-1] + 0.5*(Q_points[l] - Q_points[l-1])*(dW_points[l] + dW_points[l-1]);
		}
        inv_l_points[i] = e_sharp[jmax - 1];
        for (int l = 1; l<jmax; l++)
		{
			e_sharp[l] /= e_sharp[jmax - 1];
            if (e_sharp[l] > 0.8)
            {
                beta_points[i] = E_points[i] / Q_points[l] * tan(M_PI * e_sharp[l] / 2);
                break;
            }
		}
		i++;
	} while (feof(fd)==0);
	fclose(fd);
    imax = i - 1;
	eval_cubic_spline(E, beta, N, E_points, beta_points, imax);
    eval_cubic_spline(E, inv_lambda_in, N, E_points, inv_l_points, imax);
};
