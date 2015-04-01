/*!
    \file
    Файл содержит вспомогательные математические функции.
*/

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "calculations.h"

/*!
    \brief
	spe решает уравнение:
	\f[
		\frac{\partial f}{\partial t} =
		p(t) \frac{\partial^2 f}{\partial z^2} + g(t),
	\f]
	при граничных условиях:
	\f[
		a_1(t) f + b_1(t) \frac{\partial f}{\partial z} \Bigg|_{z = z_1} = c_1(t),
	\f]
	\f[
		a_2(t) f + b_2(t) \frac{\partial f}{\partial z} \Bigg|_{z = z_2} = c_2(t),
	\f]
	и произвольных начальных условиях неявным методом на одном шаге \f$[t, t+dt]\f$.

	\details
	\param[out] f Содержит новый результат решения уравнения на одном шаге
	\param[in] f_prev Содержит предыдущий результат решения уравнения на одном шаге
	\param[in]  z Массив узлов, над которыми решается уравнение
	\param[in]  N Число узлов
	\param[in] dt Шаг по времени (энергии)
    \param[in] pt Параметр \f$p(t)\f$
    \param[in] gt Параметр \f$g(t)\f$
    \param[in] a_1t Параметр \f$a_1(t)\f$
    \param[in] b_1t Параметр \f$b_1(t)\f$
    \param[in] c_1t Параметр \f$c_1(t)\f$
    \param[in] a_2t Параметр \f$a_2(t)\f$
    \param[in] b_2t Параметр \f$b_2(t)\f$
    \param[in] c_2t Параметр \f$c_2(t)\f$

	Код:
	\code
        A[0] = 0;
        B[0] = a_1t - b_1t/(z[1] - z[0]);
        C[0] = b_1t/(z[1] - z[0]);
        D[0] = c_1t;
        for (int i = 1; i < N - 1; i++)
        {
            double dz1, dz2, dza;
            dz1 = z[i] - z[i-1];
            dz2 = z[i+1] - z[i];
            dza = (dz1 + dz2)*0.5;
            A[i] = dt*pt/dz1/dza;
            B[i] = - 1 - 2*pt/dz1/dz2*dt;
            C[i] = dt*pt/dz2/dza;
            D[i] = - gt*dt - f_prev[i];
        }
        A[N-1] = - b_2t/(z[N-1] - z[N-2]);
        B[N-1] = a_2t + b_2t/(z[N-1] - z[N-2]);
        C[N-1] = 0;
        D[N-1] = c_2t;
	\endcode
	задаёт трёхдиагональную матрицу и правую часть уравнения. A, B, C -- диагонали
	матрицы нижняя, собственно диагональ, верхняя. D -- правая часть уравнения.
	Код:
	\code
        for (int i = 1; i<N; i++)
        {
            B[i] -= A[i]/B[i-1]*C[i-1];
            D[i] -= A[i]/B[i-1]*D[i-1];
        }
        f[N-1] = D[N-1]/B[N-1];
        for (int i = N-2; i>=0; i--)
        {
            //if (fabs(f_prev[i]) > 0) printf("%e %e\n", f[i], f_prev[i]);
            f[i] = (D[i] - C[i] * f[i+1])/B[i];
        }
	\endcode
	Реализует метод прогонки решения систем линейных уравнений с трёхдиагональными матрицами.
*/

void spe(double *f,
         double *f_prev,
         double *z,
         int N,
         double dt,
         double pt,
         double gt,
         double a_1t,
         double b_1t,
         double c_1t,
         double a_2t,
         double b_2t,
         double c_2t)
{
	double A[N], B[N], C[N], D[N];
	A[0] = 0;
	B[0] = a_1t - b_1t/(z[1] - z[0]);
	C[0] = b_1t/(z[1] - z[0]);
	D[0] = c_1t;
	for (int i = 1; i < N - 1; i++)
    {
        double dz1, dz2, dza;
        dz1 = z[i] - z[i-1];
		dz2 = z[i+1] - z[i];
        dza = (dz1 + dz2)*0.5;
		A[i] = dt*pt/dz1/dza;
		B[i] = - 1 - 2*pt/dz1/dz2*dt;
		C[i] = dt*pt/dz2/dza;
		D[i] = - gt*dt - f_prev[i];
	}
	A[N-1] = - b_2t/(z[N-1] - z[N-2]);
	B[N-1] = a_2t + b_2t/(z[N-1] - z[N-2]);
	C[N-1] = 0;
	D[N-1] = c_2t;
	for (int i = 1; i<N; i++)
    {
		B[i] -= A[i]/B[i-1]*C[i-1];
		D[i] -= A[i]/B[i-1]*D[i-1];
	}
	f[N-1] = D[N-1]/B[N-1];
	for (int i = N-2; i>=0; i--)
    {
        //if (fabs(f_prev[i]) > 0) printf("%e %e\n", f[i], f_prev[i]);
		f[i] = (D[i] - C[i] * f[i+1])/B[i];
	}
    //if (fabs(f_prev[1]) > 0) getchar();
}

/*!
    \brief
	cubic_spline ищет коэффициенты сплайна:
	\f[
		S_i(x) = a_i + b_i (x - x_i) + \frac{1}{2} c_i (x - x_i)^2 + \frac{1}{6} d_i (x - x_i)^3
		\quad x\in[x_i, x_{i+1}].
	\f]

	\details
	\param[in] x Точки \f$ x_i \f$, в которых нам известны значения функции \f$ y_i \f$
	\param[in] y Массив \f$ y_i \f$
	\param[in] N Размер массивов
	\param[out] a Массив коэффициентов \f$ a_i \f$
	\param[out] b Массив коэффициентов \f$ b_i \f$
	\param[out] c Массив коэффициентов \f$ c_i \f$
	\param[out] d Массив коэффициентов \f$ d_i \f$

	Решение задачи сводится к решению системы линейных уравнений с трёхдиагональной матрицей относительно вектора \f$ c_i \f$. Код:
	\code
        double A[N-2], B[N-2], C[N-2], D[N-2];
        double dx[N-1], dydx[N-1];
        c[0] = 0;
        c[N-1] = 0;
        dx[0] = x[1] - x[0];
        dydx[0] = (y[1] - y[0])/dx[0];
        for (int i = 0; i<(N-2); i++)
        {
            dx[i+1] = x[i+2] - x[i+1];
            dydx[i+1] = (y[i+2] - y[i+1])/dx[i+1];
            A[i+1] = dx[i];
            B[i] = 2*(dx[i] + dx[i+1]);
            C[i] = dx[i];
            D[i] = 6*(dydx[i+1] - dydx[i]);
        }
        for (int i = 1; i<N-2; i++)
        {
            B[i] -= A[i]/B[i-1]*C[i];
            D[i] -= A[i]/B[i-1]*D[i-1];
        }
        c[N-2] = D[N-3]/B[N-3];
        for (int i = N-4; i>=0; i--)
        {
            c[i+1] = (D[i] - C[i] * c[i+2])/B[i];
        }
	\endcode
	Остальные значения коэффициентов могут быть найдены как линейные комбинации \f$ c_i \f$, что реализует следующий код:
	\code
        for (int i = 0; i<N-1; i++)
        {
            a[i] = y[i];
            b[i] = dydx[i] - c[i]*dx[i]/3 - c[i+1]*dx[i]/6;
            d[i] = (c[i+1] - c[i])/dx[i];
        }
	\endcode
	Обратите внимание 3 коэффициента остаются произвольными и их выбор связан с тем, как мы задаём приближение на концах.
	В данном случае они выбраны следующими:
	\code
        d[N-1] = 0;
        b[N-1] = 0;
        a[N-1] = y[N-1];
	\endcode
*/

void cubic_spline(double *x, double *y, int N, double *a, double *b, double *c, double *d)
{
	//double a[N], b[N], c[N], d[N];
	double A[N-2], B[N-2], C[N-2], D[N-2];
	double dx[N-1], dydx[N-1];
	c[0] = 0;
	c[N-1] = 0;
	dx[0] = x[1] - x[0];
	dydx[0] = (y[1] - y[0])/dx[0];
	for (int i = 0; i<(N-2); i++)
    {
		dx[i+1] = x[i+2] - x[i+1];
		dydx[i+1] = (y[i+2] - y[i+1])/dx[i+1];
		A[i+1] = dx[i];
		B[i] = 2*(dx[i] + dx[i+1]);
		C[i] = dx[i];
		D[i] = 6*(dydx[i+1] - dydx[i]);
	}
	for (int i = 1; i<N-2; i++)
    {
		B[i] -= A[i]/B[i-1]*C[i];
		D[i] -= A[i]/B[i-1]*D[i-1];
	}
	c[N-2] = D[N-3]/B[N-3];
	for (int i = N-4; i>=0; i--)
    {
		c[i+1] = (D[i] - C[i] * c[i+2])/B[i];
	}
	for (int i = 0; i<N-1; i++)
    {
		a[i] = y[i];
		b[i] = dydx[i] - c[i]*dx[i]/3 - c[i+1]*dx[i]/6;
		d[i] = (c[i+1] - c[i])/dx[i];
	}
	d[N-1] = 0;
	b[N-1] = 0;
	a[N-1] = y[N-1];
}

/*!
    \brief
    eval_cubic_spline ищет значения функции \f$y(x'_i)\f$, интерполированное
    кубическим сплайном, в промежуточных точках \f$ x'_i \f$ основного массива
    \f$ x_i \f$.

    \details
    \param[in] xs Точки, в которых требуется найти значения функции \f$y(x'_i)\f$
    \param[out] ys Массив, искомых значений
    \param[in] M Число элементов массивов xs, ys
    \param[in] x Заданные точки \f$ x_i \f$
    \param[in] y Заданные точки \f$ y_i \f$
    \param[in] N Число элементов массивов \f$ x_i \f$, \f$ y_i \f$

    Сначала упорядочиваем массив элементов xs и находим коэффициенты сплайна
    \code
        int *in = new int[M];
        for (int i=0; i<M; i++) in[i] = i;
        std::sort(in, in + M, [&xs](int& a, int& b){return (xs[a] < xs[b]);});
        cubic_spline(x, y, N, a, b, c, d);
    \endcode
    Код:
    \code
        for (int i = 0; i<M; i++)
        {
            if (xs[in[i]]<x[0])
            {
                h = xs[in[i]] - x[0];
                ys[in[i]] = a[0] + b[0]*h + c[0]*h*h/2 + d[0]*h*h*h/6;
            }
            if (xs[in[i]]>x[N-1])
            {
                h = xs[in[i]] - x[N-1];
                ys[in[i]] = a[N-1] + b[N-1]*h + c[N-1]*h*h/2 + d[N-1]*h*h*h/6;
            }
            if ((xs[in[i]]>=x[0])&&(xs[in[i]]<x[N-1]))
            {
                while (!((xs[in[i]]<x[j+1])&&(xs[in[i]]>=x[j]))&&(j<N))
                {
                    j++;
                }
                h = xs[in[i]] - x[j];
                ys[in[i]] = a[j] + b[j]*h + c[j]*h*h/2 + d[j]*h*h*h/6;
            }
            if (xs[in[i]] == x[N-1])
            {
                ys[in[i]] = y[N-1];
            }
        }
    \endcode
    находит значения функции \f$ y(x'_i) \f$.

*/

void eval_cubic_spline(double *xs, double *ys, int M, double *x, double *y, int N)
{
	double a[N], b[N], c[N], d[N];
    int *in = new int[M];
    for (int i=0; i<M; i++) in[i] = i;
    std::sort(in, in + M, [&xs](int& a, int& b){return (xs[a] < xs[b]);});

	cubic_spline(x, y, N, a, b, c, d);
	//Пусть x сортирован по возрастанию
	int j = 0;
	double h;
	for (int i = 0; i<M; i++)
    {
        if (xs[in[i]]<x[0])
        {
            h = xs[in[i]] - x[0];
			ys[in[i]] = a[0] + b[0]*h + c[0]*h*h/2 + d[0]*h*h*h/6;
        }
        if (xs[in[i]]>x[N-1])
        {
            h = xs[in[i]] - x[N-1];
			ys[in[i]] = a[N-1] + b[N-1]*h + c[N-1]*h*h/2 + d[N-1]*h*h*h/6;
        }
		if ((xs[in[i]]>=x[0])&&(xs[in[i]]<x[N-1]))
        {
			while (!((xs[in[i]]<x[j+1])&&(xs[in[i]]>=x[j]))&&(j<N))
            {
				j++;
			}
			h = xs[in[i]] - x[j];
			ys[in[i]] = a[j] + b[j]*h + c[j]*h*h/2 + d[j]*h*h*h/6;
		}
		if (xs[in[i]] == x[N-1])
        {
			ys[in[i]] = y[N-1];
		}
	}
}

/*!
    \brief
    int_cubic_spline ищет интеграл от \f$ l_a \f$ до \f$ l_b\f$ от функции \f$ f(x) \f$,
    заданной на промежутке таблицей значений \f$ (x_i, y_i) \f$, приближая её кубическими
    сплайнами.

    \details
    \param[in] la Нижний предел интегрирования
    \param[in] lb Верхний предел интегрирования
    \param[in] x Заданные точки \f$ x_i \f$
    \param[in] y Заданные точки \f$ y_i \f$
    \param[in] N Число элементов массивов \f$ x_i \f$, \f$ y_i \f$

    Проверяем, что промежуток интегрирования внутри промежутка, на котором задана функция:
    \code
        if ((la>=x[0])&&(la<=x[N-1])&&(lb>=x[0])&&(lb<=x[N-1]))
    \endcode
    Интегрируем:
    \code
        double a[N], b[N], c[N], d[N], h;
		cubic_spline(x, y, N, a, b, c, d);
		double result = 0.;
		for (int j = 1; j<N-1; j++)
        {
            if (la<x[0])
            {
                h = x[0] - la;
                result+=a[0]*h + b[0]*h*h/2 + c[0]*h*h*h/6 + d[0]*h*h*h*h/24;
            }
            if (lb>x[N-1])
            {
                h = lb - x[N-1];
                result+=a[N-1]*h + b[N-1]*h*h/2 + c[N-1]*h*h*h/6 + d[N-1]*h*h*h*h/24;
            }
			if ((la<=x[j])&&(lb>=x[j+1]))
            {
				h = x[j+1] - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((la>x[j])&&(la<x[j+1]))
            {
				h = x[j+1] - la;
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((lb>x[j])&&(lb<x[j+1]))
            {
				h = lb - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
		}
    \endcode

*/

double int_cubic_spline(double la, double lb, double *x, double *y, int N)
{
	//Пусть x сортирован по возрастанию, la<lb
	if ((la>=x[0])&&(la<=x[N-1])&&(lb>=x[0])&&(lb<=x[N-1]))
    {
		double a[N], b[N], c[N], d[N], h;
		cubic_spline(x, y, N, a, b, c, d);
		double result = 0.;
		for (int j = 1; j<N-1; j++)
        {
            if (la<x[0])
            {
                h = x[0] - la;
                result+=a[0]*h + b[0]*h*h/2 + c[0]*h*h*h/6 + d[0]*h*h*h*h/24;
            }
            if (lb>x[N-1])
            {
                h = lb - x[N-1];
                result+=a[N-1]*h + b[N-1]*h*h/2 + c[N-1]*h*h*h/6 + d[N-1]*h*h*h*h/24;
            }
			if ((la<=x[j])&&(lb>=x[j+1]))
            {
				h = x[j+1] - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((la>x[j])&&(la<x[j+1]))
            {
				h = x[j+1] - la;
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((lb>x[j])&&(lb<x[j+1]))
            {
				h = lb - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
		}
		return result;
	} else {
		return NAN;
	}
}

/*!
    \brief random возвращает случайное число с равномерным распределением на полуинтервале
    \f$ (0, l] \f$.

    \details
    \param[in] l Длина промежутка

    0 исключается.
*/

double random(double l)
{
    return l * (rand() + 1) / (RAND_MAX + 1.);
}

/*!
    \brief linterp возвращает значение функции в случае линейной интерполяции.

    \details
    \param[in] x точка в которой ищем значении функции
    \param[in] x1 первая точка, в которой известно значение функции
    \param[in] y1 значение функции в точке \f$ x_1 \f$
    \param[in] x2 вторая точка, в которой известно значение функции
    \param[in] y2 значение функции в точке \f$ x_2 \f$
*/

double linterp(double x, double x1, double y1, double x2, double y2){
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
}
