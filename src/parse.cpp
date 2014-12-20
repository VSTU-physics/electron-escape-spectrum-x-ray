#include "parse.h"

void load_auger(int Z, const char* ch, auger_t &aug)
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
	aug.E = new double[1];
	aug.P = new double[1];
    do
	{
        fscanf(fd, "%d ", &z);
		fscanf(fd, "%d", &aug.N);
		delete aug.E;
		delete aug.P;
		aug.E = new double[aug.N];
		aug.P = new double[aug.N];
		fscanf(fd, " [");
		for (int i = 0; i<aug.N; i++)
		{
			fscanf(fd, "%lf", &aug.E[i]);
            aug.E[i] *= 1000;
		}
		fscanf(fd, "]");
		fscanf(fd, " [");
		for (int i = 0; i<aug.N; i++)
		{
			fscanf(fd, "%lf", &aug.P[i]);
		}
		fscanf(fd, "]");
		fscanf(fd, "%s%s\n", aug.atom, aug.shell);
		if ((feof(fd)!=0)&&(z!=Z))
		{
			printf("Element is absent at list\n");
			return;
		}
	} while (z!=Z);
	fclose(fd);
}

void load_subst(int Z, const char* ch, subst_t &subs){
	char chrm[50];
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL)
	{
        printf("Can't open file %s. Check that file exists\n", ch);
        return;
    };
	fscanf(fd, "%s\n", chrm);
    do
	{
        fscanf(fd, "%d %lf %lf %lf %s\n", &subs.Z, &subs.M, &subs.rho, &subs.U0, chrm);
		if ((feof(fd)!=0)&&(subs.Z!=Z))
		{
			printf("Element is absent at list\n");
			return;
		}
	} while (subs.Z!=Z);
	fclose(fd);
};

void load_ltr(double *ltr, double *E, int N, const char* ch, subst_t s)
{
	char chrm[500], filename[50];
	sprintf(filename, "%s%d_el.pl", ch, s.Z);
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
    FILE *fd2;
    fd2 = fopen("testltr.dat", "w");
    for (int i = 0; i<e_l; i++)
    {
        fscanf(fd, "%lE", &E_points[e_l - 1 - i]);
        printf("%e %e\n", E_points[e_l - 1 - i], eta(s, E_points[e_l - 1 - i]));
        for (int j = 0; j<theta_l; j++)
        {
            fscanf(fd, "%lE", &dltr[j]);
            if ((i == 0)||(i == e_l/2)||(i == e_l - 1)) fprintf(fd2, "%e %e %e\n", theta[j], dltr[j], crsec(theta[j], s, E_points[e_l - 1 - i]));
            dltr[j]*=s.rho*Na/s.M*(1 - cos(theta[j]))*sin(theta[j]);
        }
        fprintf(fd2, "\n");
        fscanf(fd, "%s\n", chrm);
        ltr_points[e_l - 1 - i] = 1./(2*M_PI*int_cubic_spline(theta[0], theta[theta_l-1], theta, dltr, theta_l));
    }
    eval_cubic_spline(E, ltr, N, E_points, ltr_points, e_l);
    fclose(fd);
    fclose(fd2);
    fd2 = fopen("testltr.gp", "w");
    fprintf(fd2, "set terminal wxt 4\n");
    fprintf(fd2, "set size square\n");
    fprintf(fd2, "plot 'testltr.dat' using 1:2 with lines title 'cross_section_T(E)',\\\n");
    fprintf(fd2, "'testltr.dat' using 1:3 with lines title 'cross_section_A(E)' \n");
    fclose(fd2);
};

void load_esharp(double *esharp, double *E, int N, const char* ch, subst_t s)
{
	char chr, chrm[500], filename[50];
	sprintf(filename, "%s%d_in.pl", ch, s.Z);
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
	int i = 0, j, jmax, imax;
	do
	{
		fscanf(fd, "%lE", &E_points[i]);
		j = 0;
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
		jmax = j - 1;
		e_sharp[i] = 0.;

		for (int l = 1; l<jmax; l++)
		{
			e_sharp[i] += 0.5*(Q_points[l] - Q_points[l-1])*(Q_points[l]*dW_points[l] + Q_points[l-1]*dW_points[l-1]);
		}
		i++;
	} while (feof(fd)==0);
	fclose(fd);
	imax = i - 1;
	eval_cubic_spline(E, esharp, N, E_points, e_sharp, imax);
};

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

void test_parse()
{
	auger_t aug;
	int Z = 14;
	load_auger(Z, "data/aug.pl", aug);
	printf("%d \n[", aug.N);
	for (int i = 0; i<aug.N; i++)
	{
		printf("%f ", aug.E[i]);
	}
	printf("]\n[");
	for (int i = 0; i<aug.N; i++)
	{
		printf("%f ", aug.P[i]);
	}
	printf("]\n");
	printf("%s \n", aug.atom);
	printf("%s \n", aug.shell);
	int N = 10;
	double ltr[N], E[N], esharp[N];
	subst_t s;
	load_subst(32, "data/subst.pl", s);
	load_ltr(ltr, E, N, "data/", s);
	load_esharp(esharp, E, N, "data/", s);
	for (int i = 0; i<N; i++)
	{
		printf("%e %e %e\n", E[i], ltr[i], esharp[i]);
	}

};

void load_mc_elastic(double *alpha, double *E, int N, double &inv_lambda_el, const char* ch, subst_t s)
{
    char chrm[500], filename[50];
	sprintf(filename, "%s%d_el.pl", ch, s.Z);
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
    int k;
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

void load_mc_inelastic(double *beta, double *E, int N, double *inv_lambda_in, const char* ch, subst_t s)
{
	char chr, chrm[500], filename[50];
	sprintf(filename, "%s%d_in.pl", ch, s.Z);
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
	int i = 0, j, jmax, imax;
	do
	{
		fscanf(fd, "%lE", &E_points[i]);
		j = 0;
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
		jmax = j - 1;
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
