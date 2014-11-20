#include "parse.h"

void load_auger(int Z, const char* ch, auger_t &aug)
{
	char chr, chrm[50];
	int z;
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL)
	{
        printf("Не могу открыть файл %s. Проверьте, что файл существует\n", ch);
        return;
    };
	fscanf(fd, "%s\n", chrm);
    while (z!=Z)
	{
        fscanf(fd, "%d ", &z);
		fscanf(fd, "%d", &aug.N);
		aug.E = new double[aug.N];
		aug.P = new double[aug.N];
		fscanf(fd, " [");
		for (int i = 0; i<aug.N; i++)
		{
			fscanf(fd, "%lf", &aug.E[i]);
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
			printf("Элемент отсутствует в списке\n", ch);
			return;
		}
	}
	fclose(fd);
}

void load_subst(int Z, const char* ch, subst_t &subs){
	char chr, chrm[50];
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL)
	{
        printf("Не могу открыть файл %s. Проверьте, что файл существует\n", ch);
        return;
    };
	fscanf(fd, "%s\n", chrm);
    while (subs.Z!=Z)
	{
        fscanf(fd, "%d %lf %lf %lf %s\n", &subs.Z, &subs.M, &subs.rho, &subs.U0);
		if ((feof(fd)!=0)&&(subs.Z!=Z)) 
		{
			printf("Элемент отсутствует в списке\n", ch);
			return;
		}
	}
	fclose(fd);	
};

void load_ltr(double *ltr, double *E, int N, const char* ch, int Z)
{
	char chr, chrm[500], filename[50];
	sprintf(filename, "%s%d_el.pl", ch, Z);
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL)
	{
        printf("Не могу открыть файл %s. Проверьте, что файл существует\n", filename);
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
			dltr[j]*=(1 - cos(theta[j]))*sin(theta[j]);
		}
		fscanf(fd, "%s\n", chrm);
		ltr_points[e_l - 1 - i] = 2*M_PI*int_cubic_spline(theta[0], theta[theta_l-1], theta, dltr, theta_l);
	}
	FILE *ftest;
	ftest = fopen("test.gp", "w");
	fprintf(ftest, "plot '-' with lines\n");
	for (int i = 0; i<e_l; i++)
	{
		fprintf(ftest, "%e %e\n", E_points[i], ltr_points[i]);
	}
	fclose(ftest);
	ftest = popen("test.gp", "w");
	pclose(ftest);
	for (int i = 0; i<N; i++)
	{
			E[i] = E_points[0] + (E_points[e_l - 1] - E_points[0])/(N-1)*i;
	}
	eval_cubic_spline(E, ltr, N, E_points, ltr_points, e_l);
}

void test_parse(){
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
	double ltr[N], E[N];
	load_ltr(ltr, E, N, "data/", 32);
	for (int i = 0; i<N; i++)
	{
		printf("%e %e\n", E[i], ltr[i]);
	}
	
}