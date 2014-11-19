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
	char chr, chrm[50], filename[50];
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
	printf("%d %d", theta_l, e_l);
	double theta[theta_l], dltr[theta_l], E_points[e_l], ltr_points[e_l];
	for (int j = 0; j<theta_l; j++)
	{
		fscanf(fd, "%lf", &theta[j]);
		printf("%f ", theta[j]);
	}
	printf("\n");
	//fscanf(fd, "%s\n", chrm);
	fscanf(fd, "%164c\n", chr);
	for (int i = 0; i<e_l; i++)
	{
		fscanf(fd, "%lf", &E_points[i]);
		printf("%f ", E_points[i]);
		for (int j = 0; j<theta_l; j++)
		{
			fscanf(fd, "%lf", &dltr[j]);
			printf("%f ", dltr[j]);
		}
		fscanf(fd, "%s\n", chrm);
		printf("\n");
	}
	
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
	
}