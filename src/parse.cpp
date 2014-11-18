#include "parse.h"

void load(int Z, const char* ch, auger_t &aug)
{
	char chr;
	int z;
	FILE *fd;
	if ((fd = fopen(ch, "r")) == NULL){
        printf("Не могу открыть файл %s. Проверьте, что файл существует\n", ch);
        return;
    };
    while (z!=Z)
	{
        fscanf(fd, "%d ", &z);
		fscanf(fd, "%d", &aug.N);
		aug.E = new double[aug.N];
		aug.P = new double[aug.N];
		for (int i = 0; i<aug.N; i++)
		{
			fscanf(fd, "%lf", &aug.E[i]);
		}
		for (int i = 0; i<aug.N; i++)
		{
			fscanf(fd, "%lf", &aug.P[i]);
		}
		fscanf(fd, "%s%s\n", aug.atom, aug.shell);
		if (feof(fd)!=0) {
			printf("Элемент отсутствует в списке\n", ch);
			return;
		}
	}
	fclose(fd);
}

void test_parse(){
	auger_t aug;
	int Z = 32;
	load(Z, "data/aug.pl", aug);
	printf("%d \n[", aug.N);
	for (int i = 0; i<aug.N; i++){
		printf("%f ", aug.E[i]);
	}
	printf("]\n[");
	for (int i = 0; i<aug.N; i++){
		printf("%f ", aug.P[i]);
	}
	printf("]\n");
	printf("%s \n", aug.atom);
	printf("%s \n", aug.shell);
}