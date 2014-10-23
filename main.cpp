#include <stdio.h>
#include <math.h>

double p(double z, double t){
	return 1;
};
double g(double z, double t){
	return 0;
};
double gamma_1(double z, double t){
	return 0;
};
double gamma_2(double z, double t){
	return 0;
};

void spe(double *f, double *f_prev, double *z, int N, double t, double dt, double alpha_1, double alpha_2, double beta_1, double beta_2){
	double A[N], B[N], C[N], D[N], dz, dq;
	A[0] = 0;
	B[0] = alpha_1 - alpha_2/(z[1] - z[0]);
	C[0] = alpha_2/(z[1] - z[0]);
	D[0] = gamma_1(z[0], t);
	for (int i = 1; i<N-1; i++){
		dz = z[i+1] - z[i];
		dq = dt/dz/dz;
		A[i] = p(z[i], t)*dq;
		B[i] = -1 - 2*p(z[i], t)*dq;
		C[i] = p(z[i], t)*dq;
		D[i] = - g(z[i], t)*dt - f_prev[i];
	}
	A[N-1] = beta_2/(z[N-1] - z[N-2]);
	B[N-1] = beta_1 - beta_2/(z[N-1] - z[N-2]);
	C[N-1] = 0;
	D[N-1] = gamma_2(z[N-1], t);
	for (int i = 1; i<N; i++){
		B[i] -= A[i]/B[i-1]*C[i];
		D[i] -= A[i]/B[i-1]*D[i-1];
	}
	f[N-1] = D[N-1]/B[N-1];
	for (int i = N-2; i>=0; i--){
		f[i] = (D[i] - C[i] * f[i+1])/B[i]; 
	}
}


int main(){
	int N = 1000, M = 100;
	double t = 0., dt = 1e-2, *phi, zmax = 3.1416;
	double dz = zmax/(N-1), z[N];
	phi = new double[N*M];
	for (int i = 0; i<N; i++){
		z[i] = (i==0) ? 0 : z[i-1]+dz;
		phi[i] = sin(z[i]);
	}
	for (int i = 1; i<M; i++){
		t+=dt;
		spe(phi+i*N, phi+(i-1)*N, z, N, t, dt, 1, 0, 1, 0);
	}
	FILE *file;
	file = fopen("results.txt", "w");
	for (int i = 0; i<N; i++){
		fprintf(file, "%f %f %f %f\n", z[i], phi[i], *(phi+N*50+i), *(phi+N*(M-1)+i));
	}
	fclose(file);
	return 0;
}