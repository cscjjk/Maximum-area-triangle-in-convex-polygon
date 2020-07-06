#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
using namespace std;

FILE * fout;
const int N = 3000000;

double x[N + 2];
double y[N + 2];

int testMyConvex(int n, double x[], double y[]){
	x[0] = x[n]; y[0] = y[n];
	x[n + 1] = x[1]; y[n + 1] = y[1];
	for (int i = 1; i <= n; i++){
		long double x1, x2, y1, y2;
		x1 = x[i] - x[i - 1];
		y1 = y[i] - y[i - 1];
		x2 = x[i + 1] - x[i];
		y2 = y[i + 1] - y[i];
		long double s = x1 * y2 - x2 * y1;
		if (s > -1e-7) return 1;
		long double r = x1 * x1 + y1 * y1;
		if (r < 1e-7) return 2;
	}
	return 0;
}


void checkData(const char * filename){
	int n, r;
	printf(filename);
	printf(" is under checking\n");
	FILE * fin = fopen(filename, "r");
	fscanf(fin, "%d", &n);
	while (n > 0){
		for (int i = 1; i <= n; i++)
			fscanf(fin, "%lf %lf", &x[i], &y[i]);
		r = testMyConvex(n, x, y);
		if (r > 0){
			printf("Bad %d\n", r);
			while (true) r = n;
		}
		fscanf(fin, "%d", &n);
	}
	fclose(fin);

	printf(filename);
	printf(" checked\n");
}

void swap(double & a, double & b){
	double c = a; a = b; b = c;
}

const double pi = 3.14159265358979323846;

void genPolygonE(int n){
	double theta = 2 * pi / n;
	double offset = (rand() % 32768) * 2 * pi / 32768.0;

	fprintf(fout, "%d\n", n);
	for (int i = n; i >= 1; i--)
		fprintf(fout, "%.4lf %.4lf\n", cos(theta * i + offset) * 1e4, sin(theta * i + offset) * 1e4);
	fprintf(fout, "\n");
}

void genByElipse(const char * filename, int instances, int n){
	fout = fopen(filename, "w");
	for (int i = 0; i < instances; i++) 
		genPolygonE(n);
	fprintf(fout, "0");
	fclose(fout);

	checkData(filename);
}

inline int smallRand(int n){ // For n<=10^4, gen a number in 0..n-1 (with equal probability)
	int A = (32768 / n) * n;  
	int x;
	do {x = rand();} while (x >= A);  // x is a number in 0..A-1, whereas A is a multiply of n
	return x % n;
}

inline double getRandDoub(){ // gen a double in [0.0,9999.9999]
	return (smallRand(10000) * 1e4 + smallRand(10000)) / 1e4;
}

double A[10002];

void genList(int n, double x[]){
	bool collision = true;
	while (collision){
		collision = false;

		A[1] = 0; A[n] = 1e4;
		for (int i = 2; i <= n - 1; i++)
			A[i] = getRandDoub();
		sort(A + 2, A + n);

		int a = 1, b = 1;
		for (int i = 2; i <= n - 1; i++){
			if ((rand() % 2) == 0){
				x[i - 1] = A[i] - A[a]; a = i;			
			}
			else{
				x[i - 1] = A[b] - A[i]; b = i;
			}
			if (x[i - 1] < 1e-5 && x[i - 1] > -1e-5){ 
				collision = true;
			}
		}
		x[n - 1] = A[n] - A[a];
		x[n] = A[b] - A[n];       //Note: the sum of x[1]..x[n] should be equal to (A[1]-A[n])+(A[n]-A[1])=0.
	}
}

void mysort(double x[], double y[], int q[], double X[], double Y[], int l, int r){  // sort by angle from large to small.
	if (l >= r) return;
	int pq = q[(l + r) / 2];
	long double pX = X[(l + r) / 2];
	long double pY = Y[(l + r) / 2];
	int i = l, j = r;
	while (i <= j){
		while (q[i] > pq || q[i] == pq && Y[i] * pX - X[i] * pY >  5e-9) i++;
		while (q[j] < pq || q[j] == pq && Y[j] * pX - X[j] * pY < -5e-9) j--;
		if (i <= j){
			swap(x[i], x[j]); swap(y[i], y[j]);
			swap(X[i], X[j]); swap(Y[i], Y[j]);
			int temp = q[i]; q[i] = q[j]; q[j] = temp;
			i++; j--;
		}
	}
	mysort(x, y, q, X, Y, l, j); mysort(x, y, q, X, Y, i, r);
}

int     q[10002];
double  X[10002];
double  Y[10002];

void genPolygonM(int n){
	bool nonconvex = true;

	while (nonconvex){
		nonconvex = false;
		genList(n, x);
		genList(n, y);
		for (int i = n; i >= 2; i--) swap(y[i], y[smallRand(i) + 1]);

		for (int i = 1; i <= n; i++)
			if (y[i] > 0)
				if (x[i] > 0){
					q[i] = 0; X[i] = x[i]; Y[i] = y[i];
				}
				else{
					q[i] = 1; X[i] = y[i]; Y[i] = -x[i];
				}
			else
				if (x[i] < 0){
					q[i] = 2; X[i] = -x[i]; Y[i] = -y[i];
				}
				else{
					q[i] = 3; X[i] = -y[i]; Y[i] = x[i];
				}

		mysort(x, y, q, X, Y, 1, n);
		x[0] = y[0] = 0;
		for (int i = 1; i <= n; i++){
			x[i] += x[i - 1]; y[i] += y[i - 1];
		}
		if (testMyConvex(n, x, y) > 0) 
			nonconvex = true;
	}

	fprintf(fout, "%d\n", n);
	for (int i = 1; i <= n; i++)
		fprintf(fout, "%.4lf %.4lf\n", x[i], y[i]);
	fprintf(fout, "\n");
}

void genByMatch(const char * filename, int instances, int n){  //n <= 10000
	fout = fopen(filename, "w");
	for (int i = 0; i < instances; i++){
		genPolygonM(n);
		if (i % 1000 == 0) printf("%d\n", i);
		}
	fprintf(fout, "0");
	fclose(fout);

	checkData(filename);
}

void genByMixed(const char * filename, int instances, int n){  //n <= 10000
	fout = fopen(filename, "w");
	for (int i = 0; i < instances; i++) 
		if (rand() % 2 == 0)
			genPolygonM(n);
		else
			genPolygonE(n);
	fprintf(fout, "0");
	fclose(fout);

	checkData(filename);
}

bool small_ang(const double & x0, const double & y0, const double & x1, const double & y1, const double & x2, const double &y2){
	double p = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
	return (p < -1e-9);
}

void XYsort(double x[], double y[], int l, int r){
	if (l >= r) return;
	double X = x[(l + r) / 2], Y = y[(l + r) / 2];
	int i = l, j = r;
	while (i <= j){
		while (small_ang(x[1], y[1], x[i], y[i], X, Y)) i++;
		while (small_ang(x[1], y[1], X, Y, x[j], y[j])) j--;
		if (i <= j){
			swap(x[i], x[j]);
			swap(y[i], y[j]);
			i++; j--;
		}
	}
	XYsort(x, y, l, j); XYsort(x, y, i, r);
}

void computeCH(const int N, int & n, double x[], double y[]){
	for (int i = 2; i <= N; i++)
		if (y[i] < y[1] || y[i] == y[1] && x[i] < x[1]){
			swap(x[1], x[i]);
			swap(y[1], y[i]);
		}
	XYsort(x, y, 2, N);

	n = 2;
	for (int i = 3; i <= N; i++){
		n ++;
		x[n] = x[i];
		y[n] = y[i];
		if (!small_ang(x[1], y[1], x[n - 1], y[n - 1], x[n], y[n])){
			double r1 = (x[n - 1] - x[1]) * (x[n - 1] - x[1]) + (y[n - 1] - y[1]) * (y[n - 1] - y[1]);
			double r2 = (x[n] - x[1]) * (x[n] - x[1]) + (y[n] - y[1]) * (y[n] - y[1]);
			if (r1 < r2){
				x[n - 1] = x[n]; y[n - 1] = y[n];
			}
			n--;
		}
		while (n >= 3){
			double p = (x[n - 1] - x[n - 2]) * (y[n] - y[n - 2]) - (x[n] - x[n - 2]) * (y[n - 1] - y[n - 2]);
			if (p < -1e-9) break;
			x[n - 1] = x[n]; y[n - 1] = y[n]; n--;
		}
	}
}

void genRandPoint(double & x, double & y, const int type){
	while (true) {
		x = ((rand() % 10000) * 1e4 + (rand() % 10000)) / 1e4; 
		y = ((rand() % 10000) * 1e4 + (rand() % 10000)) / 1e4;
		if (x * x + y * y < 1e8 || type == 0) return;
	}
}

int genPolygonH(int N, int type){	
	int n;
	
	for (int i = 1; i <= N; i++) 
		genRandPoint(x[i], y[i], type);
	computeCH(N, n, x, y);

	fprintf(fout, "%d\n", n);
	for (int i = 1; i <= n; i++)
		fprintf(fout, "%.4lf %.4lf\n", x[i], y[i]);
	fprintf(fout, "\n");

	return n;
}

void genByHull(const char * filename, int type){
	//Set type = 0 means that N poinss are drawn from a squre.  
	//Set type = 1 means that N points are drawn from a disk.  We recommand it.

	fout = fopen(filename, "w");
	int S = 0, n;
	while (S < 100000){
		n = genPolygonH(100000, type);
		printf("%d %d\n", n, S);
		S += n;
	}	
	fprintf(fout, "0");
	fclose(fout);

	checkData(filename);
}

void main(){
	genByMixed("7.in", 100, 10000);
	genByHull("1.in", 1);
	genByMatch("2.in", 5000, 200);  
	genByElipse("3.in", 5000, 200);
	genByMatch("4.in", 3333, 300); 
	genByElipse("5.in", 3333, 300); 
	genByMixed("6.in", 1000, 1000); 
}

