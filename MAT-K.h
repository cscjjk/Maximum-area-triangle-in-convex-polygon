inline double abs(const double &x){
	if (x < 0) return -x; else return x;
}

bool coincide(const double & x1, const double & y1, const double & x2, const double & y2){
	return abs(x1 - x2) < 1e-6 && abs(y1 - y2) < 1e-6;
}

inline double prod(const double & x1, const double & y1, const double & x2, const double & y2){
	return x1 * y2 - x2 * y1;
}


bool aheadorparallel(const double x[], const double y[], int j, int k){  //determine the relation between e[j] and e[k]
	double delta = (x[j + 1] - x[j]) * (y[k + 1] - y[k]) - (x[k + 1] - x[k]) * (y[j + 1] - y[j]);
	return (delta < 1e-7);
}

int chooseAngle(const double & a1, const double & b1, const double & a2, const double & b2, const double & a3, const double & b3, double & a, double & b){
	a = a1; b = b1; int r = 1;
	if (prod(a2, b2, a, b) < 0){
		a = a2; b = b2; r = 2;
	}
	if (prod(a3, b3, a, b) < 0){
		a = a3; b = b3; r = 3;
	}
	return r;
}

inline void getLine(const double & x1, const double & y1, const double & x2, const double & y2, double & a, double & b, double & c){
	a = y1 - y2;
	b = x2 - x1;
	c = x1 * y2 - x2 * y1;
}

inline void getIntersection(const double & a1, const double & b1, const double & c1, const double & a2, const double & b2, const double & c2, double &xx ,double &yy){ 
	double delta = a1 * b2 - a2 * b1;
	xx = (c2 * b1 - c1 * b2) / delta;
	yy = (c1 * a2 - c2 * a1) / delta;
}


void computeIntersection(const double x[], const double y[], const int j, const int k, double & xx, double & yy){
	double a1, b1, c1, a2, b2, c2;
	getLine(x[j], y[j], x[j + 1], y[j + 1], a1, b1, c1);
	getLine(x[k], y[k], x[k + 1], y[k + 1], a2, b2, c2);
	getIntersection(a1, b1, c1, a2, b2, c2, xx, yy);
}

void computeIntersectionII(const double x[], const double y[], const int j, const double y0, double & xx, double & yy){
	double a1 = y[j] - y[j + 1];
	double b1 = x[j + 1] - x[j];
	double c1 = x[j] * y[j + 1] - x[j + 1] * y[j];
	yy = y0;
	xx = -(b1 * yy + c1)/ a1;
}

double b_x, b_y, c_x, c_y;

void initial_abc(const int n, const double x[], const double y[], int & a, int & b, int & c){
	a = 1;
	for (int i = 2; i <= n; i++)
		if (y[i] > y[a] + 1e-6 || (y[i] > y[a] - 1e-6 && x[i] > x[a])) a = i;
	int j = a - 1, k = a + 1;
	if (j < 1) j = n;
	if (k > n) k = 1;
	if (y[j] > y[a] - 1e-6) j--;
	if (j < 1) j = n;

	double xx, yy = 0;
	while (true){
		if (y[j] > y[k]){
			if (aheadorparallel(x, y, j, k-1)) j--; 
			else{
				computeIntersection(x, y, j, k-1, xx, yy);
				if (y[a] + yy > 2 * y[j] - 1e-7) break;
				j--;
			}	
		}
		else{
			if (aheadorparallel(x, y, j, k-1)) k++;
			else{
				computeIntersection(x, y, j, k-1, xx, yy);
				if (y[a] + yy > 2 * y[k] - 1e-7) break;
				k++;
			}
		}
		if (j < 1) j = n;
		if (k > n) k = 1;
	}
	yy = (y[a] + yy)/ 2;
	if (yy > y[j + 1]) yy = y[j + 1];
	if (yy > y[k - 1]) yy = y[k - 1];
	if (yy >= y[j + 1] - 1e-6){
		c = j + 1; if (c > n) c = 1; c_x = x[c]; c_y = y[c];
	}
	else{
		c = j; computeIntersectionII(x, y, c, yy, c_x, c_y);
	}
	if (yy <= y[k] + 1e-6){
		b = k; b_x = x[b]; b_y = y[b];
	}
	else{
		b= k - 1; computeIntersectionII(x, y, b, yy, b_x, b_y);
	}
}

double rotateProcess(const int n, const double x[], const double y[], int & a, int & b, int & c, int & T_iteration){
	int ansa = 1, ansb = 2, ansc = 3;
	double maxarea2 = 0;
	T_iteration = 0;
	double xx, yy;
	double aa, bb, cc;
	int onecycle = 0;
	while (onecycle < 2){
		if (abs(prod(x[a + 1] - x[a], y[a + 1] - y[a], b_x - c_x, b_y - c_y)) < 1e-7){
			a = a % n + 1; T_iteration++;
		}
		
		double area2 = abs(prod(b_x - x[a], b_y - y[a], c_x - x[a], c_y - y[a]));
		if (area2 > maxarea2){
			maxarea2 = area2;
			ansa = a; ansb = b; ansc = c;
		}
		int r;
		if (aheadorparallel(x, y, c, b)){
			r = chooseAngle(x[a + 1] - x[a], y[a + 1] - y[a], x[b + 1] - c_x, y[b + 1] - c_y, x[a + 1] - x[a], y[a + 1] - y[a], aa, bb); // fixing c
			xx = c_x; yy = c_y;
		}
		else{
			double Dx, Dy, Ex, Ey; 
			computeIntersection(x, y, b, c, Dx, Dy);
			getLine(b_x, b_y, c_x, c_y, aa, bb, cc);
			double distA = abs(aa * x[a] + bb * y[a] + cc);
			double distI = abs(aa * Dx +   bb * Dy   + cc);
			
			if (distA < distI - 1e-5){
				Ex = 2 * c_x - Dx; Ey = 2 * c_y - Dy;
				r = chooseAngle(x[a + 1] - x[a], y[a + 1] - y[a], x[b + 1] - c_x, y[b + 1] - c_y, x[a] - Ex, y[a] - Ey, aa, bb);
				xx = c_x; yy = c_y;
			}
			else if (distA > distI + 1e-5){
				Ex = 2 * b_x - Dx; Ey = 2 * b_y - Dy;
				r = chooseAngle(x[a + 1] - x[a], y[a + 1] - y[a], b_x - x[c + 1], b_y - y[c + 1], Ex - x[a], Ey - y[a], aa, bb);
				xx = b_x, yy = b_y;
			}
			else{
				double Mx = (x[a] + Dx) / 2.0; double My = (y[a] + Dy) / 2.0;
				r = chooseAngle(x[a + 1] - x[a], y[a + 1] - y[a], x[b + 1] - Mx, y[b + 1] - My, Mx - x[c + 1], My - y[c + 1], aa, bb);
				xx = Mx; yy = My;
			}
		}

		double a1, b1, c1, a2, b2, c2;
		getLine(x[b], y[b], x[b + 1], y[b + 1], a1, b1, c1);
		getIntersection(a1, b1, c1, bb, -aa, aa * yy - bb * xx, b_x, b_y);
		if (coincide(b_x, b_y, x[b + 1], y[b + 1])){
			b = b % n + 1; b_x = x[b]; b_y = y[b];
		}
		getLine(x[c], y[c], x[c + 1], y[c + 1], a2, b2, c2);
		getIntersection(a2, b2, c2, bb, -aa, aa * yy - bb * xx, c_x, c_y);
		if (coincide(c_x, c_y, x[c + 1], y[c + 1])){
			c = c % n + 1; c_x = x[c]; c_y = y[c];
		}
		T_iteration++;
		if (r == 1){
			a = a % n + 1; T_iteration++;
		}

		//printf("%.5lf %.5lf %.5lf %.5lf \n", b_x / 10, b_y / 10, c_x / 10, c_y / 10);
		if (prod(b_x - c_x, b_y - c_y, 1, 0) < -1e-6) onecycle = 1;
		if (onecycle == 1 && prod(b_x - c_x, b_y - c_y, 1, 0) > 1e-6) onecycle = 2;
	}
	T_iteration--;
	a = ansa; b = ansb; c = ansc;
	return maxarea2;
}

double computeMAT_Kal(const int n, const double x[], const double y[], int & a, int & b, int & c, int & T_iteration){
	//printf("debug\n");
	initial_abc(n, x, y, a, b, c);
	return rotateProcess(n, x, y, a, b, c, T_iteration);
}

	 

	

