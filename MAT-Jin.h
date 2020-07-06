inline double dist(const double xb, const double yb, const double xc, const double yc, const double x0, const double y0){	   // (xb,yb), (xc,yc), (x0,y0)) lie in clockwise order
	return (yc - yb) * x0 + (xb - xc) * y0;
}

inline double triangleArea2(const double x[], const double y[], const int a, const int b, const int c){  // Pa, Pb, Pc in clockwise order
	return (x[c] - x[a]) * (y[b] - y[a]) - (x[b] - x[a]) * (y[c] - y[a]);
}

void findOne3stable(const int n, const double x[], const double y[], int & a, int & b, int & c){
	double max_root_area2 = 0;
	
	int C = 3;
	for (int B = 2; B < n; B++){
		if (C == B) C++;
		while ((C < n) && (dist(x[1], y[1], x[B], y[B], x[C + 1], y[C + 1]) > dist(x[1], y[1], x[B], y[B], x[C], y[C]))) C++;
		double nowarea2 = triangleArea2(x, y, 1, B, C);
		if (nowarea2 > max_root_area2){
			max_root_area2 = nowarea2; a = 1; b = B; c = C;
		}
	}

	while (dist(x[b], y[b], x[c], y[c], x[a + 1], y[a + 1]) > dist(x[b], y[b], x[c], y[c], x[a], y[a]) + 1e-6){
		a = a % n + 1; 
		int check = 1;
		while (check) {
			check = 0;
			while (dist(x[c], y[c], x[a], y[a], x[b + 1], y[b + 1]) > dist(x[c], y[c], x[a], y[a], x[b], y[b]) + 1e-6){
				b = b % n + 1; check = 1;
			}
			while (dist(x[a], y[a], x[b], y[b], x[c + 1], y[c + 1]) > dist(x[a], y[a], x[b], y[b], x[c], y[c]) + 1e-6){
				c = c % n + 1; check = 1;
			}
		}
	}
	while (dist(x[b], y[b], x[c], y[c], x[a - 1], y[a - 1]) > dist(x[b], y[b], x[c], y[c], x[a], y[a]) + 1e-6){
		a--; if (a == 0) a = n; 
		int check = 1;
		while (check) {
			check = 0;
			while (dist(x[c], y[c], x[a], y[a], x[b - 1], y[b - 1]) > dist(x[c], y[c], x[a], y[a], x[b], y[b]) + 1e-6){
				b--; if (b == 0) b = n; check = 1;
			}
			while (dist(x[a], y[a], x[b], y[b], x[c - 1], y[c - 1]) > dist(x[a], y[a], x[b], y[b], x[c], y[c]) + 1e-6){
				c--; if (c == 0) c = n; check = 1;
			}
		}
	}
}

double RotateAndKill(const int n, const double x[], const double y[], int & a, int & b, int & c, int & T_iteration){  
// input (a,b,c) must be 3-stable.  return (a,b,c) is the maximum area triangle.
	int ansa = 1, ansb = 2, ansc = 3;
	double maxarea2 = 0;
	T_iteration = 0;
	int r = a, t = c;
	while (b != t || c != r){		
		while (dist(x[b], y[b], x[c], y[c], x[a + 1], y[a + 1]) > dist(x[b], y[b], x[c], y[c], x[a], y[a])){
			a = a % n + 1;
			T_iteration++;
		}
		double nowarea2 = triangleArea2(x, y, a, b, c);
		if (nowarea2 > maxarea2){
			maxarea2 = nowarea2;
			ansa = a; ansb = b; ansc = c;
		}
		double A1 = y[b+1] - y[b], B1 = x[b] - x[b+1];
		double A2 = y[c+1] - y[c], B2 = x[c] - x[c+1];
		double delta = A1 * B2 - A2 * B1;
		if (delta >= -1e-9) 
			b = b % n + 1;	// in this case point I is undefined.
		else{
			double U = A1 * x[c] + B1 * y[c], V = A2 * x[b] + B2 * y[b];
			double Ix = (U * B2 - V * B1) / delta;
			double Iy = (V * A1 - U * A2) / delta;
			if (dist(x[b], y[b], x[c], y[c], x[a], y[a]) > dist(x[b], y[b], x[c], y[c], Ix, Iy)) c = c % n + 1;
			else b = b % n + 1;
		}
		T_iteration++;
	}
	a = ansa; b = ansb; c = ansc; return maxarea2;
}

double computeMAT_Jin(const int n, const double x[], const double y[], int & a, int & b, int & c, int & T_iteration){
	findOne3stable(n, x, y, a, b, c);
	return RotateAndKill(n, x, y, a, b, c, T_iteration);
}
