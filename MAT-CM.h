#include <math.h> // need sqrt

inline double ABS(const double &x){
	if (x < 0) return -x; else return x;
}

inline bool identical(const double & x1, const double & y1, const double & x2, const double & y2){
	return (ABS(x1 - x2) < 1e-4 && ABS(y1 - y2) < 1e-4);
}

inline double dist_p2p(const double & x1, const double & y1, const double & x2, const double & y2){
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

inline double dist_p2l(const double l[3], const double xx, const double yy){
	return (l[0] * xx + l[1] * yy + l[2]) / sqrt(l[0] * l[0] + l[1] * l[1]);
}

inline double crossProduct(const double x0, const double y0, const double x1, const double y1, const double x2, const double y2){
	return (x2 - x0) * (y1 - y0) - (x1 - x0) * (y2 - y0);
}

void computeLine(const double & x1, const double & y1, const double & x2, const double & y2, double l[3]){
	l[0] = y2 - y1; l[1] = x1 - x2; l[2] = x2 * y1 - x1 * y2;
}

void computeLineP(const double & x1, const double & y1, const double & x2, const double & y2, const double & x0, const double & y0, double l[3]){  
	// return a l parallel to (p1--p2) and pass p0.
	computeLine(x1, y1, x2, y2, l); l[2] = - l[0] * x0 - l[1] * y0;
}

void computeLineReflect(const double & x1, const double & y1, const double & x2, const double & y2, const double & x0, const double & y0, double l[3]){  // return l.
	// return a l that is the reflection of edge (p1-p2) around p0
	computeLine(x1, y1, x2, y2, l); l[2] = - 2 * l[0] * x0 - 2 * l[1] * y0 - l[2];
}

bool computeIntersection(const double l1[3], const double l2[3], double &xx, double &yy){
	double delta = l1[0] * l2[1] - l2[0] * l1[1];
	if (ABS(delta) < 1e-7){
		xx = yy = 1e40; return false;
	}
	else{
		xx = (l1[1] * l2[2] - l1[2] * l2[1]) / delta;
		yy = (l1[2] * l2[0] - l1[0] * l2[2]) / delta;
		return true;
	}
}

bool segmentIntersect(const double x1, const double y1, const double x2, const double y2, const double x3, const double y3, const double x4, const double y4){
	return ((crossProduct(x1, y1, x2, y2, x3, y3) * crossProduct(x1, y1, x2, y2, x4, y4) <= 1e-7)
		 && (crossProduct(x3, y3, x4, y4, x1, y1) * crossProduct(x3, y3, x4, y4, x2, y2) <= 1e-7));
}

int ea, eb, ec;
bool flusha, flushb;
double ax,ay, bx,by, cx,cy;
double Ax,Ay, Bx,By, Cx,Cy;
double Xx,Xy, Yx,Yy, Zx,Zy;
double l0[3], l1[3], l2[3], l3[3], g1[3], g2[3];

inline void update_a(const int n, const double x[], const double y[], bool proceed){
	if (proceed){
		ax = x[ea + 1]; ay = y[ea + 1]; ea = ea % n + 1; flusha = false;
	}
	else{
		ax = (Cx + Bx) / 2.0; ay = (Cy + By) / 2.0; 
	}
	computeLine(x[ea], y[ea], x[ea + 1], y[ea + 1], l1); 
}

inline void update_b(const int n, const double x[], const double y[], bool proceed){
	if (proceed){
		bx = x[eb + 1]; by = y[eb + 1]; eb = eb % n + 1; flushb = false;
	}
	else{
		bx = (Cx + Ax) / 2.0; by = (Cy + Ay) / 2.0;
	}
	computeLine(x[eb], y[eb], x[eb + 1], y[eb + 1], l2);
}

bool findRightmostCommon(const double & Lx, const double & Ly, const double & x1, const double y1, const double & x2, const double & y2, const double & Rx, const double & Ry, double & x, double & y){
	// L on the left of p1.  p2 on the left of R.    L on the left of R.
	double len = dist_p2p(Lx, Ly, Rx, Ry);
	if (dist_p2p(Lx, Ly, x1, y1) + dist_p2p(x2, y2, Rx, Ry) < len - 1e-4) return false; // no common point
	else if (dist_p2p(Lx, Ly, x1, y1) < len){
		x = x1; y = y1; return true;
	} else{
		x = Rx; y = Ry; return true;
	}
}

void toBegin(int n, const double x[], const double y[]){  //O(n)  compute a,b,c ea,eb,ec, A, B, C.
	computeLine(x[1], y[1], x[2], y[2], l0);
	cx = x[2]; cy = y[2]; ec = 2;
	computeLine(cx, cy, x[3], y[3], l3); 

	int p = 3, q = n;
	while (true){
		computeLine(x[p - 1], y[p - 1], x[p], y[p], l1); 	
		computeIntersection(l1, l0, Bx, By);
		computeLine(x[q], y[q], x[q + 1], y[q + 1], l2); 	
		computeIntersection(l2, l0, Ax, Ay);

		double JLx, JLy, JHx, JHy, KLx, KLy, KHx, KHy;
		JLx = 2 * x[p - 1] - Bx; JLy = 2 * y[p - 1] - By;
		JHx = 2 * x[p] - Bx; JHy = 2 * y[p] - By;
		KLx = 2 * x[q + 1] - Ax; KLy = 2 * y[q + 1] - Ay;
		KHx = 2 * x[q] - Ax; KHy = 2 * y[q] - Ay;

		if (segmentIntersect(JLx, JLy, JHx, JHy, KLx, KLy, KHx, KHy)){
			computeIntersection(l1, l2, Cx, Cy); 
			ea = p - 1; update_a(n, x, y, false); 
			eb = q; update_b(n, x, y, false);
			break;
		}

		double distp = dist_p2l(l0, x[p], y[p]);
		double distq = dist_p2l(l0, x[q], y[q]);
		
		if (distp < distq - 1e-5){
			computeLine(x[p], y[p], x[p + 1], y[p + 1], l1); computeIntersection(l1, l0, Bx, By);
			JLx = 2 * x[p] - Bx; JLy = 2 * y[p] - By;
			if (segmentIntersect(JLx, JLy, JHx, JHy, KLx, KLy, KHx, KHy)){
				computeLineReflect(x[1], y[1], x[2], y[2], x[p], y[p], l0);	
				computeIntersection(l0, l2, Cx, Cy); 
				Bx = 2 * x[p] - Cx; By = 2 * y[p] - Cy; 
				ea = p - 1; update_a(n, x, y, true); 
				eb = q; update_b(n, x, y, false);
				break;
			}
			else p++;
		}
		else if (distq < distp - 1e-5){
			computeLine(x[q - 1], y[q - 1], x[q], y[q], l2); computeIntersection(l2, l0, Ax, Ay);
			KLx = 2 * x[q] - Ax; KLy = 2 * y[q] - Ay;
			if (segmentIntersect(JLx, JLy, JHx, JHy, KLx, KLy, KHx, KHy)){
				computeLineReflect(x[1], y[1], x[2], y[2], x[q], y[q], l0);
				computeIntersection(l0, l1, Cx, Cy); 
				Ax = 2 * x[q] - Cx; Ay = 2 * y[q] - Cy;
				ea = p - 1; update_a(n, x, y, false); 
				eb = q - 1; update_b(n, x, y, true);
				break;
			}
			else q--;
		}
		else// distp == distq.
			if (p + 1 != q){
				computeLine(x[p], y[p], x[p + 1], y[p + 1], l1);		
				computeIntersection(l1, l0, Bx, By);
				JLx = 2 * x[p] - Bx; JLy = 2 * y[p] - By;
				
				computeLine(x[q - 1], y[q - 1], x[q], y[q], l2);		
				computeIntersection(l2, l0, Ax, Ay);
				KLx = 2 * x[q] - Ax; KLy = 2 * y[q] - Ay;
				
				if (!findRightmostCommon(JHx, JHy, JLx, JLy, KLx, KLy, KHx, KHy, Cx, Cy)) {p++; q--;}
				else{
					Ax = 2 * x[q] - Cx; Ay = 2 * y[q] - Cy; 
					Bx = 2 * x[p] - Cx; By = 2 * y[p] - Cy;
					ea = p - 1; update_a(n, x, y, true); 
					eb = q - 1; update_b(n, x, y, true);
					break;
				}
			}
			else{ // special case.
				Cx = KHx; Cy = KHy;
				Ax = 2 * x[q] - Cx; Ay = 2 * y[q] - Cy; 
				Bx = 2 * x[p] - Cx; By = 2 * y[p] - Cy;
				ea = p - 1; update_a(n, x, y, true); 
				eb = q - 1; update_b(n, x, y, true);
				break;
			}
	}
}

inline int closestPoint(const double & x, const double & y){
	double distX = (Xx - x) * (Xx - x) + (Xy - y) * (Xy - y);
	double distY = (Yx - x) * (Yx - x) + (Yy - y) * (Yy - y);
	double distZ = (Zx - x) * (Zx - x) + (Zy - y) * (Zy - y);

	if (distZ < distX + 1e-6 && distZ < distY + 1e-6) return 3;
	else if (distX < distY) return 1;
	else return 2;
}

int twoFlushLegs(const int n, const double x[], const double y[]){
	double Mx = (Cx + cx) / 2.0; 
	double My = (Cy + cy) / 2.0;
	Xx = 2 * x[ea + 1] - Cx; 
	Xy = 2 * y[ea + 1] - Cy;
	computeLineP(Mx, My, x[ea + 1], y[ea + 1], cx, cy, g1);
	computeLineP(Mx, My, x[eb + 1], y[eb + 1], cx, cy, g2);
	computeIntersection(g2, l1, Yx, Yy);
	computeIntersection(l3, l1, Zx, Zy);
	
	int d = closestPoint(Bx, By);
	switch (d){
		case 1: update_a(n, x, y, true); Bx = Xx; By = Xy; 
			computeIntersection(g1, l2, Ax, Ay); update_b(n, x, y, false); break;
		case 2: update_b(n, x, y, true); Ax = 2 * bx - Cx; Ay = 2 * by - Cy;
			Bx = Yx; By = Yy; update_a(n, x, y, false); break;
		case 3: Bx = Zx; By = Zy; computeIntersection(l3, l2, Ax, Ay); 
			update_a(n, x, y, false); update_b(n, x, y, false); break;
	}
	return d;
}

int oneFlushLeg_a(const int n, const double x[], const double y[]){
	bx = x[eb]; by = y[eb]; // For precision reason
	computeLineP(x[ea], y[ea], x[ea + 1], y[ea + 1], Ax, Ay, l0); 
	computeLineP(x[ea + 1], y[ea + 1], bx, by, cx, cy, g1);
	computeIntersection(l0, g1, Xx, Xy);
	if (computeIntersection(l0, l2, Yx, Yy)){
		if (crossProduct(Yx, Yy, Bx, By, Ax, Ay) < 0){
			Yx = Yy = 1e40;
		}
	}
	computeIntersection(l0, l3, Zx, Zy);

	int d = closestPoint(Ax, Ay);
	switch (d){
		case 1: Ax = Xx; Ay = Xy; computeIntersection(l1, g1, Bx, By);
			update_a(n, x, y, true); Cx = 2 * ax - Bx; Cy = 2 * ay - By; break;
		case 2: Ax = Yx; Ay = Yy; computeLine(Yx, Yy, cx, cy, g2); 
			computeIntersection(l1, l2, Cx, Cy);
			computeIntersection(l1, g2, Bx, By);
			update_a(n, x, y, false); flushb = true; break;
		case 3: Ax = Zx; Ay = Zy; computeIntersection(l1, l3, Bx, By);
			Cx = bx * 2 - Ax; Cy = by * 2 - Ay; update_a(n, x, y, false); break;
	}
	return d;
}

int oneFlushLeg_b(const int n, const double x[], const double y[]){
	ax = x[ea]; ay = y[ea]; // For precision reason
	if (identical(cx, cy, Bx, By)){
		update_b(n, x, y, true); Ax = 2 * bx - Cx; Ay = 2 * by - Cy; return 2;
	}
	computeLineP(x[eb], y[eb], x[eb + 1], y[eb + 1], Bx, By, l0);
	computeLineP(ax, ay, x[eb + 1], y[eb + 1], cx, cy, g2);
	computeIntersection(l0, l1, Xx, Xy);
	computeIntersection(l0, g2, Yx, Yy);
	if (computeIntersection(l0, l3, Zx, Zy))
		if (crossProduct(Zx, Zy, Ax, Ay, Bx, By) < 0){
			Zx = Zy = 1e40;
		}
	int d	= closestPoint(Bx, By);
	switch (d){
		case 1: Bx = Xx; By = Xy; computeLine(Xx, Xy, cx, cy, g1); 
			computeIntersection(g1, l2, Ax, Ay);
			computeIntersection(l1, l2, Cx, Cy);
			update_b(n, x, y, false); flusha = true; break;
		case 2: Bx = Yx; By = Yy; computeIntersection(g2, l2, Ax, Ay);
			update_b(n, x, y, true); Cx = bx * 2 - Ax; Cy = by * 2 - Ay; break;
		case 3: Bx = Zx; By = Zy; computeIntersection(l3, l2, Ax, Ay);
			Cx = ax * 2 - Bx; Cy = ay * 2 - By;  update_b(n, x, y, false); break;
	}
	return d;
}

void computeFlush_plush_MoveCanonical(const int n, const double x[], const double y[]){
	computeLine(Ax, Ay, Bx, By, l0);

	double JLx, JLy, JHx, JHy, KLx, KLy, KHx, KHy;
	computeIntersection(l1, l0, JLx, JLy);
	JHx = ax * 2 - JLx; JHy = ay * 2 - JLy;
	computeIntersection(l2, l0, KLx, KLy);
	KHx = bx * 2 - KLx; KHy = by * 2 - KLy;
	
	double dist1 = dist_p2p(Cx, Cy, JHx, JHy);
	double dist2 = dist_p2p(Cx, Cy, KHx, KHy);
	if (identical(x[ea + 1], y[ea + 1], bx, by) || (dist2 < dist1 - 1e-4)){
		Cx = KHx; Cy = KHy; 
		Ax = KLx; Ay = KLy; Bx = ax * 2 - Cx; By = ay * 2 - Cy; 
		flushb = true;
	}
	else if (dist1 < dist2 - 1e-4){
		Cx = JHx; Cy = JHy; 
		Bx = JLx; By = JLy; Ax = bx * 2 - Cx; Ay = by * 2 - Cy;
		flusha = true;
	}
	else{
		Cx = (JHx + KHx) / 2; Cy = (JHy + KHy) / 2;
		Ax = KLx; Ay = KLy; Bx = JLx; By = JLy; 
		flusha = flushb = true;
	}
}

double rotate(const int n, const double x[], const double y[], int & T_iteration){
	T_iteration = 0; 
	double maxarea2 = 0;
	double ax0 = ax, ay0 = ay, bx0 = bx, by0 = by;
	flusha = flushb = false;
	do{
		double nowarea2 = (cx - ax) * (by - ay) - (bx - ax) * (cy - ay);
		if (nowarea2 > maxarea2){
			maxarea2 = nowarea2;
		}
		if (!flusha && !flushb) computeFlush_plush_MoveCanonical(n, x, y);
		int r = 0;
		if (flusha && flushb) r = twoFlushLegs(n, x, y);
		else if (flusha) r = oneFlushLeg_a(n, x, y);
		else if (flushb) r = oneFlushLeg_b(n, x, y);
		T_iteration++;
		if (r == 3){
			T_iteration++; 
			ec = ec % n + 1; cx = x[ec]; cy = y[ec]; 
			computeLine(cx, cy, x[ec + 1], y[ec + 1], l3); 
		}		
	}while (dist_p2p(ax, ay, ax0, ay0) + dist_p2p(bx, by, bx0, by0) > 1e-3 || T_iteration < 3);

	return maxarea2;
}

double computeMAT_CM(const int n, const double x[], const double y[], int & T_iteration){
	toBegin(n, x, y);
	return rotate(n, x, y, T_iteration);
}
