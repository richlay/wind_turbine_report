//
// "Given n, r ; find a, b, eps, eps_a, eps_b"
//
//求a,b
//

#include <iostream>
#include <stdio.h>
#include <math.h>

#define PI atan(1)*4

int main()
{
	int B = 3;  //葉片數量
	double CD, CL, R, D, alpha, U1;
	CD = 0.04089;
	CL = 0.9683;
	R = 0.08;	//旋轉半徑
	D = 0.16;	//葉片旋轉的直徑
	alpha = 11/180*PI;
	U1 = 2;  //初速度
	
	double J, n, U2, Utheta, V, thi, r, omega, c, theta, sigma;
	
	
	std::cout << "Enter n and r: \n";
	std::cin >> n >> r;  //e.g. n=8 , r=0.04
	
	
	double a1, b1, a2, b2, a, b, eps, eps_a, eps_b;
	a1 = -10;
	b1 = -10;
	eps = 2;
	while (a1<=10) {
		while (b1<=10) {
			
			J = U1 / n / D;  //advance ratio
			omega = 2*PI*n;  //角速度
			//c = (2*PI*r/B)*(8.0/9.0)*(R*R)/(r*r)*(U1*U1)/(R*R*omega*omega) ;  //葉片寬度
			c=0.02;
			U2 = U1 * (1 - a1);
			Utheta = r * omega * (1 + b1);
			theta = atan(U2/Utheta)-alpha;  //pitch angle
			sigma = B * c / 2.0 / PI / r;
//			V = sqrt(U2 * U2 + Utheta * Utheta); //似乎用不到
			thi = atan(U2/Utheta);
			
			
			a2 = sigma * (1 - a1) / 4.0 / sin(thi) / sin(thi) * (CL * cos(thi) + CD * sin(thi));
			b2 = sigma * (1 - a1) / 4.0 / sin(thi) / sin(thi) * J/ (2 * PI * r /D )* (CL * sin(thi) - CD * cos(thi));
			
			if (fabs(a1 - a2)+fabs(b1 - b2) < eps){
				eps = fabs(a1 - a2)+fabs(b1 - b2);
				eps_a = fabs(a1 - a2);
				eps_b = fabs(b1 - b2);
				a = 0.5*(a1+a2);
				b = 0.5*(b1+b2);
			}
			b1 += 0.01;
		}
		b1 = -10;
		a1 += 0.01;
	}
	
	printf("\n a= %G , b= %G , eps= %G , eps_a = %G , eps_b = %G \n", a, b, eps, eps_a, eps_b);
	
	return 0;
}
