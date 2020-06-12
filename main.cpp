//
// "Given n, find J, CT, CQ (with a,b)"
//
//  main.cpp
//  find_t_q_final
//
//#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

#define PI atan(1)*4


void find_a_b(double f_r, double f_n, double *f_a, double *f_b);

int main()
{
	int B = 3;  //葉片數量
	double CD, CL, R, D, alpha, U1;
	CD = 0.04089;
	CL = 0.9683;
	R = 0.08;  //葉片旋轉的半徑
	D = 0.16;  //旋轉直徑
	alpha = 11/180*PI;
	U1 = 2;  //初速度
	
	double J, n, U2, Utheta, V, thi, r, omega, c, theta, sigma;
	
	double a1;
	double *p_a1 = &a1;
	double b1;
	double *p_b1 = &b1;
	
	double CT, CQ, sum_CT=0, sum_CQ=0;
	
	std::cout << "Enter n: ";
	std::cin >> n;
	
	r=0.004;					//r=R/(k*2)
	for(int i=1; i<=10; i++){					//i<='k' , k is accuracy , above:r=R/(k*2) , below(3個地方):r=r+R/k
		J = U1 / n / D;  //advance ratio
		omega = 2*PI*n;  //角速度
		//c = (2*PI*r/B)*(8.0/9.0)*(R*R)/(r*r)*(U1*U1)/(R*R*omega*omega) ;  //葉片寬度
		c=0.02;
		
		find_a_b(r, n, p_a1, p_b1);
		
		U2 = U1 * (1 - a1);
		Utheta = r * omega * (1 + b1);
		theta = atan(U2/Utheta)-alpha;  //pitch angle
		sigma = B * c / 2.0 / PI / r;
		V = sqrt(U2 * U2 + Utheta * Utheta);
		thi = atan(U2/Utheta);
		
		
		CT = (c/2)*pow((J*(1-a1)/(D*sin(thi))),2)*(CL*cos(thi)+CD*sin(thi));
		CQ = (c/2)*pow((J*(1-a1)/(D*sin(thi))),2)*(CL*sin(thi)-CD*cos(thi))*(r/D);
		sum_CT = sum_CT + CT;
		sum_CQ = sum_CQ + CQ;
		r=r+0.008;			//r=r+'R/k'
	}
	printf("\n J= %G , CT= %G , CQ= %G ", J, 0.008*sum_CT, 0.008*sum_CQ);					//r=r+'R/k'
	return 0;
}




// 找a,b
void find_a_b(double f_r, double f_n, double *f_a, double *f_b)
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
	
	
//	std::cout << "Enter n and r: \n";
//	std::cin >> n >> r;  //e.g. n=8 , r=0.04
	r = f_r;
	
	
	double a1, b1, a2, b2, a, b, eps, eps_a, eps_b;
	a1 = -10;
	b1 = -10;
	eps = 2;
	while (a1<=10) {
		while (b1<=10) {
			
			J = U1 / f_n / D;  //advance ratio
			omega = 2*PI*f_n;  //角速度
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
	
//	printf("\n a= %G , b= %G , eps= %G , eps_a = %G , eps_b = %G \n", a, b, eps, eps_a, eps_b);
	
	*f_a = a;
	*f_b = b;
}
