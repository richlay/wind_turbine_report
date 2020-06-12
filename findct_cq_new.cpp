//
// "Given n, find CT, CQ(without a,b)"
//
// findCTandCQ.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
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
	R = 0.08;  //葉片旋轉的半徑
	D = 0.16;  //旋轉直徑
	alpha = 11/180*PI;
	U1 = 2;  //初速度
	
	double J, n, U2, Utheta, V, thi, r, omega, c, theta, sigma;
	double a1=0;
	double b1=0;
	
	double CT, CQ, sum_CT=0, sum_CQ=0;
	
	std::cout << "Enter n: ";
	std::cin >> n;
	
	r=0.01;
	for(int i=1; i<=4; i++){
		J = U1 / n / R;  //advance ratio
		omega = 2*PI*n;  //角速度
		c = (2*PI*r/B)*(8.0/9.0)*(R*R)/(r*r)*(U1*U1)/(R*R*omega*omega) ;  //葉片寬度
		//c=0.02;
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
		r=r+0.02;
	}
	printf("\n J= %G , CT= %G , CQ= %G ", J, 0.02*sum_CT, 0.02*sum_CQ);
	return 0;
}

