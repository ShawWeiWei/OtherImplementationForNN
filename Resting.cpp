// Resting.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "work.h"
#include <ctime>
#include <cstdlib>
int _tmain(int argc, _TCHAR* argv[])
{
	long t_begin,t_end;
	t_begin=clock();	

	double aI[]={0.0};
	double aGc[]={0.25};
	int aML1[]={50};//{1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99};
	Heter(aGc,sizeof(aGc)/sizeof(double),aI,sizeof(aI)/sizeof(double),aML1,sizeof(aML1)/sizeof(int),true,false,false);

	t_end=clock();
	int all_secs=(t_end-t_begin)/(1000);
	int hours=all_secs/3600;
	int mins=(all_secs%3600)/60;
	int secs=(all_secs%3600)%60;

	printf("Process exited after %d hours %d mins %d secs\n",hours,mins,secs);    
//	instance.OutputCoupleIdx();
//	instance.OutputI();
	system("PAUSE");
	return 0;
}

