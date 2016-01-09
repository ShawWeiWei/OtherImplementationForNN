#include "stdafx.h"
#include "ML.h"
#include "SquareNetwork.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
class SquareNetwork;
const double C=20.0;
const double gK=8.0;
const double gL=2.0;
const double VK=-84.0;
const double VCa=120.0;
const double VL=-60.0;
const double V1=-1.2;
const double V2=18.0;

void MorrisLecar::SetI(double _I){
	I=_I;
}
MorrisLecar::MorrisLecar():No(0){
}

MorrisLecar::MorrisLecar(const int _No):No(_No){
}

MorrisLecar::~MorrisLecar(){}

void MorrisLecar::equation(){
		DV=(I-gL*(V-VL)-gCa*0.5*(1+tanh((V-V1)/V2))*(V-VCa)-gK*n*(V-VK))/C;
	    Dn=(0.5*(1+tanh((V-V3)/V4))-n)*phi*cosh((V-V3)/(2*V4));
}

/*void MorrisLecar::equation(){
	DV=(I-gL*(V-VL)-gCa*0.5*(1+tanh((V-V1)/V2))*(V-VCa)-gK*n*(V-VK)+Network::CoupleVector[No])/C+Network::NoiseVector[No];
	Dn=(0.5*(1+tanh((V-V3)/V4))-n)*phi*cosh((V-V3)/(2*V4));
}*/

void MorrisLecar::euler_s(){
	equation();
	V+=DV*single_dt;
	n+=Dn*single_dt;
}

void MorrisLecar::euler(){
	equation();
	V+=(DV+SquareNetwork::sm_gCouple[No]/C)*SquareNetwork::sm_dt;
	n+=Dn*SquareNetwork::sm_dt;
}


void MorrisLecar::rungekutta_s(){
		V0=V;
		n0=n;
		equation();
		DV1=DV*single_dt;
		Dn1=Dn*single_dt;
		V=V0+DV1/2.0;
		n=n0+Dn1/2.0;
		equation();
		DV2=DV*single_dt;
		Dn2=Dn*single_dt;
		V=V0+DV2/2.0;
		n=n0+Dn2/2.0;
		equation();
		DV3=DV*single_dt;
		Dn3=Dn*single_dt;
		V=V0+DV3;
		n=n0+Dn3;
		equation();
		DV4=DV*single_dt;
		Dn4=Dn*single_dt;
		V=V0+(DV1+2.0*DV2+2.0*DV3+DV4)/6.0;
		n=n0+(Dn1+2.0*Dn2+2.0*Dn3+Dn4)/6.0;
}

void MorrisLecar::rungekutta(){
	 V0=V;
	n0=n;
	double CouplePlusNoise=SquareNetwork::sm_gCouple[No]/C;//+SquareNetwork::sm_gNoise[No];
	equation();
	DV1=(DV+CouplePlusNoise)*SquareNetwork::sm_dt;
	Dn1=Dn*SquareNetwork::sm_dt;
	V=V0+DV1/2.0;
	n=n0+Dn1/2.0;
	equation();
	DV2=(DV+CouplePlusNoise)*SquareNetwork::sm_dt;
	Dn2=Dn*SquareNetwork::sm_dt;
	V=V0+DV2/2.0;
	n=n0+Dn2/2.0;
	equation();
	DV3=(DV+CouplePlusNoise)*SquareNetwork::sm_dt;
	Dn3=Dn*SquareNetwork::sm_dt;
	V=V0+DV3;
	n=n0+Dn3;
	equation();
	DV4=(DV+CouplePlusNoise)*SquareNetwork::sm_dt;
	Dn4=Dn*SquareNetwork::sm_dt;
	V=V0+(DV1+2.0*DV2+2.0*DV3+DV4)/6.0;
	n=n0+(Dn1+2.0*Dn2+2.0*Dn3+Dn4)/6.0;
}


//void MorrisLecar::InitVar(){
//	V=4.583345;//-40.0;
//	n=0.29891089;//0.0;
//}

void MorrisLecar::RandVar(){
	V=-10+30.0*rand()/(double)RAND_MAX;
	n=0.1+0.5*rand()/(double)RAND_MAX;
}

void MorrisLecar::set_class1(){
	gCa=4.0;
	V3=12.0;
	V4=17.40;
	phi=1.0/15.0;
	strcpy_s(neuron_model,"MorrisLecar1");
}

void MorrisLecar::set_class2(){
	 gCa=4.4;
	 V3=2;
	 V4=30;
	 phi=0.04;
	 strcpy_s(neuron_model,"MorrisLecar2");
}


