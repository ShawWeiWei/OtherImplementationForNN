#include "stdafx.h"
#include "work.h"
#include "SquareNetwork.h"
#include <stdio.h>

void Config(int *aML1,int sizeML1){
	for(int iML1=0;iML1<sizeML1;iML1++){
		SquareNetwork instance(16384,aML1[iML1]);
		instance.OutputNoForOneAndTwo();
	}
}

void Heter(double *aGc,int sizeGc,double *aI,int sizeI,int *aML1,int sizeML1,bool isSpiral,bool isSpiking,bool isCoupling){

	for(int iML1=0;iML1<sizeML1;iML1++){
		SquareNetwork instance(16384,aML1[iML1]);

		for(int iI=0;iI<sizeI;++iI){
			instance.ConfigureI(aI[iI],aI[iI]);


			for(int iGc=0;iGc<sizeGc;++iGc){
				instance.SetSynaptic(aGc[iGc],-25,36);

				printf("Heter_gc=%.5lf I_span=%.5f ML1=%d\n",aGc[iGc],aI[iI],aML1[iML1]);
				if(isSpiral){
					//instance.OutputExcitabilityType();
					instance.SpiralWave(4000,21000,200);
				//	instance.OutputTimeSeries();
				}
				if(isSpiking){
					instance.OutputSpikingIndex();
					instance.OutputPhaseAmplitude();
					instance.OutputAverISI();
					instance.OutputAverISIForOneAndTwo();
					instance.OutputCV();
					instance.OutputCVForOneAndTwo();
//					instance.OutputTimeAndCoupleSeries();
				}
//				if(isCoupling){
//					instance.OutputSpikeTrainAndCoupleSeries();
////					instance.OutputAllCouple();
////					instance.OutputAllCoupleAverForOneAndTwo();
//				}

			}
		}
	}
}

void HeterCouple(double *aGc,int sizeGc,double *aI,int sizeI,int *aML1,int sizeML1){
	for(int iML1=0;iML1<sizeML1;iML1++){
		SquareNetwork instance(16384,aML1[iML1]);

		for(int iI=0;iI<sizeI;++iI){
			instance.ConfigureI(aI[iI],aI[iI]);


			for(int iGc=0;iGc<sizeGc;++iGc){
				instance.SetSynaptic(aGc[iGc],-25,30);

				printf("Heter_gc=%.5lf I_span=%.5f ML1=%d\n",aGc[iGc],aI[iI],aML1[iML1]);

				instance.OutputCouple();


			}
		}
	}
}
void TypeI(double gc_1,double d_gc,int total,double *aI,int sizeI,bool isSpiral,bool isSpiking)
{
	SquareNetwork instance(16384,0);

	for(int iI=0;iI<sizeI;++iI){
		instance.ConfigureI(aI[iI],aI[iI]);

		for(int iGc=0;iGc<total;++iGc){
			instance.SetSynaptic(gc_1+d_gc*iGc,-25,30);
			printf("HomoI_gc=%.5lf I_span=%.5f\n",gc_1+d_gc*iGc,aI[iI]);
			if(isSpiral)
				instance.SpiralWaveAndTimeSeries();
			if(isSpiking){
				instance.OutputSpikingIndex();
				instance.OutputPhaseAmplitude();
				instance.OutputAverISI();
			}
	    }
	}
}