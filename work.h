#pragma once


void Heter(double *aGc,int sizeGc,double *aI,int sizeI,int *aML1,int sizeML1,bool isSpiral,bool isSpiking,bool isCoupling);

void HeterCouple(double *aGc,int sizeGc,double *aI,int sizeI,int *aML1,int sizeML1);

void TypeI(double gc_1,double d_gc,int total,double *aI,int sizeI,bool isSpiral,bool isSpiking);

void Config(int *aML1,int sizeML1);