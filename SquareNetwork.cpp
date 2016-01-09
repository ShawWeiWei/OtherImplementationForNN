#include "stdafx.h"
#include "SquareNetwork.h"
#include "ML.h"
#include <direct.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

//静态成员变量初始化
double SquareNetwork::sm_dt=0.1;
double *SquareNetwork::sm_gCouple=NULL;
//double *SquareNetwork::sm_gNoise=NULL;


SquareNetwork::SquareNetwork(const int _nNeuron,const int _percentageOfML1):m_nNeuron(_nNeuron),percentageOfML1(_percentageOfML1),m_nML1(16384.0*_percentageOfML1/100.0+0.5),m_nML2(m_nNeuron-m_nML1){
//	printf("CoupleVector is %p\n",CoupleVector);
	/*dynamically allocate memory for coupling matrix*/
	m_pCoupleIdx=new int*[m_nNeuron];
	m_pCoupleNum=new int[m_nNeuron];

	if(sm_gCouple==NULL){
		sm_gCouple=new double[m_nNeuron];
		memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	}

	m_pMorrisLecar=new MorrisLecar*[m_nNeuron];
	int part1=m_nML1;

	int *IsFilled=new int[m_nNeuron];
	memset(IsFilled,0,sizeof(int)*m_nNeuron);
	int index;
	srand((unsigned int)100);
//	srand((unsigned int)1000);

	
	for(int i=0;i<part1;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_pMorrisLecar[index]=new MorrisLecar(index);
		m_pMorrisLecar[index]->set_class1();
	}
	for(int i=0;i<m_nNeuron;++i){
		if(!IsFilled[i]){
			m_pMorrisLecar[i]=new MorrisLecar(i);
			m_pMorrisLecar[i]->set_class2();			
		}
	}
	delete[] IsFilled;
	BuildSquare();
}

void SquareNetwork::BuildSquare(){
//	sprintf_s(m_sName,40,"Square");
	int iRow,iColumn;
	int iDimension=sqrt(m_nNeuron)+0.5;
	int aTempIdx[4];
	int nTotal;

	for(int i=0;i<m_nNeuron;++i){
		iRow=i/iDimension;
		iColumn=i%iDimension;
        //up
		nTotal=0;
		if((iRow-1)>=0){
			aTempIdx[nTotal]=(iRow-1)*iDimension+iColumn;
			nTotal++;
		}

		//left
		if((iColumn-1)>=0){
			aTempIdx[nTotal]=iRow*iDimension+iColumn-1;
			nTotal++;
		}

		//down
		if((iRow+1)<iDimension){
			aTempIdx[nTotal]=(iRow+1)*iDimension+iColumn;
			nTotal++;
		}

		//right
		if((iColumn+1)<iDimension){
			aTempIdx[nTotal]=iRow*iDimension+iColumn+1;
			nTotal++;
		}
		
		m_pCoupleIdx[i]=new int[nTotal];

		for(int j=0;j<nTotal;++j){
			m_pCoupleIdx[i][j]=aTempIdx[j];
		}

		m_pCoupleNum[i]=nTotal;
	}
}

SquareNetwork::~SquareNetwork(){
	/*deallocate memory for coupling matrix*/
	
	printf("NetworkMatrix is deleted\n");
	delete [] m_pMorrisLecar;
	printf("m_gProperty is deleted\n");

	delete []sm_gCouple;
	sm_gCouple=NULL;
	printf("CoupleVector is deleted\n");
//	delete []sm_gNoise;
//	sm_gNoise=NULL;
//	printf("NoiseVector is deleted\n");
	for(int i=0;i<m_nNeuron;++i){
		delete[] m_pCoupleIdx[i];
	}
	delete []m_pCoupleNum;
	printf("m_pMorrisLecar is deleted\n");
}

//Preprocess 3:Specify neural currents

void SquareNetwork::ConfigureI(double _ML1_I_Span,double _ML2_I_Span){
	m_ML1_I_span=_ML1_I_Span;
	m_ML2_I_span=_ML2_I_Span;
	double ML1_I_End=39.7;
	double ML2_I_End=88.1;

	srand((unsigned int)100);
	double ran;
	for(int i=0;i<m_nNeuron;++i){
		ran=__Uniform_01();
		if(__IsML1(i)){
			m_pMorrisLecar[i]->SetI(ML1_I_End-ran*m_ML1_I_span);
		}
		if(__IsML2(i)){
			m_pMorrisLecar[i]->SetI(ML2_I_End-ran*m_ML2_I_span);
		}				
	}
}

void SquareNetwork::SetSynaptic(double _gc,double _threshold,double _Vsyn){
	m_coupleintensity=_gc;
	m_threshold=_threshold;
	m_Vsyn=_Vsyn;
}
//Preprocess 4:Set the initial values of variables to sychronization...
void SquareNetwork::__AlSyncVar(){
	m_t=0.0;
	printf("Initial value is randomized.\n");
	srand((unsigned int)100);               //random generator seed is 100
	for(int i=0;i<m_nNeuron;++i)
		m_pMorrisLecar[i]->RandVar();
	srand(100);
}

void SquareNetwork::SetDt(const double&_dt){
	sm_dt=_dt;
}

void SquareNetwork::__UpdateCouple(){
	memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	double sum;
	int nExc,nTotal,iIndex;
	for(int i=0;i<m_nNeuron;++i){
			sum=0.0;
			nExc=m_pCoupleNum[i];
			for(int j=0;j<nExc;++j){
				iIndex=m_pCoupleIdx[i][j];
				sum+=1.0/(1.0+exp(-(m_pMorrisLecar[iIndex]->V-m_threshold)));
			}

			sm_gCouple[i]=m_coupleintensity*(m_Vsyn-m_pMorrisLecar[i]->V)*sum;
	}
}

void SquareNetwork::__EulerIterate(){
	for(int i=0;i<m_nNeuron;++i){
		m_pMorrisLecar[i]->euler();
	}
	m_t+=sm_dt;
}

double SquareNetwork::__Uniform_01(){
	return (rand()/(double)RAND_MAX);
}
bool SquareNetwork::__IsML1(int No){
	
	if(0==strcmp(m_pMorrisLecar[No]->neuron_model,"MorrisLecar1")){
		return true;
	}
	return false;
}

bool SquareNetwork::__IsML2(int No){
	if(0==strcmp(m_pMorrisLecar[No]->neuron_model,"MorrisLecar2")){
		return true;
	}
	return false;
}

bool SquareNetwork::__IsInner(int No){
	int dimension=sqrt(m_nNeuron)+0.5;
	int row=No/dimension;
	int column=No%dimension;
	if(row>=16&&row<=111&&column>=16&&column<=111){
		return true;
	}
	return false;
}

void SquareNetwork::OutputI(){
	//Create directory
	char pre[50];
	sprintf_s(pre,"E:\\Sigmoidal_Square\\ML1_%dML2_%d",m_nML1,m_nML2);
//	_mkdir(pre);

	//Create a file pointer
	FILE *pOutputI;
	char filename[100];
	sprintf_s(filename,"%s_RandI(%.5lf,%.5lf)_I.dat",pre,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputI,filename,"w");
//	fprintf(pOutputI,"NeuronID,ML1_I,ML2_I\n");

	for(int i=0;i<m_nNeuron;++i){
//		if(__IsML1(i)){
				fprintf(pOutputI,"%d,%lf,\n",i,m_pMorrisLecar[i]->I);
//		}
//		if(__IsML2(i)){
//				fprintf(pOutputI,"%d,,%lf\n",i,m_pMorrisLecar[i]->I);
//		}
	}
}

//void SquareNetwork::OutputCouple(){
//	//Create directory
//	char *direct=new char[100];
//	char *specification=new char[40];
//	char sCouple[200];
//	__MakeDirectory(direct);
//	__MakeFilename(specification);
//	
//	sprintf_s(sCouple,"%s\\%s_RandI(%.5lf,%.5lf)_CoupleCurrent.dat",m_nML1,m_nML2);
//
//	
//	/*Create file pointer...*/
//	FILE *pOutputCouple;
//
//	fopen_s(&pOutputCouple,sCouple,"w");
//
//	__AlSyncVar();
//
//	int iter_end=1000/sm_dt+0.5;
//
//	fprintf(pOutputCouple,"%lf %lf %lf %lf %lf\n",m_t,m_pMorrisLecar[0]->V,sm_gCouple[0],m_pMorrisLecar[1]->V,sm_gCouple[1]);
//
//	for(int iter=0;iter<iter_end;++iter){
//		__UpdateCouple();
//		__EulerIterate();
//		fprintf(pOutputCouple,"%lf %lf %lf %lf %lf\n",m_t,m_pMorrisLecar[0]->V,sm_gCouple[0],m_pMorrisLecar[1]->V,sm_gCouple[1]);
//	}
//
//}

void SquareNetwork::OutputCoupleIdx(){
	FILE *pSquare;
	fopen_s(&pSquare,"E:\\Sigmoidal_Square\\ML1_%dML2_%d_CoupleIdx.dat","w");
	int iDimension=(int)(sqrt(m_nNeuron)+0.5);
	int nCount;
	for(int i=0;i<m_nNeuron;++i){
//		fprintf(pSquare,"neuron(%d)(Exc %d,Inh %d) ",i,aExcitation[i],aInhibition[i]);
		fprintf(pSquare,"neuron(%d) ",i);
		nCount=m_pCoupleNum[i];//aExcitation[i]+aInhibition[i];
		for(int j=0;j<nCount;++j){
				fprintf(pSquare,"%d ",m_pCoupleIdx[i][j]);
		}

		fprintf(pSquare,"\n");
	}
	fclose(pSquare);
}

void SquareNetwork::OutputExcitabilityType(){
	FILE *pExcitabilityType;
	char filename[200];
	sprintf_s(filename,"E:\\config\\pML1=%d_ExcitabilityType.dat",m_nML1);
	fopen_s(&pExcitabilityType,filename,"w");

	int nCount=(int)(sqrt(m_nNeuron)+0.5);
	for(int i=0;i<m_nNeuron;++i){
		for(int j=0;j<nCount;++j){
			if(__IsML1(i))
				fprintf(pExcitabilityType,"%d ",1);
			else if(__IsML2(i))
				fprintf(pExcitabilityType,"%d ",2);
			else{
				printf("ERROR EXCITABILITY TYPE!\n");
				exit(1);
			}
		}
		fprintf(pExcitabilityType,"\n");
	}
	fclose(pExcitabilityType);
}
void SquareNetwork::OutputNoForOneAndTwo(){
	srand(100);
	FILE *fp;
	char filename[100];
	sprintf_s(filename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	fopen_s(&fp,filename,"w");
	int one[5],two[5];
	int indexForOne=0,indexForTwo=0;
	int randNo;
	const int num=5;
	if(m_nML2==0){
		while(indexForOne<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
			}	
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",one[i]);
		}
		fprintf(fp,"\n");
	}
	else if(m_nML1==0){
		while(indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}		
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",two[i]);
		}
		fprintf(fp,"\n");
	}
	else{
		while(indexForOne<num&&indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}
		}
		while(indexForOne<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
			}	
		}
		while(indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}		
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",one[i]);
		}
		fprintf(fp,"\n");
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",two[i]);
		}
		fprintf(fp,"\n");
		
	}
	fclose(fp);
}





void SquareNetwork::__MakeDirectory(char *direct){
	sprintf_s(direct,100,"E:\\Sigmoidal_Square\\ML1_%d",percentageOfML1);
	_mkdir("E:\\Sigmoidal_Square");
	_mkdir(direct);
}

void SquareNetwork::__MakeFilename(char *filename){
	sprintf_s(filename,100,"gc_exc=%.5lf_Vsyn=%.5lf_threshold=%.5lf",m_coupleintensity,m_Vsyn,m_threshold);
}



//process
void SquareNetwork::OutputTimeSeries(){
	//Initialize
	__AlSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char timeSeriesFileName[200];

	__MakeDirectory(direct);
	__MakeFilename(specification);

	sprintf_s(timeSeriesFileName,"%s\\%s_RandI(%.5lf,%.5lf)_TimeSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
//	sprintf_s(coupleFileName,"%s\\%s_noise=%.5lf_Couple.dat",direct,specification,m_noiseintensity);
	/*Create file pointer...*/
	FILE* pOutputTimeSeries;
	fopen_s(&pOutputTimeSeries,timeSeriesFileName,"w");
//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans,iter_end;
	if(m_nML1==16384){
	    iter_trans=(int)(2000/sm_dt+0.5);
	    iter_end=(int)(10000/sm_dt+0.5);
	}
	else{
		iter_trans=(int)(2000/sm_dt+0.5);
		iter_end=(int)(5000/sm_dt+0.5);
	}
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputTimeSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputTimeSeries,"%lf ",m_pMorrisLecar[vec[j]]->V);
		}
		fprintf(pOutputTimeSeries,"\n");
		//printf("%d\n",i);
	}
	fclose(pOutputTimeSeries);
}

void SquareNetwork::OutputCouple(){
	//Initialize
	__AlSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char coupleSeriesFileName[200];

	__MakeDirectory(direct);
	__MakeFilename(specification);

	sprintf_s(coupleSeriesFileName,"%s\\%s_RandI(%.5lf,%.5lf)_CoupleSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
//	sprintf_s(coupleFileName,"%s\\%s_noise=%.5lf_Couple.dat",direct,specification,m_noiseintensity);
	/*Create file pointer...*/
	FILE* pOutputCoupleSeries;
	fopen_s(&pOutputCoupleSeries,coupleSeriesFileName,"w");
//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans,iter_end;
	if(percentageOfML1>95){
	    iter_trans=(int)(2000/sm_dt+0.5);
	    iter_end=(int)(5000/sm_dt+0.5);
	}
	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
	}
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputCoupleSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputCoupleSeries,"%lf ",sm_gCouple[vec[j]]);
		}
		fprintf(pOutputCoupleSeries,"\n");
		//printf("%d\n",i);
	}
	fclose(pOutputCoupleSeries);
}



void SquareNetwork::OutputAverISIForOneAndTwo(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sAverISIType1[200],sAverISIType2[200],sAverISI[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pAverISIType1,*pAverISIType2,*pAverISI;
	sprintf_s(sAverISI,"%s\\%s_RandI(%.5lf,%.5lf)_AverISI.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	sprintf_s(sAverISIType1,"%s\\%s_RandI(%.5lf,%.5lf)_AverISIType1.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sAverISIType2,"%s\\%s_RandI(%.5lf,%.5lf)_AverISIType2.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pAverISI,sAverISI,"r");
	fopen_s(&pAverISIType1,sAverISIType1,"w");
	fopen_s(&pAverISIType2,sAverISIType2,"w");
	double averISI;
	for(int i=0;i<m_nNeuron;++i){
		fscanf_s(pAverISI,"%lf \n",&averISI);
		if(0==strcmp(m_pMorrisLecar[i]->neuron_model,"MorrisLecar1")){
			fprintf(pAverISIType1,"%lf \n",averISI);
		}
		else{
			fprintf(pAverISIType2,"%lf \n",averISI);
		}
	}
	fclose(pAverISIType1);
	fclose(pAverISIType2);
	fclose(pAverISI);
}

void SquareNetwork::OutputCVForOneAndTwo(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sCVType1[200],sCVType2[200],sCV[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pCVType1,*pCVType2,*pCV;
	sprintf_s(sCV,"%s\\%s_RandI(%.5lf,%.5lf)_CV.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	sprintf_s(sCVType1,"%s\\%s_RandI(%.5lf,%.5lf)_CVType1.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sCVType2,"%s\\%s_RandI(%.5lf,%.5lf)_CVType2.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pCV,sCV,"r");
	fopen_s(&pCVType1,sCVType1,"w");
	fopen_s(&pCVType2,sCVType2,"w");
	double CV;
	for(int i=0;i<m_nNeuron;++i){
		fscanf_s(pCV,"%lf \n",&CV);
		if(0==strcmp(m_pMorrisLecar[i]->neuron_model,"MorrisLecar1")){
			fprintf(pCVType1,"%lf \n",CV);
		}
		else{
			fprintf(pCVType2,"%lf \n",CV);
		}
	}
	fclose(pCVType1);
	fclose(pCVType2);
	fclose(pCV);
}
void SquareNetwork::SpiralWave(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);
	int iter_trans,iter_end;
//	if(m_nML1==16384){
	    iter_trans=(int)(5000/sm_dt+0.5);
	    iter_end=(int)(8000/sm_dt+0.5);
//	}
//	else{
//		iter_trans=(int)(2000/sm_dt+0.5);
//		iter_end=(int)(5000/sm_dt+0.5);
//	}
//	int length=strlen(specification);
	int iColumn=int(sqrt(m_nNeuron)+0.5);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave;

	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__UpdateCouple();
		//__UpdateNoise();
		__EulerIterate();
	

		if(i%200==0){
			sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_ML1_I_span,m_ML2_I_span,m_t);
		    fopen_s(&pOutputSpiralWave,filename,"w");	

			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_pMorrisLecar[j]->V);
				if(j%iDimension==iDimension-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	}
	
	delete []direct;
	delete []specification;
	
}

void SquareNetwork::SpiralWave(double _begin,double _end,double _dt){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

//	int length=strlen(specification);
	int iColumn=int(sqrt(m_nNeuron)+0.5);


	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iBegin=int(_begin/sm_dt+0.5);
	int iEnd=int(_end/sm_dt+0.5);
	int iDt=int(_dt/sm_dt+0.5);
	for(int i=0;i<=iEnd;++i){
		if(i>=iBegin&&(i-iBegin)%iDt==0){
			/*Create file pointer...*/
			FILE* pOutputSpiralWave;
			sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_ML1_I_span,m_ML2_I_span,m_t);
			//sprintf_s(filename,"%s\\%s_noise=%.5lf_t=%.5lf.dat",direct,specification,m_t);
			fopen_s(&pOutputSpiralWave,filename,"w");
			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_pMorrisLecar[j]->V);
				if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	
		//Calculate...
		__UpdateCouple();
     	//__UpdateNoise();
		__EulerIterate();
	}
	delete []direct;
	delete []specification;
	
}

void SquareNetwork::SpiralWaveAndTimeSeries(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave,*pOutputTimeSeries;
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_TimeSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);

	fopen_s(&pOutputTimeSeries,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans,iter_end;

	if(percentageOfML1>95){
	    iter_trans=(int)(2000/sm_dt+0.5);
	    iter_end=(int)(5000/sm_dt+0.5);
	}
	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
	}

	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputTimeSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputTimeSeries,"%lf ",m_pMorrisLecar[vec[j]]->V);
		}
		fprintf(pOutputTimeSeries,"\n");
		if(i%200==0){
			sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_ML1_I_span,m_ML2_I_span,m_t);
		    fopen_s(&pOutputSpiralWave,filename,"w");	

			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_pMorrisLecar[j]->V);
				if(j%iDimension==iDimension-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	}
	
	fclose(pOutputTimeSeries);
	delete []direct;
	delete []specification;
}

void SquareNetwork::OutputSpikingIndex(){
	//Initialize
	__AlSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sSpikingIndex[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	sprintf_s(sSpikingIndex,"%s\\%s_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	/*Create file pointer...*/
	FILE* pOutputSpikingIndex;
	fopen_s(&pOutputSpikingIndex,sSpikingIndex,"w");
	//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans,iter_end;
//	if(percentageOfML1>95){
//	    iter_trans=(int)(2000/sm_dt+0.5);
//	    iter_end=(int)(5000/sm_dt+0.5);
//	}
//	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
//	}
	//Calculate
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
		__EulerIterate();
	}
	double *V1=new double[m_nNeuron];
	double *V2=new double[m_nNeuron];
	vector<int> spikingindex;
	
	vector<vector<int> > vecforspikingindex(m_nNeuron);
	for(int i=0;i<m_nNeuron;i++){
		vecforspikingindex[i].push_back(iter_trans+1);
	}
	for(int i=0;i<m_nNeuron;i++){
		V1[i]=m_pMorrisLecar[i]->V;
		V2[i]=m_pMorrisLecar[i]->V;
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__UpdateCouple();
		__EulerIterate();
		
		for(int j=0;j<m_nNeuron;++j){
			if(V2[j]>25&&V2[j]>m_pMorrisLecar[j]->V&&V2[j]>V1[j]){
				vecforspikingindex[j].push_back(i);
			}
			V1[j]=V2[j];
			V2[j]=m_pMorrisLecar[j]->V;
		}
	}
	
	for(int i=0;i<m_nNeuron;i++){
		vecforspikingindex[i].push_back(iter_end);
	}

//	int dimension=sqrt(m_nNeuron)+0.5;
//	int row,column;
	for(int i=0;i<m_nNeuron;++i){
//		row=m_nNeuron/dimension;
//		column=m_nNeuron%dimension;
//		if((row!=0)&&(row!=dimension-1)&&(column!=0)&&(column!=dimension-1)){
			for(auto it=vecforspikingindex[i].begin();it!=vecforspikingindex[i].end();it++){
				fprintf(pOutputSpikingIndex,"%d ",*it);
			}
			fprintf(pOutputSpikingIndex,"\n");
//		}
	}
	fclose(pOutputSpikingIndex);



}

void SquareNetwork::OutputPhaseAmplitude(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sSpikingIndex[200],sPhaseAmplitude[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	sprintf_s(sSpikingIndex,"%s\\%s_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sPhaseAmplitude,"%s\\%s_RandI(%.5lf,%.5lf)_PhaseAmplitude.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	getline(inputSpikingIndex,line);
	vector<int> maxiList; 
	int maxiListCount;
	istringstream findbeginandend(line);
	while(!findbeginandend.eof()){
		findbeginandend>>iSpikingIndex;
		if(!findbeginandend.fail()){
			maxiList.push_back(iSpikingIndex);
		}
	}
	maxiListCount=maxiList.size();
	int begin=maxiList[0];
	int end=maxiList[maxiListCount-1];
	const int count=end-begin+1;
	int nodecount=0;
	double *cosphase=new double[count];
	double *sinphase=new double[count];
	double *order=new double[count];
	for(int i=0;i<count;++i){
		cosphase[i]=0.0;
		sinphase[i]=0.0;
		order[i]=0.0;
	}

	double deltat,aver_deltat;
	int leftindex,rightindex;
	double phaseangle;
	
	inputSpikingIndex.seekg(0);
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			nodecount++;
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				for(int i=0;i<count;++i){
					cosphase[i]=cosphase[i]+1;
				}
//				return;
				//exit(-1);
			}
			else{
				aver_deltat=double(maxiList[maxiListCount-2]-maxiList[1])/(maxiListCount-3);
				leftindex=maxiList[0];
				rightindex=maxiList[1];
				for(int i=leftindex;i!=rightindex;i++){
					phaseangle=2.0*PI*(aver_deltat+i-rightindex)/aver_deltat;
					cosphase[i-begin]+=cos(phaseangle);
					sinphase[i-begin]+=sin(phaseangle);
				}
 
				for(int m=1;m!=maxiListCount-2;++m){
					leftindex=maxiList[m];
					rightindex=maxiList[m+1];
					deltat=double(rightindex-leftindex);
					for(int i=leftindex;i!=rightindex;++i){
						phaseangle=2.0*PI*(i-leftindex)/deltat;
						cosphase[i-begin]+=cos(phaseangle);
						sinphase[i-begin]+=sin(phaseangle);
					}
				}


                
				leftindex=maxiList[maxiListCount-2];
				rightindex=maxiList[maxiListCount-1];
				for(int i=leftindex;i<=rightindex;++i){
					phaseangle=2.0*PI*(i-leftindex)/aver_deltat;
					cosphase[i-begin]+=cos(phaseangle);
					sinphase[i-begin]+=sin(phaseangle);
				}
			}
            
		}
	}
	for(int i=0;i<count;++i){
		order[i]=sqrt(pow(cosphase[i],2)+pow(sinphase[i],2))/nodecount;
	}

	FILE *pPhaseAmplitude;
	fopen_s(&pPhaseAmplitude,sPhaseAmplitude,"w");
	for(int i=0;i<count;++i){
		fprintf(pPhaseAmplitude,"%lf\n",order[i]);
	}
	fclose(pPhaseAmplitude);
}

void SquareNetwork::OutputAverISI(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sSpikingIndex[200],sAverISI[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pAverISI;
	sprintf_s(sSpikingIndex,"%s\\%s_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sAverISI,"%s\\%s_RandI(%.5lf,%.5lf)_AverISI.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pAverISI,sAverISI,"w");
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	int iNodeIndex=0;
	double fAverISI;
	vector<int> maxiList; 
	int maxiListCount;
	vector<double> vAverISI;	
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				vAverISI.push_back(DBL_MAX);
				continue;
				//exit(-1);
			}
			fAverISI=(double)(maxiList[maxiListCount-2]-maxiList[1])*sm_dt/(maxiListCount-3);
			vAverISI.push_back(fAverISI);      
		}
	}
	for(auto it=vAverISI.begin();it!=vAverISI.end();++it){
		fprintf(pAverISI,"%lf\n",*it);
	}
	fclose(pAverISI);
}



void SquareNetwork::OutputCV(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sSpikingIndex[200],sCV[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pCV;
	sprintf_s(sSpikingIndex,"%s\\%s_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sCV,"%s\\%s_RandI(%.5lf,%.5lf)_CV.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pCV,sCV,"w");
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	int iNodeIndex=0;
	double fCV;
	vector<int> maxiList; 
	int maxiListCount;
	vector<double> vAverOfSum,vAverOfSquareSum,vCV;

	
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				vAverOfSum.push_back(0);
				vAverOfSquareSum.push_back(0);
				vCV.push_back(0.0);
				continue;
				//exit(-1);
			}
			else{
				int nISI=maxiListCount-3;
				double squareSum=0.0,sum=0.0;
				for(int i=1;i<maxiListCount-2;++i){
					double isi=(double)(maxiList[i+1]-maxiList[i])*sm_dt;
					squareSum+=isi*isi;
					sum+=isi;
				}
				double averofSum=sum/nISI;
				double averofSquareSum=squareSum/nISI;

				if((averofSquareSum-averofSum*averofSum)<0.00001){
					fCV=0.0;
				}
				else{
					fCV=sqrt(averofSquareSum-averofSum*averofSum)/averofSum;
				}
				vAverOfSum.push_back(averofSum);
				vAverOfSquareSum.push_back(averofSquareSum);
				vCV.push_back(fCV); 	
			}
		}
	}
	for(int i=0;i<m_nNeuron;++i){
		fprintf(pCV,"%lf\n",vCV[i]);
	}
	fclose(pCV);
}

void SquareNetwork::OutputTimeAndCoupleSeries(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[300];

	__MakeDirectory(direct);
	__MakeFilename(specification);
	/*Create file pointer...*/
	FILE *pOutputTimeSeries,* pOutputCoupleSeries;
	
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_TimeSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputTimeSeries,filename,"w");
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_CoupleSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputCoupleSeries,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans,iter_end;
//	if(percentageOfML1>95){
//	    iter_trans=(int)(2000/sm_dt+0.5);
//	    iter_end=(int)(8000/sm_dt+0.5);
//	}
//	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
//	}
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputTimeSeries,"%lf ",m_t);
		fprintf(pOutputCoupleSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputTimeSeries,"%lf ",m_pMorrisLecar[vec[j]]->V);
			fprintf(pOutputCoupleSeries,"%lf ",sm_gCouple[vec[j]]);
		}
		fprintf(pOutputTimeSeries,"\n");
		fprintf(pOutputCoupleSeries,"\n");
	}
	
	fclose(pOutputCoupleSeries);
	fclose(pOutputTimeSeries);
	delete []direct;
	delete []specification;
}

void SquareNetwork::OutputAllTimeSeries(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[300];

	__MakeDirectory(direct);
	__MakeFilename(specification);
	/*Create file pointer...*/
	FILE * pOutputAllTimeSeries;
	
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_AllTimeSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputAllTimeSeries,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans,iter_end;
//	if(percentageOfML1>95){
//	    iter_trans=(int)(2000/sm_dt+0.5);
//	    iter_end=(int)(8000/sm_dt+0.5);
//	}
//	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
//	}
	//Calculate
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
		if(i%10==0){
			fprintf(pOutputAllTimeSeries,"%lf ",m_t);
			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputAllTimeSeries,"%lf ",m_pMorrisLecar[j]->V);
			}
			fprintf(pOutputAllTimeSeries,"\n");
		}
	}
	
	fclose(pOutputAllTimeSeries);
	delete []direct;
	delete []specification;
}

void SquareNetwork::OutputAllCouple(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[300];

	__MakeDirectory(direct);
	__MakeFilename(specification);
	/*Create file pointer...*/
	FILE * pOutputAllCoupleSeries;
	
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_AllCoupleSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputAllCoupleSeries,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans,iter_end;
//	if(percentageOfML1>95){
//	    iter_trans=(int)(2000/sm_dt+0.5);
//	    iter_end=(int)(8000/sm_dt+0.5);
//	}
//	else{
		iter_trans=(int)(5000/sm_dt+0.5);
		iter_end=(int)(8000/sm_dt+0.5);
//	}
	//Calculate
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
		if(i%10==0){
			fprintf(pOutputAllCoupleSeries,"%lf ",m_t);
			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputAllCoupleSeries,"%lf ",sm_gCouple[j]);
			}
			fprintf(pOutputAllCoupleSeries,"\n");
		}
	}
	
	fclose(pOutputAllCoupleSeries);
	delete []direct;
	delete []specification;
}


void SquareNetwork::OutputAllCoupleAverForOneAndTwo(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sCoupleAverType1[200],sCoupleAverType2[200],sAllCoupleSeries[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pCoupleAverType1,*pCoupleAverType2,*pAllCoupleSeries;
	sprintf_s(sAllCoupleSeries,"%s\\%s_RandI(%.5lf,%.5lf)_AllCoupleSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	sprintf_s(sCoupleAverType1,"%s\\%s_RandI(%.5lf,%.5lf)_CoupleAverType1.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sCoupleAverType2,"%s\\%s_RandI(%.5lf,%.5lf)_CoupleAverType2.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pAllCoupleSeries,sAllCoupleSeries,"r");
	fopen_s(&pCoupleAverType1,sCoupleAverType1,"w");
	fopen_s(&pCoupleAverType2,sCoupleAverType2,"w");
	double *addCouple=new double[m_nNeuron];
	int nCount=0;
	double inValue;
	double inTime;
	for(int i=0;i<m_nNeuron;++i){
		addCouple[i]=0.0;
	}
	while(fscanf_s(pAllCoupleSeries,"%lf ",&inTime)!=EOF){
		nCount++;
		for(int i=0;i<m_nNeuron;++i){
			fscanf_s(pAllCoupleSeries,"%lf",&inValue);
			addCouple[i]+=inValue;

		}
		fscanf_s(pAllCoupleSeries,"\n");

	}

	for(int i=0;i<m_nNeuron;++i){
		if(0==strcmp(m_pMorrisLecar[i]->neuron_model,"MorrisLecar1")){
			fprintf(pCoupleAverType1,"%lf \n",addCouple[i]/nCount);
		}
		else{
			fprintf(pCoupleAverType2,"%lf \n",addCouple[i]/nCount);
		}
	}
	fclose(pCoupleAverType1);
	fclose(pCoupleAverType2);
	fclose(pAllCoupleSeries);
}

void SquareNetwork::OutputCouplePerSpike(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sAllCouple[200],sSpikingIndex[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	sprintf_s(sAllCouple,"%s\\%s_RandI(%.5lf,%.5lf)_AllCoupleSeries.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sSpikingIndex,"%s\\%s_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	FILE *pAllCouple;
	fopen_s(&pAllCouple,sAllCouple,"r");
	double *pSum=new double[m_nNeuron];
	for(int i=0;i<m_nNeuron;++i){
		pSum[i]=0.0;
	}
	double inDbl;
	//FILE *test;
	//fopen_s(&test,"C:\\users\\shaw\\desktop\\1.dat","w");
	int n=0;
	while(fscanf_s(pAllCouple,"%lf ",&inDbl)!=EOF){
		for(int i=0;i<m_nNeuron;++i){
			fscanf_s(pAllCouple,"%lf ",&inDbl);
			pSum[i]+=inDbl;
		}
		fscanf_s(pAllCouple,"\n");
		//n++;
		//if(n%1000==0){
		//	for(int i=0;i<m_nNeuron;++i){
		//		fprintf(test,"%lf ",pSum[i]);

		//	}
		//	fprintf(test,"\n");
		//}
	}
//	fclose(test);
	fclose(pAllCouple);
	double *pAver=new double[m_nNeuron];
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	vector<int> maxiList; 
	int maxiListCount;
	int index=0;
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount==2){        
				pAver[index]=DBL_MAX;
			}
			else{
				pAver[index]=pSum[index]/(maxiListCount-2);
			}
		}
		index++;
	}

	inputSpikingIndex.close();
	
	FILE *pCouplePerSpike;
	char sCouplePerSpike[250];
	sprintf_s(sCouplePerSpike,"%s\\%s_RandI(%.5lf,%.5lf)_CouplePerSpike.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	fopen_s(&pCouplePerSpike,sCouplePerSpike,"w");
	for(int i=0;i<m_nNeuron;++i){
		fprintf(pCouplePerSpike,"%lf\n",pAver[i]);
	}
	fclose(pCouplePerSpike);
	delete []pAver;
	delete []pSum;
}


void SquareNetwork::OutputCouplePerSpikeForOneAndTwo(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[100];
	char sCouplePerSpikeType1[200],sCouplePerSpikeType2[200],sCouplePerSpike[200];
	__MakeDirectory(direct);
	__MakeFilename(specification);

	FILE *pCouplePerSpikeType1,*pCouplePerSpikeType2,*pCouplePerSpike;
	sprintf_s(sCouplePerSpike,"%s\\%s_RandI(%.5lf,%.5lf)_CouplePerSpike.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	sprintf_s(sCouplePerSpikeType1,"%s\\%s_RandI(%.5lf,%.5lf)_CouplePerSpikeType1.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);	
	sprintf_s(sCouplePerSpikeType2,"%s\\%s_RandI(%.5lf,%.5lf)_CouplePerSpikeType2.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pCouplePerSpike,sCouplePerSpike,"r");
	fopen_s(&pCouplePerSpikeType1,sCouplePerSpikeType1,"w");
	fopen_s(&pCouplePerSpikeType2,sCouplePerSpikeType2,"w");
	double *addCouple=new double[m_nNeuron];
	double inValue;

	for(int i=0;i<m_nNeuron;++i){
		fscanf_s(pCouplePerSpike,"%lf",&inValue);
		if(0==strcmp(m_pMorrisLecar[i]->neuron_model,"MorrisLecar1")){
			fprintf(pCouplePerSpikeType1,"%lf\n",inValue);
		}
		else{
			fprintf(pCouplePerSpikeType2,"%lf\n",inValue);
		}
	}
	fclose(pCouplePerSpike);
	fclose(pCouplePerSpikeType2);
	fclose(pCouplePerSpikeType1);
}
void SquareNetwork::CalculateSigmoidal(int row,int column,double *up,double *left,double *right,double *down){
	int dimension=sqrt(m_nNeuron)+0.5;
	int upIndex=(row-1)*dimension+column;
	int leftIndex=row*dimension+column-1;
	int rightIndex=row*dimension+column+1;
	int downIndex=(row+1)*dimension+column;
	*up=1.0/(1.0+exp(-(m_pMorrisLecar[upIndex]->V-m_threshold)));
	*left=1.0/(1.0+exp(-(m_pMorrisLecar[leftIndex]->V-m_threshold)));
	*right=1.0/(1.0+exp(-(m_pMorrisLecar[rightIndex]->V-m_threshold)));
	*down=1.0/(1.0+exp(-(m_pMorrisLecar[downIndex]->V-m_threshold)));
}
void SquareNetwork::OutputSpikeTrainAndCoupleSeries(){
	/*Initialize...*/
	__AlSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[100];
	char filename[300];

	__MakeDirectory(direct);
	__MakeFilename(specification);

	char inFilename[100];
	sprintf_s(inFilename,"E:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	std::vector<int> aType1;
	std::vector<int> aType2;
	for(int i=0;i<n/2;++i){
		aType1.push_back(vec[i]);
	}
	int nType1=aType1.size();
	for(int i=n/2;i<n;++i){
		aType2.push_back(vec[i]);
	}
	int nType2=aType2.size();
	int dimension=sqrt(m_nNeuron)+0.5;
	//Calculated points

	/*Create file pointer...*/
	FILE *pOutputSynapticSeriesType1;
	FILE *pOutputSynapticSeriesType2;
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_SynapticSeriesType1.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputSynapticSeriesType1,filename,"w");
	sprintf_s(filename,"%s\\%s_RandI(%.5lf,%.5lf)_SynapticSeriesType2.dat",direct,specification,m_ML1_I_span,m_ML2_I_span);
	fopen_s(&pOutputSynapticSeriesType2,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans,iter_end;
	iter_trans=(int)(5000/sm_dt+0.5);
	iter_end=(int)(8000/sm_dt+0.5);

	//Calculate
	int header=sizeof(double *)*nType1;
	int data=sizeof(double)*nType1*4;

	double **SigmoidalType1=(double **)malloc(header+data);
	double *buf=(double *)(SigmoidalType1+nType1);
	for(int i=0;i<nType1;++i){
		SigmoidalType1[i]=buf+i*4;
	}

	header=sizeof(double *)*nType2;
	data=sizeof(double)*nType2*4;

	double **SigmoidalType2=(double **)malloc(header+data);
	buf=(double *)(SigmoidalType2+nType2);
	for(int i=0;i<nType2;++i){
		SigmoidalType2[i]=buf+i*4;
	}
	
	for(int i=1;i<=iter_trans;++i){
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
	}
	

	for(int i=iter_trans+1;i<iter_end;++i){
		for(int j=0;j<nType1;j++){
			CalculateSigmoidal(aType1[j]/dimension,aType1[j]%dimension,SigmoidalType1[j],SigmoidalType1[j]+1,SigmoidalType1[j]+2,SigmoidalType1[j]+3);
		}
		for(int j=0;j<nType2;j++){
			CalculateSigmoidal(aType2[j]/dimension,aType2[j]%dimension,SigmoidalType2[j],SigmoidalType2[j]+1,SigmoidalType2[j]+2,SigmoidalType2[j]+3);
		}
		__UpdateCouple();
//		//__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputSynapticSeriesType1,"%lf ",m_t);
		fprintf(pOutputSynapticSeriesType2,"%lf ",m_t);
		for(int j=0;j<nType1;++j){			
			fprintf(pOutputSynapticSeriesType1,"%lf ",m_pMorrisLecar[aType1[j]]->V);
			fprintf(pOutputSynapticSeriesType1,"%lf ",sm_gCouple[aType1[j]]);
			fprintf(pOutputSynapticSeriesType1,"%lf ",m_Vsyn-m_pMorrisLecar[aType1[j]]->V);
			fprintf(pOutputSynapticSeriesType1,"%lf %lf %lf %lf ",SigmoidalType1[j][0],SigmoidalType1[j][1],SigmoidalType1[j][2],SigmoidalType1[j][3]);
			
		}		
		fprintf(pOutputSynapticSeriesType1,"\n");
		for(int j=0;j<nType2;++j){			
			fprintf(pOutputSynapticSeriesType2,"%lf ",m_pMorrisLecar[aType2[j]]->V);
			fprintf(pOutputSynapticSeriesType2,"%lf ",sm_gCouple[aType2[j]]);
			fprintf(pOutputSynapticSeriesType2,"%lf ",m_Vsyn-m_pMorrisLecar[aType2[j]]->V);
			fprintf(pOutputSynapticSeriesType2,"%lf %lf %lf %lf ",SigmoidalType2[j][0],SigmoidalType2[j][1],SigmoidalType2[j][2],SigmoidalType2[j][3]);
		}	
		fprintf(pOutputSynapticSeriesType2,"\n");
	}

	fclose(pOutputSynapticSeriesType1);
	fclose(pOutputSynapticSeriesType2);
	free(SigmoidalType1);
	free(SigmoidalType2);
	delete []direct;
	delete []specification;
}
//double SquareNetwork::SyncFactor(){
//	double *sum_VikSquare=new double[m_nNeuron];
//	double *sum_Vik=new double[m_nNeuron];
//	double sum_averVkSquare=0.0,sum_averVk=0.0;
//
//	double sumVecField=0.0,averVecField=0.0;
//	double numerator=0.0,denominator=0.0,R=0.0;
//	const int iter_end=40001;
//	const int eval_begin=20000;
//	const int nTime=iter_end-eval_begin-1;
//	memset(sum_VikSquare,0,sizeof(double)*m_nNeuron);
//	memset(sum_Vik,0,sizeof(double)*m_nNeuron);
//	__AlSyncVar();
//	//Calculate
//	for(int iter=1;iter<iter_end;++iter){
//	 	__UpdateCouple();
//	//	//__UpdateNoise();
//
//		__EulerIterate();
//		if(iter>eval_begin){
//			sumVecField=0.0;
//			for(int index=0;index<m_nNeuron;++index){
//				sumVecField+=m_pMorrisLecar[index]->V;
//				sum_Vik[index]+=m_pMorrisLecar[index]->V;
//				sum_VikSquare[index]+=(m_pMorrisLecar[index]->V)*(m_pMorrisLecar[index]->V);
//			}
//			averVecField=sumVecField/m_nNeuron;
//			sum_averVk+=(averVecField);
//			sum_averVkSquare+=(averVecField*averVecField);
//		}
//	}
//
//	numerator=((sum_averVkSquare/nTime)-(sum_averVk/nTime)*(sum_averVk/nTime));
//	for(int index=0;index<m_nNeuron;++index){
//			denominator+=((sum_VikSquare[index]/nTime)-(sum_Vik[index]/nTime)*(sum_Vik[index]/nTime));
//	}
//	denominator/=m_nNeuron;
//	R=numerator/denominator;
//
//	delete []sum_Vik;
//	delete []sum_VikSquare;
//
//	return R;
//
//}

//void SquareNetwork::OutputPopulationFirings(){
//	//Create directory
//	char *direct=new char[100];
//	char *specification=new char[40];
//	char sSpikingIndex[200],sPopulationFirings[200];
//	__MakeDirectory(direct);
//	__MakeFilename(specification);
//
//	FILE *pPopulationFirings;
//	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_span,m_ML2_I_span);	
//	sprintf_s(sPopulationFirings,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_PopulationFirings.dat",direct,specification,m_noiseintensity,m_ML1_I_span,m_ML2_I_span);
//	fopen_s(&pPopulationFirings,sPopulationFirings,"w");
//	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
//	string line;
//	int iSpikingIndex;
//	int iNodeIndex=0,iNodeRow,iNodeColumn;
//	int iDimension=(int)(sqrt(m_nNeuron)+0.5);
//
//	vector<int> maxiList; 
//	int maxiListCount;
//
//	while(!inputSpikingIndex.eof()){
//		getline(inputSpikingIndex,line);
//		maxiList.clear();
//		if(!inputSpikingIndex.fail()){
//			iNodeRow=iNodeIndex/iDimension;
//			iNodeColumn=iNodeIndex%iDimension;
//			iNodeIndex++;
//			istringstream instring(line);
//			while(!instring.eof()){
//				instring>>iSpikingIndex;
//				if(!instring.fail()){
//					maxiList.push_back(iSpikingIndex);
//				}
//			}
//			maxiListCount=maxiList.size();
//			if(maxiListCount<4){
//				printf("There are not spikes!\n");
//				continue;
//				//exit(-1);
//			}
//			double index2t;
//			for(int it=1;it!=maxiListCount-1;++it){
//				index2t=(double)maxiList[it]*sm_dt;
//				fprintf(pPopulationFirings,"%d %d %lf\n",iNodeRow,iNodeColumn,index2t);
//			}    
//		}
//	}
//	
//	fclose(pPopulationFirings);
//}

//void SquareNetwork::OutputPopulationFiringsOnce(){
//	//Create directory
//	char *direct=new char[100];
//	char *specification=new char[40];
//	char sSpikingIndex[200],sPopulationFirings[200];
//	__MakeDirectory(direct);
//	__MakeFilename(specification);
//
//	FILE *pPopulationFirings;
//	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_span,m_ML2_I_span);	
//	sprintf_s(sPopulationFirings,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_PopulationFiringsOnce.dat",direct,specification,m_noiseintensity,m_ML1_I_span,m_ML2_I_span);
//	fopen_s(&pPopulationFirings,sPopulationFirings,"w");
//	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
//	string line;
//	int iSpikingIndex;
//	int iNodeIndex=0,iNodeRow,iNodeColumn;
//	int iDimension=(int)(sqrt(m_nNeuron)+0.5);
//
//	vector<int> maxiList; 
//	int maxiListCount;
//
//	while(!inputSpikingIndex.eof()){
//		getline(inputSpikingIndex,line);
//		maxiList.clear();
//		if(!inputSpikingIndex.fail()){
//			iNodeRow=iNodeIndex/iDimension;
//			iNodeColumn=iNodeIndex%iDimension;
//			iNodeIndex++;
//			istringstream instring(line);
//			while(!instring.eof()){
//				instring>>iSpikingIndex;
//				if(!instring.fail()){
//					maxiList.push_back(iSpikingIndex);
//				}
//			}
//			maxiListCount=maxiList.size();
//			if(maxiListCount<4){
//				printf("There are not spikes!\n");
//				continue;
//				//exit(-1);
//			}
//			double index2t=(double)maxiList[1]*sm_dt;
//			fprintf(pPopulationFirings,"%d %d %lf\n",iNodeRow,iNodeColumn,index2t);
//	 
//		}
//	}	
//	fclose(pPopulationFirings);
//}
//void SquareNetwork::MaximumPotential(){
//	/*Initialize...*/
//	__AlSyncVar();
//	/*Create directory...*/
//
//	char *direct=new char[100];
//	char *specification=new char[40];
//	char filename[140];
//	__MakeDirectory(direct);
//	__MakeFilename(specification);
//	int iColumn=int(sqrt(m_nNeuron)+0.5);
//	/*Create file pointer...*/
//	FILE* pOutputSpiralWave;
//	double *aMaximumV=new double[m_nNeuron];
//	double *aMinimumV=new double[m_nNeuron];
//	for(int i=0;i<m_nNeuron;++i){
//		aMaximumV[i]=-DBL_MAX;
//	}
//	for(int i=0;i<m_nNeuron;++i){
//		aMinimumV[i]=DBL_MAX;
//	}
//	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
//	for(int i=0;i<30000;++i){
//		if(i>20000&&i%20==0){
//			for(int j=0;j<m_nNeuron;++j){
//				if(aMaximumV[j]<m_pMorrisLecar[j]->V)
//					aMaximumV[j]=m_pMorrisLecar[j]->V;
//				if(aMinimumV[j]>m_pMorrisLecar[j]->V)
//					aMinimumV[j]=m_pMorrisLecar[j]->V;
//			}
//		}
//	
//		//Calculate...
//		__UpdateCouple();
//     	//__UpdateNoise();
//		__EulerIterate();
//	}
//	sprintf_s(filename,"%s\\%s_noise=%.5lf_MaximumV.dat",direct,specification,m_noiseintensity);
//	fopen_s(&pOutputSpiralWave,filename,"w");
//	for(int j=0;j<m_nNeuron;++j){
//		fprintf(pOutputSpiralWave,"%lf ",aMaximumV[j]);
//		if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
//	}
//	fprintf(pOutputSpiralWave,"\n");
//	fclose(pOutputSpiralWave);
//
//	sprintf_s(filename,"%s\\%s_noise=%.5lf_MinimumV.dat",direct,specification,m_noiseintensity);
//	fopen_s(&pOutputSpiralWave,filename,"w");
//	for(int j=0;j<m_nNeuron;++j){
//		fprintf(pOutputSpiralWave,"%lf ",aMinimumV[j]);
//		if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
//	}
//	fprintf(pOutputSpiralWave,"\n");
//	fclose(pOutputSpiralWave);
//
//	delete []direct;
//	delete []specification;
//	delete []aMaximumV;
//}



