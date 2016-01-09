#pragma once
#include <cstdlib>





#define PI 3.1415926535898

class neuron;
class SquareNetwork{

	friend class MorrisLecar;
public:

	explicit SquareNetwork(const int _nNeuron,const int _percentageOfML1);
	/*Destructor*/
	virtual ~SquareNetwork();    
	//Preprocess 1
	void BuildSquare();
	//Preprocess 2(optional):Specify synaptical decay time, rise time and coupling intensity
//	void SetNoiseIntensity(const double _noist_intensity);
	//Preprocess 3:Specify neural currents
	void ConfigureI(double _ML1_I_Span,double _ML2_I_Span);
	void SetSynaptic(double _gc,double _threshold,double _Vsyn);
	/*Set time step length...*/
	void SetDt(const double&_dt=0.01);




	//Output network for verify

	void OutputCoupleIdx();
//	void OutputNoise();
	void OutputI();
	void OutputNoForOneAndTwo();
	void OutputExcitabilityType();

	//process
	void OutputTimeSeries();
	void OutputCouple();
	void SpiralWave(double _begin,double _end,double _dt);
	void SpiralWave();
	void SpiralWaveAndTimeSeries();
	void OutputSpikingIndex();
	void OutputPhaseAmplitude();
	void OutputAverISI();
	void OutputAverISIForOneAndTwo();
	void OutputCV();
	void OutputCVForOneAndTwo();
	void OutputTimeAndCoupleSeries();

	void OutputAllTimeSeries();
	void OutputAllCouple();
	void OutputAllCoupleAverForOneAndTwo();

	void OutputCouplePerSpike();
	void OutputCouplePerSpikeForOneAndTwo();
	void CalculateSigmoidal(int row,int column,double *up,double *left,double *right,double *down);
	void OutputSpikeTrainAndCoupleSeries();
/*	void OutputPopulationFirings();
	void OutputPopulationFiringsOnce();
	void MaximumPotential();
	void Tau1Tau2SyncFactor(const double _tau_max,const double _interval,const double _g);
	double SyncFactor();
	double OrderParameter(const double _eval_begin,const double _iter_end);
	void Tau1Tau2OrderParameter(const double _tau_max,const double _interval,const double _g);*/
//	void OutputISI();

protected:
	const int m_nNeuron;
	const int percentageOfML1;

	const int m_nML1;
	const int m_nML2;
	

//	double m_noiseintensity;
	double m_coupleintensity;

	MorrisLecar **m_pMorrisLecar;
	static double *sm_gCouple;
//	static double *sm_gNoise;
	static double sm_dt;
	double m_t;
//	char *m_sName;

	double m_Vsyn;
	double m_threshold;


private:
	//Rand Function
	double __Uniform_01();

	//Network  iteration
	void __EulerIterate();
//	void __RungekuttaIterate();
//	void __UpdateNoise();
	void __UpdateCouple();

	//Preprocess 4:Set the initial values of variables to sychronization...
//	void __SyncVar();
	void __AlSyncVar();
	void __MakeDirectory(char *direct);
	void __MakeFilename(char *filename);
	bool __IsInner(int _No);
	bool __IsML1(int _No);
	bool __IsML2(int _No);
//	virtual void __DoMakeDirectory(char *direct)=0;
//	virtual void __DoMakeFilename(char *filename)=0;


	int ** m_pCoupleIdx;
	int *m_pCoupleNum;

	double m_ML1_I_span;
	double m_ML2_I_span;
};



