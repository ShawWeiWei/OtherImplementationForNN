#pragma once


class MorrisLecar{
	friend class SquareNetwork;
public:
	MorrisLecar(); 
	MorrisLecar(const int _No);
	~MorrisLecar();
	void equation();
	void euler();
	void rungekutta();
	void SetI(double _I);
//	void equation_s();
	void euler_s();
	void rungekutta_s();


    void InitVar();
	void RandVar();

	void set_class1();
	void set_class2();

private:

	const int No;

	
	double V;
	double n;
	/*parameters*/
	double gCa;
	double phi;
	double V3;
	double V4;
	double I;

	/*intermediate variables for equations*/
	double DV;
	double Dn;
	
	/*intermediate variables for updating V and n*/
	double V0;
	double n0;

	/*intermediate variables for rungekutta method*/
	double DV1;
	double Dn1;
	double DV2;
	double Dn2;
	double DV3;
	double Dn3;
	double DV4;
	double Dn4;

	/*variables to compare to determine spike time*/
	/*double V_pre;
	double V_pre1;*/

	double single_dt;
	char neuron_model[15];
};