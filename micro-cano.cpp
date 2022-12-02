/*
 * File leap_frog_langevin.cpp
 *
 * DESCRIPTION: The script Computes the dynamic of the transition from 
 *              a microcanonical ensemble to a canonical one.
 *
 * Simulation of physical systems group (SSF-UN)
 * Chaos and Complexity group.
 * National University of Colombia
 *
 * Authors : Leonel Fernando Ardila Peña, | Nicolás Torres Dominguez
 * Contact : leonelardilap@gmail.com,	  | ntorresd@unal.edu.co 
 * SOURCE  : ttp://manual.gromacs.org/documentation/current/reference-manual/algorithms/stochastic-dynamics.html
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
using namespace std;

class Oscillator
{
	private:
		double X, V, Vh, Vaux, F, k, m;   // Harmonic oscillator
		double KBT, Gamma, Alpha, dV, RG; // Thermal bath
	
	public:
		void   Init(double X0, double V0, double k0, double m0, double KBT0, double Gamma0, double dt);
		void   Move(Crandom & ran64, double dt);
		double GetX(void){return X;};
		double GetV(void){return Vh;};
		double GetP(void){return m*Vh;};
		double GetE(void){return 0.5*m*Vh*Vh + 0.5*k*X*X;};
};

void Oscillator::Init(double X0, double V0, double k0, double m0, double KBT0, double Gamma0, double dt){
	X=X0; V=V0; Vh=V0; k=k0; m=m0;
	KBT = KBT0; Gamma = Gamma0; Alpha = 1.0-exp(-Gamma0*dt);
}

void Oscillator::Move(Crandom & ran64, double dt){
	F   = -k*X;
	Vh += (1.0/m)*F*dt;
	//Comment the next two lines to get the original leap-frog integrator for the harmonic oscillator
	RG  =  ran64.gauss(0.0, 1.0); // mu=0.0, sigma=1.0
	dV  = -Alpha*Vh + sqrt((KBT/m)*(1.0 - Alpha*Alpha))*RG;
	//dV  = 0.0;
	X  +=  (Vh + 0.5*dV)*dt;
	Vh +=  dV;
}

// Microcanonical Ensemble Generation E0 < Ei < E0 + dE: where dE << E_0
void Sample_microcanonical_state(double m, double k, double & xi, double & Pi, double & Vi, double Xmax, double Pmax, double Ei, double dE, Crandom & ran64){
	bool out = true; double E;
	do{
		xi = (2*ran64.r()-1)*Xmax;
		Pi = (2*ran64.r()-1)*Pmax, Vi = Pi/m;
		E  = (Pi*Pi/m+k*xi*xi)/2;
		if (E>Ei && E<Ei+dE) out = false;
	} while(out);
}

int main(void){
	
	//Number N of oscillators sampled and iteration variable i
	int i, N = 1000000;
	
	//Arrays to save the oscillators states and oscillators sampling variables
	Oscillator *Oscillators;
	Oscillators = new Oscillator[N];
	double m=1.0, k=1.0, omega=sqrt(m/k), x, p, v;
	
	//Microcanonical ensemble sampling variables
	double Ei=2.0, dE=0.005, Xmax=sqrt(2*(Ei+dE)/k), Pmax=sqrt(2*m*(Ei+dE)), Emean=Ei, SumE=0.0;
	
	//Time evolution variables
	double T=2*M_PI/omega, tmax=50.0*T, t, dt=T/100;
	
	//Thermal bath variables
	Crandom ran64(0); //Seed-initialized random number generator
	double KBT0=2.0, Gamma0=log(2.0)/dt;
	
	//Arrays to save the oscillators energy of the sampled states on the microcanonical and canonical ensemble
	double *Microcanonical_states_energies, *Canonical_states_energies;
	Microcanonical_states_energies = new double[N], Canonical_states_energies = new double[N];
	
	//Sampling initial states from a microcanonical ensemble of Energy Ei
	ofstream microcanonical_state("microcanonical_states.dat");
	for(i=0; i<N; i++){
		Sample_microcanonical_state(m, k, x, p, v, Xmax, Pmax, Ei, dE, ran64);
		Oscillators[i].Init(x, v, k, m, KBT0, Gamma0, dt); Microcanonical_states_energies[i] = Oscillators[i].GetE();
		microcanonical_state << Oscillators[i].GetX() << ", " << Oscillators[i].GetP() << endl;
	} microcanonical_state.close();
	
	//Evolution of the interaction with the thermal bath and calculation of the mean energy <E>
	ofstream mean_energy("mean_energy.dat");
	for(t=0; t<tmax; t+=dt){
		mean_energy << t << ", " << Emean << endl; SumE = 0;
		for(i=0; i<N; i++){
			Oscillators[i].Move(ran64, dt);
			SumE += Oscillators[i].GetE();
		} Emean = SumE/N;
	} mean_energy.close();
	
	//Sampling states from a canonical ensemble of Energy Ef (once thermodynamic equilibrium have been reached)
	ofstream canonical_state("canonical_states.dat");
	for(i=0; i<N; i++){
		Canonical_states_energies[i] = Oscillators[i].GetE();
		canonical_state << Oscillators[i].GetX() << ", " << Oscillators[i].GetP() << endl;
	} canonical_state.close();
	
	//Transfer heat
	ofstream HeatFile("TransferHeat.dat");
	for(i=0; i<N; i++){
		HeatFile << Canonical_states_energies[i] - Microcanonical_states_energies[i] << endl;
	} HeatFile.close();
	
	delete Oscillators;
	delete Microcanonical_states_energies;
	delete Canonical_states_energies;
	
	return 0;
}