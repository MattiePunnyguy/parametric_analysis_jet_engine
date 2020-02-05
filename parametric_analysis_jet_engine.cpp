#include<math.h>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<random>
#include<iostream>
#include<fstream>
#include<tgmath.h>
#include<thread>

using namespace std;




// Uses stochastic gradient descent to solve ideal equations for best specific fuel consumption.
// Stocastic Gradient Descent is guaranteed to find a local minima, but not guaranteed to find a global minima.
// The ideal equaitons were taken from the book isbn13 	.





int main(){


// inputs
	// Flight Conditions
	/*	float mach_num = 2.4; 				// M/a
		float air_density = 1.225; 			// kg/m^3
		int sos_at_cruise = ;				// mph
	// Design Limits
		float max_temp = 2700;				// Rankine
		int hpr = 43150; 					// Kj/Kg K
		int thrust_req = 21700;	*/			// lbf
		

// optimization variable settings
	float alpha_learning_rate = .01;
	float tau_c_learning_rate = 10;
	float tau_f_learning_rate = .01;
	int epoc = 10000;
	float pi_fan = 1.7;
	float pi_comp = 14;
	string filename = "data.csv";



// Constant Variables
	// Constants
	float gas_constant = 1.004;
	float gamma = 1.4;
	float g_of_c = 9.81; //32.174;
	
	float hpr = 42800;
	float t_zero = 216.7;
	float mach_zero = .9; // not done
	float max_temp = 1670;
	float tau_c_max = 27;
	float tau_f_max = 4;



	// Pre-calculated Values
	float tau_lambda = max_temp/t_zero;
	float tau_r = 1 + ( (gamma - 1)*(.5)*(mach_zero*mach_zero) );
	float a_zero = sqrt((gamma - 1) * gas_constant * g_of_c * t_zero);
	float v_multiplier = (2)/(gamma - 1);
	float spec_thr_multiplier = a_zero / g_of_c;
	float lit_f_multiplier = gas_constant * t_zero /hpr;
	float *data = new float[epoc];



// Variables
float alpha = 3;
float tau_c = pow(pi_comp,((gamma-1)/gamma));
float tau_f = pow(pi_fan,((gamma-1)/gamma));



//optimizing for specific fuel consumption
for(int i=0;i<epoc;i++){

	//Checks if Physically Possible
	tau_c = (tau_c > tau_c_max) ? tau_c_max : tau_c;
	tau_f = (tau_f > tau_f_max) ? tau_f_max : tau_f;


	//Specific Fuel Consumption Calculations
	float inner_root = tau_lambda - tau_r*(tau_c - 1 + alpha*(tau_f - 1)) - (tau_lambda/(tau_r*tau_c));

	float v_nine = sqrt( v_multiplier * inner_root);

		if(v_nine != v_nine){cout<<"Domain Error"<<endl;exit(1);}

	float v_nineteen = sqrt( v_multiplier * ((tau_r * tau_f) - 1) );

		if(v_nineteen != v_nineteen){cout<<"Domain Error"<<endl;exit(1);}

	float little_f = lit_f_multiplier * (tau_lambda - (tau_r*tau_c));

	float specific_thrust = spec_thr_multiplier * (v_nine - mach_zero + (alpha * (v_nineteen - mach_zero)));

	float sfc = little_f / specific_thrust;



	//Partial-derivative Calculation
		//Help for tau_c
		float partial_wrt_tau_c_specific_thrust = (spec_thr_multiplier)*(-1/v_nine)*(v_multiplier)* ( ( tau_lambda / (tau_r*tau_c*tau_c) ) - tau_r );
		float partial_wrt_tau_c_little_f = -1*lit_f_multiplier*tau_r;
		//Help for tau_f
		float partial_specific_thrust_wrt_tau_f = spec_thr_multiplier*( (-1/v_nine)*(v_multiplier)*(-tau_r*alpha) + (alpha)*(-1/v_nineteen)*(v_multiplier*tau_r) );


	float partial_wrt_alpha = ( (-1 * tau_r * (tau_f - 1) * (.5))/sqrt(inner_root/(tau_r - 1)) ) + sqrt( ((tau_r * tau_f) - 1) / (tau_r - 1) ) - 1;
	float partial_wrt_tau_c = (partial_wrt_tau_c_little_f*specific_thrust - little_f*partial_wrt_tau_c_specific_thrust) / (specific_thrust*specific_thrust);
	float partial_wrt_tau_f = (-1*little_f/specific_thrust)*(partial_specific_thrust_wrt_tau_f);


	// Updates Variables
	alpha += alpha_learning_rate*partial_wrt_alpha; 
	tau_c -= tau_c_learning_rate*partial_wrt_tau_c;	// All come out weird
	tau_f += tau_f_learning_rate*partial_wrt_tau_f;



	//Creating array of convergence data and prints Progress
	data[i] = sfc;
	if(i % 1000 == 0){cout<<"Alpha: "<<alpha<<" Compressor Pressure Ratio: "<<pi_comp<<" Fan Pressure Ratio: "<<pi_fan<<endl;}

}



//Writes CSV File With Convergence Data
ofstream myfile;
myfile.open(filename);
	for(int i=0;i<epoc-1;i++){
	myfile<<data[i]<<",";
	}
myfile<<data[epoc-1];
myfile.close();



// Revert Temperature Ratios to Pressure Ratios and Prints to Console
pi_comp = pow(tau_c,(gamma/(gamma-1)));
pi_fan = pow(tau_f,(gamma/(gamma-1)));
cout<<"Alpha: "<<alpha<<" Compressor Pressure Ratio: "<<pi_comp<<" Fan Pressure Ratio: "<<pi_fan<<endl;



}
