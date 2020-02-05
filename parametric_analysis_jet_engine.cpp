#include<math.h>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<random>
#include<iostream>
#include<tgmath.h>
#include<thread>

using namespace std;



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
	float learning_rate = .01;
	int epoc = 1000;
	int variables = 3;



// Global Variables
	// Constant Inputs
	float gas_constant = 1.004;
	float gamma = 1.4;
	float g_of_c = 9.81; //32.174;
	
	float hpr = 42800;
	float t_zero = 216.7;
	float mach_zero = .9; // not done
	float max_temp = 1670;



	// Pre-calculated Values
	float tau_lambda = max_temp/t_zero;
	float tau_r = 1 + ( (gamma - 1)*(.5)*(mach_zero*mach_zero) );
	float a_zero = sqrt((gamma - 1) * gas_constant * g_of_c * t_zero);
	float v_multiplier = (2)/(gamma - 1);
	float spec_thr_multiplier = a_zero / g_of_c;
	float lit_f_multiplier = gas_constant * t_zero /hpr;



// Variables
float alpha = 5;
float tau_c = pow(25,((gamma-1)/gamma));
float tau_f = pow(2,((gamma-1)/gamma));

float var_array[3] = {alpha,tau_c,tau_f};





//optimizing for specific fuel consumption

	for(int i=0;i<epoc;i++){

		float temp = 0;
		float current_best = var_array[i] - learning_rate;



		for(int j=0;j<variables;j++){


			float inner_root = tau_lambda - tau_r*(tau_c - 1 + alpha*(tau_f - 1)) - (tau_lambda/(tau_r*tau_c));



			float v_nine = sqrt( v_multiplier * inner_root);

				if(v_nine != v_nine){cout<<"Domain Error"<<endl;}

			float v_nineteen = sqrt( v_multiplier * ((tau_r * tau_f) - 1) );

				if(v_nineteen != v_nineteen){cout<<"Domain Error"<<endl;}

			float little_f = lit_f_multiplier * (tau_lambda - (tau_r*tau_c));

			float specific_thrust = spec_thr_multiplier * (v_nine - mach_zero + (alpha * (v_nineteen - mach_zero)));

			float sfc = little_f / specific_thrust;



			float partial_wrt_tau_c_specific_thrust = (spec_thr_multiplier)*(1/inner_root)*(v_multiplier)* ( ( tau_lambda / (tau_r*tau_c*tau_c) ) - tau_r );

			float partial_wrt_tau_c_little_f = -1*lit_f_multiplier*tau_r;

			float partial_wrt_alpha = ( (-1 * tau_r * (tau_f - 1) * (.5))/sqrt(inner_root/(tau_r - 1)) ) + sqrt( ((tau_r * tau_f) - 1) / (tau_r - 1) ) - 1;
			float partial_wrt_tau_c = (partial_wrt_tau_c_little_f*specific_thrust - little_f*partial_wrt_tau_c_specific_thrust) / (specific_thrust*specific_thrust);
			float partial_wrt_tau_f = ;



			//alpha += learning_rate*alpha*partial_wrt_alpha;
			//tau_c += learning_rate*alpha*partial_wrt_tau_c;
			//tau_f += learning_rate*alpha*partial_wrt_tau_f;

		}

		var_array[i] = current_best;

	}
























// Revert to Pressure Ratios
float pi_comp = pow(tau_c,(gamma/(gamma-1)));
float pi_fan = pow(tau_f,(gamma/(gamma-1)));




}
