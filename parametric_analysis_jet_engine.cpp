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
// Uses Procedural Programing





int main(){


// inputs
	// Flight Conditions
		//float mach_num = 2.4; 				// M/a
		//float air_density = 1.225; 			// lbm/ft^3
		//int sos_at_cruise = 293;				// mph
		//int thrust_req = 21700;				// lbf
	// Design Limits
		float max_temp = 3200;				// Rankine
		int hpr = 18400; 					// Kj/Kg K
		float t_zero = 390.5;
		float gas_constant = .24;
		float mach_zero = 1.6;
		float g_of_c = 32.174;	//32.174 for Imperial units and 9.81 for Si
		float pi_comp_max = 27;
		float gamma = 1.4;
		float pi_fan_max = 4;
		float a_zero = 968.8;
		string filename = "data.csv";	//Where data is to end up. Default is source file location
		

// optimization settings
	float alpha_learning_rate = 1;
	float tau_c_learning_rate = 1;
	float tau_f_learning_rate = 1;
	int epoc = 100;
	float alpha = 5.3; //Initial Guess
	float pi_fan = 3.8; //Initial Guess
	float pi_comp = 16;	//Initial Guess
	


// Pre-calculated Values
	float tau_lambda = max_temp/t_zero;
	float tau_r = 1 + ( (gamma - 1)*(.5)*(mach_zero*mach_zero) );
	//float a_zero = sqrt((gamma - 1) * gas_constant * g_of_c * t_zero); //doesnt work for all calculations
	float v_multiplier = (2)/(gamma - 1);
	float spec_thr_multiplier = a_zero / g_of_c;
	float lit_f_multiplier = gas_constant * t_zero /hpr;
	float tau_c_max = pow(pi_comp_max,((gamma-1)/gamma));
	float tau_f_max = pow(pi_fan_max,((gamma-1)/gamma));
	float tau_c = pow(pi_comp,((gamma-1)/gamma));
	float tau_f = pow(pi_fan,((gamma-1)/gamma));

// Data Array Allocation
	float *sfc_data = new float[epoc];
	float *tau_c_data = new float[epoc];
	float *tau_f_data = new float[epoc];
	float *alpha_data = new float[epoc];



//optimizing for specific fuel consumption
for(int i=0;i<epoc;i++){

	//Checks if Physically Possible and Constraints
	tau_c = (tau_c > tau_c_max) ? tau_c_max : tau_c;
	tau_f = (tau_f > tau_f_max) ? tau_f_max : tau_f;
	tau_c = (tau_r*tau_c > tau_lambda) ? tau_c_max : tau_c;
	alpha = (alpha > 0) ? alpha : 0;


	//Specific Fuel Consumption Calculations
	float inner_root = tau_lambda - tau_r*( tau_c - 1 + alpha*(tau_f - 1) ) - (tau_lambda/(tau_r*tau_c));

	float v_nine = sqrt( v_multiplier * inner_root);

		if(v_nine != v_nine){cout<<"DOMAIN ERROR "<<"\n"<<"Alpha: "<<alpha<<" Compressor Pressure Ratio: "<<pi_comp<<" Fan Pressure Ratio: "<<pi_fan<<endl;break;}

	float v_nineteen = sqrt( v_multiplier * ((tau_r * tau_f) - 1) );

		if(v_nineteen != v_nineteen){cout<<"DOMAIN ERROR "<<"\n"<<"Alpha: "<<alpha<<" Compressor Pressure Ratio: "<<pi_comp<<" Fan Pressure Ratio: "<<pi_fan<<endl;break;}

	float little_f = lit_f_multiplier * (tau_lambda - (tau_r*tau_c));

	float specific_thrust = spec_thr_multiplier * ( v_nine - mach_zero + (alpha * (v_nineteen - mach_zero)) ); // Multiply with (1/(1+alpha)) to get actual specific thrust

	float sfc = (3600*little_f) / specific_thrust;



	//Partial-derivative Calculation
		//Help with alpha
		float partial_specific_thrust_wrt_alpha = spec_thr_multiplier * ( (-1*tau_r*(tau_f - 1)*v_multiplier)/(2*v_nine) + v_nineteen - mach_zero );
		//Help for tau_c
		float partial_wrt_tau_c_specific_thrust = (spec_thr_multiplier)*(-1/v_nine)*(v_multiplier)* ( ( tau_lambda / (tau_r*tau_c*tau_c) ) - tau_r );
		float partial_wrt_tau_c_little_f = -1*lit_f_multiplier*tau_r;
		//Help for tau_f
		float partial_specific_thrust_wrt_tau_f = spec_thr_multiplier*( (1/(2*v_nine))*(v_multiplier)*(-tau_r*alpha) + (alpha)*(1/(2*v_nineteen))*(v_multiplier*tau_r) );


	float partial_wrt_alpha = (-1*little_f*partial_specific_thrust_wrt_alpha)/(specific_thrust*specific_thrust);
	float partial_wrt_tau_c = (partial_wrt_tau_c_little_f*specific_thrust - little_f*partial_wrt_tau_c_specific_thrust) / (specific_thrust*specific_thrust);
	float partial_wrt_tau_f = (-1*little_f/(specific_thrust*specific_thrust))*(partial_specific_thrust_wrt_tau_f);


	// Updates Variables
	alpha -= alpha_learning_rate*partial_wrt_alpha; 
	tau_c -= tau_c_learning_rate*partial_wrt_tau_c;
	tau_f -= tau_f_learning_rate*partial_wrt_tau_f;



	//Creating array of convergence data and Prints Progress
	sfc_data[i] = sfc;
	tau_c_data[i] = tau_c;
	tau_f_data[i] = tau_f;
	alpha_data[i] = alpha;

	if(i % 1000 == 0){cout<<"Iteration: "<<i<<"\n	Sfc:   "<<sfc<<"\n	Alpha: "<<alpha<<"\n	Pi_c:  "<<pow(tau_c,(gamma/(gamma-1)))<<"\n	pi_f:  "<<pow(tau_f,(gamma/(gamma-1)))<<"\n"<<endl;}

}



//Writes CSV File With Convergence Data
ofstream myfile;
myfile.open(filename);
myfile<<"Tau_c Data , Tau_f Data , Alpha Data , Sfc Data \n";
for(int i=0;i<epoc;i++){
	myfile<<tau_c_data[i]<<","<<tau_f_data[i]<<","<<alpha_data[i]<<","<<sfc_data[i]<<"\n";
}
myfile.close();



// Revert Temperature Ratios to Pressure Ratios and Prints to Console
pi_comp = pow(tau_c,(gamma/(gamma-1)));
pi_fan = pow(tau_f,(gamma/(gamma-1)));
cout<<"=========== Final Output =========== \n	Alpha: "<<alpha<<"\n	Pi_c:  "<<pow(tau_c,(gamma/(gamma-1)))<<"\n	pi_f:  "<<pow(tau_f,(gamma/(gamma-1)))<<"\n"<<endl;



}
