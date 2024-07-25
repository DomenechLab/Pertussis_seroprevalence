// Implementation of the age-structured pertussis serotransmission model

#include <math.h>

//start_globs

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t)
static double pTrans(double r, double delta_t) {
	
	// r: event (instantaneous rate)
	// delta_t: time step
	// Returns: probability of transition
	double out = 1.0 - exp(-r * delta_t);
	return out;
}
//end_globs

//start_stochSim

int i, j; 
double birth_no; // No of births
double lambda_vec[nA]; // Force of infection
double q_vec[nA]; // Susceptibility to infection
double v_rate[nA]; // Effective vaccination rate
double F_mat[nA][nA]; // Contact matrix (density scale)
double from_E1[nA][2], from_I1[nA][2], from_E2[nA][2], from_I2[nA][2];
double from_Re[nA][2], from_Rp1[nA][2], from_Rp2[nA][2], from_Ve[nA][2], from_Vp[nA][2]; 
double from_S1[nA][3], from_S2[nA][3], from_R[nA][3], from_V[nA][3]; 
double rate2[2], trans2[2];
double rate3[3], trans3[3]; 

//Pointers to state variables
double *S1_vec = (double *) &S1_1; // Susceptible to primary infection
double *E1_vec = (double *) &E1_1; // Exposed, primary infection
double *I1_vec = (double *) &I1_1; // Infected, primary infection
double *S2_vec = (double *) &S2_1; // Susceptible to secondary infection
double *E2_vec = (double *) &E2_1; // Exposed, secondary infection
double *I2_vec = (double *) &I2_1; // Infected, secondary infection
double *R_vec = (double *) &R_1; // Recovered
double *Re_vec = (double *) &Re_1; // Exposed, from recovered state 
double *Rp1_vec = (double *) &Rp1_1; // Seropositive, after infection
double *Rp2_vec = (double *) &Rp2_1; // Seropositive, after exposure from recovered state
double *V_vec = (double *) &V_1; // Vaccinated
double *Ve_vec = (double *) &Ve_1; // Exposed, from vaccinated state 
double *Vp_vec = (double *) &Vp_1; // Seropositive, after exposure from vaccinated state
double *Ci1_vec = (double *) &Ci1_1; // Incidence of primary infections
double *Ci2_vec = (double *) &Ci2_1; // Incidence of secondary infections
double *Cs_vec = (double *) &Cs_1; // Sero-incidence

// Vectors of parameters
double *N_vec = (double *) &N_1; // age-specific population size
double *delta_vec = (double *) &delta_1; // age-specific aging rate
double *p_V_vec = (double *) &p_V_1; // age-specific vaccination coverage
double *f_vec = (double *) &f_1; // Contact matrix (density scale) stacked as vector (line by line), size nA*nA

// Extra parameters

// Susceptibility vector and vaccination rate
for(i = 0; i < nA; i++) {
	q_vec[i] = (i <= 10) ? q1 : ((i <= 20) ? (q1 * q21) : (q1 * q21 * q32)); 
	v_rate[i] = (t < t_V) ? 0.0 : (p_V_vec[i] / (1.0 - p_V_vec[i]) * delta_vec[i]);
}

// Re-create contact matrix
for(i = 0; i < nA; i++) {
	for(j = 0; j < nA; j++) {
		F_mat[i][j] = f_vec[nA * i + j];
	}
}

// Check values of parameters, for debugging
if(debug) {
	/*for(i = 0; i < nA; i++) {
		for(j = 0; j < nA; j++) {
			Rprintf("%.2f ", m_mat[i][j]);
		}
		Rprintf("\n");
	}*/
	Rprintf("\n\n");
//	for(i = 0; i < nA; i++) {
//		Rprintf("gamma_%d = %.2f, beta_%d = %.3f, delta_%d = %.4f\n\n", i, gamma_vec[i], i, beta_vec[i], i, delta_vec[i]);
//	}
//	Rprintf("t = %.2f, tV = %.1f, pV = %.1f, v_t = %.1f\n\n", t, tV, pV, v_t); 
}

// Calculate forces of infections
for(i = 0; i < nA; i++) {
	lambda_vec[i] = 0.0; 

	for(j = 0; j < nA; j++) {
		lambda_vec[i] += F_mat[i][j] * (I1_vec[j] + theta * I2_vec[j] + iota); 
	}
	lambda_vec[i] *= q_vec[i]; 
}

// Sample no of stochastic transitions 
birth_no = rpois(b_rate * N_tot * dt); 
for(i = 0; i < nA; i++) {
	
	// Sample stochastic transitions
	
	// From S1
	rate3[0] = delta_vec[i]; // Aging
	rate3[1] = lambda_vec[i]; // S1->E1, S2->E2
	rate3[2] = v_rate[i]; // S1->V, S2->V 
	reulermultinom(3, S1_vec[i], &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) from_S1[i][j] = trans3[j]; 
	
	// From S2
	reulermultinom(3, S2_vec[i], &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) from_S2[i][j] = trans3[j]; 
	
	// From E1
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = sigma; // E1->I1, E2->I2 
	reulermultinom(2, E1_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_E1[i][j] = trans2[j]; 
	
	// From E2
	reulermultinom(2, E2_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_E2[i][j] = trans2[j]; 
	
	// From I1
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = gamma; // E1->I1, E2->I2 
	reulermultinom(2, I1_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_I1[i][j] = trans2[j]; 
	
	// From I2
	reulermultinom(2, I2_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_I2[i][j] = trans2[j]; 

	// From R
	rate3[0] = delta_vec[i]; // Aging
	rate3[1] = alpha_R; // R->S2
	rate3[2] = rho_R * lambda_vec[i]; // R->Re 
	reulermultinom(3, R_vec[i], &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) from_R[i][j] = trans3[j]; 
	
	// From Re
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = 1.0 / t_p_R; // Re->Rp2
	reulermultinom(2, Re_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_Re[i][j] = trans2[j]; 
	
	// From Rp1
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = 1.0 / t_n_I; // Rp1->R
	reulermultinom(2, Rp1_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_Rp1[i][j] = trans2[j]; 
	
	// From Rp2
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = 1.0 / t_n_R; // Rp2->R
	reulermultinom(2, Rp2_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_Rp2[i][j] = trans2[j]; 
	
	// From V 
	rate3[0] = delta_vec[i]; // Aging
	rate3[1] = alpha_V; // V->S2
	rate3[2] = rho_V * lambda_vec[i]; // V->Ve 
	reulermultinom(3, V_vec[i], &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) from_V[i][j] = trans3[j];
	
	// From Ve
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = 1.0 / t_p_V; // Ve->Vp
	reulermultinom(2, Ve_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_Ve[i][j] = trans2[j]; 
	
	// From Vp
	rate2[0] = delta_vec[i]; // Aging
	rate2[1] = 1.0 / t_n_V; // Vp->V
	reulermultinom(2, Vp_vec[i], &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) from_Vp[i][j] = trans2[j]; 	
}

// Update state variables
for(i = 0; i < nA; i++) {
	
	// State variables
	S1_vec[i] += ((i == 0) ? birth_no : from_S1[i - 1][0]) - from_S1[i][0] - from_S1[i][1] - from_S1[i][2]; 
	E1_vec[i] += ((i == 0) ? 0 : from_E1[i - 1][0]) + from_S1[i][1] - from_E1[i][0] - from_E1[i][1];
	I1_vec[i] += ((i == 0) ? 0 : from_I1[i - 1][0]) + from_E1[i][1] - from_I1[i][0] - from_I1[i][1];
	S2_vec[i] += ((i == 0) ? 0 : from_S2[i - 1][0]) + from_R[i][1] + from_V[i][1] - from_S2[i][0] - from_S2[i][1] - from_S2[i][2]; 
	E2_vec[i] += ((i == 0) ? 0 : from_E2[i - 1][0]) + from_S2[i][1] - from_E2[i][0] - from_E2[i][1];
	I2_vec[i] += ((i == 0) ? 0 : from_I2[i - 1][0]) + from_E2[i][1] - from_I2[i][0] - from_I2[i][1];
	Rp1_vec[i] += ((i == 0) ? 0 : from_Rp1[i - 1][0]) + from_I1[i][1] + from_I2[i][1] - from_Rp1[i][0] - from_Rp1[i][1]; 
	R_vec[i] += ((i == 0) ? 0 : from_R[i - 1][0]) + from_Rp1[i][1] + from_Rp2[i][1] - from_R[i][0] - from_R[i][1] - from_R[i][2]; 
	Re_vec[i] += ((i == 0) ? 0 : from_Re[i - 1][0]) + from_R[i][2] - from_Re[i][0] - from_Re[i][1]; 
	Rp2_vec[i] += ((i == 0) ? 0 : from_Rp2[i - 1][0]) + from_Re[i][1] - from_Rp2[i][0] - from_Rp2[i][1]; 
	V_vec[i] += ((i == 0) ? 0 : from_V[i - 1][0]) + from_S1[i][2] + from_S2[i][2] + from_Vp[i][1] - from_V[i][0] - from_V[i][1] - from_V[i][2];
	Ve_vec[i] += ((i == 0) ? 0 : from_Ve[i - 1][0]) + from_V[i][2] - from_Ve[i][0] - from_Ve[i][1]; 
	Vp_vec[i] += ((i == 0) ? 0 : from_Vp[i - 1][0]) + from_Ve[i][1] - from_Vp[i][0] - from_Vp[i][1]; 
	
	// Incidence rates
	Ci1_vec[i] += from_I1[i][1];
	Ci2_vec[i] += from_I2[i][1];
	Cs_vec[i] += from_I1[i][1] + from_I2[i][1] + from_Re[i][1] + from_Ve[i][1];	
}

//end_stochSim


