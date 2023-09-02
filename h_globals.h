#ifndef H_GLOBALS_H_
#define H_GLOBALS_H_
// Hosna, Nov. 1st, 2011 changed the parameter values based on HOtKNots v2
//Luke Aug 2023, we need the gas constant for penalty calculation
//Luke Sep 2023 adding gas constant for scaling penalties
pf_t oneoverRT = -10.0/(1.98717*310.15);
double PS_penalty = exp(-138 * oneoverRT); //960; 		//exterior pseudoloop initiation penalty (9.6 Kcal/mol)
double PSM_penalty = exp(1007 * oneoverRT); //1500;		//penalty for introducing pseudoknot inside a multiloop (15 Kcal/mol)
double PSP_penalty = exp(1500 * oneoverRT);		//penalty for introducing pseudoknot inside a pseudoloop (15 Kcal/mol)
double PB_penalty = exp(246 * oneoverRT); //20;		//band penalty (0.2 Kcal/mol)
double PUP_penalty	= exp(6 * oneoverRT); //10;		//penalty for an un-paired base in a pseudoloop or a band (0.1 Kcal/mol)
double PPS_penalty = exp(96 * oneoverRT);//10;		//penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(0.1 Kcal/mol)

double e_stP_penalty = exp(0.89* oneoverRT); //0.83		// e_stP = 0.83 * e_s
double e_intP_penalty = exp(0.74* oneoverRT); //0.83;		// e_intP = 0.83 * e_int

// Hosna, Nov. 2nd, 2011
// changed the multiloop penalties to be the same as the simfold's values
double a_penalty = exp(339* oneoverRT); //The newest value is from DP09 //misc.multi_offset;//340;		//penalty for introducing a multiloop (3.4 Kcal/mol)
double b_penalty = exp(3* oneoverRT);  //The newest value is from DP09 //misc.multi_helix_penalty; //40;			//penalty for base pair in a multiloop (0.4 Kcal/mol)
double c_penalty = exp(2* oneoverRT);  //The newest value is from DP09 //misc.multi_free_base_penalty; //0;			//penalty for un-paired base in a multi-loop

double ap_penalty = exp(341* oneoverRT); //340;			//penalty for introducing a multiloop that spans a band (3.4 Kcal/mol)
double bp_penalty = exp(56* oneoverRT); //40;			//base pair penalty for a multiloop that spans a band (0.4 Kcal/mol)
double cp_penalty = exp(12* oneoverRT); //0;			//penalty for unpaired base in a multiloop that spans a band

// Hosna November 16, 2015
double start_hybrid_penalty = exp(310* oneoverRT);

#endif /*H_GLOBALS_H_*/
