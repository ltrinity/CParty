#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "VM_final.h"
#include "V_final.h"
#include "s_specific_functions.h"

// Ian Wark July 19 2017
// constant that defines what fres[i].pair will be compared against (>=) for impossible cases
// set to -1 because >= 0 means there is already a base pair there,
// and -1 means restricted struture says there is no base pair there.
#define FRES_RESTRICTED_MIN -1

pseudo_loop::pseudo_loop(char *seq, char* restricted, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM)
{
	this->sequence = seq;
	this->restricted = restricted;
	this->V = V;
	this->H = H;
	this->S = S;
	this->VBI = VBI;
	this->VM = VM;
    allocate_space();
    if (debug){
    	printf("an object of pseudo_loop was successfully created! \n");
    }
}

void pseudo_loop::allocate_space()
{
    int i;
    nb_nucleotides = strlen(sequence);
    needs_computation = 0; // Hosna, March 14, 2012 I need to remove this variable!! to make everything a lot faster

    index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;

    WI = new pf_t [total_length];
    if (WI == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WI[i] = 0; // if i == j -> p_up

    weakly_closed = new int[total_length];
    if (weakly_closed == NULL) giveup ("Cannot allocate memory", "weakly_closed");
    for (i=0; i < total_length; i++) weakly_closed[i] = 0;

    not_paired_all = new int[total_length];
    if (not_paired_all == NULL) giveup ("Cannot allocate memory", "not_paired_all");
    for (i=0; i < total_length; i++) not_paired_all[i] = 0;

    VP = new pf_t[total_length];
    if (VP == NULL) giveup ("Cannot allocate memory", "VP");
	//Luke Aug 2023 changing init to 0 from INF for partition function 
    for (i=0; i < total_length; i++) VP[i] = 0;

    WMB = new pf_t[total_length];
    if (WMB == NULL) giveup ("Cannot allocate memory", "WMB");
	//Luke Aug 2023 changing init to 0 from INF for partition function 
    for (i=0; i < total_length; i++) WMB[i] = 0;

    WMBP = new pf_t[total_length];
	if (WMBP == NULL) giveup("Cannot allocate memory","WMBP");
	//Luke Aug 2023 changing init to 0 from INF for partition function 
	for (i=0; i < total_length; i++) WMBP[i] = 0;

	PGPW = new pf_t[total_length];
	if (PGPW == NULL) giveup("Cannot allocate memory","PGPW");
	//Luke Aug 2023 changing init to 0 from INF for partition function 
	for (i=0; i < total_length; i++) PGPW[i] = 0;

    WIP = new pf_t[total_length];
    if (WIP == NULL) giveup ("Cannot allocate memory", "WIP");
	//Luke Aug 2023 changing init to 0 from INF for partition function 
    for (i=0; i < total_length; i++) WIP[i] = 0;


	// Luke added Aug 2023
	VPR = new pf_t[total_length];
    if (VPR == NULL) giveup ("Cannot allocate memory", "VPR");
    for (i=0; i < total_length; i++) VPR[i] = 0;
	
	VPL = new pf_t[total_length];
    if (VPL == NULL) giveup ("Cannot allocate memory", "VPL");
    for (i=0; i < total_length; i++) VPL[i] = 0;

    BE = new pf_t[total_length];
    if (BE == NULL) giveup ("Cannot allocate memory", "BE");
    for (i=0; i < total_length; i++) BE[i] = 0; //check

    border_bs = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bs[i] = new int[nb_nucleotides];

    border_bps = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bps[i] = new int[nb_nucleotides];


    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

pseudo_loop::~pseudo_loop()
{
    delete [] WI;
    delete [] WIP;
    delete [] VP;
	delete [] VPR;
	delete [] VPL;
    delete [] WMB;
    delete [] WMBP;
	delete [] PGPW;

    delete [] BE;
    delete [] weakly_closed;
    delete [] not_paired_all;


    // Ian Wark July 21 2017
    // border_bs is array of arrays
    // need to delete sub arrays as well
    for(int i = 0; i < nb_nucleotides; i++) {
        delete [] border_bs[i];
        delete [] border_bps[i];
    }

    delete [] border_bs;
    delete [] border_bps;

    delete [] index;
    delete [] int_sequence;
}

void pseudo_loop::set_features(h_str_features *f){
	fres = f;
}

int pseudo_loop::is_weakly_closed(int i, int j){
	// base case: if i > j then the region is weakly closed
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (weakly_closed[ij] == 1)
		return 1;
	return 0;
}


int pseudo_loop::is_empty_region(int i, int j){
	//base case: if i> j then the region is empty
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (not_paired_all[ij] == 1){
		return 1;
	}
	return 0;
}

void pseudo_loop::initialize(){

	int i, j;

    //Hosna: before going further, we should fill up the weakly closed array
    detect_weakly_closed(fres, weakly_closed, nb_nucleotides, index);
    detect_not_paired_all(fres, not_paired_all, nb_nucleotides, index);
    detect_border_bs(fres,border_bs, nb_nucleotides);
    detect_border_bps(fres,border_bps, nb_nucleotides);

    //testing to see if all the structures were found correctly
//    if (debug){
//	    for (i=0; i < nb_nucleotides; i++){
	        //if (fres[i].pair != -1 && fres[i].pair != -2)
//	        if(debug)
//	        printf ("%d pairs %d, is in arc %d, and has type %c\n", i, fres[i].pair, fres[i].arc, fres[i].type);

//	    }
//	    for (i = 0; i < nb_nucleotides; i++){
//	    	int j;
//	    	for (j = i; j < nb_nucleotides; j++){
//	    		int ij = index[i]+j -i;
//	    		if (weakly_closed[ij] == 1){
//		    		printf("region [%d,%d] is weakly closed. \n", i, j);
//	    		}
//	    		else{
//	    			printf("region [%d,%d] is NOT weakly closed. \n", i, j);
//	    		}
//				if (not_paired_all[ij] == 1){::compute_WI
//					printf("region [%d,%d] is EMPTY \n",i,j);
//				}
				//checking WI, VP and WMB to see if they have been initialized correctly
//				if (debug)
//				if(WI[ij] != 0){
//					printf("WI[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//				if(VP[ij] != 0){
//					printf("VP[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//				if(WMB[ij] != 0){
//					printf("WMB[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//	    	}
//	    }

//    }
//    printf("WMB was initialized successfully! \n");

}

void pseudo_loop::compute_energies(int i, int j)
{

	// Hosna, April 18th, 2007
	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)

	//	if(debug){
	//		printf("calculating VP(%d,%d) \n",i,j);
	//	}
	compute_VP(i,j,fres); // Hosna, March 14, 2012, changed the position of computing VP from after BE to befor WMBP

	compute_PGPW(i,j,fres);
//	if(debug){
//		printf("calculating WMBP(%d,%d) \n",i,j);
//	}::compute_WI
	compute_WMBP(i,j,fres);
//	if(debug){
//		printf("calculating WMB(%d,%d) \n",i,j);
//	}
    compute_WMB(i,j,fres);
//	if(debug){
//		printf("calculating WI(%d,%d) \n",i,j);
//	}
    compute_WI(i,j,fres);
//    if(debug){
//		printf("calculating WIP(%d,%d) \n",i,j);
//	}
    compute_WIP(i,j,fres);
//    if(debug){
//		printf("calculating VPP(%d,%d) \n",i,j);
//	}

    compute_VPR(i,j,fres);
//    if(debug){
//		printf("calculating VPL(%d,%d) \n",i,j);
//	}
    compute_VPL(i,j,fres);
//	if(debug){
//		printf("calculating BE(%d,%d) \n",i,j);
//	}

	compute_BE(fres[j].pair,j,fres[i].pair,i,fres);

//    if (debug){
//    	printf("WI(%d,%d) = %d \n", i,j, get_WI(i,j));
//    	printf("WIP(%d,%d) = %d \n", i,j, get_WIP(i,j));
//    	printf("VPP(%d,%d) = %d \n", i,j, get_VPP(i,j));
//    	printf("BE(%d,%d,%d,%d) = %d \n", i,fres[i].pair,j,fres[j].pair, get_BE(i,fres[i].pair,j,fres[j].pair));
//    	printf("VP(%d,%d) = %d \n", i,j, get_VP(i,j));
//    	printf("WMBP(%d,%d) = %d \n", i,j, get_WMBP(i,j));
//    	printf("WMB(%d,%d) = %d \n", i,j, get_WMB(i,j));
//    }
}

void pseudo_loop::compute_WI(int i, int j , h_str_features *fres){
	pf_t min = INF, m1 = INF, m2= INF, m3= INF;

	// Luke init partition function Aug 2023
	pf_t d2_energy_wi = 0;
	int ij = index[i]+j-i;
	int ijminus1 = index[i] + (j-1)-i;
	if (WI[ij] != 0){ //calculated before
//		if (debug)
//		{
//			printf("WI(%d,%d) was calculated before ==> WI(%d,%d) = %d \n",i,j,i,j,get_WI(i,j));
//		}
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (is_weakly_closed(i,j) == 0){
		//Luke changing to 0 from INF for part func Aug 2023
		WI[ij] = 0;
//		if (debug)
//		{
//			printf("[%d,%d] is not weakly closed ==> WI(%d,%d) = %d \n",i,j,i,j,get_WI(i,j));
//		}
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
//		if (debug){
//			printf("i ==j => WI[%d,%d]= %d \n",i,j,WI[ij]);
//		}
		return;
	}
	// Hosna: Feb 12, 2007
	// changed this part to see if it works better

	// Hosna: Feb 16, 2007:
	// we don't need to check to see if i and j are inside an arc
	// because they are not in an arc in G but they will be in an arc in G'
	//Luke changing to 0 from INF for part func Aug 2023
	if (fres[i].arc != fres[j].arc){
		WI[ij] = 0;
//		if (debug){
//			printf("i and j not in the same arc => WI[%d,%d]= %d \n",i,j,WI[ij]);
//		}
		return;
	}

// Hosna: July 2nd, 2007
// in branch 1 of WI, we can have a case like
// ((..))((...))
// such that both i and j are paired but we can chop them

	// branch 1:
//	if (fres[i].pair < 0 && fres[j].pair < 0)
//	{
	// Luke 6/27/2023: ambiguity here
	int t;
	for (t = i; t< j; t++){
		pf_t wi_1 = get_WI(i,t-1);
		//new case 1, check if V possible and add to WI(i, t-1)
		if ((fres[t].pair == j && fres[j].pair == t)
		||(fres[t].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN)){
			pf_t v_ener = V->get_energy(t,j);
			pf_t energy = wi_1 * v_ener * PPS_penalty;
			d2_energy_wi += energy;
			//m1 = (m1 > energy)? energy : m1;
		}
		//new case 3 
		pf_t energy2 = wi_1 * get_WMB(t,j) * PSP_penalty * PPS_penalty;
		d2_energy_wi += energy2;
		//m3 = (m3 > energy2)? energy2 : m3;
//		if (debug_WI){
//			printf("WI branch 1: WI[%d,%d] = %d and WI[%d,%d] = %d => energy = %d and m1 = %d \n",i,t,wi_1,(t+1),j,wi_2,energy, m1);
//		}
	}
	//new case 2
	d2_energy_wi += (WI[ijminus1] * PUP_penalty);
	//partition function
	WI[ij] = d2_energy_wi;
	//printf("Z_WI(%d,%d) = %Lf \n",i,j,d2_energy_wi);
}

// Luke modifying for CParty, cases 1-5 unchanged, new cases 6-9
void pseudo_loop::compute_VP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	//Luke Aug 2023 changing to 0 from INF
	if (VP[ij] != 0){//has been calculated before
//		if (debug)
//		{
//			printf("VP(%d,%d) was calculated before ==> VP(%d,%d) = %d \n",i,j,i,j,VP[ij]);
//		}
		return;
	}
	// base cases:
	// a) i == j => VP[ij] = INF
	// b) [i,j] is a weakly_closed region => VP[ij] = INF
	// c) i or j is paired in original structure => VP[ij] = INF

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says it should not be modified
	if (i == j || weakly_closed[ij] == 1 || fres[i].pair >= FRES_RESTRICTED_MIN || fres[j].pair >= FRES_RESTRICTED_MIN || can_pair(int_sequence[i],int_sequence[j]) != 1)	{
		//Luke changing to 0
		VP[ij] = 0;
//		if (debug){
//			printf("VP[%d,%d] = %d and can_pair(%d,%d) = %d\n", i,j, VP[ij],int_sequence[i],int_sequence[j],can_pair(int_sequence[i],int_sequence[j]));
//		}

		return;
	}
	else{
		pf_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, m7 = INF, m8= INF, m9 = INF; //different branches
		// Luke init partition function Aug 2023
		pf_t d2_energy_vp = 0;
		//branchs:
		// 1) inArc(i) and NOT_inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// need to check the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			pf_t WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1);
			pf_t WI_Bplus_jminus = get_WI(B_i + 1,j-1);
			// Luke Aug 2023 changing for part func.
			m1 =   WI_ipus1_BPminus * WI_Bplus_jminus;
			d2_energy_vp +=m1;
//			if(debug){
//				printf("VP[%d,%d] branch 1: WI(%d+1)(BP(%d)-1) = %d and WI(B(%d)+1)(%d-1) = %d => m1 = %d \n",i,j,i,i,WI_ipus1_BPminus,i,j,WI_Bplus_jminus, m1);
//			}
		}

		// 2) NOT_inArc(i) and inArc(j)
		// WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
			int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
			// Luke Aug 2023 changing for part func.
			m2 = WI_i_plus_b_minus * WI_bp_plus_j_minus;
			d2_energy_vp +=m2;
//			if(debug){
//				printf("VP[%d,%d] branch 2: WI(%d+1)(b(%d)-1) = %d and WI(bp(%d)+1)(%d-1) = %d => m2 = %d \n",i,j,i,i,WI_i_plus_b_minus,i,j,WI_bp_plus_j_minus, m2);
//			}
		}

		// 3) inArc(i) and inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			pf_t WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
			pf_t WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
			pf_t WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
			// Luke Aug 2023 changing for part func.
			m3 = WI_i_plus_Bp_minus * WI_B_plus_b_minus * WI_bp_plus_j_minus;
			d2_energy_vp +=m3;
//			if(debug){
//				printf("VP[%d,%d] branch 3: WI(%d+1)(B'(%d)-1) = %d, WI(B(%d)+1)(b(%d)-1) = %d and WI(b'(%d)+1)(%d-1) = %d => m3 = %d \n",i,j,i,i, WI_i_plus_Bp_minus,i,i,WI_B_plus_b_minus,i,j,WI_bp_plus_j_minus, m3);
//			}
		}

        // Ian Wark July 19 2017
        // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs


		// 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
		// e_stP(i,i+1,j-1,j) + VP(i+1)(j-1)
		if(fres[i+1].pair < FRES_RESTRICTED_MIN && fres[j-1].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i+1],int_sequence[j-1])){
			// Luke Aug 2023 changing for part func.
			m4 = get_e_stP(i,j) * get_VP(i+1,j-1);
			d2_energy_vp += m4;
//			if (debug){
//				printf("VP[%d,%d] branch 4: S->get_energy(%d,%d) = %d and VP[%d,%d] = %d  so m4 = %d\n", i,j,i,j, S->get_energy(i,j,int_sequence), i+1, j-1, get_VP(i+1,j-1), m4);
//			}
		}

		// 5) NOT_paired(r) and NOT_paired(rp)
		//  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
		int ip, jp;
		int max_borders;
		// Hosna, April 6th, 2007
		// whenever we use get_borders we have to check for the correct values
		int min_borders = 0; // what if both are negative
		if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
			min_borders = MIN(get_Bp(i,j),get_b(i,j));
		}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
			min_borders = get_b(i,j);
		}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
			min_borders = get_Bp(i,j);
		}
//		printf("B'(%d,%d) = %d, b(%d,%d) = %d, min_borders = %d\n",i,j,get_Bp(i,j),i,j,get_b(i,j), min_borders);
		for (ip = i+1; ip < min_borders; ip++){
			// Hosna: April 20, 2007
			// i and ip and j and jp should be in the same arc
			// also it should be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions

			// Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
			if (fres[ip].pair < FRES_RESTRICTED_MIN && (fres[i].arc == fres[ip].arc) && is_empty_region(i+1,ip-1) == 1){
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values
				max_borders= 0;
				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
					max_borders = MAX(get_bp(i,j),get_B(i,j));
				}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
					max_borders = get_B(i,j);
				}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
					max_borders = get_bp(i,j);
				}
//				printf("b'(%d,%d) = %d, B(%d,%d) = %d, max_borders = %d\n",i,j,get_bp(i,j),i,j,get_B(i,j), max_borders);
				for (jp = max_borders+1; jp < j ; jp++){
                    // Ian Wark July 19 2017
                    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
                    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
					if (fres[jp].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1) == 1){
						// Hosna: April 20, 2007
						// i and ip and j and jp should be in the same arc
						if (fres[j].arc == fres[jp].arc ){
							// Luke Aug 2023 changing for part func.
							pf_t temp = get_e_intP(i,ip,jp,j) * get_VP(ip,jp);
							d2_energy_vp += temp;
//							if (debug){
//								printf("VP(%d,%d) branch 5: e_intP(%d,%d,%d,%d) = %d, VP(%d,%d) = %d, temp = %d \n",i,j,i,ip,jp,j,get_e_intP(i,ip,jp,j),ip,jp,get_VP(ip,jp),temp);
//							}
//							printf("m5 = %d \n", m5);
							//depracated
							//if (m5 > temp){
							//	m5 = temp;
							//}
						}
					}
				}
			}
		}
		// Luke Aug 2023 case 6
		// 6) VP(i,j) = WIP(i+1,r-1) + VP(r,j-1)
		int r;
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int min_Bp_j = j;
		if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
			min_Bp_j = get_Bp(i,j);
		}
		for (r = i+1; r < min_Bp_j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				//Luke Aug 2023 : change here to get value from VP instead of VPP for part func.
				pf_t tmp = get_WIP(i+1,r-1) * get_VP(r,j-1) * (ap_penalty + (2*bp_penalty));
				d2_energy_vp += tmp;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 6: WIP(%d,%d) = %d, VPP(%d,%d) = %d ==> tmp = %d and m6 = %d \n",i,j,i+1,r-1,get_WIP(i+1,r-1),r,j-1,get_VPP(r,j-1),tmp,m6);
//				}
				//depracated
				//if (tmp < m6){
				//	m6 = tmp;
				//}
			}
		}

		// Luke Aug 2023 Case 7
		// 7) VP(i,j) = VP(i+1,r) + WIP(r+1,j-1)
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int max_i_bp = i;
		if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
			max_i_bp = get_bp(i,j);
		}
		for (r = max_i_bp+1; r < j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				// Luke Aug 2023: change here from VPP to VP
				pf_t tmp = get_VP(i+1,r) * get_WIP(r+1,j-1) * (ap_penalty + (2* bp_penalty));
				d2_energy_vp += tmp;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 7: VPP(%d,%d) = %d, WIP(%d,%d) = %d ==> tmp = %d and m7 = %d \n",i,j,i+1,r,get_VPP(i+1,r),r+1,j-1,get_WIP(r+1,j-1),tmp,m7);
//				}
				//depracated
				//if (tmp < m7){
				//	m7 = tmp;
				//}
			}
		}

		// Luke Aug 2023 case 8
		// 8) VP(i,j) = WIP(i+1,r-1) + VPR(r,j-1)
		for (r = i+1; r < min_Bp_j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				//Luke Aug 2023: change here to get value from VPR instead of VP
				pf_t tmp = get_WIP(i+1,r-1) * get_VPR(r,j-1) * (ap_penalty + (2*bp_penalty));
				d2_energy_vp += tmp;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 6: WIP(%d,%d) = %d, VPP(%d,%d) = %d ==> tmp = %d and m6 = %d \n",i,j,i+1,r-1,get_WIP(i+1,r-1),r,j-1,get_VPP(r,j-1),tmp,m6);
//				}
				//depracated
				//if (tmp < m8){
				//	m8 = tmp;
				//}
			}
		}

		// Luke Aug 2023 Case 9
		// 7) VP(i,j) = VPL(i+1,r) + WIP(r+1,j-1)
		for (r = max_i_bp+1; r < j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				//Luke Aug 2023: change here from VP to VPL
				pf_t tmp = get_VPL(i+1,r) * get_WIP(r+1,j-1) * (ap_penalty + (2* bp_penalty));
				d2_energy_vp += tmp;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 7: VPP(%d,%d) = %d, WIP(%d,%d) = %d ==> tmp = %d and m7 = %d \n",i,j,i+1,r,get_VPP(i+1,r),r+1,j-1,get_WIP(r+1,j-1),tmp,m7);
//				}
				//depracated
				//if (tmp < m9){
				//	m9 = tmp;
				//}
			}
		}
		VP[ij] = d2_energy_vp;
		//printf("Z_VP(%d,%d) = %Lf \n",i,j,d2_energy_vp);
	}
}

// Luke: adding PG'W
void pseudo_loop::compute_PGPW(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	// Luke Aug 2023 base case if already computed
	if (PGPW[ij] != 0){
		return;
	}
	// Luke Aug 2023 base case partition function 0
	if (i == j){
		PGPW[ij] = 0;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	||  (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		// Luke Aug 2023 changing to 0
		WMB[ij] = 0;
		return;
	}
	else{
		// Luke init partition function Aug 2023
		pf_t d2_energy_pgpw = 0;
		//changing to 0
		pf_t m1 = 0;
		// Luke: new case 1 for PGPW
		// 5) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
		// Hosna: Feb 5, 2007
		if(fres[j].pair < j){
			int l,l_min =-1;
			for(l = i+1; l<j; l++){
				// Hosna: March 14th, 2007
				// I think l cannot be paired

				// Hosna: April 18th, 2007
				// l and j should be in the same arc
				if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[j].arc == fres[l].arc){
					//Luke Aug 2023 part func
					pf_t temp = get_WMBP(i,l) * get_WI(l+1,j);
					//printf("Z_wmbp(%d,%d) = %d \n",i,j,get_WMBP(i,l));
					//printf("Z_wi(%d,%d) = %d \n",i,j,get_WI(l+1,j));
					d2_energy_pgpw += temp;
					//depracated
					//if (temp < m1){
					//	m1 = temp;
					//	l_min = l;
					//}

//					if (debug_WMB){
//						printf("***************\n");
//						printf("WMB(%d,%d) inside branch 5: l = %d, WMB = %d, WI = %d \n",i,j,l,get_WMB(i,l),get_WI(l+1,j));
//						printf("***************\n");
//					}
				}
			}
//			if (debug ){
//				printf("WMB(%d,%d) branch 5:  l = %d So m5 = %d\n",i,j,l_min, m5);
//			}
		}

		// get the min for WMB Aug 2023 Luke modified
		PGPW[ij] = d2_energy_pgpw;
		//printf("Z_PGPw(%d,%d) = %Lf \n",i,j,d2_energy_pgpw);
//		if (debug && i == 1 && j == 87){
//			printf("m1 = %d, m3 = %d, m4 = %d and m5 = %d ==> WMBP[%d,%d] = %d\n",m1,m3,m4,m5,i,j,WMBP[ij]);
//		}
	}


}

// Luke: modifying to follow CParty scheme
void pseudo_loop::compute_WMBP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WMBP[ij] != 0){
		return;
	}
	//base case Luke Aug 2023 set to 0 from INF
	if (i == j){
		WMBP[ij] = 0;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	||  (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		//base case Luke Aug 2023 set to 0 from INF
		WMB[ij] = 0;
		return;
	}
	else{
		pf_t m1 = INF, m3 = INF, m4 = INF, m5 = INF;
		// Luke init partition function Aug 2023
		pf_t d2_energy_wmbp = 0;
		// if not paired(j) and paired(i) then
		// WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
		if(fres[j].pair < 0 && fres[i].pair >= 0){
			int tmp = INF, l, l_min=-1;
			// Hosna: June 29, 2007
			// if j is inside i's arc then the l should be
			// less than j not bp(i)
			// check with Anne
//			for (l = i+1; l < MIN(fres[i].pair,j); l++){
			// Hosna: July 5th, 2007:
			// if we have bp(i)> j then we should not have come to the WMBP
			for (l = i+1; l < j; l++){
				// Hosna, March 14, 2007
				// fixing the for loop

				// Hosna, April 9th, 2007
				// checking the borders as they may be negative
//				if(fres[l].pair < 0 && get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
				// Hosna: July 5th, 2007:
				// removed bp(l)<0 as VP should handle that
				if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
					int bp_i_l = get_bp(i,l);
					pf_t BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
					pf_t WI_energy = get_WI(bp_i_l +1,l-1);
					pf_t VP_energy = get_VP(l,j);
					pf_t sum = BE_energy * WI_energy * VP_energy * (2*PB_penalty);
					//Luke Aug 2023 part func
					d2_energy_wmbp += sum;
					//depracated
					//if (tmp > sum){
					//	tmp = sum;
					//	l_min = l;
					//}
//					if (debug  && i == 6 && bp_i_l == 11 ){
//						printf("***************\n");
//						printf("WMBP(%d,%d) inside branch 1: %d is paired with %d  and b'(%d,%d) = %d  and bp(b')= %d \n",i,j, i,fres[i].pair, i,l,bp_i_l, fres[bp_i_l].pair);
//						printf("l = %d, BE(%d,%d,%d,%d) = %d, WI = %d, VP = %d, --> sum = %d \n",l,i,fres[i].pair,bp_i_l,fres[bp_i_l].pair,BE_energy,WI_energy,VP_energy, sum);
//						printf("***************\n");
//					}
				}
			}
			//depracated
			//m1 = 2*PB_penalty + tmp;
//			if (debug ){
//				printf("WMBP(%d,%d) branch 1:  l = %d  ==> m1 = %d \n",i,j,l_min, m1);
//			}

		}

		// 3)
		if (fres[j].pair < 0){
			int l, temp = INF, l_min=-1;
			for (l = i+1; l<j ; l++)	{
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values

				if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					// Hosna: April 19th, 2007
					// the chosen l should be less than border_b(i,j)
					if (get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && l < get_b(i,j)){

						// Hosna: June 29 2007
						// after going over the program with Cristina, we noticed that
						// l should be < B'(i,j)
	//					if (l < get_Bp(i,j) && l+TURN <= j){

						// Hosna: July 5th, 2007:
						// as long as we have i <= arc(l)< j we are fine
						if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
							pf_t sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))* get_WMBP(i,l-1)* get_VP(l,j) *(2*PB_penalty);
							//Luke adding for PGPW case
							pf_t sum2 = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))* get_PGPW(i,l-1)* get_VP(l,j)*(2*PB_penalty);

							//Luke Aug 2023 part func
							d2_energy_wmbp += sum;
							d2_energy_wmbp += sum2;
							//depracated
							//if (temp > sum){
							//	temp = sum;
							//	l_min = l;
							//}
	//						if (debug && fres[get_B(l,j)].pair == 6 && fres[get_Bp(l,j)].pair == 11){
	//							printf("***************\n");
	//							printf("WMBP(%d,%d) inside branch 3: l = %d, BE = %d, WMBP = %d, VP = %d ~~> sum = %d \n",i,j,l,get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j)),get_WMBP(i,l-1),get_VP(l,j), sum);
	//							printf("***************\n");
	//						}
						}
					}
				}
				// Hosna: April 5th
				// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
				// to the 3rd case ==> 2*P_b
				//m3 = 2*PB_penalty + temp;
	//			if (debug){
	//				printf("WMBP(%d,%d) branch 3:  l = %d So m3 = %d\n",i,j,l_min, m3);
	//			}
			}
		}
		//Luke Aug 2023 part func
		// 4) WMB(i,j) = VP(i,j) + P_b
		pf_t temp = get_VP(i,j) * PB_penalty;
		d2_energy_wmbp += temp;
		//depracated
		//if (temp < m4){
		//	m4 = temp;
		//}
//		if (debug){
//			printf("WMBP(%d,%d) branch 4:  VP = %d So m4 = %d\n",i,j,get_VP(i,j), m4);
//		}
		// Luke: removing, now in PGPW
		// 5) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
		// Hosna: Feb 5, 2007
		// if(fres[j].pair < j){
		// 	int l,l_min =-1;
		// 	for(l = i+1; l<j; l++){
		// 		// Hosna: March 14th, 2007
		// 		// I think l cannot be paired

		// 		// Hosna: April 18th, 2007
		// 		// l and j should be in the same arc
		// 		if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[j].arc == fres[l].arc){
		// 			int temp = get_WMBP(i,l) + get_WI(l+1,j);
		// 			if (temp < m5){
		// 				m5 = temp;
		// 				l_min = l;
		// 			}

//					if (debug_WMB){
//						printf("***************\n");
//						printf("WMB(%d,%d) inside branch 5: l = %d, WMB = %d, WI = %d \n",i,j,l,get_WMB(i,l),get_WI(l+1,j));
//						printf("***************\n");
//					}
			// 	}
			// }
//			if (debug ){
//				printf("WMB(%d,%d) branch 5:  l = %d So m5 = %d\n",i,j,l_min, m5);
//			}
		// }
		
		WMBP[ij] = d2_energy_wmbp;
		printf("Z_wmbp(%d,%d) = %Lf = %Lf \n",i,j,d2_energy_wmbp, WMBP[ij]);
		//depracated
		// get the min for WMB
		//WMBP[ij] = MIN(MIN(m1,m3),MIN(m4,m5));
//		if (debug && i == 1 && j == 87){
//			printf("m1 = %d, m3 = %d, m4 = %d and m5 = %d ==> WMBP[%d,%d] = %d\n",m1,m3,m4,m5,i,j,WMBP[ij]);
//		}
	}


}

void pseudo_loop::compute_WMB(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	//Luke base case changes Aug 2023
	if (WMB[ij] != 0){
		return;
	}
	if (i == j){
		WMB[ij] = 0;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

    // Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	 || (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		// Luke Aug 2023 base case change to 0
		WMB[ij] = 0;
		return;
	}
	else{
		pf_t m2 = INF, mWMBP = INF;
		// Luke Aug 2023 init part func
		pf_t d2_energy_p = 0;
		// 2)
		if (fres[j].pair >= 0 && j > fres[j].pair){
			int l, l_min=-1;
			int bp_j = fres[j].pair;
			//Luke change to 0
			pf_t temp = 0;
//			if (debug_WMB){
//				printf("\n INSIDE WMB BRANCH 2 \n where bp_j = %d and j = %d \n\n", bp_j,j);
//			}
			for (l = (bp_j +1); (l < j); l++){
				// Hosna: April 24, 2007
				// correct case 2 such that a multi-pseudoknotted
				// loop would not be treated as case 2

				// Hosna: July 5th, 2007
				// this restriction was removed as it is not needed here
//				if (l > fres[i].pair){
	//				printf("l = %d and bp_l = %d \n",l,fres[l].pair);
					// Hosna April 9th,
					// checking the borders as they may be negative numbers

					// Hosna: July 5th, 2007:
					// we don't need to check that l is unpaired here
//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
						//Luke Aug 2023 part func
						pf_t sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) * get_WMBP(i,l) * get_WI(l+1,get_Bp(l,j)-1) * PB_penalty;
						//printf("Z_BE(%d,%d) = %d \n",i,j,get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)));
						d2_energy_p += sum;
						//printf("Z_P(%d,%d) = %Lf \n",i,j,d2_energy_p);
						if (l == 600 && i == 522 && j == 615) {
							int t = 0;
						}
						//depracated
						//if (temp > sum){
						//	temp = sum;
						//	l_min = l;
						//}
//						if (debug && bp_j == 6 && fres[get_Bp(l,j)].pair == 11){
//							printf("***************\n");
//							printf("WMB(%d,%d) inside branch 2: l = %d, BE(%d,%d) = %d, WMBP = %d, WI = %d ==> sum = %d while temp = %d\n",i,j,l,bp_j,j,get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)),get_WMBP(i,l),get_WI(l+1,get_Bp(l,j)-1), sum, temp);
//							printf("***************\n");
//						}
					}
//				}

			}
			//depracated
			//m2 = PB_penalty + temp;
//			if (debug){
//				printf("WMB(%d,%d) branch 2:  l = %d So m2 = %d\n",i,j,l_min, m2);
//			}
		}
		// check the WMBP value
		mWMBP =  get_WMBP(i,j);
		printf("Z_PGP(%d,%d) = %Lf \n",i,j,mWMBP);
		//Luke Aug 2023 part func
		d2_energy_p += mWMBP;
		WMB[ij] = d2_energy_p;
		// depracated
		// get the min for WMB
		//WMB[ij] = MIN(m2,mWMBP);
//		if (debug && i == 1 && j == 87){
//			printf("m2 = %d, mWMBP = %d ==> WMB[%d,%d] = %d\n",m2,mWMBP,i,j,WMB[ij]);
//		}
	}
}

// Luke 6/28/2023 removing ambiguity
void pseudo_loop::compute_WIP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WIP[ij] !=0){ // was calculated before Luke Aug 2023 change to 0 prev INF/2
		return;
	}
//	if (debug){
//		if (i >= 27 && i <= 28 && j <= 68 && j >= 67){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when arc(%d) = %d and arc(%d) = %d and weakly_closed(%d,%d) = %d \n",i,j,i,fres[i].arc,j,fres[j].arc,i,j,weakly_closed[ij]);
//		}
//	}
	if (fres[i].arc != fres[j].arc || i == j || weakly_closed[ij]== 0){
		//Luke Aug 2023 base case change
		WIP[ij] = 0;
		return;
	}
	pf_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// Luke Aug 2023 init part func
	pf_t d2_energy_wip = 0;
    // Ian Wark July 19 2017
	// fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs

	// branch 1 removed :
	//if (fres[i].pair < FRES_RESTRICTED_MIN){
	//	m1 = get_WIP(i+1,j) + cp_penalty;
//		if (debug && (i == 27 || i == 28) && (j == 68 || j == 67)){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when  WIP(%d,%d) = %d and m1 = %d\n",i,j,i+1,j,get_WIP(i+1,j),m1);
//		}
	//}
	// new branch 2 (unchanged):
	if (fres[j].pair < FRES_RESTRICTED_MIN){
		m2 = get_WIP(i,j-1) * cp_penalty;
		d2_energy_wip += m2;
//		if (debug && (i == 27 || i == 28) && (j == 68 || j == 67)){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when  WIP(%d,%d) = %d and m2 = %d\n",i,j,i,j-1,get_WIP(i,j-1),m2);
//		}
	}
	//previous branch 3 removed:
	// now new branches 1 and 3
	int t;
	for (t = i; t <j; t++){
		pf_t wip_1 = get_WIP(i,t-1);
		if (fres[t].pair == j
		|| (fres[t].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[t],int_sequence[j]))){
			//Aug 2023 change Luke
			pf_t v_ener = V->get_energy(i,j)	* bp_penalty;
			pf_t energy = wip_1*v_ener;
			d2_energy_wip += energy;
			//depracated
			//m1 = (m1 > energy)? energy : m1;
		}
		pf_t energy2 = wip_1 * get_WMB(t,j) * PSM_penalty * bp_penalty;
		d2_energy_wip += energy2;
		//m3 = (m3 > energy2)? energy2 : m3;
//		if (debug && i ==15 && j == 20 ){
//			printf("*************************************\n");
//			printf("WIP(%d,%d) branch 4: V(%d,%d) = %d => m4 = %d \n",i,j,i,j,V->get_energy(i,j),m4);
//			printf("*************************************\n");
//		}

	}
	//	if (tmp < m3){
	//		m3 = tmp;
	//	}
	//}

	// branch 4 (removed):
	//if (fres[i].pair == j
	//|| (fres[i].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i],int_sequence[j]))){

		//m4 = V->get_energy(i,j)	+ bp_penalty;
//		if (debug && i ==15 && j == 20 ){
//			printf("*************************************\n");
//			printf("WIP(%d,%d) branch 4: V(%d,%d) = %d => m4 = %d \n",i,j,i,j,V->get_energy(i,j),m4);
//			printf("*************************************\n");
//		}

	//}

	// branch 5:
	//m5 = get_WMB(i,j) + PSM_penalty + bp_penalty;
//	if (debug && i == 6 && j == 15){
//		printf("WIP(6,15) is calling WMB \n");
//	}
	//Luke Aug 2023 part func
	WIP[ij] = d2_energy_wip; // depracatedMIN(m1,MIN(m2,m3));

//	if (debug){
//		printf("WIP(%d,%d): m1 = %d, m2 = %d, m3 = %d, m4 = %d, m5 = %d ==> min = %d \n",i,j,m1,m2,m3,m4,m5,WIP[ij]);
//	}
}

// Luke adding VPR, case 1 VP + WI' and case 2 VP (unpaired bases 3')
void pseudo_loop::compute_VPR(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (VPR[ij] != 0){ // computed before
//		if (debug){
//			printf("VPR(%d,%d) was calculated before ==> VPR(%d,%d)= %d \n",i,j,i,j,get_VPR(i,j));
//		}
		return;
	}
	if (i == j  || this->is_weakly_closed(i,j)){
		VPR[ij] = 0;
//		if (debug){
//			printf("VPR(%d,%d): i == j ==> VPR = %d \n",i,j,VPR[ij]);
//		}
		return;
	}
	pf_t m1 = INF, m2 = INF;
	// Luke Aug 2023 init part func
	pf_t d2_energy_vpr = 0;

	//case 1:
	int r=-1 ;
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int max_i_bp = i;
	if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
		max_i_bp = get_bp(i,j);
	}
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN){
			pf_t tmp = get_VP(i,r) * get_WIP(r+1,j);
			d2_energy_vpr+=tmp;
//			if (debug){
//				printf("VPP(%d,%d) branch 1: VP(%d,%d) = %d, WIP(%d,%d)= %d ==> tmp = %d  and m1 = %d\n",i,j,i,r,get_VP(i,r),r+1,j,get_WIP(r+1,j),tmp, m1);
//			}
		}
	}
	// case 2:
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(r+1,j)){
			pf_t tmp = get_VP(i,r) * (cp_penalty *(j-r)); // check the (j-r) part
			d2_energy_vpr+=tmp;
//			if (debug){
//				printf("VPR(%d,%d) branch 2: VP(%d,%d) = %d, %d *(%d-%d)= %d ==> tmp = %d  and m3 = %d\n",i,j,i,r,get_VP(i,r),cp_penalty,j,r,cp_penalty *(j-r),tmp, m3);
//			}
			//if (tmp < m2){
			//	m2 = tmp;
			//}
		}
	}

	//int min_branches = m1;
	//if (m2 < min_branches){
	//	min_branches = m2;
	//}
	VPR[ij] = d2_energy_vpr; //MIN(MIN(m1,m2));
//	if (debug){
//		printf("VPR(%d,%d): m1 = %d, m2 = %d, m3 = %d and m4 = %d ==> min = %d \n", i,j,m1,m2,m3,m4,VPP[ij]);
//	}
	return;
}

// Luke adding VPL, case 1 VP (unpaired bases 5')
void pseudo_loop::compute_VPL(int i, int j, h_str_features *fres){
		int ij = index[i]+j-i;
	if (VPL[ij] != 0){ // computed before
//		if (debug){
//			printf("VPL(%d,%d) was calculated before ==> VPP(%d,%d)= %d \n",i,j,i,j,get_VPP(i,j));
//		}
		return;
	}
	if (i == j  || this->is_weakly_closed(i,j)){
		VPL[ij] = 0;
//		if (debug){
//			printf("VPP(%d,%d): i == j ==> VPP = %d \n",i,j,VPP[ij]);
//		}
		return;
	}
	pf_t m1 = INF;
	// Luke Aug 2023 init part func
	pf_t d2_energy_vpl = 0;

	//branch 1:
	int r=-1 ;
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int min_Bp_j = j;
	if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
		min_Bp_j = get_Bp(i,j);
	}
	for (r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(i,r-1)){
			pf_t tmp = (cp_penalty * (r-i)) * get_VP(r,j);
			d2_energy_vpl += tmp;
//			if (debug){
//				printf("VPP(%d,%d) branch 4: %d *(%d-%d) = %d, VP(%d,%d)= %d ==> tmp = %d  and m4 = %d\n",i,j,cp_penalty,r,i,cp_penalty * (r-i),r,j,get_VP(r,j),tmp, m4);
//			}
			//if (tmp < m1){
			//	m1 = tmp;
			//}
		}
	}
	
	VPL[ij] = d2_energy_vpl; //MIN(MIN(m1,m2),MIN(m3,m4));
//	if (debug){
//		printf("VPL(%d,%d): m1 = %d, m2 = %d, m3 = %d and m4 = %d ==> min = %d \n", i,j,m1,m2,m3,m4,VPP[ij]);
//	}
	return;
}
void pseudo_loop::compute_BE(int i, int j, int ip, int jp, h_str_features * fres){

//	if (debug && i == 6 && ip == 11){
//		printf("coming to BE to calculate BE(6,69,11,24) \n");
//	}

    // Ian Wark July 19 2017
    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >= FRES_RESTRICTED_MIN && fres[j].pair >= FRES_RESTRICTED_MIN && fres[ip].pair >= FRES_RESTRICTED_MIN && fres[jp].pair >= FRES_RESTRICTED_MIN && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
//		if (debug && i == 1 && ip == 11){
//			printf("BE(%d,%d,%d,%d): Impossible case!! \n",i,j,ip,jp);
//		}
		return;
	}
	int iip = index[i]+ip-i;
	if (BE[iip] != 0){ // computed before
//		if (debug && i == 6 && ip == 11){
//			printf("BE(%d,%d,%d,%d) was calculated before ==> BE=%d\n",i,j,ip,jp,BE[iip]);
//		}
		return;
	}
	// base case: i.j and ip.jp must be in G
	if (fres[i].pair != j || fres[ip].pair != jp){
//		if (debug ){
//			printf("BE(%d,%d,%d,%d) = INF \n",i,j,ip,jp);
//		}
//	Luke Aug 2023 change to 0
		BE[iip] = 0;
		return;
	}

	// base case:
	if(i == ip && j == jp && i<j){
//		if (debug ){
//			printf("BE(%d,%d,%d,%d) = 0 \n",i,j,ip,jp);
//		}
		BE[iip] = 0;
		return;
	}

	pf_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// Luke Aug 2023 init part func
	pf_t d2_energy_be = 0;
	// 1) bp(i+1) == j-1
	if (fres[i+1].pair == j-1){
		m1 = get_e_stP(i,j) * get_BE(i+1,j-1,ip,jp);
		d2_energy_be += m1;
//		if(debug ){
//			printf("BE(%d,%d,%d,%d) Case 1: e_stP(%d,%d) = %d and BE(%d,%d,%d,%d) = %d ==> m1 = %d \n",i,j,ip,jp,i,j,get_e_stP(i,j),i+1,j-1,ip,jp,get_BE(i+1,j-1,ip,jp),m1);
//		}
	}

	// cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
	int l;
	for (l = i+1; l<= ip ; l++){
        // Ian Wark July 19 2017
        // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

		// Hosna: March 14th, 2007
		if (fres[l].pair >= FRES_RESTRICTED_MIN && jp <= fres[l].pair && fres[l].pair < j){
			// Hosna, March 15, 2007
			// since not_paired_all[i,l] includes i and l themselves
			// and in BE energy calculation we are looking for the oepn region (i,l)
			// we have to look at not_paired_all[i+1,l-1]
			int lp = fres[l].pair;
			int il = index[i]+l-i;
			int lpj = index[lp]+j-lp;
			// 2)
			// Hosna June 29, 2007
			// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
			// so I am checking explicitely that we won't have stems instead of internal loop
			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){//&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
				pf_t temp = get_e_intP(i,l,lp,j)* get_BE(l,lp,ip,jp);
				d2_energy_be += temp;
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 2: e_intP(%d,%d,%d,%d) = %d, BE(%d,%d,%d,%d)= %d ==> temp = %d  and m2 = %d\n",i,j,ip,jp,i,l,lp,j,get_e_intP(i,l,lp,j),l,lp,ip,jp,get_BE(l,lp,ip,jp));
//				}
				//if (m2 > temp){
				//	m2 = temp;
				//}
			}

			// 3)
			if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				pf_t temp = get_WIP(i+1,l-1) * get_BE(l,lp,ip,jp) * get_WIP(lp+1,j-1) * (ap_penalty + (2* bp_penalty));
				d2_energy_be += temp;
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 3: WIP(%d,%d) = %d, BE(%d,%d,%d,%d)= %d, WIP(%d,%d)= %d ==> temp = %d  and m3 = %d\n",i,j,ip,jp,i+1,l-1,get_WIP(i+1,l-1),l,lp,ip,jp,get_BE(l,lp,ip,jp),lp+1,j-1,get_WIP(lp+1,j-1),temp,m3);
//				}
				//if (m3 > temp){
				//	m3 = temp;
				//}
			}

			// 4)
			if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				pf_t temp = get_WIP(i+1,l-1) * get_BE(l,lp,ip,jp) * ((cp_penalty * (j-lp+1)) + ap_penalty + (2*bp_penalty));
				d2_energy_be += temp;
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 4: WIP(%d,%d) = %d, BE(%d,%d,%d,%d)= %d ==> temp = %d  and m4 = %d\n",i,j,ip,jp,i+1,l-1,get_WIP(i+1,l-1),l,lp,ip,jp,get_BE(l,lp,ip,jp),temp,m4);
//				}
				//if (m4 > temp){
				//	m4 = temp;
				//}
			}

			// 5)
			if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				pf_t temp = (ap_penalty + (2*bp_penalty) + (cp_penalty * (l-i+1))) * get_BE(l,lp,ip,jp) * get_WIP(lp+1,j-1);
				d2_energy_be += temp;
//				if (debug && l == 9 && ip == 11){
//					printf("BE(%d,%d,%d,%d) branch 5: BE(%d,%d,%d,%d)= %d and WIP(%d,%d) = %d ==> temp = %d \n", i,j,ip,jp,l,lp,ip,jp,get_BE(l,lp,ip,jp),lp+1,j-1,get_WIP(lp+1,j-1),temp);
//				}
				//if (m5 > temp){
				//	m5 = temp;
				//}
			}
		}
	}

	// finding the min and putting it in BE[iip]
	BE[iip] = d2_energy_be; //MIN(m1,MIN(MIN(m2,m3),MIN(m4,m5)));
//	if (debug && i == 6 && ip == 11){
//		printf("BE[%d,%d,%d,%d]: m1 = %d, m2 = %d, m3 = %d, m4 = %d, m5 = %d ==> min = %d \n",i,j,ip,jp,m1,m2,m3,m4,m5,BE[iip]);
//	}
}

pf_t pseudo_loop::get_WI(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
		compute_WI(i,j,fres);
	}
	 */
	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}

// Hosna, May 1st, 2012
// I don't think we need specific getter function for pkonly case
/*
int pseudo_loop::get_WI_pkonly(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	//if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
	//	compute_WI_pkonly(i,j,fres);
	//}

	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}
*/

pf_t pseudo_loop::get_VP(int i, int j){
	// Hosna, March 16, 2012
	// two bases should be at least 3 bases apart

	if (j-i < TURN || i >= j || fres[i].pair >= 0 || fres[j].pair >= 0 || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VP[ij] == INF){
		//printf("get_VP(%d,%d), and we need to compute VP (i.e. it's INF)!\n",i,j);
		compute_VP(i,j,fres);
	}
	 */
	//printf("get_VP(%d,%d), after computation its value = %d!\n",i,j, VP[ij]);
	return VP[ij];

}
pf_t pseudo_loop::get_WMB(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN ||(fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMB[ij] == INF){
		//printf("get_WMB(%d,%d), and we need to compute WMB (i.e. it's INF)!\n",i,j);
		compute_WMB(i,j,fres);
	}
	 */
	//printf("get_WMB(%d,%d), after computation its value = %d!\n",i,j, WMB[ij]);
	return WMB[ij];
	
}

// Luke Aug 2023
pf_t pseudo_loop::get_PGPW(int i, int j){
	// i and j should be at least 3 bases apart
	if (j-i< TURN || (fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMBP[ij] == INF){
		//printf("get_WMBP(%d,%d), and we need to compute WMBP (i.e. it's INF)!\n",i,j);
		compute_WMBP(i,j,fres);
	}
	 */
	//printf("get_WMBP(%d,%d), after computation its value = %d!\n",i,j, WMBP[ij]);
	return PGPW[ij];
}

// Hosna: April 18th, 2007
// changed WMB to case 2 and WMBP
pf_t pseudo_loop::get_WMBP(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i< TURN || (fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMBP[ij] == INF){
		//printf("get_WMBP(%d,%d), and we need to compute WMBP (i.e. it's INF)!\n",i,j);
		compute_WMBP(i,j,fres);
	}
	 */
	//printf("get_WMBP(%d,%d), after computation its value = %d!\n",i,j, WMBP[ij]);
	return WMBP[ij];
}

pf_t pseudo_loop::get_BE(int i, int j, int ip, int jp){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i>= TURN && i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >=0 && fres[j].pair >= 0 && fres[ip].pair >= 0 && fres[jp].pair >= 0 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip){
		if(i == ip && j == jp && i<j){
			return 0;
		}
		int iip = index[i]+ip-i;
		// Hosna, May 1st , 2012
		// these parts are not needed any more
		/*
		if (needs_computation == 1 && BE[iip] == 0){
//			if (debug ){
//				printf("computing BE(%d,%d,%d,%d) \n",i,j,ip,jp);
//			}
			//printf("get_BE(%d,%d), and we need to compute BE (i.e. it's 0)!\n",i,ip);
			compute_BE(i,j,ip,jp,fres);
		}
		 */
//		if (debug ){
//			printf("BE[%d,%d]=%d \n",i,ip,BE[iip]);
//		}
		//printf("get_BE(%d,%d), after computation its value = %d!\n",i,ip, BE[iip]);
		return BE[iip];
	}else{
//		if (debug && i == 6 && ip == 11){
//			printf("returning INF from BE(6,69,11,24) \n");
//		}
		return INF;
	}
}

pf_t pseudo_loop::get_WIP(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) != 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WIP[ij]== INF){
		//printf("get_WIP(%d,%d), and we need to compute WIP (i.e. it's INF)!\n",i,j);
		compute_WIP(i,j,fres);
	}
	 */
	//printf("get_WIP(%d,%d), after computation its value = %d!\n",i,j, WIP[ij]);
	return WIP[ij];
}

// Hosna, May 1st, 2012
// I don't think we need specific getter function for pkonly case
/*
int pseudo_loop::get_WIP_pkonly(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) != 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	//if (needs_computation == 1 && WIP[ij]== INF){
		//printf("get_WIP(%d,%d), and we need to compute WIP (i.e. it's INF)!\n",i,j);
	//	compute_WIP_pkonly(i,j,fres);
	//}

	//printf("get_WIP(%d,%d), after computation its value = %d!\n",i,j, WIP[ij]);
	return WIP[ij];
}
*/


pf_t pseudo_loop::get_VPR(int i, int j){
	// Luke Aug 2023
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VPR[ij] == INF)
	{
		//printf("get_VPR(%d,%d), and we need to compute VPR (i.e. it's INF)!\n",i,j);
		compute_VPR(i,j,fres);
	}
	 */
	//printf("get_VPR(%d,%d), after computation its value = %d!\n",i,j, VPR[ij]);
	return VPR[ij];

}

pf_t pseudo_loop::get_VPL(int i, int j){
	// Luke Aug 2023
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VPL[ij] == INF)
	{
		//printf("get_VPL(%d,%d), and we need to compute VPL (i.e. it's INF)!\n",i,j);
		compute_VPL(i,j,fres);
	}
	 */
	//printf("get_VPL(%d,%d), after computation its value = %d!\n",i,j, VPL[ij]);
	return VPL[ij];

}

// PRE: i< j
int pseudo_loop::get_b(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j){
		return INF;
	}
	int border = MIN(border_bs[j][i],INF);
//	if(debug){
//		printf("b(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}

// PRE: i<j
int pseudo_loop::get_bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bps[j][i],-1);
//	if(debug){
//		printf("bp(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}
//PRE: i<j
int pseudo_loop::get_B(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bs[i][j],-1);
//	if(debug){
//		printf("B(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}
//PRE: i<j
int pseudo_loop::get_Bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return INF;
	}
	int border = MIN(border_bps[i][j],INF);
//	if(debug){
//		printf("Bp(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}

int pseudo_loop::get_e_stP(int i, int j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	int ss = S->get_energy(i,j,int_sequence);
//	if (debug){
//		printf("stack energy got from simfold is %d \n", ss);
//	}
	return (int)round(e_stP_penalty * (double)ss);
}

int pseudo_loop::get_e_intP(int i, int ip, int jp, int j){
	// Hosna Feb 12th, 2007:
	// this function is only being called in branch 5 of VP
	// and branch 2 of BE
	// in both cases regions [i,ip] and [jp,j] are closed regions
	int e_int = VBI->get_energy(i,j,ip,jp,int_sequence);
//	int uPup = ((ip-i-1)+(j-jp-1)) * PUP_penalty;
//	return MIN(VBI_energy,uPup);

	// Hosna April 3rd, 2007
	// based on the discussion with Anne, we decided to have
	// e_intP = 0.83 * e_int
//	printf("test: e_int(5,30,7,29) = %d \n",VBI->get_energy(5,30,7,29,int_sequence));
//	printf("e_int(%d,%d,%d,%d) = %d \n",i,j,ip,jp,e_int);
	int energy = (int)round(e_intP_penalty * (double)e_int);
//	printf("e_intP(%d,%d,%d,%d) = %d \n",i,ip,jp,j,energy);
	return energy;
}

pf_t pseudo_loop::get_energy(int i, int j){
	return get_WMB(i,j);
}

void pseudo_loop::insert_node(int i, int j, char type)
{
//	if (debug)
//	{
//		printf("\n*************************\n");
//		printf("WMB insert_node: i = %d j = %d and type = %c \n",i,j,type);
//		printf("*************************\n");
//
//	}
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
//    if (debug)
//    {
//    	printf("###################\n");
//    	printf("WMB node inserted: stack_interval.i = %d, stack_interval.j = %d and atck_interval.type = %c \n",stack_interval->i,stack_interval->j,stack_interval->type);
//    	printf("###################\n\n");
//    }

}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
