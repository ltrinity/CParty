#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include <stdio.h>
#include <string.h>
#include "h_struct.h"
#include "h_common.h"
#include "V_final.h"
#include "VM_final.h"

class VM_final;
class V_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(char *seq, char* restricted, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM);

	// destructor
	~pseudo_loop();


	void set_features(h_str_features *f);
	// changed the code and added the initialized function
	// such that I can call WMB right after V
	void initialize();
    void compute_energies(int i, int j);

    pf_t get_energy(int i, int j);
	// in order to be able to check the border values consistantly
	// I am adding these get functions
	pf_t get_WI(int i, int j);
	// Hosna, May 1st, 2012
	// I don't think we need specific getter function for pkonly case
	//int get_WI_pkonly(int i, int j); // April 3, 2012

	pf_t get_VP(int i, int j);
	pf_t get_WMB(int i, int j);
	pf_t get_BE(int i, int j, int ip, int jp);
	pf_t get_WIP(int i, int j);
	// Hosna, May 1st, 2012
	// I don't think we need specific getter function for pkonly case
	//int get_WIP_pkonly(int i, int j); // April 3, 2012


	//Luke Aug 2023
	pf_t get_VPR(int i, int j);
	pf_t get_VPL(int i, int j);

	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	pf_t get_WMBP(int i, int j);
	// Luke adding for PGPW
	pf_t get_PGPW(int i, int j);

    int is_weakly_closed(int i, int j);
    int is_empty_region(int i, int j);


    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    char *get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}

private:

	int nb_nucleotides;
	int *int_sequence;
	char *sequence;
	char *restricted;

    s_hairpin_loop *H;      // hairpin loop object
    s_stacked_pair *S;      // stack pair object
    s_internal_loop *VBI;   // internal loop object
    VM_final *VM;	        // multi loop object
    V_final *V;		        // the V object

	h_str_features *fres;
	seq_interval *stack_interval;
	char *structure;
	minimum_fold *f;

	int needs_computation; // This global variable is used so that we don't compute energies in backtracking

	//Hosna
    pf_t *WI;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
    pf_t *VP;				// the loop corresponding to the pseudoknotted region of WMB
    pf_t *WMB;				// the main loop for pseudoloops and bands
    int *weakly_closed;		// the array which is keeping track of which regions are weakly closed
    int *not_paired_all;	// the array which keeps track of empty regions
    int *index;				// the array to keep the index of two dimensional arrays like WI and weakly_closed
    int **border_bs;		// keeps track of border_b and border_B
    int **border_bps;		// keeps track of border_bp and border_Bp
    pf_t *WIP;				// the loop corresponding to WI'
	//Luke Aug 2023
	pf_t *VPR;				// the loop corresponding to VPR
	pf_t *VPL;				// the loop corresponding to VPL
    pf_t *BE;				// the loop corresponding to BE

    // Hosna, April 18th, 2007
	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	pf_t *WMBP; 				// the main loop to calculate WMB
	// Luke adding
	pf_t *PGPW;

    // function to allocate space for the arrays
    void allocate_space();

    void compute_WI(int i, int j, h_str_features *fres);
	// Hosna: This function is supposed to fill in the WI array

	void compute_WI_pkonly(int i, int j, h_str_features *fres);
	// April 3, 2012
	// Hosna: This function is supposed to fill in the WI array with pk only base pairs

	void compute_VP(int i, int j, h_str_features *fres);
	// Hosna: this function is supposed to fill the VP array

	void compute_WMB(int i, int j, h_str_features *fres);
	// Hosna: this function is supposed to fill the WMB array

	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	void compute_WMBP(int i, int j, h_str_features *fres);
	// this is the helper recurrence to fill the WMB array

	// Luke adding
	void compute_PGPW(int i, int j, h_str_features *fres);
	// this is the helper HELPER recurrence to fill the WMB array

	void compute_BE(int i, int j, int ip, int jp, h_str_features *fres);
	// Hosna: this function is supposed to fill the BE array

	void compute_WIP(int i, int j, h_str_features *fres);
	// Hosna: this function is supposed to fill the WIP array

	void compute_WIP_pkonly(int i, int j, h_str_features *fres);
	// April 3, 2012
	// Hosna: this function is supposed to fill the WIP array with pk only base pairs


	void compute_VPR(int i, int j, h_str_features *fres);
	// Luke: this function is supposed to fill the VPR array

	void compute_VPL(int i, int j, h_str_features *fres);
	// Luke: this function is supposed to fill the VPL array

	// Hosna Feb 8th, 2007:
	// To get the band borders easily
	// I am adding the following functions
	int get_b(int i,int j);
	int get_bp(int i,int j);
	int get_B(int i,int j);
	int get_Bp(int i,int j);

	// Hosna Feb 8th, 2007:
	// I have to calculate the e_stP in a separate function
	pf_t get_e_stP(int i, int j);
	pf_t get_e_intP(int i,int ip, int jp, int j);

  	// Hosna: Feb 19th 2007
  	// used for backtracking
  	void insert_node (int i, int j, char type);//, seq_interval *stack_interval);

};
#endif /*PSEUDO_LOOP_H_*/
