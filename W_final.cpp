
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pseudo_loop.h"
#include "V_final.h"
#include "W_final.h"
#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"




// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(char *seq, char *res):s_min_folding(seq,res)
{
	this->nb_nucleotides = strlen(seq);
	space_allocation();
}


W_final::~W_final()
{
	delete vm;
	delete v;
	delete WMB;

}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// Hosna June 20th, 2007
	vm = new VM_final(this->int_sequence,this->nb_nucleotides);
	if (vm == NULL) giveup ("Cannot allocate memory", "W_final");
	if (debug){
		printf("nb_nucleotides = %d \n",this->nb_nucleotides);
	}

	// Hosna June 20th, 2007
	v = new V_final(nb_nucleotides);
	if (v == NULL) giveup ("Cannot allocate memory", "W_final");
	//s_min_folding::V, s_min_folding::H, s_min_folding::S, s_min_folding::VBI, vm);
	v->setloops(this->V,vm);

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (sequence,restricted,v,this->H,this->S,this->VBI,vm);
    if (WMB == NULL) giveup ("Cannot allocate memory", "W_final");

    // Hosna: June 20th 2007
    vm->set_V_matrix(v);
    vm->set_WMB_matrix(WMB);

}

double W_final::hfold(){

	pf_t energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);

    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);


	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=0; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
            V->compute_energy_restricted (i, j, fres);

        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted (j, fres);


		// test V values
		/*
		 for (i=0; i<j; i++)
		 {
		 if (fres[i].pair ==j && fres[j].pair ==i){
		 printf("---->> V(%d,%d) = %d \n",i,j, V->get_energy(i,j));

		 }
		 }
		 */
    }

	for (j=0; j < nb_nucleotides; j++)
    {
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			
			vm->compute_energy (i,j, fres);
			WMB->compute_energies(i,j);
			
			//add in the pk cases later
			vm->WM_compute_energy(i,j);
			vm->WM1_compute_energy(i,j);
			//        	if (debug){
			//        		printf("WM_final(%d,%d) = %d \n",i,j,vm->get_energy_WM(i,j));
			//        	}
        }

	}

	// end of addition at March 8, 2012, Hosna
	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted(j,fres);
    }
    energy = this->W[nb_nucleotides-1];
	//    printf("energy = %f \n", energy);

    if (debug)
    {
        print_result ();
    }
    delete [] h_fres;
    delete [] fres;
    return energy;

}

void W_final::return_structure(char *structure){
	strcpy (structure, this->structure);
	//s_min_folding::return_structure(structure);
}

/// Luke Sep 2023 Z_W exterior structure class
/// Luke init partition function
void W_final::compute_W_restricted (int j, str_features *fres)
// compute W(j)
{
    pf_t m1, m2, m3;
    int must_choose_this_branch;
	pf_t d2_energy = 0;
    m1 = W[j-1];
	//if(debug){
	//printf("W_[0,%d] = %Lf\n",j-1,W[j-1]);
	//}
	d2_energy+= m1;
    m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
	//if(debug){
	//printf("V[0,%d] = %Lf\n",j,m2);
	//}
	d2_energy+= m2;
    m3 = compute_W_br3_restricted (j, fres);
	//printf("P[0,%d] = %Lf\n",j,m3);
	d2_energy+= m3;
	
    //if (WMB->is_weakly_closed(0,j) < 0){
		// Luke changing to 0 for part func
    	//W[j] = 0;
    	//return;
    //}
	W[j] = d2_energy;
	//if(debug){
	//printf("W[j] d2 = %Lf = %Lf\n",d2_energy, W[j]);
	//printf("Z_W(%d,%d) = %Lf \n",0,j,d2_energy);
	//}
}

/// Luke modifying to multiply structure classes
/// Return sum of energy scaled energies Luke Sep 2023
pf_t W_final::compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch)
{
	pf_t energy_ij = 0, acc;
    int i;

	pf_t d2_energy_v = 0;

    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        //acc = (i-1>0) ? W[i-1]: 0;
		pf_t Wsubstruc_en = W[i-1];
		if(i-1 == -1){
			Wsubstruc_en = 1;
		}
        energy_ij = Wsubstruc_en*v->get_energy(i,j)* bw_int(AU_penalty (int_sequence[i],int_sequence[j]));
		if(energy_ij> 0){
			//printf("br2: W_i-1(%d,%d) = %Lf\n", i,j,Wsubstruc_en );
			//printf("br2: v(%d,%d) = %Lf\n", i,j,v->get_energy(i,j) );
			//printf("br2: total(%d,%d) = %Lf\n", i,j,energy_ij );
		}
		d2_energy_v+=energy_ij;
		//if(debug){
		//printf("Z_V(%d,%d) = %Lf \n",i,j,d2_energy_v);
		//}
    }
	return d2_energy_v;
}

/// Luke modifying to multiply structure classes
/// Return sum of energy scaled energies Luke Sep 2023
pf_t W_final::compute_W_br3_restricted(int j, str_features *fres){
	// Hosna June 30, 2007
	// The following would not take care of when
	// we have some unpaired bases before the start of the WMB
	//return WMB->get_energy(0,j) + PS_penalty;
	pf_t tmp, energy_ij = 0, acc;
    int i;

	pf_t d2_energy_p = 0;


    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		// Hosna: July 9, 2007
		// We only chop W to W + WMB when the bases before WMB are free
		if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

	        // We don't need to make sure i and j don't have to pair with something else,
	        //  because that would be INF - done in fold_sequence_restricted
	        //acc = (i-1>0) ? W[i-1]: 0;
			//printf("acc= %Lf \n",acc);
			
			pf_t Wsubstruc_en = W[i-1];
			if(i-1 == -1){
				Wsubstruc_en = 1;
			}
	        energy_ij = Wsubstruc_en* WMB->get_energy(i,j) * PS_penalty;
			
			//if(debug){
			//printf("Z_P(%d,%d) = %Lf; %Lf; %f\n",i,j,energy_ij,WMB->get_energy(i,j), PS_penalty);
			//}
			d2_energy_p += energy_ij;
			//if(debug){
			//printf("Z_P(%d,%d) = %Lf \n",i,j,d2_energy_p);
			//}
		}
	}
	return d2_energy_p;
}

//depracated 
void W_final::print_result ()
// PRE:  The matrix V has been calculated and the results written in f
// POST: Prints details of each elementary structure
{
    int i;
    int energy = INF, sum;

    printf ("Minimum energy: %Lf\n", W[nb_nucleotides-1]);
    sum = 0;

    for (i=0; i< nb_nucleotides; i++)
    {
        if (f[i].pair > i)
        {
			//Hosna March 8, 2012
			// changing nested ifs to switch for optimality
			switch (f[i].type){
				case HAIRP:
				//if (f[i].type == HAIRP)
					energy = V->get_energy(i, f[i].pair);
					break;
				case STACK:
				//else if (f[i].type == STACK)
					energy = V->get_energy(i, f[i].pair) - V->get_energy(i+1, f[i+1].pair);
					break;
				case P_VP:
				// Hosna: June 28th, 2007
				//else if (f[i].type == P_VP){
					energy = WMB->get_VP(i,f[i].pair);
					break;
				//case P_VPP:
				//}else if(f[i].type == P_VPP){
				//	energy = WMB->get_VPP(i,f[i].pair);
				//}
				//	break;
			}
            printf ("Pair (%d,%d), type %c,\tenergy %6d\n", i, f[i].pair, f[i].type, energy);
            sum += energy;
        }
    }
    printf ("0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n");
    printf ("%s\n", sequence);
    printf ("%s\n", structure);

}
