#include "VM_final.h"
#include "externs.h"
#include "h_externs.h"
#include "math.h"

VM_final::VM_final(int *seq, int len)
{
	length = len;
	sequence = seq;
	this->v = NULL;
	this->wmb = NULL;

	index = new int[length];    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    int i;
    for (i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    /// Luke Aug 2023 init to 0 instead of INF for part func
    /// Luke created new structure classes for inside of multiloop

    WM = new pf_t [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM[i] = 0;

    WM1 = new pf_t [total_length];
    if (WM1 == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM1[i] = 0;

    WMP = new pf_t [total_length];
    if (WMP == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WMP[i] = 0;

    VM = new pf_t [total_length];
    if (VM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) VM[i] = 0;

//    printf("an object of VM_final was successfully created! \n");
}

VM_final::~VM_final()
{
	delete [] index;
    delete [] WM;
    delete [] WM1;
    delete [] WMP;
    delete [] VM;
}

void VM_final::compute_energy(int i, int j, str_features *fres){
    /// Luke July 28, 2023
    /// Need to remove the ambiguity and add the correct recursions and structure classes
    
	// Hosna June 26, 2007
	// I have to figure out how to calculate the energy here

	// here comes the copied part from simfold with all dangling energies
	int min = INF, tmp, k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
    
    k=i+1;
    int kjminus1 = index[k] + j-1-k;

    pf_t d2_energy_vm = 0;
    //d2_energy_vm += WMP[kjminus1] * exp(misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j])*oneoverRT);
    /// Luke adding new case 3 WMP (base case)
    ///modifying for loop exit conditions sans dangles -> prev it was j-3
    for (k = i+2; k <= j-1; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;
        /// Luke adding new case 1 WM + WM1 below cases 2 and 3
        /// penalties scaled by boltzmann weights
        d2_energy_vm += (WM[iplus1k] * WM1[kplus1jminus1] * bw_int(misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j])));;
        /// Luke adding new case 2 WM + WMP
        d2_energy_vm += WM[iplus1k] * WMP[kplus1jminus1] * bw_int(misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));
        //Luke adding new case 3 WMP exit early by one iter wrt cases 1 & 2
        if(k!=j-1){
            int kjminus1 = index[k] + j-1-k;
            d2_energy_vm += WMP[kjminus1] * bw_int(((k-i-1)*misc.multi_free_base_penalty) + misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));
        }
        
    }

    // previous case 3
	//min += misc.multi_helix_penalty + misc.multi_offset +
    //       AU_penalty (sequence[i], sequence[j]);
    //
	//int wmb_energy = this->wmb->get_energy(i,j) + a_penalty + PSM_penalty;
	int ij = index[i]+j-i;
    //if(debug){
    //printf("VM[%d,%d] = %Lf \n",i,j,d2_energy_vm);
    //}
    VM[ij] = d2_energy_vm;
    //if (VM[ij] > 0.0){
	//    printf("VM[%d,%d] = %Lf \n",i,j,VM[ij]);
    //}
}

//Luke Aug 2023 base case 0 for partition function
pf_t VM_final::get_energy(int i, int j){
	int ij = index[i]+j-i;
	if (i >= j){
		return 0;
	}
	if (wmb->is_weakly_closed(i,j) == 1 || wmb->is_empty_region(i,j) == 1){
		return VM[ij];
	}
	return 0;
}

/********************************************//**
 *  PRE: simfold's WM matrix has been filled for i and j
 *  and now we need to fill in the WM matrix that hfold needs
 *  LT July 2023 adding code from simfold s_multi_loop for new cases
 *  WM 2 and 4 allowing unpaired + WMB, WM + WMB, respectively
 *  Luke Aug 2023 modifying for part func
 ***********************************************/
void VM_final::WM_compute_energy(int i, int j){
	pf_t s_wm = s_vm->get_energy_WM(i,j);
    //PARAMTYPE tmp;
    pf_t d2_energy_wm = 0;
    ///use the for loop for splits modified for CParty
    /// multiply penalties Luke Sep 2023
    for (int k=i; k < j; k++)
    {
        // wmb energy for split
        // Hosna: July 5th, 2007
        // add a b_penalty to this case to match the previous cases
        pf_t wmb_energy = wmb->get_energy(k,j)*PSM_penalty*b_penalty;
        //unpaired
        pf_t unpaired_energy =  bw_int(misc.multi_free_base_penalty*(k-i)) ;
        // new case 2 (leftmost branch pseudoknotted)
        d2_energy_wm += (unpaired_energy * wmb_energy);
        // new case 4 checking both WM for now (intermediate branch pseudoknotted)
        if (k > i && k < (j-1)){
            //4a simfold wm
            d2_energy_wm += s_vm->get_energy_WM(i, k-1) * wmb_energy;
            //4b hfold wm
            //tmp = get_energy_WM(i, k-1) + wmb_energy;
        }
    }
	int ij = index[i]+j-i;
    //pk free for now
    if(s_wm == INF){
        s_wm = 0;
    }
	this->WM[ij] = s_wm;
	//printf("hfold's WM min = %Lf \n",s_wm);
}



pf_t VM_final::get_energy_WM(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return 0;
	}
	int ij = index[i]+j-i;
//	printf("hfold's WM(%d,%d) = %d \n", i,j,WM[ij]);
	return this->WM[ij];

}

/********************************************//**
 *  LT Aug 2023 adding 
 *  WM1 should min V and WM1i,j-1 
 */
void VM_final::WM1_compute_energy(int i, int j){
    pf_t v_energy;
    pf_t d2_energy_wm1 = 0;
    //case 2 j unpaired
    pf_t unpaired_energy =  bw_int(misc.multi_free_base_penalty);
    pf_t wm1substruc_en = get_energy_WM1(i, j-1);
    d2_energy_wm1 +=  wm1substruc_en * unpaired_energy;
    //printf("CParty Z_WM1 case 2 = %Lf \n",get_energy_WM1(i, j-1) * unpaired_energy);
    //case 1
    v_energy = v->get_energy(i,j) * bw_int(
                   AU_penalty (sequence[i], sequence[j]) +
                   misc.multi_helix_penalty);
    //printf("CParty Z_WM1 case 1 = %Lf \n",v_energy);
    d2_energy_wm1 += v_energy;
	int ij = index[i]+j-i;
	this->WM1[ij] = d2_energy_wm1;
    //printf("CParty Z_WM1 = %Lf \n",d2_energy_wm1);
}

pf_t VM_final::get_energy_WM1(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return 0;
	}
	int ij = index[i]+j-i;
	return this->WM1[ij];

}

/********************************************//**
 *  LT Aug 2023 adding 
 *  WMP should min P and WMPi,j-1 
 */
void VM_final::WMP_compute_energy(int i, int j){
    //PARAMTYPE tmp = INF;
    pf_t d2_energy_wmp = 0;
    //case 2 j unpaired
    pf_t unpaired_energy =  misc.multi_free_base_penalty;
    d2_energy_wmp += get_energy_WMP(i, j-1) * unpaired_energy;
    //case 1
    int wmb_energy = this->wmb->get_energy(i,j) *PSM_penalty*b_penalty;
    d2_energy_wmp += wmb_energy;

	int ij = index[i]+j-i;
	this->WMP[ij] = d2_energy_wmp;
}

pf_t VM_final::get_energy_WMP(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return 0;
	}
	int ij = index[i]+j-i;
	return this->WMP[ij];

}
