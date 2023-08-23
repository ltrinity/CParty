#include "VM_final.h"
#include "externs.h"
#include "h_externs.h"

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

    // Luke Aug 2023 init to 0 instead of INF for part func

    WM = new int [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM[i] = 0;

    WM1 = new int [total_length];
    if (WM1 == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM1[i] = 0;

    WMP = new int [total_length];
    if (WMP == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WMP[i] = 0;

    VM = new int [total_length];
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
    // Luke July 28, 2023
    // Need to remove the ambiguity and add the correct recursions and structure classes
	// Hosna June 26, 2007
	// I have to figure out how to calculate the energy here

	// here comes the copied part from simfold with all dangling energies
	int min = INF, tmp, k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
    // Luke adding new case 3 WMP (base case)
    k=i+1;
    int kjminus1 = index[k] + j-1-k;

    pf_t d2_energy_vm;
    d2_energy_vm += WMP[kjminus1] * (misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));

    //modifying for loop exit conditions sans dangles -> prev it was j-3
    for (k = i+2; k <= j-1; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;
        // Luke adding new case 1 WM + WM1
        d2_energy_vm += WM[iplus1k] * WM1[kplus1jminus1] * (misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));;
        // Luke adding new case 2 WM + WMP
        d2_energy_vm += WM[iplus1k] * WMP[kplus1jminus1] * (misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));
        //Luke adding new case 3 WMP 
        //exit early by one iter wrt cases 1 & 2
        if(k!=j-1){
            int kjminus1 = index[k] + j-1-k;
            d2_energy_vm += WMP[kjminus1] * (((k-i-1)*misc.multi_free_base_penalty) + misc.multi_helix_penalty + misc.multi_offset + AU_penalty (sequence[i], sequence[j]));
        }
        //prev case 1: WM i+1,k + WM k+1,j-1
        //tmp = WM[iplus1k] + WM[kplus1jminus1];

        //prev case extended WM i+2,k + WM k+1,j-1 i.e. i+1 unpaired
        //if (fres[i+1].pair <= -1)
        //{
        //    tmp = WM[iplus2k] + WM[kplus1jminus1] +
        //        //dangle_top [sequence [i]]
        //        //[sequence [j]]
        //        //[sequence [i+1]] +
        //        misc.multi_free_base_penalty;
        //    if (tmp < min)
        //        min = tmp;
        //}
        //prev case extended WM i+1,k + WM k+1,j-2 i.e. j-1 unpaired
        //if (fres[j-1].pair <= -1)
        //{
        //    tmp = WM[iplus1k] + WM[kplus1jminus2] +
        //        dangle_bot [sequence[i]]
        //        [sequence[j]]
        //        [sequence[j-1]] +
        //        misc.multi_free_base_penalty;
        //    if (tmp < min)
        //        min = tmp;
        //}
        //prev case extended WM i+2,k + WM k+1,j-2 i.e. i+1 and j-1 unpaired
        //if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        //{
        //    tmp = WM[iplus2k] + WM[kplus1jminus2] +
        //        dangle_top [sequence [i]]
        //        [sequence [j]]
        //        [sequence [i+1]] +
        //        dangle_bot [sequence[i]]
        //        [sequence[j]]
        //        [sequence[j-1]] +
        //        2 * misc.multi_free_base_penalty;
        //    if (tmp < min)
        //    min = tmp;
        //}
    }

    // previous case 3
	//min += misc.multi_helix_penalty + misc.multi_offset +
    //       AU_penalty (sequence[i], sequence[j]);
    //
	//int wmb_energy = this->wmb->get_energy(i,j) + a_penalty + PSM_penalty;
	int ij = index[i]+j-i;
    if(d2_energy_vm < 0){
        VM[ij] = d2_energy_vm;
    }
	//printf("VM[%d,%d] = %d \n",i,j,VM[ij]);

}

int VM_final::get_energy(int i, int j){
	int ij = index[i]+j-i;
	if (i >= j){
		return INF;
	}
	if (wmb->is_weakly_closed(i,j) == 1 || wmb->is_empty_region(i,j) == 1){
		return VM[ij];
	}
	return INF;
}

/*
int VM_final::get_energy_pk_only(int i, int j, str_features *fres){
	int ij = index[i]+j-i;
	if (i >= j || !(fres[i].pair== j && fres[j].pair==i)){
		return INF;
	}
	if (wmb->is_weakly_closed(i,j) == 1 || wmb->is_empty_region(i,j) == 1){
		return VM[ij];
	}
	return INF;
}
*/
 
/**
 *  PRE: simfold's WM matrix has been filled for i and j
 *  and now we need to fill in the WM matrix that hfold needs
 *  LT July 2023 adding code from simfold s_multi_loop for new cases
 *  WM 2 and 4 allowing unpaired + WMB, WM + WMB, respectively
 */
void VM_final::WM_compute_energy(int i, int j){
	int s_wm = s_vm->get_energy_WM(i,j);
    PARAMTYPE tmp;
    //use the for loop for splits modified for CParty
    for (int k=i; k < j; k++)
    {
        // wmb energy for split
        // Hosna: July 5th, 2007
        // add a b_penalty to this case to match the previous cases
        int wmb_energy = wmb->get_energy(k,j)+PSM_penalty+b_penalty;
        //unpaired
        int unpaired_energy =  misc.multi_free_base_penalty *(k-i);
        // new case 2 (leftmost branch pseudoknotted)
        tmp = unpaired_energy + wmb_energy;
        if (tmp < s_wm)
            {
            s_wm = tmp;
            }
        // new case 4 checking both WM for now (intermediate branch pseudoknotted)
        if (k > i && k < (j-1)){
            //4a simfold wm
            tmp = s_vm->get_energy_WM(i, k-1) + wmb_energy;
            if (tmp < s_wm)
            {
                s_wm = tmp;
            }
            //4b hfold wm
            tmp = get_energy_WM(i, k-1) + wmb_energy;
            if (tmp < s_wm)
            {
                s_wm = tmp;
            }
        }
    }
	int ij = index[i]+j-i;
	this->WM[ij] = s_wm;
//	printf("hfold's WM min = %d \n",min);
}



int VM_final::get_energy_WM(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return INF;
	}
	int ij = index[i]+j-i;
//	printf("hfold's WM(%d,%d) = %d \n", i,j,WM[ij]);
	return this->WM[ij];

}

/**
 *  LT Aug 2023 adding 
 *  WM1 should min V and WM1i,j-1 
 */
void VM_final::WM1_compute_energy(int i, int j){
    PARAMTYPE tmp = INF;
    PARAMTYPE v_energy;
    //case 2 j unpaired
    int unpaired_energy =  misc.multi_free_base_penalty;
    tmp = get_energy_WM1(i, j-1) + unpaired_energy;
    //case 1
    v_energy = v->get_energy(i,j) +
                   AU_penalty (sequence[i], sequence[j]) +
                   misc.multi_helix_penalty;
    if (v_energy < tmp)
    {
        tmp = v_energy;
    }

	int ij = index[i]+j-i;
	this->WM1[ij] = tmp;
}

int VM_final::get_energy_WM1(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return INF;
	}
	int ij = index[i]+j-i;
	return this->WM1[ij];

}

/**
 *  LT Aug 2023 adding 
 *  WMP should min P and WMPi,j-1 
 */
void VM_final::WMP_compute_energy(int i, int j){
    PARAMTYPE tmp = INF;
    //case 2 j unpaired
    int unpaired_energy =  misc.multi_free_base_penalty;
    tmp = get_energy_WMP(i, j-1) + unpaired_energy;
    //case 1
    int wmb_energy = this->wmb->get_energy(i,j) +
                   +PSM_penalty+b_penalty;
    if (wmb_energy < tmp)
    {
        tmp = wmb_energy;
    }

	int ij = index[i]+j-i;
	this->WMP[ij] = tmp;
}

int VM_final::get_energy_WMP(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return INF;
	}
	int ij = index[i]+j-i;
	return this->WMP[ij];

}
