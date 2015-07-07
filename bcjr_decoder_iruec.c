#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/* note, the second half of the table must only be the second transition to each state */
/* last column is transition probablity */


#define MAX_CODEBITS 4


#define TRAN_LEN (MAX_CODEBITS+3)

#ifndef max
 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif
__inline  double maxstar(double a, double b)
{    
    return  max(a,b) + log(1.0+exp(-fabs(a-b)));    
}
__inline int round(double d)
{
    return (int)(d + 0.5);
}


typedef int trans_t[][TRAN_LEN];

int calculate_uec_edge_transitions(int *tran_cur, int states_cur, int trans_cur, int *tran_nxt, int states_nxt, int trans_nxt,  int *trans_out)
{
	int c,max_even,max_odd,outtrans,c1;

	if (states_cur & 1)
		mexErrMsgTxt("states needs to be even");
	if (states_nxt & 1)
		mexErrMsgTxt("states needs to be even");
	

	outtrans = 0;
	if (states_cur <= states_nxt){     /* if number of states increasing, just output that larger trellis (but skip some states */
		for (c = 0; c < (trans_nxt); c++){
			if (tran_nxt[c*TRAN_LEN + 0] < states_cur){
				for (c1 = 0; c1 < TRAN_LEN; c1++)
					trans_out[outtrans*TRAN_LEN + c1] = tran_nxt[c*TRAN_LEN + c1];	
				outtrans++;
			}
		}
		return outtrans;
	}
	
		
	/* copy current trellis as a starting point */
	for (c = 0; c < (trans_cur*TRAN_LEN); c++)
		trans_out[c] = tran_cur[c];		
		
	max_even = states_nxt-1;
	max_odd = states_nxt-1-1;
		
	/* go through each transition and see if it needs modifying */
	for (c = 0; c < trans_cur; c++){
		if (trans_out[c*TRAN_LEN + 1] >= states_nxt){
			if ((trans_out[c*TRAN_LEN + 1]+1) & 1)
				trans_out[c*TRAN_LEN + 1] = max_odd;
			else
				trans_out[c*TRAN_LEN + 1] = max_even;
		}
	}	
	return trans_cur;
	//return trans_out;
}

void bcjr_decoder_iruec(double *uncoded_in, double *coded_in1, double *coded_in2, double *coded_in3,double *coded_in4,
        double *uncoded_out, double *coded_out1, double *coded_out2,double *coded_out3,double *coded_out4,
        int len, int* trans[], int *state_count, int *trans_count, int codebit_count, int codebook_count,
        double *trans_probs[], double *fractions, int last_state)
//void bcjr_decoder_iruec(double *uncoded_in, double *coded_in1, double *coded_in2, double *coded_in3,double *coded_in4,
//        double *uncoded_out, double *coded_out1, double *coded_out2,double *coded_out3,double *coded_out4,
//        int len, trans_t trans, int state_count, int trans_count, int codebit_count,
//        int reg_trellis, double *trans_prob, int last_state)
{
    
    int i,j;
    //int transc = trans_count;
    int set0,set1;
    double p1,p0;
    int a,b,c,sb,eb,tc;
    double **gammas;
    double **alphas;
    double **betas;
	int *startbit;
	int *endbit;
	double fraccum;
    //int states = state_count;
	int max_states = 0;
	int max_trans = 0;
	int *tran;
	double *trans_prob;
	int *tran_edge;

	
	/* calculate the start/end index of each codebook */
	startbit = malloc(codebook_count*sizeof(int));
	endbit = malloc(codebook_count*sizeof(int));
	fraccum = 0;
	for (c = 0; c < codebook_count; c++){
		if (c == 0)
			startbit[0] = 0;
		else
			startbit[c] = endbit[c-1];
		
		fraccum += fractions[c];
		if (fraccum > 1)
			mexErrMsgIdAndTxt("MyToolbox:inputError","sum(fractions) is greater than 1");
		//if (c == (codebook_count-1))
		//	endbit[codebook_count-1] = len-1;
		//else
			endbit[c] = (int)round((double)len*fraccum);
	}
	/*
	for (c = 0; c < codebook_count; c++){
		mexPrintf("%i  %i\n",startbit[c],endbit[c]);
	}*/

	/*work out max states and transitions */
	for (c = 0; c < codebook_count; c++){
		if (max_states < state_count[c])
			max_states = state_count[c];
		if (max_trans < trans_count[c])
			max_trans = trans_count[c];	
	}
	tran_edge = malloc(max_trans*sizeof(int)*TRAN_LEN);
	
	
	gammas = (double**) malloc (max_trans*sizeof(double *));
		for (i = 0; i < max_trans; i++)
			gammas[i] = (double*) malloc(len*sizeof(double));
	/* calculate gammas */ 		
	for (c = 0; c < codebook_count; c++){				
		sb = startbit[c];
		eb = endbit[c];
		tran = trans[c];
		trans_prob = trans_probs[c];
		for (i = 0; i < trans_count[c]; i++)
		{
			if (tran[i*TRAN_LEN + 2] && tran[i*TRAN_LEN + 3]){
				for (j = sb; j < eb; j++)
					gammas[i][j] = uncoded_in[j] + coded_in1[j];
			}
			else if (tran[i*TRAN_LEN + 2]){
				for (j = sb; j < eb; j++)
					gammas[i][j] = uncoded_in[j];
			}
			else if (tran[i*TRAN_LEN + 3]){
				for (j = sb; j < eb; j++)
					gammas[i][j] = coded_in1[j];
			}
			else{
				for (j = sb; j < eb; j++)
					gammas[i][j] = 0;
			}
		}
		
		if (codebit_count >= 2){
			for (i = 0; i < trans_count[c]; i++){
				if (tran[i*TRAN_LEN + 4]){
					for (j = sb; j < eb; j++)
						gammas[i][j] += coded_in2[j];
				}
			}
		}
		if (codebit_count >= 3){
			for (i = 0; i < trans_count[c]; i++){
				if (tran[i*TRAN_LEN + 5]){
					for (j = sb; j < eb; j++)
						gammas[i][j] += coded_in3[j];
				}
			}
		}
		if (codebit_count >= 4){
			for (i = 0; i < trans_count[c]; i++){
				if (tran[i*TRAN_LEN + 6]){
					for (j = sb; j < eb; j++)
						gammas[i][j] += coded_in4[j];
				}
			}
		}
	

		for (i = 0; i < trans_count[c]; i++){
			for (j = sb; j < eb; j++)
				gammas[i][j] += trans_prob[i];
		}
	}
    
    
    
    /* set and initialise memory */
    alphas = (double**) malloc (max_states*sizeof(double *));
    for (i = 0; i < max_states; i++)
        alphas[i] = (double*) malloc(len*sizeof(double));
    betas = (double**) malloc (max_states*sizeof(double *));
    for (i = 0; i < max_states; i++)
        betas[i] = (double*) malloc(len*sizeof(double));
    
	for (i = 0; i < max_states; i++) {
		for (j=0; j< len; j++)
			alphas[i][j] = -9000;
	}
	for (i = 0; i < max_states; i++){
		for (j=0; j< len; j++)
			betas[i][j] = -9000;
	}

    
    
    /* forward recursion (alphas) */
    alphas[0][0] = 0;        /* first state */
    for (i = 1; i < max_states; i++)
        alphas[i][0] = -9000;
	for (c = 0; c < codebook_count; c++){				
		sb = startbit[c];
		if (sb == 0)
			sb = 1;	
		eb = endbit[c];
		if ((c > 0) && (sb > 1)){/* if  there is another trellis after this, the last bit is special */						
			i = sb;
			tc = calculate_uec_edge_transitions(trans[c-1], state_count[c-1], trans_count[c-1],trans[c], state_count[c],trans_count[c], tran_edge);
			tran = tran_edge;
			for (j = 0; j < tc; j++)
				alphas[tran[j*TRAN_LEN + 1]][i] = maxstar(  alphas[tran[j*TRAN_LEN + 1]][i],  alphas[tran[j*TRAN_LEN + 0]][i-1] + gammas[j][i-1]  );
			sb++;
		}
		tran = trans[c];
		for (i = sb; i < eb; i++){
			for (j = 0; j < trans_count[c]; j++)
				alphas[tran[j*TRAN_LEN + 1]][i] = maxstar(  alphas[tran[j*TRAN_LEN + 1]][i],  alphas[tran[j*TRAN_LEN + 0]][i-1] + gammas[j][i-1]  );
		}
		
	}
    
    
    
    /* backwards recursion (betas) */
    /* double betas[8][len]; */
    if (last_state < 1 || last_state > max_states){
        for (i = 0; i < max_states; i++)
            betas[i][len-1] = 0;     /* end state unknown */
    }
    else
    {
        for (i = 0; i < max_states; i++)
            betas[i][len-1] = -9000;
        betas[last_state-1][len-1] = 0;
    }   
	
	for (c = codebook_count-1; c >= 0; c--){				
		sb = startbit[c]-1;
		eb = endbit[c];
		if (eb == len)   /* last state is already initialised */
			eb = len-1;
		else{
			if ((c+1) >= codebook_count)
				mexErrMsgTxt("Fractions do not add to 1");
			eb-=2;  /* if  there is another trellis after this, the last bit is special */
			i = eb;
			tc = calculate_uec_edge_transitions(trans[c], state_count[c], trans_count[c],trans[c+1], state_count[c+1],trans_count[c+1], tran_edge);
			tran = tran_edge;
			for (j = 0; j < tc; j++)
				betas[tran[j*TRAN_LEN + 0]][i] = maxstar( betas[tran[j*TRAN_LEN + 0]][i],  betas[tran[j*TRAN_LEN + 1]][i+1] + gammas[j][i+1]  );
		}			
		tran = trans[c];
		for (i = eb-1; (i >= sb) && (i>=0); i--){
			for (j = 0; j < trans_count[c]; j++)
				betas[tran[j*TRAN_LEN + 0]][i] = maxstar( betas[tran[j*TRAN_LEN + 0]][i],  betas[tran[j*TRAN_LEN + 1]][i+1] + gammas[j][i+1]  );
		}
	}

    
    
    /*
    mexPrintf("trans probs\n");
    for (i=0;i<max_trans;i++)
    {
		mexPrintf("%f ",trans_prob[i]);    
    }

    mexPrintf("GAMMAS\n");
    for (i=0;i<max_trans;i++)
    {
    for (j=0;j<len;j++)
    mexPrintf("%f ",gammas[i][j]);
    mexPrintf("\n");
    
    } 
    mexPrintf("ALPHAS\n");
    for (i=0;i<max_states;i++)
    {
    for (j=0;j<len;j++)
    mexPrintf("%f ",alphas[i][j]);
    mexPrintf("\n");
    }
    mexPrintf("\n");
    mexPrintf("\n");
    mexPrintf("BETAS\n");
    for (i=0;i<max_states;i++)
    {
    for (j=0;j<len;j++)
    mexPrintf("%f ",betas[i][j]);
    mexPrintf("\n");
    
    }
     mexPrintf("\n");
     */
	 
    /* deltas */
    /* reuse gammas memory */
	for (c = codebook_count-1; c >= 0; c--){				
		sb = startbit[c];
		eb = endbit[c];
		tran = trans[c];
		if (c < codebook_count - 1)/* if  there is another trellis after this, the last bit is special */						
			eb = endbit[c] - 1;
		for (i = 0; i < trans_count[c]; i++){
			a = tran[i*TRAN_LEN + 0];
			b = tran[i*TRAN_LEN + 1];
			for (j = sb; j < eb; j++)
				gammas[i][j] += alphas[a][j] + betas[b][j];
		}
		if (c < codebook_count - 1){ /* if  there is another trellis after this, the last bit is special */	
			i = endbit[c] - 1;
			tc = calculate_uec_edge_transitions(trans[c], state_count[c], trans_count[c],trans[c+1], state_count[c+1],trans_count[c+1], tran_edge);
			tran = tran_edge;
			for (i = 0; i < tc; i++){
				a = tran[i*TRAN_LEN + 0];
				b = tran[i*TRAN_LEN + 1];
				gammas[i][j] += alphas[a][j] + betas[b][j];			
			}
		}
	}
    
    
    /* extrinisic uncoded llr */
	for (c = codebook_count-1; c >= 0; c--){				
		sb = startbit[c];
		eb = endbit[c];
		tran = trans[c];
		for (j = sb; j < eb; j++)
		{
			set0 = 0;
			set1 = 0;
			
			for (i = 0; i < trans_count[c]; i++)
			{
				if (tran[i*TRAN_LEN + 2]){
					if (set1)
						p1 = maxstar(p1,gammas[i][j]);
					else
					{
						set1 = 1;
						p1 = gammas[i][j];
					}
				}
				else{
					if (set0)
						p0 = maxstar(p0,gammas[i][j]);
					else
					{
						set0 = 1;
						p0 = gammas[i][j];
					}
				}
			}
			uncoded_out[j] = p1 - p0 - uncoded_in[j];
		}
    
    
		/* coded llr out */
		for (j = sb; j < eb; j++)
		{
			set0 = 0;
			set1 = 0;
			
			for (i = 0; i < trans_count[c]; i++)
			{
				if (tran[i*TRAN_LEN + 3]){
					if (set1)
						p1 = maxstar(p1,gammas[i][j]);
					else
					{
						set1 = 1;
						p1 = gammas[i][j];
					}
				}
				else{
					if (set0)
						p0 = maxstar(p0,gammas[i][j]);
					else
					{
						set0 = 1;
						p0 = gammas[i][j];
					}
				}
			}
			coded_out1[j] = p1 - p0 - coded_in1[j];
		}
		
		/* coded llr out2 */
		if (codebit_count >= 2)
		{
			for (j = sb; j < eb; j++)
			{
				set0 = 0;
				set1 = 0;
				
				for (i = 0; i < trans_count[c]; i++)
				{
					if (tran[i*TRAN_LEN + 4]){
						if (set1)
							p1 = maxstar(p1,gammas[i][j]);
						else
						{
							set1 = 1;
							p1 = gammas[i][j];
						}
					}
					else{
						if (set0)
							p0 = maxstar(p0,gammas[i][j]);
						else
						{
							set0 = 1;
							p0 = gammas[i][j];
						}
					}
				}
				coded_out2[j] = p1 - p0 - coded_in2[j];
			}
		}
		
		/* coded llr out3 */
		if (codebit_count >= 3)
		{
			for (j = sb; j < eb; j++)
			{
				set0 = 0;
				set1 = 0;
				
				for (i = 0; i < trans_count[c]; i++)
				{
					if (tran[i*TRAN_LEN + 5]){
						if (set1)
							p1 = maxstar(p1,gammas[i][j]);
						else
						{
							set1 = 1;
							p1 = gammas[i][j];
						}
					}
					else{
						if (set0)
							p0 = maxstar(p0,gammas[i][j]);
						else
						{
							set0 = 1;
							p0 = gammas[i][j];
						}
					}
				}
				coded_out3[j] = p1 - p0 - coded_in3[j];
			}
		}
		/* coded llr out4 */
		if (codebit_count >= 4)
		{
			for (j = sb; j < eb; j++)
			{
				set0 = 0;
				set1 = 0;
				
				for (i = 0; i < trans_count[c]; i++)
				{
					if (tran[i*TRAN_LEN + 6]){
						if (set1)
							p1 = maxstar(p1,gammas[i][j]);
						else
						{
							set1 = 1;
							p1 = gammas[i][j];
						}
					}
					else{
						if (set0)
							p0 = maxstar(p0,gammas[i][j]);
						else
						{
							set0 = 1;
							p0 = gammas[i][j];
						}
					}
				}
				coded_out4[j] = p1 - p0 - coded_in3[j];
			}
		}
	}
    
    for (i = 0; i < max_trans; i++){
        free(gammas[i]);
    }
    free(gammas);
    for (i = 0; i < max_states; i++){
        free(alphas[i]);
    }
    free(alphas);
    for (i = 0; i < max_states; i++){
        free(betas[i]);
    }
    free(betas);
	free(startbit);
	free(endbit);
	free(tran_edge);
    
    
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    double *uncoded_in;
    double *coded_in1;
	double *coded_in2;
	double *coded_in3;
	double *coded_in4;
    double *uncoded_out;
    double *coded_out1;
    double *coded_out2;
	double *coded_out3;
	double *coded_out4;
	double *fractions;
    int len,i,j,last_state,ck_i,ck_mem;
    char *input_buf;
//    int blank_transitions[600][TRAN_LEN];
    double *trans_probs[20];
	int transc_list[20];
	int statec_list[20];
	double *current_prob;
    double *p;
//	int *codelibrary;
	int *codelibrary_pointers[20];
	int *current_tran;
	int codebit_count = 1;
	int codebook_count = 1;
	const mxArray *cell_element_ptr;
	int this_cb;
	
    last_state = -1;
    
    /* check for proper number of arguments */
    if(nrhs!=6 && nrhs!=5)
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
		"Four inputs required. [uncoded_out, coded_out] = bcjr_decoder(uncoded_in, coded_in, 'iruec', transitions_cell, fractions, {, last_state}). States start at 1, set last_state < 0 for unknown");
    
    
    /* make sure the 1st input argument is type double */
    if( !mxIsDouble(prhs[0]) ||
            mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    /* make sure the 2nd input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
            mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
   
    
    /* check that number of rows in 1st input argument is 1 */
    
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector 1.");
    }
    
    len = mxGetN(prhs[0]);
    /*if (len < 3)
    *    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    *if (mxGetN(prhs[1]) < 3)
    *    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector."); */
    if (mxGetN(prhs[0]) != mxGetN(prhs[1]))
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be same dimentions.");
    
    
    input_buf = mxArrayToString(prhs[2]);
    
    if ((strncmp(input_buf,"iruec",5)!=0) || (strncmp(input_buf,"new",5)==0))
    {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be 'iruec'.");
	}
        
	/*
	if( !mxIsDouble(prhs[3]) ||
			mxIsComplex(prhs[3])) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
	}	
	if(mxGetM(prhs[3])==1 || mxGetM(prhs[3])>599) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword), Tx6 (two codewords) or Tx7 (two codewords).  Max transitions = 600");
	}
	if(mxGetN(prhs[3])!=5 && mxGetN(prhs[3])!=6 && mxGetN(prhs[3])!=7 && mxGetN(prhs[3])!=8) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword), Tx6 (two codewords) or Tx7 (two codewords).  Max transitions = 600");
	}
	if(mxGetN(prhs[3])==6){
		codebit_count = 2;
	}
	if(mxGetN(prhs[3])==7){
		codebit_count = 3;
	}
	if(mxGetN(prhs[3])==8){
		codebit_count = 4;
	}
	
	if (codebit_count == 2)
	{		
		if(mxGetM(prhs[1])!=2) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs two rows.");
		}
	}
	else if (codebit_count == 3)
	{
		if(mxGetM(prhs[1])!=3) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs three rows.");
		}
	}
	else if (codebit_count == 4)
	{
		if(mxGetM(prhs[1])!=4) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs four rows.");
		}
	}
	else
	{
		if(mxGetM(prhs[1])!=1) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector 2.");
		}
	}
	*/
	if(nrhs==6){
		if( !mxIsDouble(prhs[5]) ||
				mxIsComplex(prhs[5])) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input last_state must be type double.");
		}
		last_state = (int)mxGetScalar(prhs[5]);
	}
	
	
	/*  transition probabilites cell handling   */
	if(!mxIsCell(prhs[3]))
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Transitions input must be a cell.");
	if (mxGetM(prhs[3]) != 1)
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Transitions cell must be a row cell.");
	codebook_count = mxGetN(prhs[3]);
	
	ck_mem = 0;
	if (codebook_count > 20)
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Max number of codebooks is 20");
	for (ck_i=0; ck_i < codebook_count; ck_i++){	
		cell_element_ptr = mxGetCell(prhs[3], ck_i);
		/*mexPrintf("%i\n",ck_i);*/
		if (cell_element_ptr == NULL)
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:emptyCell","\tEmpty Cell\n");
		if (mxGetNumberOfDimensions(cell_element_ptr) != 2)
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Each transitions matrix needs to be 2D");
		if (mxGetM(cell_element_ptr) == 1)
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Each transitions matrix needs to be a matrix");
		/*mexPrintf("%ix%i\n",mxGetM(cell_element_ptr),mxGetN(cell_element_ptr));*/
		if ((mxGetN(cell_element_ptr)-4)>4)
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Max number of codebits is 4");
		if (codebit_count < (mxGetN(cell_element_ptr)-4))
			codebit_count = (mxGetN(cell_element_ptr)-4);
		
		/* work out how much memory we will need to malloc */
		ck_mem += mxGetM(cell_element_ptr);
		transc_list[ck_i] = mxGetM(cell_element_ptr);
	}
	
	/* check fractions vector matches with codebook count */
	if( !mxIsDouble(prhs[4]) ||
			mxIsComplex(prhs[4]) || (mxGetM(prhs[4]) != 1)) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input 'fractions' must be type double row vector.");
	}
	if  (mxGetN(prhs[4]) != codebook_count)
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input 'fractions' must be same length as number of codebooks.");
	if(mxGetM(prhs[1])!=codebit_count) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs to match max codebit count.");
	}
	
	/* codelibrary_pointers */
	for (ck_i=0; ck_i < codebook_count; ck_i++){
	
		cell_element_ptr = mxGetCell(prhs[3], ck_i);	
		codelibrary_pointers[ck_i] = malloc(sizeof(int) * mxGetM(cell_element_ptr) * TRAN_LEN);
		memset(codelibrary_pointers[ck_i],0,sizeof(int) * mxGetM(cell_element_ptr) * TRAN_LEN); /*set this to zero to make sure unused codebits are set to zero.*/
		trans_probs[ck_i] = malloc(sizeof(double) * mxGetM(cell_element_ptr));
		current_tran = codelibrary_pointers[ck_i];
		current_prob = trans_probs[ck_i];
		p = mxGetPr(cell_element_ptr);
		this_cb = mxGetN(cell_element_ptr)-4;

		for (i=0; i < 2; i++){
			for (j=0; j <  mxGetM(cell_element_ptr); j++)
				current_tran[j*TRAN_LEN + i]= ((int)*p++)-1;
		}
		for (i=2; i < 4; i++){
			for (j=0; j <  mxGetM(cell_element_ptr); j++)
				current_tran[j*TRAN_LEN + i] = (int)*p++;
		}
		if (this_cb >= 2){		
			for (i=4; i < 5; i++){
				for (j=0; j <  mxGetM(cell_element_ptr); j++)
					current_tran[j*TRAN_LEN + i] = (int)*p++;
			}			
		}
		if (this_cb >= 3){		
			for (i=5; i < 6; i++){
				for (j=0; j <  mxGetM(cell_element_ptr); j++)
					current_tran[j*TRAN_LEN + i] = (int)*p++;
			}			
		}
		if (this_cb >= 4){		
			for (i=6; i < 7; i++){
				for (j=0; j <  mxGetM(cell_element_ptr); j++)
					current_tran[j*TRAN_LEN + i] = (int)*p++;
			}			
		}
		for (i=0; i < mxGetM(cell_element_ptr); i++)
			current_prob[i] = *p++;
	}	

    
    /* get pointers */
	uncoded_in = mxGetPr(prhs[0]);
	fractions = mxGetPr(prhs[4]);
	if (codebit_count == 1)
		coded_in1 = mxGetPr(prhs[1]);
	else if (codebit_count == 2)
	{	
		coded_in1 = malloc(sizeof(double) * len);
		coded_in2 = malloc(sizeof(double) * len);
		p = mxGetPr(prhs[1]);
		for (i=0; i < len; i++){
			coded_in1[i] = *p++;
			coded_in2[i] = *p++;
		}
	}
	else if (codebit_count == 3)
	{	
		coded_in1 = malloc(sizeof(double) * len);
		coded_in2 = malloc(sizeof(double) * len);
		coded_in3 = malloc(sizeof(double) * len);
		p = mxGetPr(prhs[1]);
		for (i=0; i < len; i++){
			coded_in1[i] = *p++;
			coded_in2[i] = *p++;
			coded_in3[i] = *p++;
		}
	}
	else if (codebit_count == 4)
	{	
		coded_in1 = malloc(sizeof(double) * len);
		coded_in2 = malloc(sizeof(double) * len);
		coded_in3 = malloc(sizeof(double) * len);
		coded_in4 = malloc(sizeof(double) * len);
		p = mxGetPr(prhs[1]);
		for (i=0; i < len; i++){
			coded_in1[i] = *p++;
			coded_in2[i] = *p++;
			coded_in3[i] = *p++;
			coded_in4[i] = *p++;
		}
	}
	else
		mexErrMsgIdAndTxt("Error","Error code: 48573869730");
		
    
    
    /* create the output matrixs */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)len,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    uncoded_out = mxGetPr(plhs[0]);    
	if (codebit_count == 1){
		plhs[1] = mxCreateDoubleMatrix(1,(mwSize)len,mxREAL);
		coded_out1 = mxGetPr(plhs[1]);
	}
	else if (codebit_count == 2)
	{	
		plhs[1] = mxCreateDoubleMatrix(2,(mwSize)len,mxREAL);
		coded_out1 = malloc(sizeof(double) * len);
		coded_out2 = malloc(sizeof(double) * len);
	}
	else if (codebit_count == 3)
	{
		plhs[1] = mxCreateDoubleMatrix(3,(mwSize)len,mxREAL);
		coded_out1 = malloc(sizeof(double) * len);
		coded_out2 = malloc(sizeof(double) * len);
		coded_out3 = malloc(sizeof(double) * len);
	}
	else if (codebit_count == 4)
	{
		plhs[1] = mxCreateDoubleMatrix(4,(mwSize)len,mxREAL);
		coded_out1 = malloc(sizeof(double) * len);
		coded_out2 = malloc(sizeof(double) * len);
		coded_out3 = malloc(sizeof(double) * len);
		coded_out4 = malloc(sizeof(double) * len);
	}
	else
		mexErrMsgIdAndTxt("Error","Error code: 5812645");
    
    /* call the computational routine */
	for (ck_i = 0; ck_i < codebook_count; ck_i++){
		statec_list[ck_i] = 0;
		for (i=0;i<transc_list[ck_i];i++){
			statec_list[ck_i] = max(statec_list[ck_i], (codelibrary_pointers[ck_i])[i*TRAN_LEN + 0]+1);
			statec_list[ck_i] = max(statec_list[ck_i], (codelibrary_pointers[ck_i])[i*TRAN_LEN + 1]+1);
		}
	}
	
	bcjr_decoder_iruec(uncoded_in, coded_in1, coded_in2, coded_in3,coded_in4,
			uncoded_out, coded_out1, coded_out2, coded_out3,coded_out4, len, codelibrary_pointers,
			statec_list,transc_list,codebit_count, codebook_count, trans_probs, fractions, last_state);
	
           
    mxFree(input_buf);


	if (codebit_count == 2)
	{	
		p = mxGetPr(plhs[1]);
		for (i=0; i < len; i++){
			*p++ = coded_out1[i];
			*p++ = coded_out2[i];
		}
		free(coded_in1);
		free(coded_in2);
		free(coded_out1);
		free(coded_out2);
	}
	if (codebit_count == 3)
	{	
		p = mxGetPr(plhs[1]);
		for (i=0; i < len; i++){
			*p++ = coded_out1[i];
			*p++ = coded_out2[i];
			*p++ = coded_out3[i];
		}
		free(coded_in1);
		free(coded_in2);
		free(coded_in3);
		free(coded_out1);
		free(coded_out2);
		free(coded_out3);
	}
	if (codebit_count == 4)
	{	
		p = mxGetPr(plhs[1]);
		for (i=0; i < len; i++){
			*p++ = coded_out1[i];
			*p++ = coded_out2[i];
			*p++ = coded_out3[i];
			*p++ = coded_out4[i];
		}
		free(coded_in1);
		free(coded_in2);
		free(coded_in3);
		free(coded_in4);
		free(coded_out1);
		free(coded_out2);
		free(coded_out3);
		free(coded_out4);
	}
	for (ck_i=0; ck_i < codebook_count; ck_i++){		
		free(codelibrary_pointers[ck_i]);
		free(trans_probs[ck_i]);
	}

    /* mexPrintf("%f %f %f %f %f",maxstar(2.4,2.4),maxstar(6.7,8.1), maxstar(3,200), maxstar(-5,-8999), maxstar(-6.7,8.1));*/ 
}


