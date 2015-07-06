#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/* note, the second half of the table must only be the second transition to each state */
/* last column is transition probablity */
static int transurc8[16][6] = {
    {0,          0,          0,          0,          0,          0},
    {1,          4,          0,          1,          0,          0},
    {2,          5,          0,          1,          0,          0},
    {3,          1,          0,          0,          0,          0},
    {4,          6,          0,          1,          0,          0},
    {5,          2,          0,          0,          0,          0},
    {6,          3,          0,          0,          0,          0},
    {7,          7,          0,          1,          0,          0},
    {0,          4,          1,          1,          0,          0},
    {1,          0,          1,          0,          0,          0},
    {2,          1,          1,          0,          0,          0},
    {3,          5,          1,          1,          0,          0},
    {4,          2,          1,          0,          0,          0},
    {5,          6,          1,          1,          0,          0},
    {6,          7,          1,          1,          0,          0},
    {7,          3,          1,          0,          0,          0}
};

static int translte8[16][6] =  {
    {1-1,        1-1,        0,          0,          0,          0},
    {2-1,        5-1,        0,          0,          0,          0},
    {3-1,        6-1,        0,          1,          0,          0},
    {4-1,        2-1,        0,          1,          0,          0},
    {5-1,        3-1,        0,          1,          0,          0},
    {6-1,        7-1,        0,          1,          0,          0},
    {7-1,        8-1,        0,          0,          0,          0},
    {8-1,        4-1,        0,          0,          0,          0},
    {1-1,        5-1,        1,          1,          0,          0},
    {2-1,        1-1,        1,          1,          0,          0},
    {3-1,        2-1,        1,          0,          0,          0},
    {4-1,        6-1,        1,          0,          0,          0},
    {5-1,        7-1,        1,          0,          0,          0},
    {6-1,        3-1,        1,          0,          0,          0},
    {7-1,        4-1,        1,          1,          0,          0},
    {8-1,        8-1,        1,          1,          0,          0}
};

static int transurc4[16][6]   =  {
    {1-1,          1-1,          0,          0,          0,          0},
    {1-1,          3-1,          1,          1,          0,          0},
    {2-1,          3-1,          0,          1,          0,          0},
    {2-1,          1-1,          1,          0,          0,          0},
    {3-1,          4-1,          0,          1,          0,          0},
    {3-1,          2-1,          1,          0,          0,          0},
    {4-1,          2-1,          0,          0,          0,          0},
    {4-1,          4-1,          1,          1,          0,          0}
};
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

/* reg_trellis - does each state have 2 transitions ending there? */
void bcjr_decoder(double *uncoded_in, double *coded_in1, double *coded_in2, double *coded_in3,
        double *uncoded_out, double *coded_out1, double *coded_out2,double *coded_out3,
        int len, int trans[][6], int state_count, int trans_count, int codeword_count,
        int reg_trellis, double *trans_prob, int last_state)
{
    
    int i,j;
    int transc = trans_count;
    int set0,set1;
    double p1,p0;
    int a,b;
    double **gammas;
    double **alphas;
    double **betas;
    double temp[16];
    int states = state_count;

    /* calculate gammas */ 
    gammas = (double**) malloc (transc*sizeof(double *));
    for (i = 0; i < transc; i++)
        gammas[i] = (double*) malloc(len*sizeof(double));
    
    for (i = 0; i < transc; i++)
    {
        if (trans[i][2] && trans[i][3]){
            for (j = 0; j < len; j++){
                gammas[i][j] = uncoded_in[j] + coded_in1[j];
            }
        }
        else if (trans[i][2]){
            for (j = 0; j < len; j++)
                gammas[i][j] = uncoded_in[j];
        }
        else if (trans[i][3]){
            for (j = 0; j < len; j++)
                gammas[i][j] = coded_in1[j];
        }
        else{
            for (j = 0; j < len; j++)
                gammas[i][j] = 0;
        }
    }
	
	if (codeword_count >= 2)
	{
		for (i = 0; i < transc; i++)
		{
			if (trans[i][4]){
				for (j = 0; j < len; j++)
					gammas[i][j] += coded_in2[j];
			}
		}
	}
	if (codeword_count >= 3)
	{
		for (i = 0; i < transc; i++)
		{
			if (trans[i][5]){
				for (j = 0; j < len; j++)
					gammas[i][j] += coded_in3[j];
			}
		}
	}
	

    if (trans_prob != (double*)0)
    {
        for (i = 0; i < transc; i++)
        {
            for (j=0; j < len; j++)
                gammas[i][j] += trans_prob[i];
        }
    }
    
    
    
    /* set and initialise memory */
    alphas = (double**) malloc (states*sizeof(double *));
    for (i = 0; i < states; i++)
        alphas[i] = (double*) malloc(len*sizeof(double));
    betas = (double**) malloc (states*sizeof(double *));
    for (i = 0; i < states; i++)
        betas[i] = (double*) malloc(len*sizeof(double));
    
    if (reg_trellis == 0)
    {
        for (i = 0; i < states; i++) {
            for (j=0; j< len; j++)
                alphas[i][j] = -9000;
        }
        for (i = 0; i < states; i++){
            for (j=0; j< len; j++)
                betas[i][j] = -9000;
        }
    }
    
    
    
    /* forward recursion (alphas) */
    alphas[0][0] = 0;        /* first state */
    for (i = 1; i < states; i++)
        alphas[i][0] = -9000;
    if (reg_trellis == 0){
        for (i = 1; i < len; i++)
        {
            for (j = 0; j < transc; j++)
                alphas[trans[j][1]][i] = maxstar(  alphas[trans[j][1]][i],  alphas[trans[j][0]][i-1] + gammas[j][i-1]  );
        }
    }
    else
    {
        for (i = 1; i < len; i++)
        {
            for (j = 0; j < states; j++)
                temp[trans[j][1]] = alphas[trans[j][0]][i-1] + gammas[j][i-1];
            for (j = states; j < transc; j++)
                alphas[trans[j][1]][i] = maxstar( temp[trans[j][1]],  alphas[trans[j][0]][i-1] + gammas[j][i-1]  );
        }
    }
    
    
    /* backwards recursion (betas) */
    /* double betas[8][len]; */
    if (last_state < 1 || last_state > states){
        for (i = 0; i < states; i++)
            betas[i][len-1] = 0;     /* end state unknown */
    }
    else
    {
        for (i = 0; i < states; i++)
            betas[i][len-1] = -9000;
        betas[last_state-1][len-1] = 0;
    }
    
    if (reg_trellis == 0){
        for (i = len-2; i >= 0; i--)
        {
            for (j = 0; j < transc; j++)
                betas[trans[j][0]][i] = maxstar( betas[trans[j][0]][i],  betas[trans[j][1]][i+1] + gammas[j][i+1]  );
        }
    }
    else
    {
        for (i = len-2; i >= 0; i--)
        {
            for (j = 0; j < states; j++)
                temp[trans[j][0]] = betas[trans[j][1]][i+1] + gammas[j][i+1];
            for (j = states; j < transc; j++)
                betas[trans[j][0]][i] = maxstar( temp[trans[j][0]],  betas[trans[j][1]][i+1] + gammas[j][i+1]  );
        }
    }
    
    /*
     * if (trans_prob != (double*)0)
     * {
     * mexPrintf("trans probs\n");
     * for (i=0;i<transc;i++)
     * {
     * mexPrintf("%f ",trans_prob[i]);
     *
     *
     * }
     * }
     * mexPrintf("GAMMAS\n");
     * for (i=0;i<transc;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",gammas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("ALPHAS\n");
     * for (i=0;i<states;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",alphas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("\n");
     * mexPrintf("\n");
     * mexPrintf("BETAS\n");
     * for (i=0;i<states;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",betas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("\n");
     */
	 
    /* deltas */
    /* reuse gammas memory */
    for (i = 0; i < transc; i++)
    {
        a = trans[i][0];
        b = trans[i][1];
        for (j = 0; j < len; j++)
            gammas[i][j] += alphas[a][j] + betas[b][j];
    }
    
    
    /* extrinisic uncoded llr */
    for (j = 0; j < len; j++)
    {
        set0 = 0;
        set1 = 0;
        
        for (i = 0; i < transc; i++)
        {
            if (trans[i][2]){
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
    for (j = 0; j < len; j++)
    {
        set0 = 0;
        set1 = 0;
        
        for (i = 0; i < transc; i++)
        {
            if (trans[i][3]){
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
	if (codeword_count >= 2)
	{
		for (j = 0; j < len; j++)
		{
			set0 = 0;
			set1 = 0;
			
			for (i = 0; i < transc; i++)
			{
				if (trans[i][4]){
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
	if (codeword_count >= 3)
	{
		for (j = 0; j < len; j++)
		{
			set0 = 0;
			set1 = 0;
			
			for (i = 0; i < transc; i++)
			{
				if (trans[i][5]){
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
    
    for (i = 0; i < transc; i++){
        free(gammas[i]);
    }
    free(gammas);
    for (i = 0; i < states; i++){
        free(alphas[i]);
    }
    free(alphas);
    for (i = 0; i < states; i++){
        free(betas[i]);
    }
    free(betas);
    
    
}





/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    double *uncoded_in;
    double *coded_in1;
	double *coded_in2;
	double *coded_in3;
    double *uncoded_out;
    double *coded_out1;
    double *coded_out2;
	double *coded_out3;
    int len,i,j,states,last_state;
    char *input_buf;
    int blank_transitions[600][6];
    double *trans_probs;
    double *p;
	int codeword_count = 1;

    last_state = -1;
    
    /* check for proper number of arguments */
    if(nrhs<3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",">=Three inputs required. [uncoded_out, coded_out] = bcjr_decoder(uncoded_in, coded_in, {'urc8','urc4','lte8','uec'})");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required. [uncoded_out, coded_out] = bcjr_decoder(uncoded_in, coded_in, {'urc8','urc4','lte8','uec'})");
    }
    
    
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
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
    
    len = mxGetN(prhs[0]);
    /*if (len < 3)
    *    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    *if (mxGetN(prhs[1]) < 3)
    *    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector."); */
    if (mxGetN(prhs[0]) != mxGetN(prhs[1]))
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be same dimentions.");
    
    
    input_buf = mxArrayToString(prhs[2]);
    
    if ((strncmp(input_buf,"uec",5)==0) || (strncmp(input_buf,"new",5)==0))
    {
        
        if(nrhs!=5 && nrhs!=4)
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Four inputs required. [uncoded_out, coded_out] = bcjr_decoder(uncoded_in, coded_in, 'uec', transitions{, last_state}). States start at 1, set last_state < 0 for unknown");
        if( !mxIsDouble(prhs[3]) ||
                mxIsComplex(prhs[3])) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
        }
        
        if(mxGetM(prhs[3])==1 || mxGetM(prhs[3])>599) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword), Tx6 (two codewords) or Tx7 (two codewords).  Max transitions = 600");
        }
        if(mxGetN(prhs[3])!=5 && mxGetN(prhs[3])!=6 && mxGetN(prhs[3])!=7) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword), Tx6 (two codewords) or Tx7 (two codewords).  Max transitions = 600");
        }
		if(mxGetN(prhs[3])==6){
			codeword_count = 2;
		}
		if(mxGetN(prhs[3])==7){
			codeword_count = 3;
		}
		if (codeword_count == 2)
		{		
			if(mxGetM(prhs[1])!=2) {
				mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs two rows.");
			}
		}
		else if (codeword_count == 3)
		{
			if(mxGetM(prhs[1])!=3) {
				mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","Input coded LLRs needs threee rows.");
			}
		}
		else
		{
			if(mxGetM(prhs[1])!=1) {
				mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
			}
		}
        
        if(nrhs==5){
            if( !mxIsDouble(prhs[4]) ||
                    mxIsComplex(prhs[4])) {
                mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input last_state must be type double.");
            }
            last_state = mxGetScalar(prhs[4]);
        }
        
        trans_probs = malloc(mxGetM(prhs[3])*sizeof(double *));
        p = mxGetPr(prhs[3]);
         
        for (i=0; i < 2; i++){
            for (j=0; j <  mxGetM(prhs[3]); j++)
                blank_transitions[j][i] = ((int)*p++)-1;
        }
        for (i=2; i < 4; i++){
            for (j=0; j <  mxGetM(prhs[3]); j++)
                blank_transitions[j][i] = (int)*p++;
        }
		if (codeword_count >= 2){		
			for (i=4; i < 5; i++){
				for (j=0; j <  mxGetM(prhs[3]); j++)
					blank_transitions[j][i] = (int)*p++;
			}			
		}
		if (codeword_count >= 3){		
			for (i=5; i < 6; i++){
				for (j=0; j <  mxGetM(prhs[3]); j++)
					blank_transitions[j][i] = (int)*p++;
			}			
		}
        for (i=0; i < mxGetM(prhs[3]); i++)
            trans_probs[i] = *p++;
        
    }
	else
	{
		if(mxGetM(prhs[1])!=1) {
			mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
		}
	}
    
    
    /* get pointers */
	uncoded_in = mxGetPr(prhs[0]);
	if (codeword_count == 1)
		coded_in1 = mxGetPr(prhs[1]);
	else if (codeword_count == 2)
	{	
		coded_in1 = malloc(sizeof(double) * len);
		coded_in2 = malloc(sizeof(double) * len);
		p = mxGetPr(prhs[1]);
		for (i=0; i < len; i++){
			coded_in1[i] = *p++;
			coded_in2[i] = *p++;
		}
	}
	else if (codeword_count == 3)
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
	else
		mexErrMsgIdAndTxt("Error","Error code: 48573869730");
		
    
    
    
    /* create the output matrixs */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)len,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    uncoded_out = mxGetPr(plhs[0]);    
	if (codeword_count == 1){
		plhs[1] = mxCreateDoubleMatrix(1,(mwSize)len,mxREAL);
		coded_out1 = mxGetPr(plhs[1]);
	}
	else if (codeword_count == 2)
	{	
		plhs[1] = mxCreateDoubleMatrix(2,(mwSize)len,mxREAL);
		coded_out1 = malloc(sizeof(double) * len);
		coded_out2 = malloc(sizeof(double) * len);
	}
	else if (codeword_count == 3)
	{
		plhs[1] = mxCreateDoubleMatrix(3,(mwSize)len,mxREAL);
		coded_out1 = malloc(sizeof(double) * len);
		coded_out2 = malloc(sizeof(double) * len);
		coded_out3 = malloc(sizeof(double) * len);
	}
	else
		mexErrMsgIdAndTxt("Error","Error code: 5812645");
    
    /* call the computational routine */
    if (strncmp(input_buf,"urc8",5)==0)
    {
        bcjr_decoder(uncoded_in, coded_in1, 0, 0,
                uncoded_out, coded_out1, 0, 0, len,transurc8,8,16,1,1,(double*)0 , last_state);
    }
    else if (strncmp(input_buf,"urc4",5)==0)
    {
        bcjr_decoder(uncoded_in, coded_in1, 0, 0,
                uncoded_out, coded_out1, 0, 0, len,transurc4,4,8,1,1,(double*)0 , last_state);
    }
    else if (strncmp(input_buf,"lte8",5)==0)
    {
        bcjr_decoder(uncoded_in, coded_in1, 0, 0,
                uncoded_out, coded_out1, 0, 0, len,translte8,8,16,1,1,(double*)0 , last_state);
    }
    else if ((strncmp(input_buf,"uec",5)==0)  || (strncmp(input_buf,"new",5)==0))
    {
        states = 0;
        for (i=0;i<mxGetM(prhs[3]);i++)
        {
            states = max(states, blank_transitions[i][0]);
            states = max(states, blank_transitions[i][1]);
        }
        
        bcjr_decoder(uncoded_in, coded_in1, coded_in2, coded_in3,
                uncoded_out, coded_out1, coded_out2, coded_out3, len,blank_transitions,
                states+1,mxGetM(prhs[3]),codeword_count,0,trans_probs, last_state);
        
        free(trans_probs);
        
        
    }
    else
    {
        mexErrMsgIdAndTxt("MyToolbox:inputError","Select from 'urc8','urc4' or 'lte8'");
    }
    
    mxFree(input_buf);

	if (codeword_count == 2)
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
	if (codeword_count == 3)
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

    /* mexPrintf("%f %f %f %f %f",maxstar(2.4,2.4),maxstar(6.7,8.1), maxstar(3,200), maxstar(-5,-8999), maxstar(-6.7,8.1));*/ 
}


