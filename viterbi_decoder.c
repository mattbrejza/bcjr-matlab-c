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
    {1,          1,          0,          0,          0,          0},
    {1,          3,          1,          1,          0,          0},
    {2,          3,          0,          1,          0,          0},
    {2,          1,          1,          0,          0,          0},
    {3,          4,          0,          1,          0,          0},
    {3,          2,          1,          0,          0,          0},
    {4,          2,          0,          0,          0,          0},
    {4,          4,          1,          1,          0,          0}
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
void viterbi_decoder(double *uncoded_in, 
        double *uncoded_out,
        int len, int trans[][6], int state_count, int trans_count,
        double *trans_prob)
{
    
    int i,j;
    int transc = trans_count;
    int set0,set1;
    double p1,p0;
    int a,b;
    double **gammas;
    double **alphas;
    int states = state_count;
	int state,best_state;
	double best_alpha;

    /* calculate gammas */ 
    gammas = (double**) malloc (transc*sizeof(double *));
    for (i = 0; i < transc; i++)
        gammas[i] = (double*) malloc(len*sizeof(double));
    
    for (i = 0; i < transc; i++)
    {        
        if (trans[i][2]){
            for (j = 0; j < len; j++)
                gammas[i][j] = uncoded_in[j];
        }else{
            for (j = 0; j < len; j++)
                gammas[i][j] = 0;
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
    

	for (i = 0; i < states; i++) {
		for (j=0; j< len; j++)
			alphas[i][j] = -9000;
	}
    
    
    
    /* forward recursion (alphas) */
    alphas[0][0] = 0;        /* first state */
    for (i = 1; i < states; i++)
        alphas[i][0] = -9000;
	for (i = 1; i < len; i++)
	{
		for (j = 0; j < transc; j++)
			alphas[trans[j][1]][i] = maxstar(  alphas[trans[j][1]][i],  alphas[trans[j][0]][i-1] + gammas[j][i-1]  );
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
	 
	 
	/* viterbi time */
	state = 0;	
	for (i = len-1; i>=0; i--){
		best_alpha = -900000000;
		for (j = 0; j < transc; j++){
			if ((trans[j][1] == state) && (alphas[trans[j][0]][i] >  best_alpha)){
				best_alpha = alphas[trans[j][0]][i];
				uncoded_out[i] = trans[j][2];
				best_state = trans[j][0];
			}			
		}
		state = best_state;
	}
	    
    for (i = 0; i < transc; i++){
        free(gammas[i]);
    }
    free(gammas);
    for (i = 0; i < states; i++){
        free(alphas[i]);
    }
    free(alphas);
}





/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    double *uncoded_in;
    double *uncoded_out;
    int len,i,j,states;
    int blank_transitions[600][6];
    double *trans_probs;
    double *p;
	int codeword_count = 1;

    
    /* check for proper number of arguments */
    if(nrhs<2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",">=Two inputs required. [uncoded_out] = viterbi_decoder(uncoded_in, transitions)");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required. [uncoded_out] = viterbi_decoder(uncoded_in, transitions)");
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
    if (len < 3)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");


    
    
    {    
        if(nrhs!=2)
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required. [uncoded_out, coded_out] = viterbi_decoder(uncoded_in, transitions).");
        if( !mxIsDouble(prhs[1]) ||
                mxIsComplex(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
        }
        
        if(mxGetM(prhs[1])==1 || mxGetM(prhs[1])>599) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword).  Max transitions = 600");
        }
        if(mxGetN(prhs[1])!=5 ) {
            mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notMatrix","transitions must be Tx5 (one codeword).  Max transitions = 600");
        }
		
        
        trans_probs = malloc(mxGetM(prhs[1])*sizeof(double *));
        p = mxGetPr(prhs[1]);
         
        for (i=0; i < 2; i++){
            for (j=0; j <  mxGetM(prhs[1]); j++)
                blank_transitions[j][i] = ((int)*p++)-1;
        }
        for (i=2; i < 4; i++){
            for (j=0; j <  mxGetM(prhs[1]); j++)
                blank_transitions[j][i] = (int)*p++;
        }		
        for (i=0; i < mxGetM(prhs[1]); i++)
            trans_probs[i] = *p++;
        
    }
    
    
    /* get pointers */
	uncoded_in = mxGetPr(prhs[0]);
		
    
    
    
    /* create the output matrixs */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)len,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    uncoded_out = mxGetPr(plhs[0]);    
	
    
    /* call the computational routine */
	states = 0;
	for (i=0;i<mxGetM(prhs[1]);i++)
	{
		states = max(states, blank_transitions[i][0]);
		states = max(states, blank_transitions[i][1]);
	}
	
	viterbi_decoder(uncoded_in, 
			uncoded_out, len,blank_transitions,
			states+1,mxGetM(prhs[1]),trans_probs);
	
	free(trans_probs);


    /* mexPrintf("%f %f %f %f %f",maxstar(2.4,2.4),maxstar(6.7,8.1), maxstar(3,200), maxstar(-5,-8999), maxstar(-6.7,8.1));*/ 
}


