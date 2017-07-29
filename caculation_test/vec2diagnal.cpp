#include "mex.h"  // 使用MEX 文件必须包含的头文件

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m,n;
	if(nlhs>1)
		mexErrMsgTxt("Number of outputs must be 1 or 0!");		
	if(nrhs!=1)
		mexErrMsgTxt("One inputs required.");
		
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	if (! (mxIsDouble(prhs[0]) && (m==1 || n==1)) )
		mexErrMsgTxt("Input must be a vector!");
		
	n = m>=n ? m:n;
	plhs[0]= mxCreateDoubleMatrix(n, n, mxREAL);
	
	double *M,*v;
	M = mxGetPr(plhs[0]);
	v = mxGetPr(prhs[0]);
	
	for (m=0; m<n*n; ++m)
	{
	if (!(m %(n+1)))
	{
	*M= *v;
	M++;
	v++;
	}
	else
	{
	*M=0;
	M++;
	}
	}
}