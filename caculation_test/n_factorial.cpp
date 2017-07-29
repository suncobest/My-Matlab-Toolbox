# include "mex.h"  // 使用MEX 文件必须包含的头文件


// 执行具体工作的C 函数
int factorial(int n)
{
     if(n==0)
	    return 1;
     return n*factorial(n-1);
}

// MEX 文件接口函数
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *a;
	int b;
//	if(nlhs!=1)
//		mexErrMsgTxt("One outputs required.");
		
	if(nrhs!=1)
		mexErrMsgTxt("One inputs required.");
	if (!mxIsDouble(prhs[0])||mxGetN(prhs[0])*mxGetM(prhs[0])!=1)
		mexErrMsgTxt("Input must be scalars.");
	plhs[0]= mxCreateDoubleMatrix(1, 1, mxREAL);
	a=mxGetPr(plhs[0]);
	b=*(mxGetPr(prhs[0]));
	*a=factorial(b);
}