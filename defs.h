//
//Defined macros
#define MakeVector(x, y, z, v) (v)[0]=(x),(v)[1]=(y),(v)[2]=(z)
#define VecNegate(a) (a)[0]=(-(a)[0]); (a)[1]=(-(a)[1]); (a)[2]=(-(a)[2])
#define VecDot(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define VecLen(a) ((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])
#define VecLenSq(a) sqrtf((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])
#define VecCopy(a,b) (b)[0]=(a)[0];(b)[1]=(a)[1];(b)[2]=(a)[2]
#define VecAdd(a,b,c) \
	(c)[0]=(a)[0]+(b)[0]; \
	(c)[1]=(a)[1]+(b)[1]; \
	(c)[2]=(a)[2]+(b)[2]
#define VecSub(a,b,c) \
	(c)[0]=(a)[0]-(b)[0]; \
	(c)[1]=(a)[1]-(b)[1]; \
	(c)[2]=(a)[2]-(b)[2]
#define VecS(A,a,b) (b)[0]=(A)*(a)[0]; (b)[1]=(A)*(a)[1]; (b)[2]=(A)*(a)[2]
#define VecAddS(A,a,b,c) \
	(c)[0]=(A)*(a)[0]+(b)[0]; \
	(c)[1]=(A)*(a)[1]+(b)[1]; \
	(c)[2]=(A)*(a)[2]+(b)[2]
#define VecCross(a,b,c) \
	(c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; \
	(c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; \
	(c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]
#define VecZero(v) (v)[0]=0.0;(v)[1]=0.0;v[2]=0.0
#define ABS(a) (((a)<0)?(-(a)):(a))
#define SGN(a) (((a)<0)?(-1):(1))
#define SGND(a) (((a)<0.0)?(-1.0):(1.0))
#define rnd() (float)rand() / (float)RAND_MAX
#define PI 3.14159265358979323846264338327950
#define MY_INFINITY 1000.0


