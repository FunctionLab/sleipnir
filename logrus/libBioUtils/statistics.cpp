/*
 * Numerical Recipes in C
 * RANLIB
 * http://www.csit.fsu.edu/~burkardt
 */

#include "stdafx.h"
#include "statistics.h"
#include "measure.h"

namespace libBioUtils {

double ranf( ) {

	return ( (double)rand( ) / RAND_MAX ); }

double fsign( double dNum, double dSign ) {

	return ( ( ( ( dSign > 0 ) && ( dNum < 0 ) ) ||
		( ( dSign < 0 ) && ( dNum > 0 ) ) ) ? -dNum : dNum ); }

double CStatistics::LjungBox( const float* adX, size_t iN ) {

	return LjungBox( adX, iN, iN - 1 ); }

double CStatistics::LjungBox( const float* adX, size_t iN, size_t iH ) {
	float*	adCor;
	size_t	i, iShift;
	double	dRet;

	adCor = new float[ iN ];
	for( iShift = 0; iShift < iN; ++iShift ) {
		adCor[ iShift ] = 0;
		for( i = 0; i <= iShift; ++i )
			adCor[ iShift ] += adX[ i + iN - iShift - 1 ] * adX[ i ]; }
	for( i = 0; i < iN; ++i )
		if( adCor[ iN - 1 ] )
			adCor[ i ] /= adCor[ iN - 1 ];
		else
			adCor[ i ] = 0;
	adCor[ iN - 1 ] = 1;

	dRet = 0;
	for( i = 0; i < iH; ++i )
		dRet += adCor[ i ] * adCor[ i ] / ( iN - i - 1 );
	delete[] adCor;

	return 1 - Chi2CDF( sqrt( iN * ( iN + 2 ) * dRet ), 0, 1, iH ); }

double CStatisticsImpl::Chi2CDF( double dX, double dA, double dB, double dC ) {
	double	dRet, dX2, dY, dP2;

	if( dX <= dA )
		dRet = 0;
	else {
		dY = ( dX - dA ) / dB;
		dX2 = dY * dY / 2;
		dP2 = dC / 2;
		dRet = IncompleteGamma( dP2, dX2 ); }

	return dRet; }

double CStatisticsImpl::IncompleteGamma( double p, double x ) {
	double	a, arg, b, c, pn1, pn2, pn3, pn4, pn5, pn6, rn, value;
	double	exp_arg_min	= -88.0E+00;
	double	overflow	= 1.0E+37;
	double	plimit		= 1000.0E+00;
	double	tol			= 1.0E-07;
	double	xbig		= 1.0E+08;

	value = 0;

	if( p <= 0 )
		return 0;
	if( x <= 0 )
		return 0;

//
//  Use a normal approximation if PLIMIT < P.
//
	if( plimit < p ) {
		pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 ) + 1.0 / ( 9.0 * p ) - 1.0 );
		return Normal01CDF( pn1 ); }

//
//  Is X extremely large compared to P?
//
	if( xbig < x )
		return 1;

//
//  Use Pearson's series expansion.
//  (P is not large enough to force overflow in the log of Gamma.
//
	if( x <= 1.0 || x < p ) {
		arg = p * log ( x ) - x - GammaLog( p + 1.0 );
		c = 1.0;
		value = 1.0;
		a = p;

		for( ; ; ) {
			a = a + 1.0;
			c = c * x / a;
			value = value + c;

			if( c <= tol )
				break; }
		arg = arg + log ( value );
		value = ( exp_arg_min <= arg ) ? exp( arg ) : 0; }
	else {
//
//  Use a continued fraction expansion.
//
		arg = p * log ( x ) - x - GammaLog( p );
		a = 1.0 - p;
		b = a + x + 1.0;
		c = 0.0;
		pn1 = 1.0;
		pn2 = x;
		pn3 = x + 1.0;
		pn4 = x * b;
		value = pn3 / pn4;

		for( ; ; ) {
			a = a + 1.0;
			b = b + 2.0;
			c = c + 1.0;
			pn5 = b * pn3 - a * c * pn1;
			pn6 = b * pn4 - a * c * pn2;

			if( 0 < fabs( pn6 ) ) {
				rn = pn5 / pn6;
				if( fabs( value - rn ) <= CMath::Min( tol, tol * rn ) ) {
					arg = arg + log( value );
					return ( ( exp_arg_min <= arg ) ? ( 1 - exp( arg ) ) : 1 ); }

				value = rn; }
			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
//
//  Rescale terms in continued fraction if terms are large.
//
			if( overflow <= fabs( pn5 ) ) {
				pn1 = pn1 / overflow;
				pn2 = pn2 / overflow;
				pn3 = pn3 / overflow;
				pn4 = pn4 / overflow; } } }

	return value; }

double CStatisticsImpl::Normal01CDF( double x ) {
	double	a1	= 0.398942280444E+00;
	double	a2	= 0.399903438504E+00;
	double	a3	= 5.75885480458E+00;
	double	a4	= 29.8213557808E+00;
	double	a5	= 2.62433121679E+00;
	double	a6	= 48.6959930692E+00;
	double	a7	= 5.92885724438E+00;
	double	b0	= 0.398942280385E+00;
	double	b1	= 3.8052E-08;
	double	b2	= 1.00000615302E+00;
	double	b3	= 3.98064794E-04;
	double	b4	= 1.98615381364E+00;
	double	b5	= 0.151679116635E+00;
	double	b6	= 5.29330324926E+00;
	double	b7	= 4.8385912808E+00;
	double	b8	= 15.1508972451E+00;
	double	b9	= 0.742380924027E+00;
	double	b10	= 30.789933034E+00;
	double	b11	= 3.99019417011E+00;
	double	cdf;
	double	q;
	double	y;

//
//  |X| <= 1.28.
//
	if( fabs( x ) <= 1.28 ) {
		y = 0.5 * x * x;
		q = 0.5 - fabs( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) ); }
//
//  1.28 < |X| <= 12.7
//
	else if( fabs( x ) <= 12.7 ) {
		y = 0.5 * x * x;
		q = exp ( - y ) * b0 / ( fabs ( x ) - b1 +
			b2 / ( fabs ( x ) + b3 +
			b4 / ( fabs ( x ) - b5 +
			b6 / ( fabs ( x ) + b7 -
			b8 / ( fabs ( x ) + b9 +
			b10 / ( fabs ( x ) + b11 ) ) ) ) ) ); }
//
//  12.7 < |X|
//
	else
		q = 0;

//
//  Take account of negative X.
//
	cdf = ( x < 0 ) ? q : ( 1 - q );

	return cdf; }

double CStatisticsImpl::GammaLog( double x ) {
	double	c[ 7 ]	= {
		-1.910444077728E-03, 
		 8.4171387781295E-04, 
		-5.952379913043012E-04, 
		 7.93650793500350248E-04, 
		-2.777777777777681622553E-03, 
		 8.333333333333333331554247E-02, 
		 5.7083835261E-03 };
	double	corr;
	double	d1		= -5.772156649015328605195174E-01;
	double	d2		= 4.227843350984671393993777E-01;
	double	d4		= 1.791759469228055000094023E+00;
	double	frtbig	= 1.42E+09;
	int i;
	double	p1[ 8 ]	= {
		4.945235359296727046734888E+00, 
		2.018112620856775083915565E+02, 
		2.290838373831346393026739E+03, 
		1.131967205903380828685045E+04, 
		2.855724635671635335736389E+04, 
		3.848496228443793359990269E+04, 
		2.637748787624195437963534E+04, 
		7.225813979700288197698961E+03 };
	double	p2[ 8 ]	= {
		4.974607845568932035012064E+00, 
		5.424138599891070494101986E+02, 
		1.550693864978364947665077E+04, 
		1.847932904445632425417223E+05, 
		1.088204769468828767498470E+06, 
		3.338152967987029735917223E+06, 
		5.106661678927352456275255E+06, 
		3.074109054850539556250927E+06 };
	double	p4[ 8 ]	= {
		1.474502166059939948905062E+04, 
		2.426813369486704502836312E+06, 
		1.214755574045093227939592E+08, 
		2.663432449630976949898078E+09, 
		2.940378956634553899906876E+010,
		1.702665737765398868392998E+011,
		4.926125793377430887588120E+011, 
		5.606251856223951465078242E+011 };
	double	pnt68	= 0.6796875E+00;
	double	q1[ 8 ]	= {
		6.748212550303777196073036E+01, 
		1.113332393857199323513008E+03, 
		7.738757056935398733233834E+03, 
		2.763987074403340708898585E+04, 
		5.499310206226157329794414E+04, 
		6.161122180066002127833352E+04, 
		3.635127591501940507276287E+04, 
		8.785536302431013170870835E+03 };
	double	q2[ 8 ]	= {
		1.830328399370592604055942E+02, 
		7.765049321445005871323047E+03, 
		1.331903827966074194402448E+05, 
		1.136705821321969608938755E+06, 
		5.267964117437946917577538E+06, 
		1.346701454311101692290052E+07, 
		1.782736530353274213975932E+07, 
		9.533095591844353613395747E+06 };
	double	q4[ 8 ]	= {
		2.690530175870899333379843E+03, 
		6.393885654300092398984238E+05, 
		4.135599930241388052042842E+07, 
		1.120872109616147941376570E+09, 
		1.488613728678813811542398E+010, 
		1.016803586272438228077304E+011, 
		3.417476345507377132798597E+011, 
		4.463158187419713286462081E+011 };
	double	res;
	double	sqrtpi	= 0.9189385332046727417803297E+00;
	double	xbig	= 4.08E+36;
	double	xden;
	double	xm1;
	double	xm2;
	double	xm4;
	double	xnum;
	double	xsq;

//
//	Return immediately if the argument is out of range.
//
	if( x <= 0 || xbig < x )
		return HUGE_VAL;

	if( x <= EpsilonDouble( ) )
		res = -log( x );
	else if( x <= 1.5 ) {
		if ( x < pnt68 ) {
			corr = -log( x );
			xm1 = x; }
		else {
			corr = 0;
			xm1 = ( x - 0.5 ) - 0.5; }

		if ( x <= 0.5 || pnt68 <= x ) {
			xden = 1.0;
			xnum = 0.0;

			for( i = 0; i < 8; i++ ) {
				xnum = xnum * xm1 + p1[i];
				xden = xden * xm1 + q1[i]; }
			res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) ); }
		else {
			xm2 = ( x - 0.5 ) - 0.5;
			xden = 1.0;
			xnum = 0.0;
			for( i = 0; i < 8; i++ ) {
				xnum = xnum * xm2 + p2[i];
				xden = xden * xm2 + q2[i]; }
			res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) ); } }
	else if( x <= 4.0 ) {
		xm2 = x - 2.0;
		xden = 1.0;
		xnum = 0.0;
		for( i = 0; i < 8; i++ ) {
			xnum = xnum * xm2 + p2[i];
			xden = xden * xm2 + q2[i]; }

		res = xm2 * ( d2 + xm2 * ( xnum / xden ) ); }
	else if( x <= 12.0 ) {
		xm4 = x - 4.0;
		xden = - 1.0;
		xnum = 0.0;
		for( i = 0; i < 8; i++ ) {
			xnum = xnum * xm4 + p4[i];
			xden = xden * xm4 + q4[i]; }
		res = d4 + xm4 * ( xnum / xden ); }
	else {
		res = 0.0;

		if( x <= frtbig ) {
			res = c[6];
			xsq = x * x;

			for( i = 0; i < 6; i++ )
				res = res / xsq + c[i]; }

		res = res / x;
		corr = log ( x );
		res = res + sqrtpi - 0.5 * corr;
		res = res + x * ( corr - 1.0 ); }

	return res; }

double CStatisticsImpl::EpsilonDouble( ) {
	double	dRet;

	dRet = 1;
	while( 1 < (double)( 1 + dRet ) )
		dRet = dRet / 2;

	return ( 2 * dRet ); }

double CStatistics::SampleChi2( size_t iDF ) {

	return ( 2 * SampleGamma( 1, iDF / 2 ) ); }

double CStatistics::SampleGammaLogStandard( double dXX ) {
	static const double	c_adCof[]	= { 76.18009172947146, -86.50532032941677, 24.01409824083091,
		-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
	double	dX, dY, dTmp, dSer;
	size_t	j;

	dX = dY = dXX;
	dTmp = dX + 5.5;
	dTmp -= ( dX + 0.5 ) * log( dTmp );
	dSer = 1.000000000190015;
	for( j = 0; j <= 5; ++j )
		dSer += c_adCof[ j ] / ++dY;

	return ( -dTmp + log( 2.5066282746310005 * dSer / dX ) ); }

double CStatistics::SampleGamma( double dLocation, double dShape ) {

	return ( SampleGammaStandard( dShape ) / dLocation ); }

double CStatistics::SampleGammaStandard( double dShape ) {
static double q1 = 4.166669E-2;
static double q2 = 2.083148E-2;
static double q3 = 8.01191E-3;
static double q4 = 1.44121E-3;
static double q5 = -7.388E-5;
static double q6 = 2.4511E-4;
static double q7 = 2.424E-4;
static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 4.07753E-2;
static double e5 = 1.0293E-2;
static double aa = 0.0;
static double aaa = 0.0;
static double sqrt32 = 5.656854;
static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
double a = dShape;
    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = SampleNormalStandard();
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = ranf();
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s+0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = SampleExponentialStandard();
    u = ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0+0.3678794*a;
S130:
    p = b*ranf();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(SampleExponentialStandard() < sgamma) goto S130;
    return sgamma;
S140:
    sgamma = -log((b-p)/ a);
    if(SampleExponentialStandard() < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;
}

double CStatistics::SampleNormalStandard( ) {
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static double snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

double CStatistics::SampleExponentialStandard( ) {
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
};
static long i;
static double sexpo,a,u,ustar,umin;
static double *q1 = q;
    a = 0.0;
    u = ranf();
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
    if(u <= 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
S70:
    ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    sexpo = a+umin**q1;
    return sexpo;
}

double CStatistics::PValueLognormal( double dSample, double dMean, double dVariance ) {
	double	dCDF;

	dCDF = LognormalCDF( dSample, dMean, dVariance );
	if( dSample > exp( dMean ) )
		dCDF = 1 - dCDF;

	return ( 2 * dCDF ); }

double CStatistics::LognormalCDF( double dX, double dMean, double dVariance ) {

	if( dX <= 0 )
		return 0;

	return NormalCDF( log( dX ), dMean, dVariance ); }

double CStatistics::NormalCDF( double dX, double dMean, double dVariance ) {

	return Normal01CDF( ( dX - dMean ) / dVariance ); }

double CStatistics::TTest( double dMeanOne, double dVarOne, size_t iNOne, double dMeanTwo,
	double dVarTwo, size_t iNTwo ) {
	size_t	iDegFree;
	double	dPoolVar, dT;

	iDegFree = iNOne + iNTwo - 2;
	dPoolVar = ( ( ( iNOne - 1 ) * dVarOne ) + ( ( iNTwo - 1 ) * dVarTwo ) ) / iDegFree;
	dT = ( dMeanOne - dMeanTwo ) / sqrt( dPoolVar * ( ( 1.0 / iNOne ) + ( 1.0 / iNTwo ) ) );

	return IncompleteBeta( 0.5 * iDegFree, 0.5, iDegFree / ( iDegFree + ( dT * dT ) ) ); }

double CStatisticsImpl::IncompleteBeta( double dA, double dB, double dX ) {
	double	dBT;

	if( ( dX < 0 ) || ( dX > 1 ) )
		return -1;
	dBT = ( ( dX == 0 ) || ( dX == 1 ) ) ? 0 :
		exp( GammaLog( dA + dB ) - GammaLog( dA ) - GammaLog( dB ) + ( dA * log( dX ) ) +
		( dB * log( 1 - dX ) ) );

	return ( ( dX < ( ( dA + 1 ) / ( dA + dB + 2 ) ) ) ?
		( dBT * IncompleteBetaCF( dA, dB, dX ) / dA ) :
		( 1 - ( dBT * IncompleteBetaCF( dB, dA, 1 - dX ) / dB ) ) ); }

double CStatisticsImpl::IncompleteBetaCF(double a, double b, double x)
{
static const double	c_dFPMin	= 1e-30;
static const double	c_iMaxIt	= 100;
static const double	c_dEPS		= 3e-7;
int m,m2;
double aa,c,d,del,h,qab,qam,qap;
qab=a+b;
qap=a+1.0;
qam=a-1.0;
c=1.0;
d=1.0-qab*x/qap;
if (fabs(d) < c_dFPMin) d=c_dFPMin;
d=1.0/d;
h=d;
for (m=1;m<=c_iMaxIt;m++) {
m2=2*m;
aa=m*(b-m)*x/((qam+m2)*(a+m2));
d=1.0+aa*d;
if (fabs(d) < c_dFPMin) d=c_dFPMin;
c=1.0+aa/c;
if (fabs(c) < c_dFPMin) c=c_dFPMin;
d=1.0/d;
h *= d*c;
aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
d=1.0+aa*d;
if (fabs(d) < c_dFPMin) d=c_dFPMin;
c=1.0+aa/c;
if (fabs(c) < c_dFPMin) c=c_dFPMin;
d=1.0/d;
del=d*c;
h *= del;
if (fabs(del-1.0) < c_dEPS) break;
}
if (m > c_iMaxIt) return -1;
return h;
}

double CStatistics::BinomialCDF( size_t iObservations, size_t iSample, double dProb ) {
	double	d;

	d = ( iObservations - ( iSample * dProb ) ) / sqrt( iSample * dProb * ( 1 - dProb ) );
	return ( 1 - NormalCDF( d, 0, 1 ) ); }

double CStatistics::HypergeometricPDF( size_t iHitsOne, size_t iSizeOne, size_t iHitsTwo,
	size_t iSizeTwo ) {

	return exp( LogFact( iSizeTwo - iHitsTwo ) + LogFact( iHitsTwo ) + LogFact( iSizeOne ) +
		LogFact( iSizeTwo - iSizeOne ) - LogFact( iHitsOne ) -
		LogFact( iHitsTwo - iHitsOne ) - LogFact( iSizeTwo - iHitsTwo + iHitsOne - iSizeOne ) -
		LogFact( iSizeOne - iHitsOne ) - LogFact( iSizeTwo ) ); }

double CStatistics::HypergeometricCDF( size_t iHitsOne, size_t iSizeOne, size_t iHitsTwo,
	size_t iSizeTwo ) {
	size_t	i, iHits;
	double	dRet;

	dRet = 0;
	iHits = ( iHitsTwo < iSizeOne ) ? iHitsTwo : iSizeOne;
	for( i = iHitsOne; i <= iHits; ++i )
		dRet += HypergeometricPDF( i, iSizeOne, iHitsTwo, iSizeTwo );

	return ( ( dRet > 1 ) ? 1 : dRet ); }

double CStatistics::Precision( size_t iTP, size_t iFP, size_t iTN, size_t iFN ) {

	return ( (double)iTP / ( iTP + iFP ) ); }

double CStatistics::Recall( size_t iTP, size_t iFP, size_t iTN, size_t iFN ) {

	return ( (double)iTP / ( iTP + iFN ) ); }

double CStatistics::FScore( size_t iTP, size_t iFP, size_t iTN, size_t iFN, double dBeta ) {
	double	dP, dR;

	dP = Precision( iTP, iFP, iTN, iFN );
	dR = Recall( iTP, iFP, iTN, iFN );
	dBeta *= dBeta;

	return ( ( dBeta + 1 ) * dP * dR / ( ( dBeta * dP ) + dR ) ); }

}
