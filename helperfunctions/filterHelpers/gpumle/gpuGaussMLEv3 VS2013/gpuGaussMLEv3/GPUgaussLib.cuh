/*!
 * \file GPUgaussLib.cuh
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief This contains the definitions for all the cuda device functions.  It is only 
 * included in the GPUgaussMLEv2.cu file to prevent multiple definitions of functions.
 */
#include "GPUgaussLib.h"
#include "definitions.h"
#include "PropLib.cuh"
//*******************************************************************************************
__device__ void kernel_CalcLLRProp(float Cr, float Itheta, float Div, float * LR) {
	/*!
	* \brief returns propabilities corresponding to LLR
	* \param Diag
	* \param thetaH1
	* \param Div
	* \param LR
	*/
	LR[0] = Div;
	LR[1] = pow(Itheta, 2) / (Cr + 1e-5);
	LR[2] = 2 * (1 - normcdf(sqrt(max(Div, 0.0)), 0, 1));
	LR[3] = (1 - normcdf(sqrt(max(Div, 0.0)), sqrt(LR[1]), 1)) + (1 - normcdf(sqrt(max(Div, 0.0)), -sqrt(LR[1]), 1));
	LR[4] = 2 * normpdf(sqrt(max(Div, 0.0)), 0, 1);
	LR[5] = normpdf(sqrt(max(Div, 0.0)), sqrt(LR[1]), 1) + normpdf(sqrt(max(Div, 0.0)), -sqrt(LR[1]), 1);
	/*LR[5]=LR[5]/(LR[4]+LR[5]+1e-4);
	LR[4]=LR[4]/(LR[4]+LR[5]+1e-4);*/
}
//*******************************************************************************************
__device__ float kernel_IntGauss1D(const int ii, const float x, const float sigma) {
	/*! 
	 * \brief /f$ \frac{1}{2} /f$
	 * \param ii ???
	 * \param x ???
	 * \param sigma sigma value of the PSF
	 * \return float
	 */
	const float norm=1.0f/2.0f/sigma/sigma;
    return 1.0f/2.0f*(erf((ii-x+0.5f)*sqrt(norm))-erf((ii-x-0.5f)*sqrt(norm)));
}

//*******************************************************************************************
__device__ void kernel_DxIntGaus1D(const int ii, const float x, const float sigma, float *dx, float *dxx) {
	/*!
	* \brief compute the derivative of the 1D gaussian
	* \param ii ???
	* \param x ???
	* \param sigma ???
	* \param N ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float a, b;
	a = exp(-1.0f / 2.0f*pow(((ii + 0.5f - x) / sigma), 2.0f));
	b = exp(-1.0f / 2.0f*pow((ii - 0.5f - x) / sigma, 2.0f));

	*dx = -1.0f / sqrt(2.0f*pi) / sigma*(a - b);

	*dxx = -1.0f / sqrt(2.0f*pi) / pow(sigma, 3)*((ii + 0.5f - x)*a - (ii - 0.5f - x)*b);
}

__device__ void kernel_DerivativeIntGauss1DSigma(const int ii, const float x,
	const float Sx, const float N, const float PSFy, float *dudt, float *d2udt2) {
	/*!
	* \brief compute the derivative of the 1D gaussian
	* \param ii ???
	* \param x ???
	* \param Sx ???
	* \param N ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float ax, bx;

	ax = exp(-1.0f / 2.0f*pow(((ii + 0.5f - x) / Sx), 2.0f));
	bx = exp(-1.0f / 2.0f*pow((ii - 0.5f - x) / Sx, 2.0f));
	*dudt = -N / sqrt(2.0f*pi) / Sx / Sx*(ax*(ii - x + 0.5f) - bx*(ii - x - 0.5f))*PSFy;

	if (d2udt2)
		*d2udt2 = -2.0f / Sx*dudt[0] - N / sqrt(2.0f*pi) / pow(Sx, 5)*(ax*pow((ii - x + 0.5f), 3) - bx*pow((ii - x - 0.5f), 3))*PSFy;
}

//

//*******************************************************************************************
__device__ void kernel_DerivativeIntGauss1D(const int ii, const float x, const float sigma, const float N,
	const float PSFy, float *dudt, float *d2udt2) {
	/*!
	* \brief compute the derivative of the 1D gaussian
	* \param ii ???
	* \param x ???
	* \param sigma ???
	* \param N ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/
	float a, b;
	a = exp(-1.0f / 2.0f*pow(((ii + 0.5f - x) / sigma), 2.0f));
	b = exp(-1.0f / 2.0f*pow((ii - 0.5f - x) / sigma, 2.0f));

	*dudt = -N / sqrt(2.0f*pi) / sigma*(a - b)*PSFy;

	if (d2udt2)
		*d2udt2 = -N / sqrt(2.0f*pi) / pow(sigma, 3)*((ii + 0.5f - x)*a - (ii - 0.5f - x)*b)*PSFy;
}


//*******************************************************************************************
__device__ void kernel_DerivativeIntGauss2DSigma(const int ii, const int jj, const float x, const float y,
	const float S, const float N, const float PSFx, const float PSFy, float *dudt, float *d2udt2) {
	/*!
	* \brief compute the derivative of the 2D gaussian
	* \param ii ???
	* \param jj ???
	* \param x ???
	* \param y ???
	* \param S ???
	* \param N ???
	* \param PSFx ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float dSx, dSy, ddSx, ddSy;

	kernel_DerivativeIntGauss1DSigma(ii, x, S, N, PSFy, &dSx, &ddSx);
	kernel_DerivativeIntGauss1DSigma(jj, y, S, N, PSFx, &dSy, &ddSy);

	*dudt = dSx + dSy;
	if (d2udt2) *d2udt2 = ddSx + ddSy; //+ 2.0f*dSx * dSy;

}

//*******************************************************************************************
__device__ void kernel_DsIntGaus1D(const int ii, const float x, const float sigma, float *ds, float *dss) {
	/*!
	* \brief compute the derivative of the 1D gaussian
	* \param ii ???
	* \param x ???
	* \param sigma ???
	* \param N ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float ax, bx;
	//TODO: check derivatives
	ax = exp(-1.0f / 2.0f*pow(((ii + 0.5f - x) / sigma), 2.0f));
	bx = exp(-1.0f / 2.0f*pow((ii - 0.5f - x) / sigma, 2.0f));
	
	*ds = -1.0f / sqrt(2.0f*pi) / sigma / sigma*(ax*(ii - x + 0.5f) - bx*(ii - x - 0.5f));

	*dss = -2.0f / sigma*ds[0] - 1.0f / sqrt(2.0f*pi) / pow(sigma, 5)*(ax*pow((ii - x + 0.5f), 3) - bx*pow((ii - x - 0.5f), 3));
}


//*******************************************************************************************
__device__ float kernel_alpha(const float z, const float Ax, const float Bx, const float d){
	/*! 
	 * \brief compute coefficient for alpha
	 * \param z ???
	 * \param Ax ???
	 * \param Bx ???
	 * \param d ???
	 * \return float alpha value
	 */
	
	return 1.0f+pow(z/d, 2)+Ax*pow(z/d, 3)+Bx*pow(z/d, 4);
}

//*******************************************************************************************
__device__ float kernel_dalphadz(const float z, const float Ax, const float Bx, const float d){
	/*! 
	 * \brief compute first derivative for alpha in relation to z
	 * \param z ???
	 * \param Ax ???
	 * \param Bx ???
	 * \param d ???
	 * \return float alpha value
	 */
    return (2.0f*z/(d*d) + 3.0f*Ax*pow(z, 2)/(d*d*d)+4.0f*Bx*pow(z, 3)/pow(d, 4));
}

//*******************************************************************************************
__device__ float kernel_d2alphadz2(const float z, const float Ax, const float Bx, const float d){
	/*! 
	 * \brief compute second derivative for alpha in relation to z
	 * \param z ???
	 * \param Ax ???
	 * \param Bx ???
	 * \param d ???
	 * \return float alpha value
	 */
    return (2.0f/(d*d) + 6.0f*Ax*z/(d*d*d)+12.0f*Bx*pow(z, 2)/pow(d, 4));
}

//*******************************************************************************************
__device__ void kernel_DerivativeIntGaussPSF1D(const int ii, const float x, const float sigma, const float N,
	const float PSFx, const float PSFy, float *dudt, float *d2udt2) {
	/*!
	* \brief compute the derivative of the 1D gaussian
	* \param ii ???
	* \param x ???
	* \param sigma ???
	* \param N ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/
	float dx= 0, dxx = 0;
	
	kernel_DxIntGaus1D(ii, x, sigma, &dx, &dxx);

	*dudt = N*dx*PSFx*PSFy;

	if (d2udt2)
		*d2udt2 = N*dxx*PSFx*PSFy;
}


//*******************************************************************************************
__device__ void kernel_DerivativeIntGaussPSF1D(const int ii, const float x, const float sigma, const float N,
        const float PSFy, float *dudt, float *d2udt2) {
	/*! 
	 * \brief compute the derivative of the 1D gaussian
	 * \param ii ???
	 * \param x ???
	 * \param sigma ???
	 * \param N ???
	 * \param PSFy ???
	 * \param dudt ???
	 * \param d2udt2 ???
	 */    
	float dx = 0, dxx = 0;

	kernel_DxIntGaus1D(ii, x, sigma, &dx, &dxx);

    
    *dudt = N*dx*PSFy;
    
    if (d2udt2)
        *d2udt2 =N*dxx*PSFy;
}

//*******************************************************************************************
__device__ void kernel_DerivativeIntGaussPSF1DSigma(const int ii, const float x,
        const float Sx, const float N, const float PSFy, float *dudt, float *d2udt2) {
	/*! 
	 * \brief compute the derivative of the 1D gaussian
	 * \param ii ???
	 * \param x ???
	 * \param Sx ???
	 * \param N ???
	 * \param PSFy ???
	 * \param dudt ???
	 * \param d2udt2 ???
	 */    
    
    float ds = 0, dss = 0;

	kernel_DsIntGaus1D(ii, x, Sx, &ds, &dss);
    

    *dudt = N*ds*PSFy;
    
    if (d2udt2)
        *d2udt2 =N*dss*PSFy;
}

//*******************************************************************************************
__device__ void kernel_DerivativeIntGaussPSF2DSigma(const int ii, const int jj, const float x, const float y,
	const float S, const float N, const float PSFx, const float PSFy, float *dudt, float *d2udt2) {
	/*!
	* \brief compute the derivative of the 2D gaussian
	* \param ii ???
	* \param jj ???
	* \param x ???
	* \param y ???
	* \param S ???
	* \param N ???
	* \param PSFx ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float dSx, dSy, ddSx, ddSy;

	kernel_DerivativeIntGaussPSF1DSigma(ii, x, S, N, PSFy, &dSx, &ddSx);
	kernel_DerivativeIntGaussPSF1DSigma(jj, y, S, N, PSFx, &dSy, &ddSy);


	*dudt = dSx + dSy;
	if (d2udt2) *d2udt2 = ddSx + ddSy;

}

__device__ void kernel_DerivativeIntGaussPSF3DSigma(const int ii, const int jj, const int hh, const float x, const float y, const float z,
	const float Sx, const float Sz, const float N, const float PSFx, const float PSFy, const float PSFz, float *dudtxy, float *dudtz, float *d2udt2xy, float *d2udt2z) {
	
	/*!
	* \brief compute the derivative of the 2D gaussian
	* \param ii ???
	* \param jj ???
	* \param x ???
	* \param y ???
	* \param S ???
	* \param N ???
	* \param PSFx ???
	* \param PSFy ???
	* \param dudt ???
	* \param d2udt2 ???
	*/

	float dSxy, ddSxy;
	kernel_DerivativeIntGaussPSF2DSigma(ii, jj, x, y, Sx, N, PSFx, PSFy, &dSxy, &ddSxy);

	*dudtxy = dSxy*PSFz;
	if (d2udt2xy) *d2udt2xy = ddSxy*PSFz;

	float ds = 0, dss = 0;

	kernel_DsIntGaus1D(hh, z, Sz, &ds, &dss);


	*dudtz = N*ds*PSFy*PSFx;

	if (d2udt2z)
		*d2udt2z = N*dss*PSFy*PSFx;
}

//*******************************************************************************************
__device__ void kernel_DerivativeIntGaussPSF2Dz(const int ii, const int jj, const float *theta,
        const float PSFSigma_x, const float PSFSigma_y, const float Ax, const float Ay, 
		const float Bx, const float By, const float gamma, const float d, float *pPSFx, float *pPSFy, float *dudt, float *d2udt2) {
	/*! 
	 * \brief compute the derivative of the 2D gaussian
	 * \param ii ???
	 * \param jj ???
	 * \param theta ???
	 * \param PSFSigma_x ???
	 * \param PSFSigma_y ???
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param pPSFx ???
	 * \param pPSFy ???
	 * \param dudt ???
	 * \param d2udt2 ???
	 */    
    
    float Sx, Sy, dSx, dSy, ddSx, ddSy, dSdzx, dSdzy,ddSddzx,ddSddzy;
    float z, PSFx, PSFy,alphax,alphay,ddx,ddy;
    float dSdalpha_x,dSdalpha_y,d2Sdalpha2_x,d2Sdalpha2_y;
    z=theta[4];
    
    alphax  = kernel_alpha(z-gamma, Ax, Bx, d);
    alphay  = kernel_alpha(z+gamma, Ay, By, d);
 
    Sx=PSFSigma_x*sqrt(alphax);
    Sy=PSFSigma_y*sqrt(alphay);
    
    PSFx=kernel_IntGauss1D(ii, theta[0], Sx);
    PSFy=kernel_IntGauss1D(jj, theta[1], Sy);
    *pPSFx=PSFx;
    *pPSFy=PSFy;
    
    kernel_DerivativeIntGaussPSF1D(ii, theta[0], Sx, theta[2], PSFy, &dudt[0], &ddx);
	kernel_DerivativeIntGaussPSF1D(jj, theta[1], Sy, theta[2], PSFx, &dudt[1], &ddy);
	kernel_DerivativeIntGaussPSF1DSigma(ii, theta[0], Sx, theta[2], PSFy, &dSx, &ddSx);
	kernel_DerivativeIntGaussPSF1DSigma(jj, theta[1], Sy, theta[2], PSFx, &dSy, &ddSy);

    dSdalpha_x=PSFSigma_x/2.0f/sqrt(alphax);
    dSdalpha_y=PSFSigma_y/2.0f/sqrt(alphay);
    
    dSdzx  = dSdalpha_x*kernel_dalphadz(z-gamma, Ax, Bx, d); 
    dSdzy  = dSdalpha_y*kernel_dalphadz(z+gamma, Ay, By, d);
    dudt[4] = dSx*dSdzx+dSy*dSdzy;
    
    if (d2udt2){
    d2udt2[0] =ddx;
    d2udt2[1] =ddy;
    
    d2Sdalpha2_x=-PSFSigma_x/4.0f/pow(alphax,1.5f);
    d2Sdalpha2_y=-PSFSigma_y/4.0f/pow(alphay,1.5f);
    
    ddSddzx  = d2Sdalpha2_x*pow(kernel_dalphadz(z-gamma, Ax, Bx, d),2)+dSdalpha_x*kernel_d2alphadz2(z-gamma, Ax, Bx, d); 
    ddSddzy  = d2Sdalpha2_y*pow(kernel_dalphadz(z+gamma, Ay, By, d),2)+dSdalpha_y*kernel_d2alphadz2(z+gamma, Ay, By, d); 
    
    d2udt2[4] =ddSx*(dSdzx*dSdzx)+dSx*ddSddzx+
            ddSy*(dSdzy*dSdzy)+dSy*ddSddzy;
    }
}


__device__ void avgDerivative(int sz, float * data, float *aX, float *aY){

	float tmpx=0.0;
    float tmpy=0.0;
	for (int ii=0;ii<sz-1;ii++) for(int jj=0;jj<sz-1;jj++) {
        tmpx+=data[sz*jj+ii+1]-data[sz*jj+ii];
        tmpy+=data[sz*(jj+1)+ii]-data[sz*jj+ii];
    }
    *aX=(float)tmpx/(sz-1.0)/(sz-1.0);
    *aY=(float)tmpy/(sz-1.0)/(sz-1.0);

}

//*******************************************************************************************
__device__ void kernel_CenterofMass2D(const int sz, const float *data, float *x, float *y,float aX, float aY) {
	/*!
	 * \brief compute the 2D center of mass of a subregion
	 * \param sz nxn size of the subregion
	 * \param data subregion to search
	 * \param x x coordinate to return
	 * \param y y coordinate to return
	 */
    float tmpx=0.0f;
    float tmpy=0.0f;
    float tmpsum=0.0f;
    
    //for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
    //    tmpx+=data[sz*jj+ii]*ii;
    //    tmpy+=data[sz*jj+ii]*jj;
    //    tmpsum+=data[sz*jj+ii];
    //}
	aX=0;
	aY=0;
    for (int ii=0;ii<sz;ii++) for(int jj=0;jj<sz;jj++) {
        tmpx+=(data[sz*jj+ii]-aX*ii-aY*jj)*ii;
        tmpy+=(data[sz*jj+ii]-aY*jj-aX*ii)*jj;
        tmpsum+=data[sz*jj+ii]-aX*ii-aY*jj;
        }
    *x=tmpx/tmpsum;
    *y=tmpy/tmpsum;
}


//*******************************************************************************************
__device__ void kernel_GaussFMaxMin2D(const int sz, const float sigma,const float * data, float *MaxN, float *MinBG) {
    /*!
	 * \brief returns filtered min and pixels of a given subregion
	 * \param sz nxn size of the subregion
	 * \param sigma used in filter calculation
	 * \param data the subregion to search
	 * \param MaxN maximum pixel value
	 * \param MinBG minimum background value
	 */
    int ii, jj, kk, ll;
    float filteredpixel=0, sum=0;
    *MaxN=0.0f;
    *MinBG=10e10f; //big
    
    float norm=1.0f/2.0f/sigma/sigma;
    //loop over all pixels
    for (kk=0;kk<sz;kk++) for (ll=0;ll<sz;ll++){
        filteredpixel=0.0f;
        sum=0.0f;
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++){
            filteredpixel+=exp(-pow((float)(ii-kk), 2)*norm)*exp(-pow((float)(ll-jj), 2)*norm)*data[ii*sz+jj];
            sum+=exp(-pow((float)(ii-kk), 2)*norm)*exp(-pow((float)(ll-jj), 2)*norm);
        }
        filteredpixel/=sum;
        
        *MaxN=max(*MaxN, filteredpixel);
        *MinBG=min(*MinBG, filteredpixel);
    }

}

//*******************************************************************************************
__device__ void kernel_CenterofMass3D(const int szXY, const int szZ, const float *data, float *x, float *y, float *z, float aX, float aY, float aZ) {
	/*!
	* \brief compute the 2D center of mass of a subregion
	* \param sz nxn size of the subregion
	* \param data subregion to search
	* \param x x coordinate to return
	* \param y y coordinate to return
	*/
	float tmpx = 0.0f;
	float tmpy = 0.0f;
	float tmpz = 0.0f;
	float tmpsum = 0.0f;

	//for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
	//    tmpx+=data[sz*jj+ii]*ii;
	//    tmpy+=data[sz*jj+ii]*jj;
	//    tmpsum+=data[sz*jj+ii];
	//}
	aX = 0;
	aY = 0;
	aZ = 0;
	for (int hh = 0; hh<szZ; hh++) for (int ii = 0; ii<szXY; ii++) for (int jj = 0; jj<szXY; jj++) {
		tmpx += (data[hh*szXY*szXY +szXY*ii+jj] - aX*ii - aY*jj-aZ*hh)*ii;
		tmpy += (data[hh*szXY*szXY +szXY*ii+jj] - aY*jj - aX*ii - aZ*hh)*jj;
		tmpz += (data[hh*szXY*szXY + szXY*ii+jj] - aY*jj - aX*ii - aZ*hh)*hh;
		tmpsum += data[hh*szXY*szXY + szXY*ii+jj] - aX*ii - aY*jj - aZ*hh;
	}
	*x = tmpx / tmpsum;
	*y = tmpy / tmpsum;
	*z = tmpz / tmpsum;
}


//*******************************************************************************************
__device__ void kernel_GaussFMaxMin3D(const int szXY, const int szZ, float sigmax, float sigmaz, const float * data, float *MaxN, float *MinBG) {
	/*!
	* \brief returns filtered min and pixels of a given subregion
	* \param sz nxn size of the subregion
	* \param sigma used in filter calculation
	* \param data the subregion to search
	* \param MaxN maximum pixel value
	* \param MinBG minimum background value
	*/
	int hh,ii, jj, kk, ll,oo;
	float filteredpixel = 0, sum = 0, gauss;
	*MaxN = 0.0f;
	*MinBG = 10e10f; //big

	float norm = 1.0f / pow(sqrt(2.0f),3) / sigmax/sigmax/sigmaz;
	//loop over all pixels
	for (kk = 0; kk<szZ; kk++) for (ll = 0; ll<szXY; ll++)for (oo = 0; oo<szXY; oo++){
		filteredpixel = 0.0f;
		sum = 0.0f;
		for (hh = 0; hh < szZ; hh++) for (ii = 0; ii < szXY; ii++) for (jj = 0; jj < szXY; jj++){
			gauss = exp(-pow((float)(hh - kk), 2)*norm)*exp(-pow((float)(ii - ll), 2)*norm)*exp(-pow((float)(oo - jj), 2)*norm);
			filteredpixel += gauss*data[hh*szXY*szXY + szXY*ii+jj];
			sum += gauss;
		}
		filteredpixel /= sum;

		*MaxN = max(*MaxN, filteredpixel);
		*MinBG = min(*MinBG, filteredpixel);
	}

}


















