/*
 *  molgroups.cc
 *  Gauss
 *
 *  Created by Frank Heinrich on 27/10/08.
 *  updated 9-May-2012
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "molgroups.h"


//------------------------------------------------------------------------------------------------------
//Parent Object Implementation

nSLDObj::nSLDObj()
{
    bWrapping=true;
    absorb=0;
};

nSLDObj::~nSLDObj(){};

double nSLDObj::fnGetAbsorb(double z){return absorb;};
void nSLDObj::fnWriteData2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    double dLowerLimit, dUpperLimit, d, dmirror, dAreaInc, dnSLDInc;
	int i;
    
    fprintf(fp, "z%s a%s nsl%s \n",cName, cName, cName);
	
	dLowerLimit=fnGetLowerLimit();
	dUpperLimit=fnGetUpperLimit();
	d=floor(dLowerLimit/stepsize+0.5)*stepsize;
	
	for (i=0; i<dimension; i++)
	{
        d=double(i)*stepsize;
        dmirror=d-float(2*i)*stepsize;
        if ((bWrapping==true) && (dmirror>=dLowerLimit))
        {
            dAreaInc=fnGetArea(d)+fnGetArea(dmirror);
            dnSLDInc=(fnGetnSLD(d)*fnGetArea(d)+fnGetnSLD(dmirror)*fnGetArea(dmirror))/(fnGetArea(d)+fnGetArea(dmirror));
        }
        else
        {
            dAreaInc=fnGetArea(d);
            dnSLDInc=fnGetnSLD(d);
            //printf("Bin %i Area %f nSLD %e nSL %e \n", i, dAreaInc, fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize);
        }
        fprintf(fp, "%lf %lf %e \n", d, dAreaInc, dnSLDInc*dAreaInc*stepsize);
	};
    fprintf(fp, "\n");
}


//Philosophy for this first method: You simply add more and more volume and nSLD to the
//volume and nSLD array. After all objects have filled up those arrays the maximal area is
//determined which is the area per molecule and unfilled volume is filled with bulk solvent.
//Hopefully the fit algorithm finds a physically meaningful solution. There has to be a global
//hydration paramter for the bilayer.
//Returns maximum area


//does a Catmull-Rom Interpolation on an equal distance grid
// 0<t<=1 is the relative position on the interval between p0 and p1
// p-1 and p2 are needed for derivative calculation
double nSLDObj::CatmullInterpolate(double t, double pm1, double p0, double p1, double p2){
    
    double m0, m1, t_2, t_3, h00, h10, h01, h11;
        
    m0=(p1-pm1)/2;
    m1=(p2-p0) /2;
            
    t_2=t*t;
    t_3=t_2*t;
    h00=   2*t_3-3*t_2+1;
    h10=     t_3-2*t_2+t;
    h01=(-2)*t_3+3*t_2;
    h11=     t_3-t_2;
    
    return h00*p0+h10*m0+h01*p1+h11*m1;
    
};

double nSLDObj::fnTriCubicCatmullInterpolate(double p[4][4][4],double t[3]){
    double dFirstStage[4][4];
    double dSecondStage[4];
    int i,j;
    
    for (i=0; i<4; i++){
        for (j=0; j<4; j++){
            dFirstStage[i][j]=CatmullInterpolate(t[0],p[0][i][j],p[1][i][j],p[2][i][j],p[3][i][j]);
        }
    }
    
    for (i=0; i<4; i++){
        dSecondStage[i]=CatmullInterpolate(t[1],dFirstStage[0][i],dFirstStage[1][i],dFirstStage[2][i],dFirstStage[3][i]);
    }
    
    return CatmullInterpolate(t[2],dSecondStage[0],dSecondStage[1],dSecondStage[2],dSecondStage[3]);        

};



double nSLDObj::fnWriteProfile(double aArea[], double anSL[], int dimension, double stepsize, double dMaxArea)
{
    double dLowerLimit, dUpperLimit, d, dAreaInc, dprefactor;
	int i;
    
	dLowerLimit=fnGetLowerLimit();
	dUpperLimit=fnGetUpperLimit();
    if (dUpperLimit==0)
    {
        dUpperLimit=double(dimension)*stepsize;
    }
	d=floor(dLowerLimit/stepsize+0.5)*stepsize;
	
	while (d<=dUpperLimit)
	{
	    i=int(d/stepsize);
		dprefactor=1;
		//printf("Here we are %i, dimension %i \n", i, dimension);
		if ((i<0) && (bWrapping==true)) {i=-1*i;};
		if ((i==0) && (bWrapping==true)) {dprefactor=2;}											//avoid too low filling when mirroring
		if ((i>=0) && (i<dimension))
		{
		    dAreaInc=fnGetArea(d);
		    aArea[i]=aArea[i]+dAreaInc*dprefactor;
			if (aArea[i]>dMaxArea) {dMaxArea=aArea[i];};
			anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
			//printf("Bin %i Area %f total %f nSL %f total %f \n", i, dAreaInc, aArea[i], fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
		}
		d=d+stepsize;  
	};
	
	return dMaxArea;
    
};
double nSLDObj::fnWriteProfile(double aArea[], double anSL[], double aAbsorb[], int dimension, double stepsize, double dMaxArea)
{
    double dLowerLimit, dUpperLimit, d, dAreaInc, dprefactor;
	int i;
	
	dLowerLimit=fnGetLowerLimit();
	dUpperLimit=fnGetUpperLimit();
    if (dUpperLimit==0)
    {
        dUpperLimit=double(dimension)*stepsize;
    }
	d=floor(dLowerLimit/stepsize+0.5)*stepsize;
	
	while (d<=dUpperLimit)
	{
	    i=int(d/stepsize);
		dprefactor=1;
		//printf("Here we are %i, dimension %i \n", i, dimension);
		if ((i<0) && (bWrapping==true)) {i=-1*i;};
		if ((i==0) && (bWrapping==true)) {dprefactor=2;}											//avoid too low filling when mirroring
		if ((i>=0) && (i<dimension))
		{
            //printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
		    dAreaInc=fnGetArea(d);
		    aArea[i]=aArea[i]+dAreaInc*dprefactor;
			if (aArea[i]>dMaxArea) {dMaxArea=aArea[i];};
			anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
            aAbsorb[i]=aAbsorb[i]+fnGetAbsorb(d)*dAreaInc*stepsize*dprefactor;
			//printf("Bin %i Area %f total %f nSL %f total %f \n", i, dAreaInc, aArea[i], fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
		}
		d=d+stepsize;  
	};
	
	return dMaxArea;
    
};
void nSLDObj::fnOverlayProfile(double aArea[], double anSL[], int dimension, double stepsize, double dMaxArea)
{
    double dLowerLimit, dUpperLimit, d, dAreaInc, dprefactor, temparea;
	int i;
	
	dLowerLimit=fnGetLowerLimit();
	dUpperLimit=fnGetUpperLimit();
    if (dUpperLimit==0)
    {
        dUpperLimit=double(dimension)*stepsize;
    }
	d=floor(dLowerLimit/stepsize+0.5)*stepsize;
	
	while (d<=dUpperLimit)
	{
	    i=int(d/stepsize);
		dprefactor=1;
		//printf("Here we are %i, dimension %i, maxarea %f \n", i, dimension, dMaxArea);
		if ((i<0) && (bWrapping==true)) {i=-1*i;};
		if ((i==0) && (bWrapping==true)) {dprefactor=2;}											//avoid too low filling when mirroring
		if ((i>=0) && (i<dimension))
		{
		    dAreaInc=fnGetArea(d);
			temparea=dAreaInc*dprefactor+aArea[i];
			if (temparea>dMaxArea) {
                //printf("Bin %i Areainc %f area now %f nSLD %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
                anSL[i]=anSL[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
                anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
                aArea[i]=dMaxArea;
			}
            else {
                aArea[i]=aArea[i]+dAreaInc*dprefactor;
                anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
            }
		}
		d=d+stepsize;  
	};
};

void nSLDObj::fnOverlayProfile(double aArea[], double anSL[], double aAbsorb[], int dimension, double stepsize, double dMaxArea)
{
    double dLowerLimit, dUpperLimit, d, dAreaInc, dprefactor, temparea;
	int i;
	
	dLowerLimit=fnGetLowerLimit();
	dUpperLimit=fnGetUpperLimit();
    if (dUpperLimit==0)
    {
        dUpperLimit=double(dimension)*stepsize;
    }
	d=floor(dLowerLimit/stepsize+0.5)*stepsize;
	
	while (d<=dUpperLimit)
	{
	    i=int(d/stepsize);
		dprefactor=1;
		//printf("Here we are %i, dimension %i, maxarea %f \n", i, dimension, dMaxArea);
		if ((i<0) && (bWrapping==true)) {i=-1*i;};
		if ((i==0) && (bWrapping==true)) {dprefactor=2;}											//avoid too low filling when mirroring
		if ((i>=0) && (i<dimension))
		{
		    dAreaInc=fnGetArea(d);
			temparea=dAreaInc*dprefactor+aArea[i];
			if (temparea>dMaxArea) {
                //printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
                anSL[i]=anSL[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
                anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
                aAbsorb[i]=aAbsorb[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
                aAbsorb[i]=aAbsorb[i]+fnGetAbsorb(d)*dAreaInc*stepsize*dprefactor;
                aArea[i]=dMaxArea;
			}
            else {
                //printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
                aArea[i]=aArea[i]+dAreaInc*dprefactor;
                anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
                aAbsorb[i]=aAbsorb[i]+fnGetAbsorb(d)*dAreaInc*stepsize*dprefactor;
            }
		}
		d=d+stepsize;  
	};
};

//------------------------------------------------------------------------------------------------------
//Function Object Implementation
//------------------------------------------------------------------------------------------------------

BoxErr::BoxErr(double dz, double dsigma, double dlength, double dvolume, double dnSL, double dnumberfraction=1)
{
    z=dz; sigma=dsigma; l=dlength, vol=dvolume, nSL=dnSL, nf=dnumberfraction;
};

BoxErr::~BoxErr(){};

//Gaussian function definition, integral is volume, return value is area at position z
double BoxErr::fnGetArea(double dz) {
    
    return (vol/l)*0.5*(erf((dz-z+0.5*l)/sqrt(2)/sigma)-erf((dz-z-0.5*l)/sqrt(2)/sigma))*nf;
};

//constant nSLD
double BoxErr::fnGetnSLD(double dz) {return nSL/vol;};

//Gaussians are cut off below and above 3 sigma
double BoxErr::fnGetLowerLimit() {return z-0.5*l-3*sigma;};
double BoxErr::fnGetUpperLimit() {return z+0.5*l+3*sigma;};

void   BoxErr::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "BoxErr %s z %lf sigma %lf l %lf vol %lf nSL %e nf %lf \n",cName, z, sigma, l, vol, nSL, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
}


//------------------------------------------------------------------------------------------------------

Box2Err::Box2Err(double dz, double dsigma1, double dsigma2, double dlength, double dvolume, double dnSL, double dnumberfraction=1)
{
    z=dz; sigma1=dsigma1; sigma2=dsigma2; l=dlength, vol=dvolume, nSL=dnSL, nf=dnumberfraction;
};

Box2Err::~Box2Err(){};

//Gaussian function definition, integral is volume, return value is area at position z
double Box2Err::fnGetArea(double dz) {
    
    return (vol/l)*0.5*(erf((dz-z+0.5*l)/sqrt(2)/sigma1)-erf((dz-z-0.5*l)/sqrt(2)/sigma2))*nf;
};

//constant nSLD
double Box2Err::fnGetnSLD(double dz) {return nSL/vol;};

//Gaussians are cut off below and above 3 sigma
double Box2Err::fnGetLowerLimit() {return z-0.5*l-3*sigma1;};
double Box2Err::fnGetUpperLimit() {return z+0.5*l+3*sigma2;};

void Box2Err::fnSetSigma(double sigma)
{
	sigma1=sigma;
	sigma2=sigma;
}
void Box2Err::fnSetSigma(double dsigma1, double dsigma2)
{
	sigma1=dsigma1;
	sigma2=dsigma2;
}

void Box2Err::fnSetZ(double dz)
{
	z=dz;
};

void   Box2Err::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Box2Err %s z %lf sigma1 %lf sigma2 %lf l %lf vol %lf nSL %e nf %lf \n",cName, z, sigma1, sigma2, l, vol, nSL, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
}


//------------------------------------------------------------------------------------------------------

Gaussian::Gaussian(double dz, double dsigma, double dvolume, double dnSL, double dnumberfraction=1)
{
    z=dz; sigma=dsigma; vol=dvolume, nSL=dnSL, nf=dnumberfraction;
};

Gaussian::~Gaussian(){};

//Gaussian function definition, integral is volume, return value is area at position z
double Gaussian::fnGetArea(double dz) {return (vol/sqrt(2*3.141592654)/sigma)*exp(-0.5*(z-dz)*(z-dz)/sigma/sigma)*nf;};

//constant nSLD
double Gaussian::fnGetnSLD(double dz) {return nSL/vol;};

//Gaussians are cut off below and above 3 sigma
double Gaussian::fnGetLowerLimit() {return z-3*sigma;};
double Gaussian::fnGetUpperLimit() {return z+3*sigma;};

void   Gaussian::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Gaussian %s z %lf sigma %lf vol %lf nSL %e nf %lf \n",cName, z, sigma, vol, nSL, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
}

//------------------------------------------------------------------------------------------------------

Parabolic::Parabolic(double dC, double dH, double dn, double dnSLD, double dnumberfraction=1)
{
    C=dC; H=dH, n=dn, nSLD=dnSLD, nf=dnumberfraction;
    bWrapping=false;
};

Parabolic::~Parabolic(){};

//Gaussian function definition, integral is volume, return value is area at position z
double Parabolic::fnGetArea(double dz) {
    if (dz<H) {return C*(1-pow(dz/H,n))*nf;}
    else {return 0;}
};

//constant nSLD
double Parabolic::fnGetnSLD(double dz) {return nSLD;};

//Gaussians are cut off below and above 3 sigma
double Parabolic::fnGetLowerLimit() {return 0;};
double Parabolic::fnGetUpperLimit() {return 0;};

void   Parabolic::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Parabolic %s C %lf H %lf n %lf  nSLD %e nf %lf \n",cName, C, H, n, nSLD, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
}
//------------------------------------------------------------------------------------------------------

StretchGaussian::StretchGaussian(double dz, double dsigma, double dlength, double dvolume, double dnSL, double dnumberfraction=1)
{
    z=dz; sigma=dsigma; l=dlength, vol=dvolume, nSL=dnSL, nf=dnumberfraction;
};

StretchGaussian::~StretchGaussian(){};

//Gaussian function definition, integral is volume, return value is area at position z
double StretchGaussian::fnGetArea(double dz) {
    
    double returnvalue;
    double temp, dvgauss;
    
    temp=sqrt(2*3.141592654)*sigma;
    dvgauss=vol/(1+l/temp);
    
    if (dz<(z-0.5*l))
    {
        returnvalue=dvgauss/temp*exp(-0.5*(z-dz-0.5*l)*(z-dz-0.5*l)/sigma*sigma)*nf;
    }
    else if ((dz>=(z-0.5*l)) && (dz<=(z+0.5*l)))
    {
        returnvalue=dvgauss/temp*nf;
    }
    else
    {
        returnvalue=dvgauss/temp*exp(-0.5*(dz-z-0.5*l)*(dz-z-0.5*l)/sigma*sigma)*nf;
    }
    
    return returnvalue;
};

//constant nSLD
double StretchGaussian::fnGetnSLD(double dz) {return nSL/vol;};

//Gaussians are cut off below and above 3 sigma
double StretchGaussian::fnGetLowerLimit() {return z-0.5*l-3*sigma;};
double StretchGaussian::fnGetUpperLimit() {return z+0.5*l+3*sigma;};

void   StretchGaussian::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Gaussian %s z %lf sigma %lf l %lf vol %lf nSL %e nf %lf \n",cName, z, sigma, l, vol, nSL, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
}

//------------------------------------------------------------------------------------------------------
// Combined Object Implementation
//------------------------------------------------------------------------------------------------------

PC::PC()
{
    
    cg = new Box2Err(0,0,0,0,0,0,1);
    phosphate = new Box2Err(0,0,0,0,0,0,1);
    choline = new Box2Err(0,0,0,0,0,0,1);
    
    cg->l=4.21; phosphate->l=3.86; choline->l=6.34;                             //from fit to Feller data
    cg->sigma1=2.53; cg->sigma2=2.29;
    phosphate->sigma1=2.29; phosphate->sigma2=2.02;
    choline->sigma1=2.02; choline->sigma2=2.26;
    //from fit to Feller data
    l=9.575;                                                                    //group cg    phosphate choline
    //z     15.00 18.44     19.30
    //l      4.21  3.86      6.34
    
    cg->vol=147; phosphate->vol=54; choline->vol=120;                           //nominal values
    cg->nSL=3.7755e-4; phosphate->nSL=2.8350e-4; choline->nSL=-6.0930e-5;
    cg->nf=1; phosphate->nf=1; choline->nf=1;
    
    fnAdjustParameters();
    
};

PC::~PC(){
    delete cg;
    delete phosphate;
    delete choline;
};

void PC::fnAdjustParameters(){
    cg->z=z-0.5*l+0.5*cg->l; phosphate->z=z-0.5*l+cg->l+0.5*phosphate->l;
    choline->z=z+0.5*l-0.5*choline->l;                                           
};

//Return value is area at position z
double PC::fnGetArea(double dz) {
    return (cg->fnGetArea(dz)+phosphate->fnGetArea(dz)+choline->fnGetArea(dz))*nf;
};

//get nSLD from molecular subgroups
double PC::fnGetnSLD(double dz) {
    double cgarea, pharea, charea, sum;
    
    cgarea=cg->fnGetArea(dz);
    pharea=phosphate->fnGetArea(dz);
    charea=choline->fnGetArea(dz);
    sum=cgarea+pharea+charea;
    
    if (sum==0) {return 0;}
    else {
        return (cg->fnGetnSLD(dz)*cgarea+
                phosphate->fnGetnSLD(dz)*pharea+
                choline->fnGetnSLD(dz)*charea)/sum;
    }
};

//Use limits of molecular subgroups
double PC::fnGetLowerLimit() {return cg->fnGetLowerLimit();};
double PC::fnGetUpperLimit() {return choline->fnGetUpperLimit();};

void PC::fnSetSigma(double sigma)
{
    cg->sigma1=sigma;
    cg->sigma2=sigma;
    phosphate->sigma1=sigma;
    phosphate->sigma2=sigma;
    choline->sigma1=sigma;
    choline->sigma2=sigma;
};


void PC::fnSetZ(double dz){
    z=dz;
    fnAdjustParameters();
};

void PC::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    //cg->fnWriteGroup2File(fp, "cg", dimension, stepsize);
    //phosphate->fnWriteGroup2File(fp, "phosphate", dimension, stepsize);
    //choline->fnWriteGroup2File(fp, "choline", dimension, stepsize);
    
    //delete []str;
    
}
//------------------------------------------------------------------------------------------------------
PCm::PCm() 
{
    cg->sigma2=2.53; cg->sigma1=2.29;										//from fit to Feller data
    phosphate->sigma2=2.29; phosphate->sigma1=2.02;
    choline->sigma2=2.02; choline->sigma1=2.26;
    
    fnAdjustParameters();
    
};

PCm::~PCm() {};

void PCm::fnAdjustParameters(){
    cg->z=z+0.5*l-0.5*cg->l; phosphate->z=z+0.5*l-cg->l-0.5*phosphate->l;
    choline->z=z-0.5*l+0.5*choline->l;                                           
};

//Use limits of molecular subgroups
double PCm::fnGetLowerLimit() {return cg->fnGetLowerLimit();};
double PCm::fnGetUpperLimit() {return choline->fnGetUpperLimit();};

void   PCm::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    fprintf(fp, "PCm %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    //cg->fnWriteGroup2File(fp, "cg_m", dimension, stepsize);
    //phosphate->fnWriteGroup2File(fp, "phosphate_m", dimension, stepsize);
    //choline->fnWriteGroup2File(fp, "choline_m", dimension, stepsize);
    
    //delete []str;
    
}
//---------------------------------------------------------------------------------------------------------------
PS::PS()
{
    
    cg = new Box2Err(0,0,0,0,0,0,1);
    phosphate = new Box2Err(0,0,0,0,0,0,1);
    serine = new Box2Err(0,0,0,0,0,0,1);
    
    cg->l=3.6; phosphate->l=2.0; serine->l=3.0;                             //from fit to Feller data
    cg->sigma1=2.53; cg->sigma2=2.29;
    phosphate->sigma1=2.29; phosphate->sigma2=2.02;
    serine->sigma1=2.02; serine->sigma2=2.26;
    //from fit to Feller data
    l=8.6;                                                                    //group cg    phosphate choline
    //z     15.00 18.44     19.30
    //l      4.21  3.86      6.34
    
    cg->vol=147; phosphate->vol=54; serine->vol=80;                           //nominal values
    cg->nSL=3.7755e-4; phosphate->nSL=2.8350e-4; serine->nSL=1.8408E-04;
    cg->nf=1; phosphate->nf=1; serine->nf=1;
    
    fnAdjustParameters();
    
};

PS::~PS(){
    delete cg;
    delete phosphate;
    delete serine;
};

void PS::fnAdjustParameters(){
    cg->z=z-0.5*l+0.5*cg->l; phosphate->z=z-0.5*l+cg->l+0.5*phosphate->l;
    serine->z=z+0.5*l-0.5*serine->l;                                           
};

//Return value is area at position z
double PS::fnGetArea(double dz) {
    return (cg->fnGetArea(dz)+phosphate->fnGetArea(dz)+serine->fnGetArea(dz))*nf;
};

//get nSLD from molecular subgroups
double PS::fnGetnSLD(double dz) {
    double cgarea, pharea, searea, sum;
    
    cgarea=cg->fnGetArea(dz);
    pharea=phosphate->fnGetArea(dz);
    searea=serine->fnGetArea(dz);
    sum=cgarea+pharea+searea;
    
    if (sum==0) {return 0;}
    else {
        return (cg->fnGetnSLD(dz)*cgarea+
                phosphate->fnGetnSLD(dz)*pharea+
                serine->fnGetnSLD(dz)*searea)/sum;
    }
};

//Use limits of molecular subgroups
double PS::fnGetLowerLimit() {return cg->fnGetLowerLimit();};
double PS::fnGetUpperLimit() {return serine->fnGetUpperLimit();};

void PS::fnSetSigma(double sigma)
{
    cg->sigma1=sigma;
    cg->sigma2=sigma;
    phosphate->sigma1=sigma;
    phosphate->sigma2=sigma;
    serine->sigma1=sigma;
    serine->sigma2=sigma;
};


void PS::fnSetZ(double dz){
    z=dz;
    fnAdjustParameters();
};

void PS::fnSetnSL(double nSL_cg, double nSL_phosphate, double nSL_serine){
    //printf("nSL cg %e nSL phosphate %e nSL serine %e \n", nSL_cg, nSL_phosphate, nSL_serine);
    cg->nSL=nSL_cg;
    phosphate->nSL=nSL_phosphate;
    serine->nSL=nSL_serine;
}


void PS::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    cg->fnWriteGroup2File(fp, "cg", dimension, stepsize);
    phosphate->fnWriteGroup2File(fp, "phosphate", dimension, stepsize);
    serine->fnWriteGroup2File(fp, "serine", dimension, stepsize);
    
    //delete []str;
    
}

//-----------------------------------------------------------------------------------------------------------
//Amino Acids
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//general implementation
//-----------------------------------------------------------------------------------------------------------
//ExchangeRatio is not yet fully implemented

AminoAcid::AminoAcid()
{
	ExchangeRatio=0;
	Deuterated=0;
}

double AminoAcid::fnGetnSLD(double z) 
{
	double temp;
    
	if (Deuterated==1) {temp=1;}
	else {temp=0.;}
    //correct for deuterated residues
	return (nSL-temp*float(nH)*(-3.741e-5)+temp*float(nH)*(6.674e-5))/vol;
}
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//Specific implementations
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

AA_Lys::AA_Lys()
{
	nSL=1.5660e-4;
	vol=171.3;
	nH=13;
	nExch=4;
}
AA_Arg::AA_Arg()
{
	nSL=3.4260e-4;
	vol=202.1;
	nH=13;
	nExch=6;
}
AA_His::AA_His()
{
	nSL=4.7406e-4;
	vol=167.3;
	nH=7;
	nExch=3;
}
AA_Asn::AA_Asn()
{
	nSL=3.4356e-4;
	vol=135.2;
	nH=6;
	nExch=3;
}
AA_Asp::AA_Asp()
{
	nSL=3.8343e-4;
	vol=124.5;
	nH=4;
	nExch=1;
}
AA_Cys::AA_Cys()
{
	nSL=1.9191e-4;
	vol=105.6;
	nH=5;
	nExch=1;
}
AA_Thr::AA_Thr()
{
	nSL=2.1315e-4;
	vol=122.1;
	nH=7;
	nExch=2;
}
AA_Ser::AA_Ser()
{
	nSL=2.2149e-4;
	vol=99.1;
	nH=5;
	nExch=2;
}
AA_Gln::AA_Gln()
{
	nSL=3.3522e-4;
	vol=161.1;
	nH=8;
	nExch=3;
}
AA_Glu::AA_Glu()
{
	nSL=3.7509e-4;
	vol=155.1;
	nH=6;
	nExch=1;
}
AA_Pro::AA_Pro()
{
	nSL=2.2158e-4;
	vol=129.3;
	nH=7;
	nExch=0;
}
AA_Gly::AA_Gly()
{
	nSL=1.7178e-4;
	vol=66.4;
	nH=3;
	nExch=1;
}
AA_Ala::AA_Ala()
{
	nSL=1.6344e-4;
	vol=91.5;
	nH=5;
	nExch=1;
}
AA_Val::AA_Val()
{
	nSL=1.4676e-4;
	vol=141.7;
	nH=9;
	nExch=1;
}
AA_Ile::AA_Ile()
{
	nSL=1.3842e-4;
	vol=168.8;
	nH=11;
	nExch=1;
}
AA_Leu::AA_Leu()
{
	nSL=1.3842e-4;
	vol=167.9;
	nH=11;
	nExch=1;
}
AA_Met::AA_Met()
{
	nSL=1.7523e-4;
	vol=170.8;
	nH=9;
	nExch=1;
}
AA_Tyr::AA_Tyr()
{
	nSL=4.7073e-4;
	vol=203.6;
	nH=9;
	nExch=2;
}
AA_Phe::AA_Phe()
{
	nSL=4.1268e-4;
	vol=203.4;
	nH=9;
	nExch=1;
}
AA_Trp::AA_Trp()
{
	nSL=6.0123e-4;
	vol=237.6;
	nH=10;
	nExch=2;
}

//------------------------------------------------------------------------------------------------------
// Monolayer - single PC lipid
//------------------------------------------------------------------------------------------------------
Monolayer::Monolayer(){
    
	substrate = new Box2Err();
	lipid     = new	Box2Err();
	methyl    = new	Box2Err();
	
	substrate->l=20;
	substrate->z=10;
	substrate->nf=1;
	substrate->sigma1=2.0;
    rho_substrate=0;
    absorb_substrate=0;
	
    volacyllipid=925;                                                       //DOPC
    nslacyllipid=-2.67e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
	
	hc_substitution=0.;

    //fnAdjustParameters();
};

Monolayer::~Monolayer(){
    delete substrate;
    delete lipid;
    delete methyl;
};

void Monolayer::fnAdjustParameters(){
    
    double l_hc;
    double V_hc;
    double nf_hc_lipid, nSL_hc, absorb_hc;
    double c_s_hc, c_A_hc, c_V_hc;
    
    double l_m;
    double V_m;
    double nf_m_lipid, nSL_m, absorb_m;
    double c_s_m, c_A_m, c_V_m;  
    
    //set all sigma
    
    fnSetSigma(sigma);
    
    //outer hydrocarbons
    l_hc=l_lipid;
    nf_hc_lipid=1;
    V_hc=nf_hc_lipid*(volacyllipid-volmethyllipid);
    nSL_hc=nf_hc_lipid*(nslacyllipid-nslmethyllipid);
    absorb_hc=nf_hc_lipid*(absorbacyllipid-absorbmethyllipid);
    
    normarea=V_hc/l_hc;
    c_s_hc=vf_bilayer;
    c_A_hc=1;
    c_V_hc=1;
    
    //printf("%e %e %e \n", normarea, l_hc, nf_hc_lipid);

    
    lipid->l=l_hc;
    lipid->vol=V_hc;
    lipid->nSL=nSL_hc;
    lipid->absorb=absorb_hc;
    lipid->nf=c_s_hc*c_A_hc*c_V_hc;
    
    //outher methyl
    nf_m_lipid=1;
    V_m=nf_m_lipid*volmethyllipid;
    l_m=l_hc*V_m/V_hc;
    nSL_m=nf_m_lipid*nslmethyllipid;
    absorb_m=nf_m_lipid*absorbmethyllipid;
    
    c_s_m=c_s_hc;
    c_A_m=1;
    c_V_m=1;
    
    methyl->l=l_m;
    methyl->vol=V_m;
    methyl->nSL=nSL_m;
    methyl->absorb=absorb_m;
    methyl->nf=c_s_m*c_A_m*c_V_m;
    
    //PC headgroup
    headgroup->nf=c_s_hc*c_A_hc*nf_hc_lipid*(1-hc_substitution);
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    substrate->absorb=absorb_substrate*substrate->vol;
    
    
    // set all lengths
    methyl->z=substrate->l+0.5*(methyl->l);
    lipid->z=methyl->z+0.5*(methyl->l+lipid->l);
    headgroup->fnSetZ(lipid->z+0.5*(lipid->l+headgroup->l));    
    
};

//Return value is area at position z
double Monolayer::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+lipid->fnGetArea(dz)+headgroup->fnGetArea(dz)
            +methyl->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double Monolayer::fnGetnSLD(double dz) {
    double substratearea, lipidarea, headgrouparea;
    double methylarea, sum;
    
    substratearea=substrate->fnGetArea(dz);
    lipidarea=lipid->fnGetArea(dz);
    headgrouparea=headgroup->fnGetArea(dz);
    methylarea=methyl->fnGetArea(dz);
    
    sum=substratearea+lipidarea+headgrouparea+methylarea;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                headgroup->fnGetnSLD(dz)*headgrouparea+
                lipid->fnGetnSLD(dz)*lipidarea+
                methyl->fnGetnSLD(dz)*methylarea
                )/sum;
    }
};

//Use limits of molecular subgroups
double Monolayer::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double Monolayer::fnGetUpperLimit() {return headgroup->fnGetUpperLimit();};

double Monolayer::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};
double Monolayer::fnWriteProfile(double aArea[], double anSLD[], double aAbsorb[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,aAbsorb,dimension,stepsize,dMaxArea);
	return normarea;
};

void Monolayer::fnSetSigma(double dsigma)
{
    //set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    headgroup->fnSetSigma(dsigma);
    lipid->sigma2=dsigma;
    lipid->sigma1=global_rough;
    methyl->fnSetSigma(global_rough);
    substrate->sigma2=global_rough;

}
void Monolayer::fnSetnSL(double nSL_methyl, double nSL_lipid, double nSL_headgroup1, double nSL_headgroup2, double nSL_headgroup3)
{
    nslmethyllipid=nSL_methyl;
    nslacyllipid=nSL_lipid;
    headgroup->fnSetnSL(nSL_headgroup1,nSL_headgroup2,nSL_headgroup3);
}


void Monolayer::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    methyl->fnWriteGroup2File(fp, "methyl", dimension, stepsize);
    lipid->fnWriteGroup2File(fp, "lipid", dimension, stepsize);
    headgroup->fnWriteGroup2File(fp, "headgroup", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
        
}

//------------------------------------------------------------------------------------------------------
// Lipid bilayer single lipid
//------------------------------------------------------------------------------------------------------
ssBLM::ssBLM(){
    
	substrate = new Box2Err();
	headgroup1= new	PCm();                                                    //mirrored PC head group
	lipid1    = new	Box2Err();
	methyl1   = new	Box2Err();
	methyl2	  = new	Box2Err();
	lipid2	  = new Box2Err();
	headgroup2 = new PC();                                                   //PC head group
	
	substrate->l=20;
	substrate->z=10;
	substrate->nf=1;
	substrate->sigma1=2.0;
	
    volacyllipid=925;
    nslacyllipid=-2.67e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
	
	l_submembrane=3.0;
    
	hc_substitution_1=0.;
	hc_substitution_2=0.;	
	
    fnAdjustParameters();
};

ssBLM::~ssBLM(){
    delete substrate;
    delete headgroup1;
    delete lipid1;
    delete methyl1;
    delete methyl2;
    delete lipid2;
    delete headgroup2;
};

void ssBLM::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    lipid1->sigma1=sigma;
    lipid1->sigma2=sigma;
    methyl1->sigma1=sigma;
    methyl1->sigma2=sigma;
    methyl2->sigma1=sigma;
    methyl2->sigma2=sigma;
    lipid2->sigma1=sigma;
    lipid2->sigma2=sigma;
    headgroup2->fnSetSigma(sigma);
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid=1;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid);
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid);
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid=1;
    V_om=nf_om_lipid*volmethyllipid;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_lipid=1;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    V_im=nf_im_lipid*volmethyllipid;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer PC headgroup
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    
    //inner PC headgroup
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    headgroup1->fnSetZ(substrate->l+l_submembrane+0.5*headgroup1->l);
    lipid1->z=headgroup1->fnGetZ()+0.5*(headgroup1->l+lipid1->l);
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    
    
    
};

//Return value is area at position z
double ssBLM::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double ssBLM::fnGetnSLD(double dz) {
    double substratearea, lipid1area, headgroup1area;
    double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    
    substratearea=substrate->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    
    sum=substratearea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area
                )/sum;
    }
};

//Use limits of molecular subgroups
double ssBLM::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double ssBLM::fnGetUpperLimit() {return headgroup2->fnGetUpperLimit();};

double ssBLM::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void ssBLM::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}



//------------------------------------------------------------------------------------------------------
// Tethered Lipid bilayer - binary system
//------------------------------------------------------------------------------------------------------
tBLM::tBLM(){
    
	substrate  = new Box2Err();
	bME        = new Box2Err();
	tether	   = new Box2Err();
	tetherg    = new	Box2Err(); 
	headgroup1 = new	PCm();                                                    //mirrored PC head group
	lipid1     = new	Box2Err();
	methyl1    = new	Box2Err();
	methyl2	   = new	Box2Err();
	lipid2	   = new Box2Err();
	headgroup2 = new PC();                                                   //PC head group
	
	substrate->l=20;
	substrate->z=10;
	substrate->nf=1;
	substrate->sigma1=2.0;
	bME->vol=110;
	bME->nSL=3.243e-5;
	bME->l=5.2;
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=925;
    nslacyllipid=-2.67e-4;
	volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=982;
    nslacyltether=-2.85e-4;
	
	hc_substitution_1=0;
	hc_substitution_2=0;
    
    fnAdjustParameters();
};

tBLM::~tBLM(){
    delete substrate;
    delete bME;
    delete tether;
	delete tetherg;
	delete headgroup1;
	delete lipid1;
	delete methyl1;
	delete methyl2;
	delete lipid2;
	delete headgroup2;
};

void tBLM::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid,nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->sigma1=sigma; 
    tetherg->sigma2=sigma;
    lipid1->sigma1=sigma;
    lipid1->sigma2=sigma;
    methyl1->sigma1=sigma;
    methyl1->sigma2=sigma;
    methyl2->sigma1=sigma;
    methyl2->sigma2=sigma;
    lipid2->sigma1=sigma;
    lipid2->sigma2=sigma;
    headgroup2->fnSetSigma(sigma);
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid=1;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid);
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid);
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid=1;
    V_om=nf_om_lipid*volmethyllipid;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    nf_ihc_lipid=1-nf_ihc_tether;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer PC headgroup
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    
    //inner PC headgroup
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    
    
    
};

//Return value is area at position z
double tBLM::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
    double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM::fnGetUpperLimit() {return headgroup2->fnGetUpperLimit();};

double tBLM::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}

//------------------------------------------------------------------------------------------------------
// Lipid bilayer - binary system
//------------------------------------------------------------------------------------------------------
tBLM_binary::tBLM_binary(){
    
	headgroup1_2=  new Box2Err();                                                    //second headgroups
	headgroup2_2 = new Box2Err();                                                    
	
	headgroup1_2->vol=330;       //was 330
	headgroup2_2->vol=330;       //was 330
	headgroup1_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    iSeparateLeafletComposition=0;
    
    fnAdjustParameters();
};

tBLM_binary::~tBLM_binary(){
    delete headgroup1_2;
    delete headgroup2_2;
};

void tBLM_binary::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nf_ohc_lipid_2, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nf_om_lipid_2, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nf_ihc_lipid_2, nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_lipid_2, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    fnSetSigma(sigma);
        
    //outer hydrocarbons
    l_ohc=l_lipid2;
    if (iSeparateLeafletComposition==0) {
        nf_ohc_lipid  =1-nf_lipid_2;
        nf_ohc_lipid_2=nf_lipid_2;
    }
    else{
        nf_ohc_lipid  =1-nf_lipid2_2;
        nf_ohc_lipid_2=nf_lipid2_2;        
    }
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2);
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2);
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    if (iSeparateLeafletComposition==0){
        nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
        nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;        
    }
    else{
        nf_ihc_lipid=(1-nf_ihc_tether)*(1-nf_lipid1_2);
        nf_ihc_lipid_2=(1-nf_ihc_tether)*(nf_lipid1_2);                
    }
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer headgroups
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    headgroup2_2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2);
    
    //inner headgroups
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    headgroup1_2->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    headgroup1_2->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_2->l);  
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    headgroup2_2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_2->l);
    
    
    
};

//Return value is area at position z
double tBLM_binary::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz)+headgroup1_2->fnGetArea(dz)+headgroup2_2->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM_binary::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
    double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    double headgroup1_2_area, headgroup2_2_area;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    headgroup1_2_area=headgroup1_2->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    headgroup2_2_area=headgroup2_2->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area+headgroup1_2_area+headgroup2_2_area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                headgroup1_2->fnGetnSLD(dz)*headgroup1_2_area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area+
                headgroup2_2->fnGetnSLD(dz)*headgroup2_2_area
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM_binary::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM_binary::fnGetUpperLimit() 
{
	double a,b;
	a=headgroup2->fnGetUpperLimit();
	b=headgroup2_2->fnGetUpperLimit();
	
	if (a>b) {return a;} else {return b;}
    
};

void tBLM_binary::fnSetSigma(double sigma)
{
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    headgroup1_2->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->fnSetSigma(sigma); 
    lipid1->fnSetSigma(sigma);
    methyl1->fnSetSigma(sigma);
    methyl2->fnSetSigma(sigma);
    lipid2->fnSetSigma(sigma);
    headgroup2->fnSetSigma(sigma);
    headgroup2_2->fnSetSigma(sigma);
    
}

double tBLM_binary::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_binary::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    headgroup1_2->fnWriteGroup2File(fp, "headgroup1_2", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    headgroup2_2->fnWriteGroup2File(fp, "headgroup2_2", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}
//------------------------------------------------------------------------------------------------------
// Lipid bilayer - ternary system
//------------------------------------------------------------------------------------------------------
tBLM_ternary::tBLM_ternary(){
    
	headgroup1_2=  new Box2Err();                                                    //second headgroups
	headgroup2_2 = new Box2Err();  
	headgroup1_3=  new Box2Err();                                                    //third headgroups
	headgroup2_3 = new Box2Err();                                                  
	
	headgroup1_2->vol=330;       //was 330
	headgroup2_2->vol=330;       //was 330
	headgroup1_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
    
	headgroup1_3->vol=330;       //was 330
	headgroup2_3->vol=330;       //was 330
	headgroup1_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_3->l=9.5;
	headgroup2_3->l=9.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    volacyllipid_3=925;
    nslacyllipid_3=-2.67e-4;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    fnAdjustParameters();
};

tBLM_ternary::~tBLM_ternary(){
    delete headgroup1_2;
    delete headgroup2_2;
    delete headgroup1_3;
    delete headgroup2_3;
};

void tBLM_ternary::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nf_ohc_lipid_2, nf_ohc_lipid_3, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nf_om_lipid_2, nf_om_lipid_3, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nf_ihc_lipid_2, nf_ihc_lipid_3, nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_lipid_2, nf_im_lipid_3, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    
    fnSetSigma(sigma);
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid  =1-nf_lipid_2-nf_lipid_3;
    nf_ohc_lipid_2=nf_lipid_2;
    nf_ohc_lipid_3=nf_lipid_3;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ohc_lipid_3*(volacyllipid_3-volmethyllipid_3);
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ohc_lipid_3*(nslacyllipid_3-nslmethyllipid_3);
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    nf_om_lipid_3=nf_ohc_lipid_3;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2+nf_om_lipid_3*volmethyllipid_3;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2+nf_om_lipid_3*nslmethyllipid_3;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    nf_ihc_lipid_3=(1-nf_ihc_tether)*nf_ohc_lipid_3;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_lipid_3=nf_ihc_lipid_3;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_lipid_3*volmethyllipid_3+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_lipid_3*nslmethyllipid_3+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer headgroups
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    headgroup2_2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2);
    headgroup2_3->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_3*(1-hc_substitution_2);
    
    //inner headgroups
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    headgroup1_2->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1);
    headgroup1_3->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_3*(1-hc_substitution_1);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    headgroup1_2->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_2->l);
    headgroup1_3->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_3->l);  
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    headgroup2_2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_2->l);
    headgroup2_3->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_3->l);
    
    
};

//Return value is area at position z
double tBLM_ternary::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz)+headgroup1_2->fnGetArea(dz)+headgroup2_2->fnGetArea(dz)
            +headgroup1_3->fnGetArea(dz)+headgroup2_3->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM_ternary::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
	double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    double headgroup1_2_area, headgroup2_2_area;
	double headgroup1_3_area, headgroup2_3_area;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    headgroup1_2_area=headgroup1_2->fnGetArea(dz);
    headgroup1_3_area=headgroup1_3->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    headgroup2_2_area=headgroup2_2->fnGetArea(dz);
    headgroup2_3_area=headgroup2_3->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area+headgroup1_2_area+headgroup2_2_area
    +headgroup1_3_area+headgroup2_3_area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                headgroup1_2->fnGetnSLD(dz)*headgroup1_2_area+
                headgroup1_3->fnGetnSLD(dz)*headgroup1_3_area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area+
                headgroup2_2->fnGetnSLD(dz)*headgroup2_2_area+
                headgroup2_3->fnGetnSLD(dz)*headgroup2_3_area
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM_ternary::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM_ternary::fnGetUpperLimit() 
{
	double a,b,c;
	a=headgroup2->fnGetUpperLimit();
	b=headgroup2_2->fnGetUpperLimit();
	c=headgroup2_3->fnGetUpperLimit();
	
	if (a>b) {
		if (a>c) return a; else return c;}
	else{
		if (b>c) return b; else return c;}
    
};

void tBLM_ternary::fnSetSigma(double sigma)
{
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    headgroup1_2->fnSetSigma(sigma);
    headgroup1_3->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->fnSetSigma(sigma); 
    lipid1->fnSetSigma(sigma);
    methyl1->fnSetSigma(sigma);
    methyl2->fnSetSigma(sigma);
    lipid2->fnSetSigma(sigma);
    headgroup2->fnSetSigma(sigma);
    headgroup2_2->fnSetSigma(sigma);
    headgroup2_3->fnSetSigma(sigma);
    
}

double tBLM_ternary::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_ternary::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    headgroup1_2->fnWriteGroup2File(fp, "headgroup1_2", dimension, stepsize);
    headgroup1_3->fnWriteGroup2File(fp, "headgroup1_3", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    headgroup2_2->fnWriteGroup2File(fp, "headgroup2_2", dimension, stepsize);
    headgroup2_3->fnWriteGroup2File(fp, "headgroup2_3", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}
//------------------------------------------------------------------------------------------------------
// Lipid bilayer - ternary system
//------------------------------------------------------------------------------------------------------
tBLM_ternary_chol::tBLM_ternary_chol(){
    
	headgroup1_2=  new Box2Err();                                                    //second headgroups
	headgroup2_2 = new Box2Err();                                                   
	
	headgroup1_2->vol=330;       //was 330
	headgroup2_2->vol=330;       //was 330
	headgroup1_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    volchol=630;
    nslchol=1.3215e-4;
    
    fnAdjustParameters();
};

tBLM_ternary_chol::~tBLM_ternary_chol(){
    delete headgroup1_2;
    delete headgroup2_2;
};

void tBLM_ternary_chol::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nf_ohc_lipid_2, nf_ohc_chol, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nf_om_lipid_2, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nf_ihc_lipid_2, nf_ihc_chol, nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_lipid_2, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    fnSetSigma(sigma);
    
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid  =1-nf_lipid_2-nf_chol;
    nf_ohc_lipid_2=nf_lipid_2;
    nf_ohc_chol=nf_chol;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ohc_chol*volchol;
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ohc_chol*nslchol;
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    nf_ihc_chol=(1-nf_ihc_tether)*nf_ohc_chol;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_chol*volchol+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_chol*nslchol+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer headgroups
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    headgroup2_2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2);
    
    //inner headgroups
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    headgroup1_2->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    headgroup1_2->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_2->l);
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    headgroup2_2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_2->l);
    
};

//Return value is area at position z
double tBLM_ternary_chol::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz)+headgroup1_2->fnGetArea(dz)+headgroup2_2->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM_ternary_chol::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
	double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    double headgroup1_2_area, headgroup2_2_area;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    headgroup1_2_area=headgroup1_2->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    headgroup2_2_area=headgroup2_2->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area+headgroup1_2_area+headgroup2_2_area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                headgroup1_2->fnGetnSLD(dz)*headgroup1_2_area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area+
                headgroup2_2->fnGetnSLD(dz)*headgroup2_2_area
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM_ternary_chol::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM_ternary_chol::fnGetUpperLimit() 
{
	double a,b;
	a=headgroup2->fnGetUpperLimit();
	b=headgroup2_2->fnGetUpperLimit();
	
	if (a>b) return a;
	else return b;
    
};

void tBLM_ternary_chol::fnSetSigma(double sigma)
{
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    headgroup1_2->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->fnSetSigma(sigma); 
    lipid1->fnSetSigma(sigma);
    methyl1->fnSetSigma(sigma);
    methyl2->fnSetSigma(sigma);
    lipid2->fnSetSigma(sigma);
    headgroup2->fnSetSigma(sigma);
    headgroup2_2->fnSetSigma(sigma);
    
}

double tBLM_ternary_chol::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_ternary_chol::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    headgroup1_2->fnWriteGroup2File(fp, "headgroup1_2", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    headgroup2_2->fnWriteGroup2File(fp, "headgroup2_2", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}
//------------------------------------------------------------------------------------------------------
// Lipid bilayer - quaternary system
//------------------------------------------------------------------------------------------------------
tBLM_quaternary_chol::tBLM_quaternary_chol(){
    
	headgroup1_2=  new Box2Err();                                                    //second headgroups
	headgroup2_2 = new Box2Err();
	headgroup1_3 = new Box2Err();
	headgroup2_3 = new Box2Err();                                                   
	
	headgroup1_2->vol=330;       //was 330
	headgroup2_2->vol=330;       //was 330
	headgroup1_3->vol=330;       //was 330
	headgroup2_3->vol=330;       //was 330
	headgroup1_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	headgroup1_3->l=9.5;
	headgroup2_3->l=9.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
	volacyllipid_3=925;
    nslacyllipid_3=-2.67e-4;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    volchol=630;
    nslchol=1.3215e-4;
    
    fnAdjustParameters();
};

tBLM_quaternary_chol::~tBLM_quaternary_chol(){
    delete headgroup1_2;
    delete headgroup2_2;
    delete headgroup1_3;
    delete headgroup2_3;
};

void tBLM_quaternary_chol::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nf_ohc_lipid_2, nf_ohc_lipid_3, nf_ohc_chol, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nf_om_lipid_2, nf_om_lipid_3, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nf_ihc_lipid_2, nf_ihc_lipid_3, nf_ihc_chol, nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_lipid_2, nf_im_lipid_3, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    fnSetSigma(sigma);
    
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid  =1-nf_lipid_2-nf_lipid_3-nf_chol;
    nf_ohc_lipid_2=nf_lipid_2;
    nf_ohc_lipid_3=nf_lipid_3;
    nf_ohc_chol=nf_chol;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ohc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ohc_chol*volchol;
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ohc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ohc_chol*nslchol;
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc;
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    nf_om_lipid_3=nf_ohc_lipid_3;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2+nf_om_lipid_3*volmethyllipid_3;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2+nf_om_lipid_3*nslmethyllipid_3;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om;
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    nf_ihc_lipid_3=(1-nf_ihc_tether)*nf_ohc_lipid_3;
    nf_ihc_chol=(1-nf_ihc_tether)*nf_ohc_chol;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ihc_chol*volchol+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ihc_chol*nslchol+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_lipid_3=nf_ihc_lipid_3;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_lipid_3*volmethyllipid_3+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_lipid_3*nslmethyllipid_3+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im;
    
    //outer headgroups
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2);
    headgroup2_2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2);
    headgroup2_3->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_3*(1-hc_substitution_2);
    
    //inner headgroups
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1);
    headgroup1_2->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1);
    headgroup1_3->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_3*(1-hc_substitution_1);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    //substrate
    substrate->vol=normarea*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    headgroup1_2->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_2->l);
    headgroup1_3->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_3->l);
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    headgroup2_2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_2->l);
    headgroup2_3->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_3->l);
    
};

//Return value is area at position z
double tBLM_quaternary_chol::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz)+headgroup1_2->fnGetArea(dz)+headgroup2_2->fnGetArea(dz)
            +headgroup1_3->fnGetArea(dz)+headgroup2_3->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM_quaternary_chol::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
	double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    double headgroup1_2_area, headgroup2_2_area, headgroup1_3_area, headgroup2_3_area;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    headgroup1_2_area=headgroup1_2->fnGetArea(dz);
    headgroup1_3_area=headgroup1_3->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    headgroup2_2_area=headgroup2_2->fnGetArea(dz);
    headgroup2_3_area=headgroup2_3->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area+headgroup1_2_area+headgroup2_2_area+headgroup1_3_area+headgroup2_3_area;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                headgroup1_2->fnGetnSLD(dz)*headgroup1_2_area+
                headgroup1_3->fnGetnSLD(dz)*headgroup1_3_area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area+
                headgroup2_2->fnGetnSLD(dz)*headgroup2_2_area+
                headgroup2_3->fnGetnSLD(dz)*headgroup2_3_area
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM_quaternary_chol::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM_quaternary_chol::fnGetUpperLimit() 
{
	double a,b,c;
	a=headgroup2->fnGetUpperLimit();
	b=headgroup2_2->fnGetUpperLimit();
	c=headgroup2_3->fnGetUpperLimit();
	
	if (a>b) {
		if (a>c) return a;
		else return c;}
	else {
		if (b>c) return b;
		else return c;}
};

void tBLM_quaternary_chol::fnSetSigma(double sigma)
{
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    headgroup1_2->fnSetSigma(sigma);
    headgroup1_3->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->fnSetSigma(sigma); 
    lipid1->fnSetSigma(sigma);
    methyl1->fnSetSigma(sigma);
    methyl2->fnSetSigma(sigma);
    lipid2->fnSetSigma(sigma);
    headgroup2->fnSetSigma(sigma);
    headgroup2_2->fnSetSigma(sigma);
    headgroup2_3->fnSetSigma(sigma);
    
}

double tBLM_quaternary_chol::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_quaternary_chol::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    headgroup1_2->fnWriteGroup2File(fp, "headgroup1_2", dimension, stepsize);
    headgroup1_3->fnWriteGroup2File(fp, "headgroup1_3", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    headgroup2_2->fnWriteGroup2File(fp, "headgroup2_2", dimension, stepsize);
    headgroup2_3->fnWriteGroup2File(fp, "headgroup2_3", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}

//------------------------------------------------------------------------------------------------------
// Lipid bilayer - quaternary system with domains
//------------------------------------------------------------------------------------------------------
tBLM_quaternary_chol_domain::tBLM_quaternary_chol_domain(){
    
	headgroup1_2        = new Box2Err();                                              //second headgroups
	headgroup2_2        = new Box2Err();
	headgroup1_3        = new Box2Err();
	headgroup2_3        = new Box2Err();                                                   
    headgroup1_domain   = new PCm();
    lipid1_domain       = new Box2Err();
    methyl1_domain      = new Box2Err();
    methyl2_domain      = new Box2Err();
    lipid2_domain       = new Box2Err();
    headgroup2_domain   = new PC();
    headgroup1_2_domain = new Box2Err();	  
    headgroup2_2_domain = new Box2Err();
    headgroup1_3_domain = new Box2Err();	  
    headgroup2_3_domain = new Box2Err();
    tetherg_domain      = new Box2Err();
    tether_domain       = new Box2Err();
	
	headgroup1_2->vol=330;       //was 330
	headgroup2_2->vol=330;       //was 330
	headgroup1_3->vol=330;       //was 330
	headgroup2_3->vol=330;       //was 330
	headgroup1_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_3->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	headgroup1_3->l=9.5;
	headgroup2_3->l=9.5;
	headgroup1_2_domain->vol=330;       //was 330
	headgroup2_2_domain->vol=330;       //was 330
	headgroup1_3_domain->vol=330;       //was 330
	headgroup2_3_domain->vol=330;       //was 330
	headgroup1_2_domain->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_2_domain->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_3_domain->nSL=6.0012e-4; // was 6.0122e-4
	headgroup2_3_domain->nSL=6.0012e-4; // was 6.0122e-4
	headgroup1_2_domain->l=9.5;
	headgroup2_2_domain->l=9.5;
	headgroup1_3_domain->l=9.5;
	headgroup2_3_domain->l=9.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
	volacyllipid_3=925;
    nslacyllipid_3=-2.67e-4;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    volchol=630;
    nslchol=1.3215e-4;
    
    fnAdjustParameters();
};

tBLM_quaternary_chol_domain::~tBLM_quaternary_chol_domain(){
    delete headgroup1_2;
    delete headgroup2_2;
    delete headgroup1_3;
    delete headgroup2_3;
    delete headgroup1_domain;
    delete lipid1_domain;
    delete methyl1_domain;
    delete methyl2_domain;
    delete lipid2_domain;
    delete headgroup2_domain;
    delete headgroup1_2_domain;	  
    delete headgroup2_2_domain;
    delete headgroup1_3_domain;	  
    delete headgroup2_3_domain;
    delete tether_domain;
    delete tetherg_domain;
};

void tBLM_quaternary_chol_domain::fnAdjustParameters(){
    
    double l_ohc;
    double V_ohc;
    double nf_ohc_lipid, nf_ohc_lipid_2, nf_ohc_lipid_3, nf_ohc_chol, nSL_ohc;
    double c_s_ohc, c_A_ohc, c_V_ohc;
    
    double l_om;
    double V_om;
    double nf_om_lipid, nf_om_lipid_2, nf_om_lipid_3, nSL_om;
    double c_s_om, c_A_om, c_V_om;  
    
    double l_ihc;
    double V_ihc;
    double nf_ihc_lipid, nf_ihc_lipid_2, nf_ihc_lipid_3, nf_ihc_chol, nf_ihc_tether, nSL_ihc;
    double c_s_ihc, c_A_ihc, c_V_ihc;
    
    double l_im;
    double V_im;
    double nf_im_lipid, nf_im_lipid_2, nf_im_lipid_3, nf_im_tether, nSL_im;
    double c_s_im, c_A_im, c_V_im;  
    
    double V_tg;
    double c_s_tg, c_A_tg, c_V_tg;
    
    double l_EO,V_EO;
    double c_s_EO, c_A_EO, c_V_EO;
    
    double l_bME,V_bME;
    
    double d1;
    
    
    fnSetSigma(sigma);
    
    //----------------first domain----------------
    
    //outer hydrocarbons
    l_ohc=l_lipid2;
    nf_ohc_lipid  =1-nf_lipid_2-nf_lipid_3-nf_chol;
    nf_ohc_lipid_2=nf_lipid_2;
    nf_ohc_lipid_3=nf_lipid_3;
    nf_ohc_chol=nf_chol;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ohc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ohc_chol*volchol;
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ohc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ohc_chol*nslchol;
    
    normarea=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2->l=l_ohc;
    lipid2->vol=V_ohc;
    lipid2->nSL=nSL_ohc;
    lipid2->nf=c_s_ohc*c_A_ohc*c_V_ohc*(1-frac_domain);
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    nf_om_lipid_3=nf_ohc_lipid_3;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2+nf_om_lipid_3*volmethyllipid_3;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2+nf_om_lipid_3*nslmethyllipid_3;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2->l=l_om;
    methyl2->vol=V_om;
    methyl2->nSL=nSL_om;
    methyl2->nf=c_s_om*c_A_om*c_V_om*(1-frac_domain);
    
    //inner hydrocarbons
    l_ihc=l_lipid1;
    
    nf_ihc_tether=nf_tether;
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    nf_ihc_lipid_3=(1-nf_ihc_tether)*nf_ohc_lipid_3;
    nf_ihc_chol=(1-nf_ihc_tether)*nf_ohc_chol;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ihc_chol*volchol+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ihc_chol*nslchol+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1->l=l_ihc;
    lipid1->vol=V_ihc;
    lipid1->nSL=nSL_ihc;
    lipid1->nf=c_s_ihc*c_A_ihc*c_V_ihc*(1-frac_domain);
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_lipid_3=nf_ihc_lipid_3;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_lipid_3*volmethyllipid_3+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_lipid_3*nslmethyllipid_3+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1->l=l_im;
    methyl1->vol=V_im;
    methyl1->nSL=nSL_im;
    methyl1->nf=c_s_im*c_A_im*c_V_im*(1-frac_domain);
    
    //outer headgroups
    headgroup2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2)*(1-frac_domain);
    headgroup2_2->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2)*(1-frac_domain);
    headgroup2_3->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_3*(1-hc_substitution_2)*(1-frac_domain);
    
    //inner headgroups
    headgroup1->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1)*(1-frac_domain);
    headgroup1_2->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1)*(1-frac_domain);
    headgroup1_3->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_3*(1-hc_substitution_1)*(1-frac_domain);
    
    //tether glycerol part
    V_tg=tetherg->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1->l)/0.9;
    tetherg->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether;
    V_EO=tether->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether->nf=c_s_EO*c_A_EO*c_V_EO;
    tether->l=l_EO;
    
    if ((tether->nf*tether->vol/tether->l)>normarea) {
        tether->l=(tether->nf*tether->vol)/normarea;
    }
    
    l_tether=tether->l;
    
    //------------second domain--------------------
    
    //outer hydrocarbons
    l_ohc=l_lipid2_domain;
    nf_ohc_lipid  =1-nf_lipid_2_domain-nf_lipid_3_domain-nf_chol_domain;
    nf_ohc_lipid_2=nf_lipid_2_domain;
    nf_ohc_lipid_3=nf_lipid_3_domain;
    nf_ohc_chol=nf_chol_domain;
    V_ohc=nf_ohc_lipid*(volacyllipid-volmethyllipid)+nf_ohc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ohc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ohc_chol*volchol;
    nSL_ohc=nf_ohc_lipid*(nslacyllipid-nslmethyllipid)+nf_ohc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ohc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ohc_chol*nslchol;
    
    normarea_domain=V_ohc/l_ohc;
    c_s_ohc=vf_bilayer;
    c_A_ohc=1;
    c_V_ohc=1;
    
    lipid2_domain->l=l_ohc;
    lipid2_domain->vol=V_ohc;
    lipid2_domain->nSL=nSL_ohc;
    lipid2_domain->nf=c_s_ohc*c_A_ohc*c_V_ohc*frac_domain;
    
    //outher methyl
    nf_om_lipid  =nf_ohc_lipid;
    nf_om_lipid_2=nf_ohc_lipid_2;
    nf_om_lipid_3=nf_ohc_lipid_3;
    V_om=nf_om_lipid*volmethyllipid+nf_om_lipid_2*volmethyllipid_2+nf_om_lipid_3*volmethyllipid_3;
    l_om=l_ohc*V_om/V_ohc;
    nSL_om=nf_om_lipid*nslmethyllipid+nf_om_lipid_2*nslmethyllipid_2+nf_om_lipid_3*nslmethyllipid_3;
    
    c_s_om=c_s_ohc;
    c_A_om=1;
    c_V_om=1;
    
    methyl2_domain->l=l_om;
    methyl2_domain->vol=V_om;
    methyl2_domain->nSL=nSL_om;
    methyl2_domain->nf=c_s_om*c_A_om*c_V_om*frac_domain;
    
    //inner hydrocarbons
    l_ihc=l_lipid1_domain;
    
    nf_ihc_tether=nf_tether*l_lipid1/l_lipid1_domain;
    nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid;
    nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2;
    nf_ihc_lipid_3=(1-nf_ihc_tether)*nf_ohc_lipid_3;
    nf_ihc_chol=(1-nf_ihc_tether)*nf_ohc_chol;
    V_ihc=nf_ihc_lipid*(volacyllipid-volmethyllipid)+nf_ihc_lipid_2*(volacyllipid_2-volmethyllipid_2)+nf_ihc_lipid_3*(volacyllipid_3-volmethyllipid_3)+nf_ihc_chol*volchol+nf_ihc_tether*(volacyltether-volmethyltether);
    nSL_ihc=nf_ihc_lipid*(nslacyllipid-nslmethyllipid)+nf_ihc_lipid_2*(nslacyllipid_2-nslmethyllipid_2)+nf_ihc_lipid_3*(nslacyllipid_3-nslmethyllipid_3)+nf_ihc_chol*nslchol+nf_ihc_tether*(nslacyltether-nslmethyltether);
    
    c_s_ihc=vf_bilayer;
    c_A_ihc=normarea_domain*l_ihc/V_ihc;
    c_V_ihc=1;
    
    lipid1_domain->l=l_ihc;
    lipid1_domain->vol=V_ihc;
    lipid1_domain->nSL=nSL_ihc;
    lipid1_domain->nf=c_s_ihc*c_A_ihc*c_V_ihc*frac_domain;
    
    //inner methyl
    nf_im_lipid=nf_ihc_lipid;
    nf_im_lipid_2=nf_ihc_lipid_2;
    nf_im_lipid_3=nf_ihc_lipid_3;
    nf_im_tether=nf_ihc_tether;
    V_im=nf_im_lipid*volmethyllipid+nf_im_lipid_2*volmethyllipid_2+nf_im_lipid_3*volmethyllipid_3+nf_im_tether*volmethyltether;
    l_im=l_ihc*V_im/V_ihc;
    nSL_im=nf_im_lipid*nslmethyllipid+nf_im_lipid_2*nslmethyllipid_2+nf_im_lipid_3*nslmethyllipid_3+nf_im_tether*nslmethyltether;
    
    c_s_im=c_s_ihc;
    c_A_im=c_A_ihc;
    c_V_im=1;
    
    methyl1_domain->l=l_im;
    methyl1_domain->vol=V_im;
    methyl1_domain->nSL=nSL_im;
    methyl1_domain->nf=c_s_im*c_A_im*c_V_im*frac_domain;
    
    //outer headgroups
    headgroup2_domain->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid*(1-hc_substitution_2)*frac_domain;
    headgroup2_2_domain->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_2*(1-hc_substitution_2)*frac_domain;
    headgroup2_3_domain->nf=c_s_ohc*c_A_ohc*nf_ohc_lipid_3*(1-hc_substitution_2)*frac_domain;
    
    //inner headgroups
    headgroup1_domain->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid*(1-hc_substitution_1)*frac_domain;
    headgroup1_2_domain->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_2*(1-hc_substitution_1)*frac_domain;
    headgroup1_3_domain->nf=c_s_ihc*c_A_ihc*nf_ihc_lipid_3*(1-hc_substitution_1)*frac_domain;
    
    //tether glycerol part
    V_tg=tetherg_domain->vol;
    
    c_s_tg=c_s_ihc;
    c_A_tg=c_A_ihc;
    c_V_tg=nf_ihc_tether;
    
    tetherg_domain->l=tetherg->vol/((volacyltether-volmethyltether)/lipid1_domain->l)/0.9;
    tetherg_domain->nf=c_s_tg*c_A_tg*c_V_tg;
    
    //tether EO part
    l_EO=l_tether_domain;
    V_EO=tether_domain->vol;
    
    c_s_EO=c_s_ihc;
    c_A_EO=c_A_ihc;
    c_V_EO=nf_ihc_tether;
    
    tether_domain->nf=c_s_EO*c_A_EO*c_V_EO;
    tether_domain->l=l_EO;
    
    if ((tether_domain->nf*tether->vol/tether_domain->l)>normarea_domain) {
        tether_domain->l=(tether_domain->nf*tether_domain->vol)/normarea_domain;
    }
    
    l_tether_domain=tether_domain->l;
    
    //-------------substrate----------------    
    
    //bME
    bME->l=5.2;
    l_bME=bME->l;
    headgroup1->l=9.575;
    V_bME=bME->vol;
    
    bME->nf=((1-frac_domain)*tether->nf+frac_domain*tether_domain->nf)*mult_tether; //2.333;

    
    d1=headgroup1->l+bME->l-tether->l-tetherg->l;
    if (d1>0) {
        bME->l=bME->l-d1/2;
        headgroup1->l=headgroup1->l-d1/2;
    }
    
    
    if ((tether->nf*tether->vol/tether->l+mult_tether*tether->nf*bME->vol/bME->l)>normarea) {
		mult_tether=((normarea-tether->nf*tether->vol/tether->l)/(bME->vol/bME->l))/tether->nf;
		if (mult_tether<0) {
			mult_tether=0;
		}
    }
    
    bME->nf=tether->nf*mult_tether; //2.333;
    
    
    
    //substrate
    substrate->vol=((1-frac_domain)*normarea+frac_domain*normarea_domain)*substrate->l;
    substrate->nSL=rho_substrate*substrate->vol;
    
    
    // set all lengths
    bME->z=0.5*bME->l+substrate->l;
    
    tether->z=0.5*tether->l+substrate->l;
    tetherg->z=tether->z+0.5*tether->l+0.5*tetherg->l;
    lipid1->z=tetherg->z+0.5*(tetherg->l+lipid1->l);
    headgroup1->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1->l);
    headgroup1_2->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_2->l);
    headgroup1_3->fnSetZ(lipid1->z-0.5*lipid1->l-0.5*headgroup1_3->l);
    methyl1->z=lipid1->z+0.5*(lipid1->l+methyl1->l);
    methyl2->z=methyl1->z+0.5*(methyl1->l+methyl2->l);
    lipid2->z=methyl2->z+0.5*(methyl2->l+lipid2->l);
    headgroup2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2->l);
    headgroup2_2->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_2->l);
    headgroup2_3->fnSetZ(lipid2->z+0.5*lipid2->l+0.5*headgroup2_3->l);

    tether_domain->z=0.5*tether_domain->l+substrate->l;
    tetherg_domain->z=tether_domain->z+0.5*tether_domain->l+0.5*tetherg_domain->l;
    lipid1_domain->z=tetherg_domain->z+0.5*(tetherg_domain->l+lipid1_domain->l);
    headgroup1_domain->fnSetZ(lipid1_domain->z-0.5*lipid1_domain->l-0.5*headgroup1_domain->l);
    headgroup1_2_domain->fnSetZ(lipid1_domain->z-0.5*lipid1_domain->l-0.5*headgroup1_2_domain->l);
    headgroup1_3_domain->fnSetZ(lipid1_domain->z-0.5*lipid1_domain->l-0.5*headgroup1_3_domain->l);
    methyl1_domain->z=lipid1_domain->z+0.5*(lipid1_domain->l+methyl1_domain->l);
    methyl2_domain->z=methyl1_domain->z+0.5*(methyl1_domain->l+methyl2_domain->l);
    lipid2_domain->z=methyl2_domain->z+0.5*(methyl2_domain->l+lipid2_domain->l);
    headgroup2_domain->fnSetZ(lipid2_domain->z+0.5*lipid2_domain->l+0.5*headgroup2_domain->l);
    headgroup2_2_domain->fnSetZ(lipid2_domain->z+0.5*lipid2_domain->l+0.5*headgroup2_2_domain->l);
    headgroup2_3_domain->fnSetZ(lipid2_domain->z+0.5*lipid2_domain->l+0.5*headgroup2_3_domain->l);

};

//Return value is area at position z
double tBLM_quaternary_chol_domain::fnGetArea(double dz) {
    return (substrate->fnGetArea(dz)+bME->fnGetArea(dz)+tether->fnGetArea(dz)
            +tetherg->fnGetArea(dz)+lipid1->fnGetArea(dz)+headgroup1->fnGetArea(dz)
            +methyl1->fnGetArea(dz)+methyl2->fnGetArea(dz)+lipid2->fnGetArea(dz)
            +headgroup2->fnGetArea(dz)+headgroup1_2->fnGetArea(dz)+headgroup2_2->fnGetArea(dz)
            +headgroup1_3->fnGetArea(dz)+headgroup2_3->fnGetArea(dz)+tether_domain->fnGetArea(dz)
            +tetherg_domain->fnGetArea(dz)+lipid1_domain->fnGetArea(dz)+headgroup1_domain->fnGetArea(dz)
            +methyl1_domain->fnGetArea(dz)+methyl2_domain->fnGetArea(dz)+lipid2_domain->fnGetArea(dz)
            +headgroup2_domain->fnGetArea(dz)+headgroup1_2_domain->fnGetArea(dz)+headgroup2_2_domain->fnGetArea(dz)
            +headgroup1_3_domain->fnGetArea(dz)+headgroup2_3_domain->fnGetArea(dz));
};

//get nSLD from molecular subgroups
double tBLM_quaternary_chol_domain::fnGetnSLD(double dz) {
    double substratearea, bMEarea, tetherarea, tethergarea, lipid1area, headgroup1area;
	double methyl1area, methyl2area, lipid2area, headgroup2area, sum;
    double headgroup1_2_area, headgroup2_2_area, headgroup1_3_area, headgroup2_3_area;
    double tetherarea_domain, tethergarea_domain, lipid1area_domain, headgroup1area_domain;
	double methyl1area_domain, methyl2area_domain, lipid2area_domain, headgroup2area_domain;
    double headgroup1_2_area_domain, headgroup2_2_area_domain, headgroup1_3_area_domain, headgroup2_3_area_domain;
    
    substratearea=substrate->fnGetArea(dz);
    bMEarea=bME->fnGetArea(dz);
    tetherarea=tether->fnGetArea(dz);
    tethergarea=tetherg->fnGetArea(dz);
    lipid1area=lipid1->fnGetArea(dz);
    headgroup1area=headgroup1->fnGetArea(dz);
    headgroup1_2_area=headgroup1_2->fnGetArea(dz);
    headgroup1_3_area=headgroup1_3->fnGetArea(dz);
    methyl1area=methyl1->fnGetArea(dz);
    methyl2area=methyl2->fnGetArea(dz);
    lipid2area=lipid2->fnGetArea(dz);
    headgroup2area=headgroup2->fnGetArea(dz);
    headgroup2_2_area=headgroup2_2->fnGetArea(dz);
    headgroup2_3_area=headgroup2_3->fnGetArea(dz);
    tetherarea_domain=tether_domain->fnGetArea(dz);
    tethergarea_domain=tetherg_domain->fnGetArea(dz);
    lipid1area_domain=lipid1_domain->fnGetArea(dz);
    headgroup1area_domain=headgroup1_domain->fnGetArea(dz);
    headgroup1_2_area_domain=headgroup1_2_domain->fnGetArea(dz);
    headgroup1_3_area_domain=headgroup1_3_domain->fnGetArea(dz);
    methyl1area_domain=methyl1_domain->fnGetArea(dz);
    methyl2area_domain=methyl2_domain->fnGetArea(dz);
    lipid2area_domain=lipid2_domain->fnGetArea(dz);
    headgroup2area_domain=headgroup2_domain->fnGetArea(dz);
    headgroup2_2_area_domain=headgroup2_2_domain->fnGetArea(dz);
    headgroup2_3_area_domain=headgroup2_3_domain->fnGetArea(dz);
    
    sum=substratearea+bMEarea+tetherarea+tethergarea+lipid1area+headgroup1area+methyl1area
    +methyl2area+lipid2area+headgroup2area+headgroup1_2_area+headgroup2_2_area+headgroup1_3_area+headgroup2_3_area+
    tetherarea_domain+tethergarea_domain+lipid1area_domain+headgroup1area_domain+methyl1area_domain
    +methyl2area_domain+lipid2area_domain+headgroup2area_domain+headgroup1_2_area_domain+headgroup2_2_area_domain
    +headgroup1_3_area_domain+headgroup2_3_area_domain;
    
    if (sum==0) {return 0;}
    else {
        return (
                substrate->fnGetnSLD(dz)*substratearea+
                bME->fnGetnSLD(dz)*bMEarea+
                tether->fnGetnSLD(dz)*tetherarea+
                tetherg->fnGetnSLD(dz)*tethergarea+
                headgroup1->fnGetnSLD(dz)*headgroup1area+
                headgroup1_2->fnGetnSLD(dz)*headgroup1_2_area+
                headgroup1_3->fnGetnSLD(dz)*headgroup1_3_area+
                lipid1->fnGetnSLD(dz)*lipid1area+
                methyl1->fnGetnSLD(dz)*methyl1area+
                methyl2->fnGetnSLD(dz)*methyl2area+
                lipid2->fnGetnSLD(dz)*lipid2area+
                headgroup2->fnGetnSLD(dz)*headgroup2area+
                headgroup2_2->fnGetnSLD(dz)*headgroup2_2_area+
                headgroup2_3->fnGetnSLD(dz)*headgroup2_3_area+
                tether_domain->fnGetnSLD(dz)*tetherarea_domain+
                tetherg_domain->fnGetnSLD(dz)*tethergarea_domain+
                headgroup1_domain->fnGetnSLD(dz)*headgroup1area_domain+
                headgroup1_2_domain->fnGetnSLD(dz)*headgroup1_2_area_domain+
                headgroup1_3_domain->fnGetnSLD(dz)*headgroup1_3_area_domain+
                lipid1_domain->fnGetnSLD(dz)*lipid1area_domain+
                methyl1_domain->fnGetnSLD(dz)*methyl1area_domain+
                methyl2_domain->fnGetnSLD(dz)*methyl2area_domain+
                lipid2_domain->fnGetnSLD(dz)*lipid2area_domain+
                headgroup2_domain->fnGetnSLD(dz)*headgroup2area_domain+
                headgroup2_2_domain->fnGetnSLD(dz)*headgroup2_2_area_domain+
                headgroup2_3_domain->fnGetnSLD(dz)*headgroup2_3_area_domain
                )/sum;
    }
};

//Use limits of molecular subgroups
double tBLM_quaternary_chol_domain::fnGetLowerLimit() {return substrate->fnGetLowerLimit();};
double tBLM_quaternary_chol_domain::fnGetUpperLimit() 
{
	double a,b,c,d,e,f,temp;
	a=headgroup2->fnGetUpperLimit();
	b=headgroup2_2->fnGetUpperLimit();
	c=headgroup2_3->fnGetUpperLimit();
	d=headgroup2_domain->fnGetUpperLimit();
	e=headgroup2_2_domain->fnGetUpperLimit();
	f=headgroup2_3_domain->fnGetUpperLimit();
	
    temp=a;
	if (b>temp) {temp=b;}
    if (c>temp) {temp=c;}
    if (d>temp) {temp=d;}
    if (e>temp) {temp=e;}
    if (f>temp) {temp=f;}
    
    return temp;
};

double tBLM_quaternary_chol_domain::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};

void tBLM_quaternary_chol_domain::fnSetSigma(double sigma)
{
    // set all sigma
    
    //sigma=sqrt(2.4*2.4 + global_rough*global_rough);
    
    substrate->sigma2=global_rough;
    bME->sigma1=global_rough;
    bME->sigma2=global_rough;
    headgroup1->fnSetSigma(sigma);
    headgroup1_2->fnSetSigma(sigma);
    headgroup1_3->fnSetSigma(sigma);
    tether->sigma1=global_rough; 
    tether->sigma2=sigma;
    tetherg->fnSetSigma(sigma); 
    lipid1->fnSetSigma(sigma);
    methyl1->fnSetSigma(sigma);
    methyl2->fnSetSigma(sigma);
    lipid2->fnSetSigma(sigma);
    headgroup2->fnSetSigma(sigma);
    headgroup2_2->fnSetSigma(sigma);
    headgroup2_3->fnSetSigma(sigma);
    headgroup1_domain->fnSetSigma(sigma);
    headgroup1_2_domain->fnSetSigma(sigma);
    headgroup1_3_domain->fnSetSigma(sigma);
    lipid1_domain->fnSetSigma(sigma);
    methyl1_domain->fnSetSigma(sigma);
    methyl2_domain->fnSetSigma(sigma);
    lipid2_domain->fnSetSigma(sigma);
    headgroup2_domain->fnSetSigma(sigma);
    headgroup2_2_domain->fnSetSigma(sigma);
    headgroup2_3_domain->fnSetSigma(sigma);
    tether_domain->sigma1=global_rough; 
    tether_domain->sigma2=sigma;
    tetherg_domain->fnSetSigma(sigma); 
    
}


void tBLM_quaternary_chol_domain::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
    //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    substrate->fnWriteGroup2File(fp, "substrate", dimension, stepsize);
    bME->fnWriteGroup2File(fp, "bME", dimension, stepsize);
    tether->fnWriteGroup2File(fp, "tether", dimension, stepsize);
    tetherg->fnWriteGroup2File(fp, "tetherg", dimension, stepsize);
    headgroup1->fnWriteGroup2File(fp, "headgroup1", dimension, stepsize);
    headgroup1_2->fnWriteGroup2File(fp, "headgroup1_2", dimension, stepsize);
    headgroup1_3->fnWriteGroup2File(fp, "headgroup1_3", dimension, stepsize);
    lipid1->fnWriteGroup2File(fp, "lipid1", dimension, stepsize);
    methyl1->fnWriteGroup2File(fp, "methyl1", dimension, stepsize);
    methyl2->fnWriteGroup2File(fp, "methyl2", dimension, stepsize);
    lipid2->fnWriteGroup2File(fp, "lipid2", dimension, stepsize);
    headgroup2->fnWriteGroup2File(fp, "headgroup2", dimension, stepsize);
    headgroup2_2->fnWriteGroup2File(fp, "headgroup2_2", dimension, stepsize);
    headgroup2_3->fnWriteGroup2File(fp, "headgroup2_3", dimension, stepsize);
    tether_domain->fnWriteGroup2File(fp, "tether_domain", dimension, stepsize);
    tetherg_domain->fnWriteGroup2File(fp, "tetherg_domain", dimension, stepsize);
    headgroup1_domain->fnWriteGroup2File(fp, "headgroup1_domain", dimension, stepsize);
    headgroup1_2_domain->fnWriteGroup2File(fp, "headgroup1_2_domain", dimension, stepsize);
    headgroup1_3_domain->fnWriteGroup2File(fp, "headgroup1_3_domain", dimension, stepsize);
    lipid1_domain->fnWriteGroup2File(fp, "lipid1_domain", dimension, stepsize);
    methyl1_domain->fnWriteGroup2File(fp, "methyl1_domain", dimension, stepsize);
    methyl2_domain->fnWriteGroup2File(fp, "methyl2_domain", dimension, stepsize);
    lipid2_domain->fnWriteGroup2File(fp, "lipid2_domain", dimension, stepsize);
    headgroup2_domain->fnWriteGroup2File(fp, "headgroup2_domain", dimension, stepsize);
    headgroup2_2_domain->fnWriteGroup2File(fp, "headgroup2_2_domain", dimension, stepsize);
    headgroup2_3_domain->fnWriteGroup2File(fp, "headgroup2_3_domain", dimension, stepsize);
    fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);
    
    
    //delete []str;
    
}



//---------------------------------------------------------------------------------------------------------
//Discrete nSL, area profile
//---------------------------------------------------------------------------------------------------------

Discrete::Discrete(double dstartposition, double dnormarea, const char *cFileName)
{
    
    double dz, darea, dnSLProt, dnSLDeut;
    int i;
    
    dStartPosition=dstartposition;
	normarea=dnormarea;
    dProtExchange=0;
    dnSLDBulkSolvent=-0.566e-6;
    
    FILE *fp;									
    fp=fopen(cFileName,"r");
    
    if(fp==NULL) {
        printf("Error: can't open file.\n");}
    
    
    iNumberOfPoints=0;
    while(!feof(fp)) {
        fscanf(fp, "%lf %lf %lf %lf", &dz, &darea, &dnSLProt, &dnSLDeut);
        //printf("%lf %lf %e %e \n", dz, darea, dnSLProt, dnSLDeut);
        iNumberOfPoints+=1;
    }
    
    rewind(fp);
    zcoord  = new double[iNumberOfPoints];
    area    = new double[iNumberOfPoints];
    nSLProt = new double[iNumberOfPoints];
    nSLDeut = new double[iNumberOfPoints];
    
    i=0;
    while(!feof(fp)) {
        fscanf(fp, "%lf %lf %lf %lf", &dz, &darea, &dnSLProt, &dnSLDeut);
        
        zcoord[i]=dz;
        area[i]=darea;
        nSLProt[i]=dnSLProt;
        nSLDeut[i]=dnSLDeut;
        
        i+=1;
    }
    
    fclose(fp);
    
    dZSpacing=zcoord[1]-zcoord[0];
    
};

Discrete::~Discrete(){
    delete [] zcoord;
    delete [] area;
    delete [] nSLProt;
    delete [] nSLDeut;
};

//Return value is area at position z
double Discrete::fnGetArea(double dz) {
    int iBinLow, iBinHigh;
    double dFraction, dtemp;
    
    dz=dz-dStartPosition;                       //internal z for profile
    dz=dz/dZSpacing;                            //floating point bin
    
    iBinHigh=int(ceil(dz));                     //get bin for dz
    dFraction=modf(dz,&dtemp);
    iBinLow=int(dtemp);
    
    if ((iBinLow>=0) && (iBinHigh<iNumberOfPoints)) {
        
        return((dFraction*area[iBinLow]+(1-dFraction)*area[iBinHigh])*nf);
        
    }
    else {
        return 0;
    }
};

//get nSLD from molecular subgroups
double Discrete::fnGetnSLD(double dz) {
    
    int iBinLow, iBinHigh;
    double dFraction, dtemp1, dtemp2, dtemp3,dtemp4;
    
    dz=dz-dStartPosition;                       //internal z for profile
    dz=dz/dZSpacing;                            //floating point bin
    
    iBinHigh=int(ceil(dz));                     //get bin for dz
    dFraction=modf(dz,&dtemp1);
    iBinLow=int(dtemp1);
    
    if((iBinLow>=0) && (iBinHigh<iNumberOfPoints)) {
        
        dtemp1=dFraction*nSLProt[iBinLow]+(1-dFraction)*nSLProt[iBinHigh];
        dtemp2=dFraction*nSLDeut[iBinLow]+(1-dFraction)*nSLDeut[iBinHigh];
        
        dtemp3=dProtExchange*(dnSLDBulkSolvent+0.566e-6)/(6.34e-6+0.566e-6);
        dtemp4=(dFraction*area[iBinLow]+(1-dFraction)*area[iBinHigh])*dZSpacing;
        
        //printf ("nf %e binlow %i binhigh %i dtemp1 %e dtemp2 %e dtemp3 %e dtemp4 %e areahigh % e protexch %e nSLDBulk %e nSLD %e \n", nf, iBinLow, iBinHigh, dtemp1, dtemp2, dtemp3, dtemp4, area[iBinHigh],dProtExchange,dnSLDBulkSolvent, (((1-dtemp3)*dtemp1+dtemp3*dtemp2)/dtemp4));
        
        return(((1-dtemp3)*dtemp1+dtemp3*dtemp2)/dtemp4);
        
    }
    else {
        return 0;
    }
    
};

//Use limits of molecular subgroups
double Discrete::fnGetLowerLimit() {return (dStartPosition);}
double Discrete::fnGetUpperLimit() {return (dStartPosition+double(iNumberOfPoints)*dZSpacing);}


void Discrete::fnSetNormarea(double dnormarea)
{
    normarea=dnormarea;
};


void Discrete::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Discrete %s StartPosition %e \n",cName, dStartPosition);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);    
}

DiscreteEuler::DiscreteEuler(double dstartposition, double dnormarea, double BetaStart, double BetaEnd, double BetaInc,
                             double GammaStart, double GammaEnd, double GammaInc, char* strFileNameRoot, 
                             char* strFileNameBeta, char* strFileNameGamma, char* strFileNameEnding)
{
    
    double dz, dzold, darea, dnSLProt, dnSLDeut, dB, dG;
    int i,j,k;
    FILE *fp;									
    char strFilename[200], buf[20];
    
    
    dStartPosition=dstartposition;
	normarea=dnormarea;
    dProtExchange=0;
    dnSLDBulkSolvent=-0.566e-6;
    
    dBetaStart=BetaStart;
    dBetaEnd=BetaEnd;
    dBetaInc=BetaInc;
    dGammaStart=GammaStart;
    dGammaEnd=GammaEnd;
    dGammaInc=GammaInc;
    
    strcpy(strFilename,"");
    strcat(strFilename,strFileNameRoot);
    strcat(strFilename,strFileNameBeta);
    sprintf(buf,"%f",dBetaStart);
    strcat(strFilename,buf);
    strcat(strFilename,strFileNameGamma);
    sprintf(buf,"%f",dGammaStart);
    strcat(strFilename,buf);
    strcat(strFilename,strFileNameEnding);
    
    fp=fopen(strFilename,"r");
    
    if(fp==NULL) {
        printf("Error: can't open file.\n");}
    iNumberOfPoints=0;
    
    j=0;
    while(!feof(fp)) {
        fscanf(fp, "%lf %lf %lf %lf", &dz, &darea, &dnSLProt, &dnSLDeut);
        //printf("%lf %lf %e %e \n", dz, darea, dnSLProt, dnSLDeut);
        if (j==1){
            dZSpacing=dzold-dz;
        }
        iNumberOfPoints+=1;
        dzold=dz;
        j++;
    }
    
    fclose(fp);
    
    iNumberOfBeta=int((dBetaEnd-dBetaStart)/dBetaInc);
    iNumberOfGamma=int((dGammaEnd-dGammaStart)/dGammaInc);
                      
    j=fn3Cto1C(iNumberOfBeta,iNumberOfGamma,iNumberOfPoints);
    zcoord  = new double[j];
    area    = new double[j];
    nSLProt = new double[j];
    nSLDeut = new double[j];
    
    j=0;
    for (dB=dBetaStart; dB<dBetaEnd; dB+=dBetaInc) {
        k=0;
        for (dG=dGammaStart; dG<dGammaEnd; dG+=dGammaInc) {
            
            strcpy(strFilename,"");
            strcat(strFilename,strFileNameRoot);
            strcat(strFilename,strFileNameBeta);
            sprintf(buf,"%f",dB);
            strcat(strFilename,buf);
            strcat(strFilename,strFileNameGamma);
            sprintf(buf,"%f",dG);
            strcat(strFilename,buf);
            strcat(strFilename,strFileNameEnding);
            
            fp=fopen(strFilename,"r");
            
            if(fp==NULL) {
                printf("Error: can't open file.\n");}
            
            
            i=0;
            while(!feof(fp)) {
                fscanf(fp, "%lf %lf %lf %lf", &dz, &darea, &dnSLProt, &dnSLDeut);
                
                zcoord [fn3Cto1C(j,k,i)]=dz;
                area   [fn3Cto1C(j,k,i)]=darea;
                nSLProt[fn3Cto1C(j,k,i)]=dnSLProt;
                nSLDeut[fn3Cto1C(j,k,i)]=dnSLDeut;
                
                i++;
            }            
            k++;
        }
        j++;
    }
};

DiscreteEuler::~DiscreteEuler(){
    delete [] zcoord;
    delete [] area;
    delete [] nSLProt;
    delete [] nSLDeut;
};

//Do coordinate conversion to go from 3D array to 1D array
int DiscreteEuler::fn3Cto1C(int c1, int c2, int c3) {
    return c1*iNumberOfBeta+c2*iNumberOfGamma+c3*iNumberOfPoints;
}

//Return value is area at position z
double DiscreteEuler::fnGetArea(double dz) {
    
    int i,j,k, ii, jj, kk;
    int iPosBinLow, iPosBinHigh;
    int iBetaBinLow, iBetaBinHigh;
    int iGammaBinLow, iGammaBinHigh;
    double dPosT, dBetaT, dGammaT;
    double returnvalue, dtemp;
    double t[3];                                // tvalues for beta, gamma, position
    double p[4][4][4];                          // all function values for tricubic spline interpolations
                                                // [beta][gamma][position]
                                                // index values: 0->-1, 1->0, 2->1, 3->2
    
    dz=dz-dStartPosition;                       //internal z for profile
    dz=dz/dZSpacing;                            //floating point bin    
    dPosT=modf(dz,&dtemp);
    t[2]=dPosT;
    iPosBinLow=int(dtemp);
    iPosBinHigh=iPosBinLow+1;
    
    dBetaT=modf((dBeta-dBetaStart)/dBetaInc,&dtemp);
    t[0]=dPosT;
    iBetaBinLow=int(dtemp);
    iBetaBinHigh=iBetaBinLow+1;
    
    dGammaT=modf((dGamma-dGammaStart)/dGammaInc,&dtemp);
    t[1]=dPosT;
    iGammaBinLow=int(dtemp);
    iGammaBinHigh=iGammaBinLow+1;
    
    if ((iPosBinLow>=0) && (iPosBinHigh<iNumberOfPoints) && (iBetaBinLow>=0) &&
        (iBetaBinHigh<iNumberOfBeta) && (iGammaBinLow>=0) && (iGammaBinHigh<iNumberOfGamma)) {
        
        for (i=0; i<4; i++){
            ii=iBetaBinLow+i-1;
            if (ii<0) {ii=0;}                                                   // beta does not wrap
            if (ii>=iNumberOfBeta) {ii=iNumberOfBeta-1;}
            for (j=0; j<4; j++){
                jj=iGammaBinLow+j-1;
                if (jj<0) {jj=iNumberOfGamma-jj;}                               // gamma does wrap
                if (jj>=iNumberOfGamma) {jj=jj-iNumberOfGamma;}                
                for (k=0; k<4; k++) {
                    kk=iPosBinLow+k-1;
                    if (kk<0) {jj=0;}                                           // position does not wrap
                    if (kk>=iNumberOfPoints) {kk=iNumberOfPoints-1;}        
                    
                    p[i][j][k]=area[fn3Cto1C(ii,jj,kk)];                         // get area for this point                    
                    
                }
            }
        }
        
        returnvalue=fnTriCubicCatmullInterpolate(p,t);
        
    }
    else {
        printf("Euler Angles or Position Out of Range in DescreteEuler::fnGetArea");
        returnvalue=0;
    }

        return returnvalue*nf;
};

//get nSLD from molecular subgroups
double DiscreteEuler::fnGetnSLD(double dz) {
    
    int i,j,k, ii, jj, kk;
    int iPosBinLow, iPosBinHigh;
    int iBetaBinLow, iBetaBinHigh;
    int iGammaBinLow, iGammaBinHigh;
    double dPosT, dBetaT, dGammaT;
    double returnvalue, dtemp;
    double dtemp1, dtemp2, dtemp3,dtemp4;
    double t[3];                                           // tvalues for beta, gamma, position
    double parea[4][4][4],pprot[4][4][4],pdeut[4][4][4];   // all function values for tricubic spline interpolations
    // [beta][gamma][position]
    // index values: 0->-1, 1->0, 2->1, 3->2
    
    dz=dz-dStartPosition;                       //internal z for profile
    dz=dz/dZSpacing;                            //floating point bin    
    dPosT=modf(dz,&dtemp);
    t[2]=dPosT;
    iPosBinLow=int(dtemp);
    iPosBinHigh=iPosBinLow+1;
    
    dBetaT=modf((dBeta-dBetaStart)/dBetaInc,&dtemp);
    t[0]=dPosT;
    iBetaBinLow=int(dtemp);
    iBetaBinHigh=iBetaBinLow+1;
    
    dGammaT=modf((dGamma-dGammaStart)/dGammaInc,&dtemp);
    t[1]=dPosT;
    iGammaBinLow=int(dtemp);
    iGammaBinHigh=iGammaBinLow+1;
    
    if ((iPosBinLow>=0) && (iPosBinHigh<iNumberOfPoints) && (iBetaBinLow>=0) &&
        (iBetaBinHigh<iNumberOfBeta) && (iGammaBinLow>=0) && (iGammaBinHigh<iNumberOfGamma)) {
        
        for (i=0; i<4; i++){
            ii=iBetaBinLow+i-1;
            if (ii<0) {ii=0;}                                                   // beta does not wrap
            if (ii>=iNumberOfBeta) {ii=iNumberOfBeta-1;}
            for (j=0; j<4; j++){
                jj=iGammaBinLow+j-1;
                if (jj<0) {jj=iNumberOfGamma-jj;}                               // gamma does wrap
                if (jj>=iNumberOfGamma) {jj=jj-iNumberOfGamma;}                
                for (k=0; k<4; k++) {
                    kk=iPosBinLow+k-1;
                    if (kk<0) {jj=0;}                                           // position does not wrap
                    if (kk>=iNumberOfPoints) {kk=iNumberOfPoints-1;}        
                    
                    parea[i][j][k]=area[fn3Cto1C(ii,jj,kk)];
                    pprot[i][j][k]=nSLProt[fn3Cto1C(ii,jj,kk)];
                    pdeut[i][j][k]=nSLDeut[fn3Cto1C(ii,jj,kk)];

                    
                }
            }
        }
        
        dtemp1=fnTriCubicCatmullInterpolate(pprot,t);
        dtemp2=fnTriCubicCatmullInterpolate(pdeut,t);        
        dtemp3=dProtExchange*(dnSLDBulkSolvent+0.566e-6)/(6.34e-6+0.566e-6);
        dtemp4=(fnTriCubicCatmullInterpolate(parea,t))*dZSpacing;
        returnvalue=(((1-dtemp3)*dtemp1+dtemp3*dtemp2)/dtemp4);
        
    }
    else {
        printf("Euler Angles or Position Out of Range in DescreteEuler::fnGetArea");
        returnvalue=0;
    }
    
    return returnvalue;

};

//Use limits of molecular subgroups
double DiscreteEuler::fnGetLowerLimit() {return (dStartPosition);}
double DiscreteEuler::fnGetUpperLimit() {return (dStartPosition+double(iNumberOfPoints)*dZSpacing);}


void DiscreteEuler::fnSetNormarea(double dnormarea)
{
    normarea=dnormarea;
};


void DiscreteEuler::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Discrete %s StartPosition %e \n",cName, dStartPosition);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);    
}


//---------------------------------------------------------------------------------------------------------
//Freeform group 4 boxes
//---------------------------------------------------------------------------------------------------------
FreeBox::FreeBox(int n, double dstartposition, double dnSLD, double dnormarea)
{
    numberofboxes = n;
	startposition = dstartposition;
	nSLD=dnSLD;
	normarea=dnormarea;
    
    if (n>0) {box1 = new Box2Err(); box1->nf=1;};
    if (n>1) {box2 = new Box2Err(); box2->nf=1;};
    if (n>2) {box3 = new Box2Err(); box3->nf=1;};
	if (n>3) {box4 = new Box2Err(); box4->nf=1;};
	if (n>4) {box5 = new Box2Err(); box5->nf=1;};
	if (n>5) {box6 = new Box2Err(); box6->nf=1;};
	if (n>6) {box7 = new Box2Err(); box7->nf=1;};
	if (n>7) {box8 = new Box2Err(); box8->nf=1;};
	if (n>8) {box9 = new Box2Err(); box9->nf=1;};
	if (n>9) {box10 = new Box2Err(); box10->nf=1;};
    
    fnAdjustParameters();
    
};

FreeBox::~FreeBox(){
    if (numberofboxes>0) {delete box1;};
    if (numberofboxes>1) {delete box2;};
    if (numberofboxes>2) {delete box3;};
    if (numberofboxes>3) {delete box4;};
    if (numberofboxes>4) {delete box5;};
    if (numberofboxes>5) {delete box6;};
    if (numberofboxes>6) {delete box7;};
    if (numberofboxes>7) {delete box8;};
    if (numberofboxes>8) {delete box9;};
    if (numberofboxes>9) {delete box10;};
};

void FreeBox::fnAdjustParameters(){
	if (numberofboxes>0) {
		box1->z=startposition+0.5*box1->l;
		box1->vol=box1->l*normarea*vf1;
		box1->nSL=nSLD*box1->vol;
		if (numberofboxes>1) {
			box2->z=box1->z+0.5*box1->l+0.5*box2->l;
			box2->vol=box2->l*normarea*vf2;
			box2->nSL=nSLD*box2->vol;
			if (numberofboxes>2) {
				box3->z=box2->z+0.5*box2->l+0.5*box3->l;
				box3->vol=box3->l*normarea*vf3;
				box3->nSL=nSLD*box3->vol;
				if (numberofboxes>3) {
					box4->z=box3->z+0.5*box3->l+0.5*box4->l;
					box4->vol=box4->l*normarea*vf4;
					box4->nSL=nSLD*box4->vol;
					if (numberofboxes>4) {
						box5->z=box4->z+0.5*box4->l+0.5*box5->l;
						box5->vol=box5->l*normarea*vf5;
						box5->nSL=nSLD*box5->vol;
						if (numberofboxes>5) {
							box6->z=box5->z+0.5*box5->l+0.5*box6->l;	
							box6->vol=box6->l*normarea*vf6;
							box6->nSL=nSLD*box6->vol;
							if (numberofboxes>6) {
								box7->z=box6->z+0.5*box6->l+0.5*box7->l;
								box7->vol=box7->l*normarea*vf7;
								box7->nSL=nSLD*box7->vol;
								if (numberofboxes>7) {
									box8->z=box7->z+0.5*box7->l+0.5*box8->l;
									box8->vol=box8->l*normarea*vf8;
									box8->nSL=nSLD*box8->vol;
									if (numberofboxes>8) {
										box9->z=box8->z+0.5*box8->l+0.5*box9->l;
										box9->vol=box9->l*normarea*vf9;
										box9->nSL=nSLD*box9->vol;
										if (numberofboxes>9) {
                                            box10->z=box9->z+0.5*box9->l+0.5*box10->l;
                                            box10->vol=box10->l*normarea*vf10;
                                            box10->nSL=nSLD*box10->vol;
										};
									};
								};
							};
						};
					};
				};
			};
		};
	};
};

//Return value is area at position z
double FreeBox::fnGetArea(double dz) {
    double sum;
    sum=0;
    if (numberofboxes>0) {sum=box1->fnGetArea(dz);};
    if (numberofboxes>1) {sum=sum+box2->fnGetArea(dz);};
    if (numberofboxes>2) {sum=sum+box3->fnGetArea(dz);};
    if (numberofboxes>3) {sum=sum+box4->fnGetArea(dz);};
    if (numberofboxes>4) {sum=sum+box5->fnGetArea(dz);};
    if (numberofboxes>5) {sum=sum+box6->fnGetArea(dz);};
    if (numberofboxes>6) {sum=sum+box7->fnGetArea(dz);};
    if (numberofboxes>7) {sum=sum+box8->fnGetArea(dz);};
    if (numberofboxes>8) {sum=sum+box9->fnGetArea(dz);};
    if (numberofboxes>9) {sum=sum+box10->fnGetArea(dz);};	 
    
    return sum;
};

//get nSLD from molecular subgroups
double FreeBox::fnGetnSLD(double dz) {
    //printf("nSLD %e \n", nSLD);
	return nSLD;
};

//Use limits of molecular subgroups
double FreeBox::fnGetLowerLimit() {return box1->fnGetLowerLimit();};
double FreeBox::fnGetUpperLimit() {
    if (numberofboxes>9) {return box10->fnGetUpperLimit();}
    else if (numberofboxes>8) {return box9->fnGetUpperLimit();}
    else if (numberofboxes>7) {return box8->fnGetUpperLimit();}
    else if (numberofboxes>6) {return box7->fnGetUpperLimit();}
    else if (numberofboxes>5) {return box6->fnGetUpperLimit();}
    else if (numberofboxes>4) {return box5->fnGetUpperLimit();}
    else if (numberofboxes>3) {return box4->fnGetUpperLimit();}
    else if (numberofboxes>2) {return box3->fnGetUpperLimit();}
    else if (numberofboxes>1) {return box2->fnGetUpperLimit();}
    else {return box1->fnGetUpperLimit();}
    
    
    
};

void FreeBox::fnSetSigma(double sigma)
{
    if (numberofboxes>0) {box1->sigma1=sigma; box1->sigma2=sigma;};
    if (numberofboxes>1) {box2->sigma1=sigma; box2->sigma2=sigma;};
    if (numberofboxes>2) {box3->sigma1=sigma; box3->sigma2=sigma;};
    if (numberofboxes>3) {box4->sigma1=sigma; box4->sigma2=sigma;};
    if (numberofboxes>4) {box5->sigma1=sigma; box5->sigma2=sigma;};
    if (numberofboxes>5) {box6->sigma1=sigma; box6->sigma2=sigma;};
    if (numberofboxes>6) {box7->sigma1=sigma; box7->sigma2=sigma;};
    if (numberofboxes>7) {box8->sigma1=sigma; box8->sigma2=sigma;};
    if (numberofboxes>8) {box9->sigma1=sigma; box9->sigma2=sigma;};
    if (numberofboxes>9) {box10->sigma1=sigma; box10->sigma2=sigma;};
};

void FreeBox::fnSetStartposition(double dz)
{
    startposition=dz;
	fnAdjustParameters();
};

void FreeBox::fnSetNormarea(double dnormarea)
{
    normarea=dnormarea;
	fnAdjustParameters();
};

void FreeBox::fnSetnSLD(double dnSLD)
{
    nSLD=dnSLD;
	fnAdjustParameters();
};


void FreeBox::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    //char *str = new char[80];
    
    fprintf(fp, "FreeBox %s numberofboxes %i \n",cName, numberofboxes);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
    //cg->fnWriteGroup2File(fp, "cg", dimension, stepsize);
    //phosphate->fnWriteGroup2File(fp, "phosphate", dimension, stepsize);
    //choline->fnWriteGroup2File(fp, "choline", dimension, stepsize);
    
    //delete []str;
    
}
//---------------------------------------------------------------------------------------------------------
//Hermite spline interpolation
//---------------------------------------------------------------------------------------------------------

Hermite::Hermite(int n, double dstartposition, double dnSLD, double dnormarea)
{
    
    numberofcontrolpoints = n;
	nSLD=dnSLD;
	normarea=dnormarea;
    dp = new double[n];
    vf = new double[n];
    
    
};

Hermite::~Hermite(){
    delete [] dp;
    delete [] vf;
};

//Return value is area at position z
double Hermite::fnGetArea(double dz) {
    double p0, p1, m0, m1, h00, h01, h10, h11, t, dd, t_2, t_3, volfrac;
    int interval;  
    
    //printf("fnGetArea z %e dp[0] %e dp[5] %e vf[0] %e vf[5] %e", dz, dp[0], dp[5], vf[0], vf[5]);
    
    
    interval=-1;
    for (int i=0; i<(numberofcontrolpoints-1); i++) {
        //printf("i %i dz %e dp[i] %e dp[i+1] %e \n", i, dz, dp[i], dp[i+1]);
        if ((dp[i]<=dz) && (dp[i+1]>dz)) {
            //printf("Found interval %i \n", i);
            interval=i;
        }
    }
    
    if (interval>=0) {
        
        if (interval==0) {
            m0=0;
            m1=(vf[2]-vf[0])/(dp[2]-dp[0]);
        }
        else if (interval==(numberofcontrolpoints-2)) {
            m0=(vf[interval+1]-vf[interval-1])/(dp[interval+1]-dp[interval-1]);
            m1=0;
        }
        else {
            m0=(vf[interval+1]-vf[interval-1])/(dp[interval+1]-dp[interval-1]);
            m1=(vf[interval+2]-vf[interval])/(dp[interval+2]-dp[interval]);
        }
        p0=vf[interval];
        p1=vf[interval+1];
        
        
        dd=dp[interval+1]-dp[interval];
        t=(dz-dp[interval])/dd;
        t_2=t*t;
        t_3=t_2*t;
        h00=   2*t_3-3*t_2+1;
        h10=     t_3-2*t_2+t;
        h01=(-2)*t_3+3*t_2;
        h11=     t_3-t_2;
        
        volfrac=h00*p0+h10*dd*m0+h01*p1+h11*dd*m1;
        //printf("m0 %e m1 %e p0 %e p1 %e dd %e dz %e t %e normarea %e", m0, m1, p0, p1, dd, dz, t, normarea);
        
        //printf("vf %e normarea %e \n", volfrac, normarea);
        
        return volfrac*normarea;
        
    }
    else {
        //printf("zero return");
        return 0; 
    }
    
    
};

//Return value is integral at position z
double Hermite::fnGetIntegral(double dz) {
    double p0, p1, m0, m1, h00, h01, h10, h11, t, dd, t_2, t_3, t_4,integral;
    int interval;  
    
    //printf("fnGetArea z %e dp[0] %e dp[5] %e vf[0] %e vf[5] %e", dz, dp[0], dp[5], vf[0], vf[5]);
    
    
    interval=-1;
    for (int i=0; i<(numberofcontrolpoints-1); i++) {
        //printf("i %i dz %e dp[i] %e dp[i+1] %e \n", i, dz, dp[i], dp[i+1]);
        if ((dp[i]<=dz) && (dp[i+1]>dz)) {
            //printf("Found interval %i \n", i);
            interval=i;
        }
    }
    
    if (interval>=0) {
        
        if (interval==0) {
            m0=0;
            m1=(vf[2]-vf[0])/(dp[2]-dp[0]);
        }
        else if (interval==(numberofcontrolpoints-2)) {
            m0=(vf[interval+1]-vf[interval-1])/(dp[interval+1]-dp[interval-1]);
            m1=0;
        }
        else {
            m0=(vf[interval+1]-vf[interval-1])/(dp[interval+1]-dp[interval-1]);
            m1=(vf[interval+2]-vf[interval])/(dp[interval+2]-dp[interval]);
        }
        p0=vf[interval];
        p1=vf[interval+1];
        
        
        dd=dp[interval+1]-dp[interval];
        t=(dz-dp[interval])/dd;
        t_2=t*t;
        t_3=t_2*t;
        t_4=t_3*t;
        h00= (1/2)*t_4-      t_3          +t;
        h01=(-1/2)*t_4+      t_3            ;
        h10= (1/4)*t_4-(2/3)*t_3+(1/2)*t_2  ;
        h11= (1/4)*t_4-(1/3)*t_3            ;
        
        integral=dd*(h00*p0+h10*dd*m0+h01*p1+h11*dd*m1);
        //printf("m0 %e m1 %e p0 %e p1 %e dd %e dz %e t %e normarea %e", m0, m1, p0, p1, dd, dz, t, normarea);
        
        //printf("vf %e normarea %e \n", volfrac, normarea);
        
        return integral*normarea;
        
    }
    else {
        //printf("zero return");
        return 0; 
    }
    
    
};

double Hermite::fnGetVolume(double dz1, double dz2) {
    
    double temp, integral;
    int i1, i2;
    
    if (dz1>dz2){
        temp=dz2;
        dz2=dz1;
        dz1=temp;
    }
    
    //check for boundaries
    if (dz1<dp[0]) {dz1=dp[0];}
    if (dz2>dp[numberofcontrolpoints-1]) {dz2=dp[numberofcontrolpoints];}
    
    //check in which intervals of the spline the two values lie
    i1=-1; i2=-1;
    for (int i=0; i<(numberofcontrolpoints-1); i++) {
        if ((dp[i]<=dz1) && (dp[i+1]>dz1)) {
            i1=i;
        }
        if ((dp[i]<=dz2) && (dp[i+1]>dz2)) {
            i2=i;
        }
    }
       
    //sum first over all complete inteval between support points
    integral=0;
    for (int i=i1+1; i<i2; i++) {
        integral+=fnGetIntegral(dp[i]-(0.00001*(dp[i]-dp[i-1])))-fnGetIntegral(dp[i-1]);
    }
    
    //add now the incomplete intervals at the beginning and the end of the integration interval
    if (i1!=0) {
        integral+=fnGetIntegral(dp[i1]-(0.00001*(dp[i1]-dp[i1-1])))-fnGetIntegral(dz1);        
    }
    integral+=fnGetIntegral(dz2)-fnGetIntegral(dp[i2]);
    
    return integral;
};

//get nSLD from molecular subgroups
double Hermite::fnGetnSLD(double dz) {
    //printf("nSLD %e \n", nSLD);
	return nSLD;
};

//Use limits of molecular subgroups
double Hermite::fnGetLowerLimit() {return dp[0];};
double Hermite::fnGetUpperLimit() {return dp[numberofcontrolpoints-1];}


void Hermite::fnSetNormarea(double dnormarea)
{
    normarea=dnormarea;
};

void Hermite::fnSetnSLD(double dnSLD)
{
    nSLD=dnSLD;
};


void Hermite::fnWriteGroup2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
    fprintf(fp, "Hermite %s numberofcontrolpoints %i normarea %e nf %e\n",cName, numberofcontrolpoints,normarea, nf);
    nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);    
}


//-----------------------------------------------------------------------------------------------------------
// Bilayer Library
//-----------------------------------------------------------------------------------------------------------
Monolayer_DOPS::Monolayer_DOPS()
{

    headgroup = new PS();
    
    headgroup->vol=260;                //PS volume and length are estimates
	headgroup->nSL=8.4513e-4;
	headgroup->l=7.5;
    
    volacyllipid=925;
    nslacyllipid=-2.6688e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    
    fnAdjustParameters();

}
Monolayer_DOPS::~Monolayer_DOPS()
{
    delete headgroup;
}
Monolayer_DOPS_xray::Monolayer_DOPS_xray()
{
    
    fnSetnSL(5.07E-4,6.81E-3,1.885e-3,1.33E-3,1.323E-3);
    fnAdjustParameters();
    
}
Monolayer_DOPS_xray::~Monolayer_DOPS_xray()
{
}

Monolayer_DPPS::Monolayer_DPPS()
{
    
    headgroup = new PS(); 
    
    headgroup->vol=260;                //PS volume and length are estimates
	headgroup->nSL=8.4513e-4;
	headgroup->l=7.5;
    
    volacyllipid=789;
    nslacyllipid=-3.2502E-04;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    
    fnAdjustParameters();
    
}
Monolayer_DPPS::~Monolayer_DPPS()
{
    delete headgroup;
}
Monolayer_DPPS_xray::Monolayer_DPPS_xray()
{
    fnSetnSL(5.07E-4,6.81E-3,1.885e-3,1.33E-3,1.323E-3);
    fnAdjustParameters();
    
}
Monolayer_DPPS_xray::~Monolayer_DPPS_xray()
{
}

ssBLM_POPC::ssBLM_POPC()
{
	
    volacyllipid=925;
    nslacyllipid=-2.6688e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    
    fnAdjustParameters();
}

tBLM_HC18_DOPC::tBLM_HC18_DOPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
    fnAdjustParameters();
    
}

tBLM_HC18_d54DMPC::tBLM_HC18_d54DMPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=735;
    nslacyllipid=5.3324E-03;
    volmethyllipid=98.8;
    nslmethyllipid=5.334e-4;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
    fnAdjustParameters();
    
}

tBLM_WC14_DOPC::tBLM_WC14_DOPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=850;
    nslacyltether=-3.5834e-4;
    
    fnAdjustParameters();
    
}

tBLM_HC18_POPC_POPA::tBLM_HC18_POPC_POPA()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=925;
    nslacyllipid=-2.67e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=174;                //was 174
	headgroup2_2->vol=174;                //was 174
	headgroup1_2->nSL=6.2364e-4;          //was 6.2364e-4
	headgroup2_2->nSL=6.2364e-4;          //was 6.2364e-4
	headgroup1_2->l=5;
	headgroup2_2->l=5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_HC18_POPC_POPG::tBLM_HC18_POPC_POPG()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=925;
    nslacyllipid=-2.67e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=270;                //PG volume and length are estimates
	headgroup2_2->vol=270;
	headgroup1_2->nSL=7.1472e-4;
	headgroup2_2->nSL=7.1472e-4;
	headgroup1_2->l=7.8;
	headgroup2_2->l=7.8;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}
tBLM_HC18_POPC_POPS::tBLM_HC18_POPC_POPS()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=925;
    nslacyllipid=-2.67e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=925;
    nslacyllipid_2=-2.67e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_DOPS::tBLM_HC18_DOPC_DOPS()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}
tBLM_HC18_DOPC_d54DMPC::tBLM_HC18_DOPC_d54DMPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=330;       
	headgroup2_2->vol=330;       
	headgroup1_2->nSL=6.0012e-4; 
	headgroup2_2->nSL=6.0012e-4; 
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	
    volacyllipid_2=735;
    nslacyllipid_2=5.3324E-03;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=5.334e-4;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DMPC_d54DMPC::tBLM_HC18_DMPC_d54DMPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=735;
    nslacyllipid=-2.9166E-04;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=330;       
	headgroup2_2->vol=330;       
	headgroup1_2->nSL=6.0012e-4; 
	headgroup2_2->nSL=6.0012e-4; 
	headgroup1_2->l=9.5;
	headgroup2_2->l=9.5;
	
    volacyllipid_2=735;
    nslacyllipid_2=5.3324E-03;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=5.334e-4;
    
    fnAdjustParameters();
    
}

tBLM_WC14_DOPC_DOPS::tBLM_WC14_DOPC_DOPS()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=850;
    nslacyltether=-3.5834e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_WC14_DOPC_PIP::tBLM_WC14_DOPC_PIP()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=850;
    nslacyltether=-3.5834e-4;
    
	headgroup1_2->vol=500;                //PIP volume and length are estimates
	headgroup2_2->vol=500;
	headgroup1_2->nSL=1.22e-3;
	headgroup2_2->nSL=1.22e-3;
	headgroup1_2->l=12.0;
	headgroup2_2->l=12.0;
	
    volacyllipid_2=1025;
    nslacyllipid_2=-7.5785e-5;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_PIP::tBLM_HC18_DOPC_PIP()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=500;                //PS volume and length are estimates
	headgroup2_2->vol=500;
	headgroup1_2->nSL=1.22e-3;
	headgroup2_2->nSL=1.22e-3;
	headgroup1_2->l=12.0;
	headgroup2_2->l=12.0;
	
    volacyllipid_2=1025;
    nslacyllipid_2=-7.5785e-5;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_DOPS_PIP::tBLM_HC18_DOPC_DOPS_PIP()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
	headgroup1_3->vol=500;                //PIP volume and length are estimates
	headgroup2_3->vol=500;
	headgroup1_3->nSL=1.22e-3;
	headgroup2_3->nSL=1.22e-3;
	headgroup1_3->l=12.0;
	headgroup2_3->l=12.0;
	
    volacyllipid_3=1025;
    nslacyllipid_3=-7.5785e-5;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_WC14_DOPC_DOPS_CHOL::tBLM_WC14_DOPC_DOPS_CHOL()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=850;
    nslacyltether=-3.5834e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
    volchol=630;
    nslchol=1.3215e-4;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_DOPS_CHOL::tBLM_HC18_DOPC_DOPS_CHOL()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
    volchol=630;
    nslchol=1.3215e-4;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_DOPS_PIP_CHOL::tBLM_HC18_DOPC_DOPS_PIP_CHOL()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
    volchol=630;
    nslchol=1.3215e-4;
	
	headgroup1_3->vol=500;                //PIP volume and length are estimates
	headgroup2_3->vol=500;
	headgroup1_3->nSL=1.22e-3;
	headgroup2_3->nSL=1.22e-3;
	headgroup1_3->l=12.0;
	headgroup2_3->l=12.0;
	
    volacyllipid_3=1025;
    nslacyllipid_3=-7.5785e-5;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    fnAdjustParameters();
    
}

tBLM_HC18_DOPC_DOPS_PIP_CHOL_domain::tBLM_HC18_DOPC_DOPS_PIP_CHOL_domain()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	tether_domain->vol=380;
	tether_domain->nSL=2.1864e-4;
	tetherg_domain->vol=110;
	tetherg_domain->nSL=1.8654e-4;
	
    volacyllipid=972;
    nslacyllipid=-2.09e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=999;
    nslacyltether=-2.25e-4;
    
	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	headgroup1_2_domain->vol=260;                //PS volume and length are estimates
	headgroup2_2_domain->vol=260;
	headgroup1_2_domain->nSL=8.4513e-4;
	headgroup2_2_domain->nSL=8.4513e-4;
	headgroup1_2_domain->l=7.5;
	headgroup2_2_domain->l=7.5;
	
    volacyllipid_2=972;
    nslacyllipid_2=-2.09e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
	
    volchol=630;
    nslchol=1.3215e-4;
	
	headgroup1_3->vol=500;                //PIP volume and length are estimates
	headgroup2_3->vol=500;
	headgroup1_3->nSL=1.22e-3;
	headgroup2_3->nSL=1.22e-3;
	headgroup1_3->l=12.0;
	headgroup2_3->l=12.0;
	headgroup1_3_domain->vol=500;                //PIP volume and length are estimates
	headgroup2_3_domain->vol=500;
	headgroup1_3_domain->nSL=1.22e-3;
	headgroup2_3_domain->nSL=1.22e-3;
	headgroup1_3_domain->l=12.0;
	headgroup2_3_domain->l=12.0;
	
    volacyllipid_3=1025;
    nslacyllipid_3=-7.5785e-5;
    volmethyllipid_3=98.8;
    nslmethyllipid_3=-9.15e-5;
    
    fnAdjustParameters();
    
}

//------------------------------------------------------------------------------------------------------
void fnWriteConstant(FILE *fp, const char *cName, double area, double nSLD, int dimension, double stepsize)
{
    int i;
    double d;
    
    fprintf(fp, "Constant %s area %lf \n",cName, area);
    fprintf(fp, "z%s a%s nsl%s \n",cName, cName, cName);
	for (i=0; i<dimension; i++)
	{
        d=double(i)*stepsize;
        fprintf(fp, "%lf %lf %e \n", d, area, nSLD*area*stepsize);
	};
	fprintf(fp,"\n");
    
}


//------------------------------------------------------------------------------------------------------
double fnClearCanvas(double aArea[], double anSL[], int dimension)
{
	int j;
	
	for (j=0; j<dimension; j++) {
		aArea[j]=0; anSL[j]=0;
	};
	
	return 0;
}
double fnClearCanvas(double aArea[], double anSL[], double aAbsorb[], int dimension)
{
	int j;
	
	for (j=0; j<dimension; j++) {
		aArea[j]=0; anSL[j]=0; aAbsorb[j]=0;
	};
	
	return 0;
}
//------------------------------------------------------------------------------------------------------
// Overlays one canvas onto another 

void fnOverlayCanvasOnCanvas(double aArea[], double anSL[], double aArea2[], double anSL2[], int dimension, double dMaxArea)
{
    double temparea;
	int i;
	
	
	for(i=0; i<dimension; i++)
	{
        temparea=aArea2[i]+aArea[i];
        if (temparea>dMaxArea) 
        {
            anSL[i]=anSL[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
            anSL[i]=anSL[i]+anSL2[i];
            aArea[i]=dMaxArea;
        }
        else 
        {
            //printf("Bin %i Areainc %f area now %f nSLinc %g nSL now %g \n", i, aArea2[i], aArea[i], anSL2[1], anSL[i]);
            aArea[i]=aArea[i]+aArea2[i];
            anSL[i]=anSL[i]+anSL2[i];
        }
	};
};
void fnOverlayCanvasOnCanvas(double aArea[], double anSL[], double aAbsorb[], double aArea2[], double anSL2[], double aAbsorb2[], int dimension, double dMaxArea)
{
    double temparea;
	int i;
	
	
	for(i=0; i<dimension; i++)
	{
        temparea=aArea2[i]+aArea[i];
        if (temparea>dMaxArea) 
        {
            anSL[i]=anSL[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
            anSL[i]=anSL[i]+anSL2[i];
            aAbsorb[i]=aAbsorb[i]*(1-((temparea-dMaxArea)/aArea[i]));						//eliminate the overfilled portion using original content
            aAbsorb[i]=aAbsorb[i]+aAbsorb2[i];
            aArea[i]=dMaxArea;
        }
        else 
        {
            //printf("Bin %i Areainc %f area now %f nSLinc %g nSL now %g \n", i, aArea2[i], aArea[i], anSL2[1], anSL[i]);
            aArea[i]=aArea[i]+aArea2[i];
            anSL[i]=anSL[i]+anSL2[i];
            aAbsorb[i]=aAbsorb[i]+aAbsorb2[i];
        }
	};
};

//------------------------------------------------------------------------------------------------------
//writes out canvas to reflectivity model taking into account bulk nSLD

void fnWriteCanvas2Model(double aArea[], double anSL[], fitinfo fit[], int gaussstart, int dimension, double stepsize, double dMaxArea, double normarea, int modelstart, int modelend)
{
    int i, j;
    if (dMaxArea!=0)  {
        for (i=modelstart; i<modelend+1; i++)
            for (j=0; j<dimension; j++)  {
                fit[i].m.rho[j+gaussstart]=(anSL[j]/(normarea*stepsize))+(1-(aArea[j]/normarea))*fit[i].m.rho[fit[i].m.n-1];
            }
    }
    else  {
        for (i=modelstart; i<modelend+1; i++)  
            for (j=0; j<dimension; j++)  {fit[i].m.rho[j+gaussstart]=fit[i].m.rho[fit[i].m.n-1];}
    }
}
void fnWriteCanvas2Model(double aArea[], double anSL[], double aAbsorb[], fitinfo fit[], int gaussstart, int dimension, double stepsize, double dMaxArea, double normarea, int modelstart, int modelend)
{
    int i, j;
    if (dMaxArea!=0)  {
        for (i=modelstart; i<modelend+1; i++)
            for (j=0; j<dimension; j++)  {
                //printf("bin %i area %e normarea %e areafraction %e bulk mu %e absorption %e result %e \n",j, aArea[j], normarea, aArea[j]/normarea,fit[i].m.mu[fit[i].m.n-1], aAbsorb[j],(aAbsorb[j]/(normarea*stepsize))+(1-(aArea[j]/normarea))*fit[i].m.mu[fit[i].m.n-1]);
                fit[i].m.rho[j+gaussstart]=(anSL[j]/(normarea*stepsize))+(1-(aArea[j]/normarea))*fit[i].m.rho[fit[i].m.n-1];
                fit[i].m.mu[j+gaussstart]=(aAbsorb[j]/(normarea*stepsize))+(1-(aArea[j]/normarea))*fit[i].m.mu[fit[i].m.n-1];
            }
    }
    else  {
        for (i=modelstart; i<modelend+1; i++)  
            for (j=0; j<dimension; j++)  {
                fit[i].m.rho[j+gaussstart]=fit[i].m.rho[fit[i].m.n-1];
                fit[i].m.mu[j+gaussstart]=fit[i].m.mu[fit[i].m.n-1];
            }
    }
}

