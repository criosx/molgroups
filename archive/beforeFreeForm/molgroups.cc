/*
 *  molgroups.cc
 *  Gauss
 *
 *  Created by Frank Heinrich on 27/10/08.
 *  updated 7-July-2011 after minor updates
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "molgroups.h"
#include "stdio.h"
#include "stdlib.h"


//------------------------------------------------------------------------------------------------------
//Parent Object Implementation

nSLDObj::nSLDObj()
{
    bWrapping=true;
};

nSLDObj::~nSLDObj(){};

double nSLDObj::fnGetLowerLimit(){return 0;};
double nSLDObj::fnGetUpperLimit(){return 0;};
double nSLDObj::fnGetArea(double z){return 0;};
double nSLDObj::fnGetnSLD(double z){return 0;};
void nSLDObj::fnWritePar2File (FILE *fp, const char *cName, int dimension, double stepsize){};
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
			  //printf("Bin %i Areainc %f area now %f nSLD %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize, anSL[i]);
			}
                        else {
			  aArea[i]=aArea[i]+dAreaInc*dprefactor;
			  anSL[i]=anSL[i]+fnGetnSLD(d)*dAreaInc*stepsize*dprefactor;
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

void   BoxErr::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
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

void Box2Err::fnSetAllSigma(double sigma)
{
	sigma1=sigma;
	sigma2=sigma;
}
                
void Box2Err::fnSetZ(double dz)
{
	z=dz;
};

void   Box2Err::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
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

void   Gaussian::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
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

void   Parabolic::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
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

void   StretchGaussian::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
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

void PC::fnSetAllSigma(double sigma)
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

void PC::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        //cg->fnWritePar2File(fp, "cg", dimension, stepsize);
        //phosphate->fnWritePar2File(fp, "phosphate", dimension, stepsize);
        //choline->fnWritePar2File(fp, "choline", dimension, stepsize);
        
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

void   PCm::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        fprintf(fp, "PCm %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        //cg->fnWritePar2File(fp, "cg_m", dimension, stepsize);
        //phosphate->fnWritePar2File(fp, "phosphate_m", dimension, stepsize);
        //choline->fnWritePar2File(fp, "choline_m", dimension, stepsize);
        
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
  headgroup1->fnSetAllSigma(sigma);
  lipid1->sigma1=sigma;
  lipid1->sigma2=sigma;
  methyl1->sigma1=sigma;
  methyl1->sigma2=sigma;
  methyl2->sigma1=sigma;
  methyl2->sigma2=sigma;
  lipid2->sigma1=sigma;
  lipid2->sigma2=sigma;
  headgroup2->fnSetAllSigma(sigma);
  
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


void ssBLM::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
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
  headgroup1->fnSetAllSigma(sigma);
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
  headgroup2->fnSetAllSigma(sigma);
  
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


void tBLM::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        bME->fnWritePar2File(fp, "bME", dimension, stepsize);
        tether->fnWritePar2File(fp, "tether", dimension, stepsize);
        tetherg->fnWritePar2File(fp, "tetherg", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
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

  substrate->sigma2=global_rough;
  bME->sigma1=global_rough;
  bME->sigma2=global_rough;
  headgroup1->fnSetAllSigma(sigma);
  headgroup1_2->fnSetAllSigma(sigma);
  tether->sigma1=global_rough; 
  tether->sigma2=sigma;
  tetherg->fnSetAllSigma(sigma); 
  lipid1->fnSetAllSigma(sigma);
  methyl1->fnSetAllSigma(sigma);
  methyl2->fnSetAllSigma(sigma);
  lipid2->fnSetAllSigma(sigma);
  headgroup2->fnSetAllSigma(sigma);
  headgroup2_2->fnSetAllSigma(sigma);
  
  //outer hydrocarbons
  l_ohc=l_lipid2;
  nf_ohc_lipid  =1-nf_lipid_2;
  nf_ohc_lipid_2=nf_lipid_2;
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

double tBLM_binary::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_binary::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        bME->fnWritePar2File(fp, "bME", dimension, stepsize);
        tether->fnWritePar2File(fp, "tether", dimension, stepsize);
        tetherg->fnWritePar2File(fp, "tetherg", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        headgroup1_2->fnWritePar2File(fp, "headgroup1_2", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
        headgroup2_2->fnWritePar2File(fp, "headgroup2_2", dimension, stepsize);
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

  // set all sigma
  
  //sigma=sqrt(2.4*2.4 + global_rough*global_rough);

  substrate->sigma2=global_rough;
  bME->sigma1=global_rough;
  bME->sigma2=global_rough;
  headgroup1->fnSetAllSigma(sigma);
  headgroup1_2->fnSetAllSigma(sigma);
  headgroup1_3->fnSetAllSigma(sigma);
  tether->sigma1=global_rough; 
  tether->sigma2=sigma;
  tetherg->fnSetAllSigma(sigma); 
  lipid1->fnSetAllSigma(sigma);
  methyl1->fnSetAllSigma(sigma);
  methyl2->fnSetAllSigma(sigma);
  lipid2->fnSetAllSigma(sigma);
  headgroup2->fnSetAllSigma(sigma);
  headgroup2_2->fnSetAllSigma(sigma);
  headgroup2_3->fnSetAllSigma(sigma);
  
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

double tBLM_ternary::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_ternary::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        bME->fnWritePar2File(fp, "bME", dimension, stepsize);
        tether->fnWritePar2File(fp, "tether", dimension, stepsize);
        tetherg->fnWritePar2File(fp, "tetherg", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        headgroup1_2->fnWritePar2File(fp, "headgroup1_2", dimension, stepsize);
	    headgroup1_3->fnWritePar2File(fp, "headgroup1_3", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
        headgroup2_2->fnWritePar2File(fp, "headgroup2_2", dimension, stepsize);
	    headgroup2_3->fnWritePar2File(fp, "headgroup2_3", dimension, stepsize);
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

  // set all sigma
  
  //sigma=sqrt(2.4*2.4 + global_rough*global_rough);

  substrate->sigma2=global_rough;
  bME->sigma1=global_rough;
  bME->sigma2=global_rough;
  headgroup1->fnSetAllSigma(sigma);
  headgroup1_2->fnSetAllSigma(sigma);
  tether->sigma1=global_rough; 
  tether->sigma2=sigma;
  tetherg->fnSetAllSigma(sigma); 
  lipid1->fnSetAllSigma(sigma);
  methyl1->fnSetAllSigma(sigma);
  methyl2->fnSetAllSigma(sigma);
  lipid2->fnSetAllSigma(sigma);
  headgroup2->fnSetAllSigma(sigma);
  headgroup2_2->fnSetAllSigma(sigma);
  
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

double tBLM_ternary_chol::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_ternary_chol::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        bME->fnWritePar2File(fp, "bME", dimension, stepsize);
        tether->fnWritePar2File(fp, "tether", dimension, stepsize);
        tetherg->fnWritePar2File(fp, "tetherg", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        headgroup1_2->fnWritePar2File(fp, "headgroup1_2", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
        headgroup2_2->fnWritePar2File(fp, "headgroup2_2", dimension, stepsize);
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

  // set all sigma
  
  //sigma=sqrt(2.4*2.4 + global_rough*global_rough);

  substrate->sigma2=global_rough;
  bME->sigma1=global_rough;
  bME->sigma2=global_rough;
  headgroup1->fnSetAllSigma(sigma);
  headgroup1_2->fnSetAllSigma(sigma);
  headgroup1_3->fnSetAllSigma(sigma);
  tether->sigma1=global_rough; 
  tether->sigma2=sigma;
  tetherg->fnSetAllSigma(sigma); 
  lipid1->fnSetAllSigma(sigma);
  methyl1->fnSetAllSigma(sigma);
  methyl2->fnSetAllSigma(sigma);
  lipid2->fnSetAllSigma(sigma);
  headgroup2->fnSetAllSigma(sigma);
  headgroup2_2->fnSetAllSigma(sigma);
  headgroup2_3->fnSetAllSigma(sigma);
  
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

double tBLM_quaternary_chol::fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea)
{
    nSLDObj::fnWriteProfile(aArea,anSLD,dimension,stepsize,dMaxArea);
	return normarea;
};


void tBLM_quaternary_chol::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        //fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf);
        //nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        substrate->fnWritePar2File(fp, "substrate", dimension, stepsize);
        bME->fnWritePar2File(fp, "bME", dimension, stepsize);
        tether->fnWritePar2File(fp, "tether", dimension, stepsize);
        tetherg->fnWritePar2File(fp, "tetherg", dimension, stepsize);
        headgroup1->fnWritePar2File(fp, "headgroup1", dimension, stepsize);
        headgroup1_2->fnWritePar2File(fp, "headgroup1_2", dimension, stepsize);
		headgroup1_3->fnWritePar2File(fp, "headgroup1_3", dimension, stepsize);
        lipid1->fnWritePar2File(fp, "lipid1", dimension, stepsize);
        methyl1->fnWritePar2File(fp, "methyl1", dimension, stepsize);
        methyl2->fnWritePar2File(fp, "methyl2", dimension, stepsize);
        lipid2->fnWritePar2File(fp, "lipid2", dimension, stepsize);
        headgroup2->fnWritePar2File(fp, "headgroup2", dimension, stepsize);
        headgroup2_2->fnWritePar2File(fp, "headgroup2_2", dimension, stepsize);
		headgroup2_3->fnWritePar2File(fp, "headgroup2_3", dimension, stepsize);
		fnWriteConstant(fp, "normarea", normarea, 0, dimension, stepsize);

        
        //delete []str;

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

void FreeBox::fnSetAllSigma(double sigma)
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


void FreeBox::fnWritePar2File(FILE *fp, const char *cName, int dimension, double stepsize)
{
        //char *str = new char[80];
        
        fprintf(fp, "FreeBox %s numberofboxes %i \n",cName, numberofboxes);
        nSLDObj::fnWriteData2File(fp, cName, dimension, stepsize);
        //cg->fnWritePar2File(fp, "cg", dimension, stepsize);
        //phosphate->fnWritePar2File(fp, "phosphate", dimension, stepsize);
        //choline->fnWritePar2File(fp, "choline", dimension, stepsize);
        
        //delete []str;

}

//-----------------------------------------------------------------------------------------------------------
// Bilayer Library
//-----------------------------------------------------------------------------------------------------------
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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=982;
    nslacyltether=-2.85e-4;
    
    fnAdjustParameters();

}

tBLM_WC14_DOPC::tBLM_WC14_DOPC()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
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
    volacyltether=982;
    nslacyltether=-2.85e-4;

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
    volacyltether=982;
    nslacyltether=-2.85e-4;

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

tBLM_HC18_DOPC_DOPS::tBLM_HC18_DOPC_DOPS()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=982;
    nslacyltether=-2.85e-4;

	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
    volmethyllipid_2=98.8;
    nslmethyllipid_2=-9.15e-5;
    
    fnAdjustParameters();

}

tBLM_WC14_DOPC_DOPS::tBLM_WC14_DOPC_DOPS()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
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
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=850;
    nslacyltether=-3.5834e-4;

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

tBLM_HC18_DOPC_PIP::tBLM_HC18_DOPC_PIP()
{
	tether->vol=380;
	tether->nSL=2.1864e-4;
	tetherg->vol=110;
	tetherg->nSL=1.8654e-4;
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=982;
    nslacyltether=-2.85e-4;

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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=982;
    nslacyltether=-2.85e-4;

	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
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
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=980;
    nslacyltether=-2.85e-4;

	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
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
	
    volacyllipid=980;
    nslacyllipid=-2.00e-4;
    volmethyllipid=98.8;
    nslmethyllipid=-9.15e-5;
    volmethyltether=98.8;
    nslmethyltether=-9.15e-5;
    volacyltether=980;
    nslacyltether=-2.85e-4;

	headgroup1_2->vol=260;                //PS volume and length are estimates
	headgroup2_2->vol=260;
	headgroup1_2->nSL=8.4513e-4;
	headgroup2_2->nSL=8.4513e-4;
	headgroup1_2->l=7.5;
	headgroup2_2->l=7.5;
	
    volacyllipid_2=980;
    nslacyllipid_2=-2.00e-4;
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

