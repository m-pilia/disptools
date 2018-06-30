/*!
 * \file VolumeMatching3D.c
 *
 * This file comes from
 *      PREDICT - Atrophy Simulation
 *      https://www.med.upenn.edu/sbia/atrophysimulation.html
 *      https://www.nitrc.org/projects/atrophysim/
 *
 * Publications:
 *  + B.Karacali and C.Davatzikos, "Estimating Topology Preserving and Smooth
 *    Displacement Fields", IEEE Transactions on Medical Imaging, Vol 23,
 *    No.7, p.868-880, July 2004.
 *  + B.Karacali and C.Davatzikos, "Simulation of Tissue Atrophy Using a
 *    Topology Preserving Transformation Model", submitted to IEEE Trans. 
 *    on Medical Imaging, 2004.
 *  + Zhong Xue, Dinggang Shen, Bilge Karacali, Joshua Stern, David Rottenberg,
 *    and Christos Davatzikos, "Simulating Deformations of MR Brain Images for
 *    Validation of Atlas-based Segmentation and Registration Algorithms",
 *    submitted, 2005.
 *
 * LICENSE
 *
 * Unless otherwise stated, users of nitrc.org are hereby granted a 
 * nonexclusive, royalty-free copyright and design patent license to
 * use nitrc.org content and material in individual software or
 * documentation.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
/** #include <dfieldio.h> */
#include <float.h>
#include "headers/field.h"
#include "headers/VolumeMatching3D.h"

/********************************************************************************************************************************/	

/* Predefined constants */
#define PI 3.14159265358979

/********************************************************************************************************************************/	

/* Global Variables */
int r=256;
int c=256;
int s=198;
int BS=4;
char Deformationfieldfilename[600];
char Volumefilename[600];
char Outputfilename[600];
char *Defaultoutputfilename="defaultoutput.dat";
char Restrictionfilename[600];
int NOOUTFILE=1;
float DX=1.0;
float DY=1.0;
float DZ=1.0;
int MAXITNUM=1000;
int NORESTRICTION=1;
int VERBOSITY=1;
int NODEFORMATION=1;
float max_voxel_error=DBL_MAX;
float eps_jacobian=0.001;
float zeta_jacobian=100.0;


/********************************************************************************************************************************/	

/* Function Declarations */
void Usage(char*);
int HandleOptions(int, char**);
float VolumeMatchingCost3D(float*, float*, float*, float*, int, int, int);
int VolumeMatchingGrads3D(float*, float*, float*, float*, int, int, int, float**, float**, float**);
int VolumeMatching3D(float**, float**, float**, float*, int, int, int, int, int, char*, const FLOATING);

/********************************************************************************************************************************/	
/* Function Definitions */
float VolumeMatchingCost3D(float *f, float *g, float *e, float *TV, int r, int c, int s)
{
    float res=0.0, cres;
    int i, j, k, rc=r*c, Vind;
    int i000, i001, i010, i011, i100, i101, i110, i111;
    float f000, f001, f010, f011, f100, f101, f110, f111;
    float fx000, fx001, fx010, fx011, fx100, fx101, fx110, fx111;
    float fy000, fy001, fy010, fy011, fy100, fy101, fy110, fy111;
    float fz000, fz001, fz010, fz011, fz100, fz101, fz110, fz111;
    float g000, g001, g010, g011, g100, g101, g110, g111;
    float gx000, gx001, gx010, gx011, gx100, gx101, gx110, gx111;
    float gy000, gy001, gy010, gy011, gy100, gy101, gy110, gy111;
    float gz000, gz001, gz010, gz011, gz100, gz101, gz110, gz111;
    float e000, e001, e010, e011, e100, e101, e110, e111;
    float ex000, ex001, ex010, ex011, ex100, ex101, ex110, ex111;
    float ey000, ey001, ey010, ey011, ey100, ey101, ey110, ey111;
    float ez000, ez001, ez010, ez011, ez100, ez101, ez110, ez111;
    float J000, J001, J010, J011, J100, J101, J110, J111;
    float cV;

    max_voxel_error = 0.0;
    
    Vind=0;
    for (k=0;k<s-1;k++)
        for (j=0;j<c-1;j++)
            for (i=0;i<r-1;i++)
            {
                i000=k*rc+j*r+i;
                i001=i000+r;
                i010=i000+1;
                i011=i001+1;
                i100=i000+rc;
                i101=i100+r;
                i110=i100+1;
                i111=i101+1;
                
                f000=f[i000]; f001=f[i001]; f010=f[i010]; f011=f[i011]; f100=f[i100]; f101=f[i101]; f110=f[i110]; f111=f[i111]; 
                g000=g[i000]; g001=g[i001]; g010=g[i010]; g011=g[i011]; g100=g[i100]; g101=g[i101]; g110=g[i110]; g111=g[i111]; 
                e000=e[i000]; e001=e[i001]; e010=e[i010]; e011=e[i011]; e100=e[i100]; e101=e[i101]; e110=e[i110]; e111=e[i111]; 
                
                fx000=f001-f000; fy000=f010-f000; fz000=f100-f000; 
		fx001=f001-f000; fy001=f011-f001; fz001=f101-f001; 
		fx010=f011-f010; fy010=f010-f000; fz010=f110-f010; 
		fx011=f011-f010; fy011=f011-f001; fz011=f111-f011; 
		fx100=f101-f100; fy100=f110-f100; fz100=f100-f000; 
		fx101=f101-f100; fy101=f111-f101; fz101=f101-f001; 
		fx110=f111-f110; fy110=f110-f100; fz110=f110-f010; 
		fx111=f111-f110; fy111=f111-f101; fz111=f111-f011; 

                gx000=g001-g000; gy000=g010-g000; gz000=g100-g000; 
		gx001=g001-g000; gy001=g011-g001; gz001=g101-g001; 
		gx010=g011-g010; gy010=g010-g000; gz010=g110-g010; 
		gx011=g011-g010; gy011=g011-g001; gz011=g111-g011; 
		gx100=g101-g100; gy100=g110-g100; gz100=g100-g000; 
		gx101=g101-g100; gy101=g111-g101; gz101=g101-g001; 
		gx110=g111-g110; gy110=g110-g100; gz110=g110-g010; 
		gx111=g111-g110; gy111=g111-g101; gz111=g111-g011; 
		
                ex000=e001-e000; ey000=e010-e000; ez000=e100-e000; 
		ex001=e001-e000; ey001=e011-e001; ez001=e101-e001; 
		ex010=e011-e010; ey010=e010-e000; ez010=e110-e010; 
		ex011=e011-e010; ey011=e011-e001; ez011=e111-e011; 
		ex100=e101-e100; ey100=e110-e100; ez100=e100-e000; 
		ex101=e101-e100; ey101=e111-e101; ez101=e101-e001; 
		ex110=e111-e110; ey110=e110-e100; ez110=e110-e010; 
		ex111=e111-e110; ey111=e111-e101; ez111=e111-e011; 

                J000=fx000*(gy000*ez000-gz000*ey000)-fy000*(gx000*ez000-gz000*ex000)+fz000*(gx000*ey000-gy000*ex000);
                J001=fx001*(gy001*ez001-gz001*ey001)-fy001*(gx001*ez001-gz001*ex001)+fz001*(gx001*ey001-gy001*ex001);
                J010=fx010*(gy010*ez010-gz010*ey010)-fy010*(gx010*ez010-gz010*ex010)+fz010*(gx010*ey010-gy010*ex010);
                J011=fx011*(gy011*ez011-gz011*ey011)-fy011*(gx011*ez011-gz011*ex011)+fz011*(gx011*ey011-gy011*ex011);
                J100=fx100*(gy100*ez100-gz100*ey100)-fy100*(gx100*ez100-gz100*ex100)+fz100*(gx100*ey100-gy100*ex100);
                J101=fx101*(gy101*ez101-gz101*ey101)-fy101*(gx101*ez101-gz101*ex101)+fz101*(gx101*ey101-gy101*ex101);
                J110=fx110*(gy110*ez110-gz110*ey110)-fy110*(gx110*ez110-gz110*ex110)+fz110*(gx110*ey110-gy110*ex110);
                J111=fx111*(gy111*ez111-gz111*ey111)-fy111*(gx111*ez111-gz111*ex111)+fz111*(gx111*ey111-gy111*ex111);
                
                cres=0.0;
                cV=TV[Vind];
                if (cV>0.0)
                    cres+=(J000+J001+J010+J011+J100+J101+J110+J111-8*cV)*(J000+J001+J010+J011+J100+J101+J110+J111-8*cV);
                if (J000<eps_jacobian)
                    cres+=zeta_jacobian*(J000-eps_jacobian)*(J000-eps_jacobian);
                if (J001<eps_jacobian)
                    cres+=zeta_jacobian*(J001-eps_jacobian)*(J001-eps_jacobian);
                if (J010<eps_jacobian)
                    cres+=zeta_jacobian*(J010-eps_jacobian)*(J010-eps_jacobian);
                if (J011<eps_jacobian)
                    cres+=zeta_jacobian*(J011-eps_jacobian)*(J011-eps_jacobian);
                if (J100<eps_jacobian)
                    cres+=zeta_jacobian*(J100-eps_jacobian)*(J100-eps_jacobian);
                if (J101<eps_jacobian)
                    cres+=zeta_jacobian*(J101-eps_jacobian)*(J101-eps_jacobian);
                if (J110<eps_jacobian)
                    cres+=zeta_jacobian*(J110-eps_jacobian)*(J110-eps_jacobian);
                if (J111<eps_jacobian)
                    cres+=zeta_jacobian*(J111-eps_jacobian)*(J111-eps_jacobian);
                res+=0.5*cres;

                float error = fabs(J000+J001+J010+J011+J100+J101+J110+J111-8*cV) / 8.0;
                if (error > max_voxel_error) {
                    max_voxel_error = error;
                }
                
                Vind++;
            }
    return res;
}

/********************************************************************************************************************************/	

int VolumeMatchingGrads3D(float *f, float *g, float *e, float *TV, int r, int c, int s, float **df, float **dg, float **de)
{
    int i, j, k, rc=r*c, Vind;
    int i000, i001, i010, i011, i100, i101, i110, i111;
    float f000, f001, f010, f011, f100, f101, f110, f111;
    float fx000, fx001, fx010, fx011, fx100, fx101, fx110, fx111;
    float fy000, fy001, fy010, fy011, fy100, fy101, fy110, fy111;
    float fz000, fz001, fz010, fz011, fz100, fz101, fz110, fz111;
    float g000, g001, g010, g011, g100, g101, g110, g111;
    float gx000, gx001, gx010, gx011, gx100, gx101, gx110, gx111;
    float gy000, gy001, gy010, gy011, gy100, gy101, gy110, gy111;
    float gz000, gz001, gz010, gz011, gz100, gz101, gz110, gz111;
    float e000, e001, e010, e011, e100, e101, e110, e111;
    float ex000, ex001, ex010, ex011, ex100, ex101, ex110, ex111;
    float ey000, ey001, ey010, ey011, ey100, ey101, ey110, ey111;
    float ez000, ez001, ez010, ez011, ez100, ez101, ez110, ez111;
    float J000, J001, J010, J011, J100, J101, J110, J111;
    float cV, cE;
    float E000, E001, E010, E011, E100, E101, E110, E111;
    
    Vind=0;
    
    for (k=0;k<s-1;k++)
        for (j=0;j<c-1;j++)
            for (i=0;i<r-1;i++)
            {
                i000=k*rc+j*r+i;
                i001=i000+r;
                i010=i000+1;
                i011=i001+1;
                i100=i000+rc;
                i101=i100+r;
                i110=i100+1;
                i111=i101+1;
                
                f000=f[i000]; f001=f[i001]; f010=f[i010]; f011=f[i011]; f100=f[i100]; f101=f[i101]; f110=f[i110]; f111=f[i111]; 
                g000=g[i000]; g001=g[i001]; g010=g[i010]; g011=g[i011]; g100=g[i100]; g101=g[i101]; g110=g[i110]; g111=g[i111]; 
                e000=e[i000]; e001=e[i001]; e010=e[i010]; e011=e[i011]; e100=e[i100]; e101=e[i101]; e110=e[i110]; e111=e[i111]; 
                
                fx000=f001-f000; fy000=f010-f000; fz000=f100-f000; fx001=f001-f000; fy001=f011-f001; fz001=f101-f001; fx010=f011-f010; fy010=f010-f000; fz010=f110-f010; fx011=f011-f010; fy011=f011-f001; fz011=f111-f011; fx100=f101-f100; fy100=f110-f100; fz100=f100-f000; fx101=f101-f100; fy101=f111-f101; fz101=f101-f001; fx110=f111-f110; fy110=f110-f100; fz110=f110-f010; fx111=f111-f110; fy111=f111-f101; fz111=f111-f011; 
                gx000=g001-g000; gy000=g010-g000; gz000=g100-g000; gx001=g001-g000; gy001=g011-g001; gz001=g101-g001; gx010=g011-g010; gy010=g010-g000; gz010=g110-g010; gx011=g011-g010; gy011=g011-g001; gz011=g111-g011; gx100=g101-g100; gy100=g110-g100; gz100=g100-g000; gx101=g101-g100; gy101=g111-g101; gz101=g101-g001; gx110=g111-g110; gy110=g110-g100; gz110=g110-g010; gx111=g111-g110; gy111=g111-g101; gz111=g111-g011; 
                ex000=e001-e000; ey000=e010-e000; ez000=e100-e000; ex001=e001-e000; ey001=e011-e001; ez001=e101-e001; ex010=e011-e010; ey010=e010-e000; ez010=e110-e010; ex011=e011-e010; ey011=e011-e001; ez011=e111-e011; ex100=e101-e100; ey100=e110-e100; ez100=e100-e000; ex101=e101-e100; ey101=e111-e101; ez101=e101-e001; ex110=e111-e110; ey110=e110-e100; ez110=e110-e010; ex111=e111-e110; ey111=e111-e101; ez111=e111-e011; 
                
                J000=fx000*(gy000*ez000-gz000*ey000)-fy000*(gx000*ez000-gz000*ex000)+fz000*(gx000*ey000-gy000*ex000);
                J001=fx001*(gy001*ez001-gz001*ey001)-fy001*(gx001*ez001-gz001*ex001)+fz001*(gx001*ey001-gy001*ex001);
                J010=fx010*(gy010*ez010-gz010*ey010)-fy010*(gx010*ez010-gz010*ex010)+fz010*(gx010*ey010-gy010*ex010);
                J011=fx011*(gy011*ez011-gz011*ey011)-fy011*(gx011*ez011-gz011*ex011)+fz011*(gx011*ey011-gy011*ex011);
                J100=fx100*(gy100*ez100-gz100*ey100)-fy100*(gx100*ez100-gz100*ex100)+fz100*(gx100*ey100-gy100*ex100);
                J101=fx101*(gy101*ez101-gz101*ey101)-fy101*(gx101*ez101-gz101*ex101)+fz101*(gx101*ey101-gy101*ex101);
                J110=fx110*(gy110*ez110-gz110*ey110)-fy110*(gx110*ez110-gz110*ex110)+fz110*(gx110*ey110-gy110*ex110);
                J111=fx111*(gy111*ez111-gz111*ey111)-fy111*(gx111*ez111-gz111*ex111)+fz111*(gx111*ey111-gy111*ex111);
                
                cV=TV[Vind];
                if (cV>0.0)
                {
                    cE=(J000+J001+J010+J011+J100+J101+J110+J111-8.0*cV);
                    (*df)[i000] +=cE*(-gy000*ez000+gz000*ey000+gx000*ez000-gz000*ex000-gx000*ey000+gy000*ex000-gy001*ez001+gz001*ey001+gx010*ez010-gz010*ex010-gx100*ey100+gy100*ex100);
                    (*df)[i001] +=cE*(gy000*ez000-gz000*ey000+gy001*ez001-gz001*ey001+gx000*ez001-gz001*ex000-gx000*ey001+gy001*ex000+gx010*ez011-gz011*ex010-gx100*ey101+gy101*ex100);
                    (*df)[i010] +=cE*(-gx000*ez000+gz000*ex000-gy000*ez010+gz010*ey000-gx010*ez010+gz010*ex010-gx010*ey000+gy000*ex010-gy001*ez011+gz011*ey001-gx110*ey100+gy100*ex110);
                    (*df)[i011] +=cE*(-gx000*ez001+gz001*ex000+gy000*ez010-gz010*ey000+gy001*ez011-gz011*ey001-gx010*ez011+gz011*ex010-gx010*ey001+gy001*ex010-gx110*ey101+gy101*ex110);
                    (*df)[i100] +=cE*(gx000*ey000-gy000*ex000-gy100*ez000+gz000*ey100+gx100*ez000-gz000*ex100+gx100*ey100-gy100*ex100-gy101*ez001+gz001*ey101+gx110*ez010-gz010*ex110);
                    (*df)[i101] +=cE*(gx000*ey001-gy001*ex000+gy100*ez000-gz000*ey100+gy101*ez001-gz001*ey101+gx100*ez001-gz001*ex100+gx100*ey101-gy101*ex100+gx110*ez011-gz011*ex110);
                    (*df)[i110] +=cE*(gx010*ey000-gy000*ex010-gx100*ez000+gz000*ex100-gy100*ez010+gz010*ey100-gx110*ez010+gz010*ex110+gx110*ey100-gy100*ex110-gy101*ez011+gz011*ey101);
                    (*df)[i111] +=cE*(gx010*ey001-gy001*ex010-gx100*ez001+gz001*ex100+gy100*ez010-gz010*ey100+gy101*ez011-gz011*ey101-gx110*ez011+gz011*ex110+gx110*ey101-gy101*ex110);
                    (*dg)[i000] +=cE*(fx000*(-e100+e010)-fy000*(-e100+e001)+fz000*(-e010+e001)-fy001*(-ez001)+fz001*(-ey001)+fx010*(-ez010)+fz010*ex010+fx100*ey100-fy100*ex100);
                    (*dg)[i001] +=cE*(-fy000*ez000+fz000*ey000+fx000*(-e101+e011)-fy001*(e101-e000)+fz001*(e011-e000)+fx010*(-ez011)+fz011*ex010+fx100*ey101-fy101*ex100);
                    (*dg)[i010] +=cE*(fx000*ez000+fz000*(-ex000)+fx010*(e110-e000)-fy000*(-e110+e011)+fz010*(e000-e011)-fy001*(-ez011)+fz011*(-ey001)+fx110*ey100-fy100*ex110);
                    (*dg)[i011] +=cE*(fx000*ez001+fz001*(-ex000)-fy000*ez010+fz010*ey000+fx010*(e111-e001)-fy001*(e111-e010)+fz011*(-e001+e010)+fx110*ey101-fy101*ex110);
                    (*dg)[i100] +=cE*(fx000*(-ey000)-fy000*(-ex000)+fx100*(e000-e110)-fy100*(e000-e101)+fz000*(-e110+e101)-fy101*(-ez001)+fz001*(-ey101)+fx110*(-ez010)+fz010*ex110);
                    (*dg)[i101] +=cE*(fx000*(-ey001)-fy001*(-ex000)-fy100*ez000+fz000*ey100+fx100*(e001-e111)-fy101*(-e001+e100)+fz001*(e111-e100)+fx110*(-ez011)+fz011*ex110);
                    (*dg)[i110] +=cE*(fx010*(-ey000)-fy000*(-ex010)+fx100*ez000+fz000*(-ex100)+fx110*(-e010+e100)-fy100*(e010-e111)+fz010*(e100-e111)-fy101*(-ez011)+fz011*(-ey101));
                    (*dg)[i111] +=cE*(fx010*(-ey001)-fy001*(-ex010)+fx100*ez001+fz001*(-ex100)-fy100*ez010+fz010*ey100+fx110*(-e011+e101)-fy101*(-e011+e110)+fz011*(-e101+e110));
                    (*de)[i000] +=cE*(fx000*(-g010+g100)-fy000*(-g001+g100)+fz000*(-g001+g010)-fy001*gz001+fz001*gy001+fx010*gz010+fz010*(-gx010)+fx100*(-gy100)-fy100*(-gx100));
                    (*de)[i001] +=cE*(-fy000*(-gz000)+fz000*(-gy000)+fx000*(-g011+g101)-fy001*(g000-g101)+fz001*(g000-g011)+fx010*gz011+fz011*(-gx010)+fx100*(-gy101)-fy101*(-gx100));
                    (*de)[i010] +=cE*(fx000*(-gz000)+fz000*gx000+fx010*(g000-g110)-fy000*(-g011+g110)+fz010*(g011-g000)-fy001*gz011+fz011*gy001+fx110*(-gy100)-fy100*(-gx110));
                    (*de)[i011] +=cE*(fx000*(-gz001)+fz001*gx000-fy000*(-gz010)+fz010*(-gy000)+fx010*(g001-g111)-fy001*(g010-g111)+fz011*(-g010+g001)+fx110*(-gy101)-fy101*(-gx110));
                    (*de)[i100] +=cE*(fx000*gy000-fy000*gx000+fx100*(g110-g000)-fy100*(g101-g000)+fz000*(-g101+g110)-fy101*gz001+fz001*gy101+fx110*gz010+fz010*(-gx110));
                    (*de)[i101] +=cE*(fx000*gy001-fy001*gx000-fy100*(-gz000)+fz000*(-gy100)+fx100*(g111-g001)-fy101*(-g100+g001)+fz001*(g100-g111)+fx110*gz011+fz011*(-gx110));
                    (*de)[i110] +=cE*(fx010*gy000-fy000*gx010+fx100*(-gz000)+fz000*gx100+fx110*(-g100+g010)-fy100*(g111-g010)+fz010*(g111-g100)-fy101*gz011+fz011*gy101);
                    (*de)[i111] +=cE*(fx010*gy001-fy001*gx010+fx100*(-gz001)+fz001*gx100-fy100*(-gz010)+fz010*(-gy100)+fx110*(-g101+g011)-fy101*(-g110+g011)+fz011*(-g110+g101));
                }
                if (J000<eps_jacobian)
                {
                    E000=zeta_jacobian*(J000-eps_jacobian);
                    (*df)[i000] +=E000*(-gy000*ez000+gz000*ey000+gx000*ez000-gz000*ex000-gx000*ey000+gy000*ex000);
                    (*df)[i001] +=E000*(gy000*ez000-gz000*ey000);
                    (*df)[i010] +=E000*(-gx000*ez000+gz000*ex000);
                    (*df)[i100] +=E000*(gx000*ey000-gy000*ex000);
                    (*dg)[i000] +=E000*(fx000*(-e100+e010)-fy000*(-e100+e001)+fz000*(-e010+e001));
                    (*dg)[i001] +=E000*(-fy000*ez000+fz000*ey000);
                    (*dg)[i010] +=E000*(fx000*ez000+fz000*(-ex000));
                    (*dg)[i100] +=E000*(fx000*(-ey000)-fy000*(-ex000));
                    (*de)[i000] +=E000*(fx000*(-g010+g100)-fy000*(-g001+g100)+fz000*(-g001+g010));
                    (*de)[i001] +=E000*(-fy000*(-gz000)+fz000*(-gy000));
                    (*de)[i010] +=E000*(fx000*(-gz000)+fz000*gx000);
                    (*de)[i100] +=E000*(fx000*gy000-fy000*gx000);
                }
                if (J001<eps_jacobian)
                {
                    E001=zeta_jacobian*(J001-eps_jacobian);
                    (*df)[i000] +=E001*(-gy001*ez001+gz001*ey001);
                    (*df)[i001] +=E001*(gy001*ez001-gz001*ey001+gx000*ez001-gz001*ex000-gx000*ey001+gy001*ex000);
                    (*df)[i011] +=E001*(-gx000*ez001+gz001*ex000);
                    (*df)[i101] +=E001*(gx000*ey001-gy001*ex000);
                    (*dg)[i000] +=E001*(-fy001*(-ez001)+fz001*(-ey001));
                    (*dg)[i001] +=E001*(fx000*(-e101+e011)-fy001*(e101-e000)+fz001*(e011-e000));
                    (*dg)[i011] +=E001*(fx000*ez001+fz001*(-ex000));
                    (*dg)[i101] +=E001*(fx000*(-ey001)-fy001*(-ex000));
                    (*de)[i000] +=E001*(-fy001*gz001+fz001*gy001);
                    (*de)[i001] +=E001*(fx000*(-g011+g101)-fy001*(g000-g101)+fz001*(g000-g011));
                    (*de)[i011] +=E001*(fx000*(-gz001)+fz001*gx000);
                    (*de)[i101] +=E001*(fx000*gy001-fy001*gx000);
                }
                if (J010<eps_jacobian)
                {
                        E010=zeta_jacobian*(J010-eps_jacobian);
                    (*df)[i000] +=E010*(gx010*ez010-gz010*ex010);
                    (*df)[i010] +=E010*(-gy000*ez010+gz010*ey000-gx010*ez010+gz010*ex010-gx010*ey000+gy000*ex010);
                    (*df)[i011] +=E010*(gy000*ez010-gz010*ey000);
                    (*df)[i110] +=E010*(gx010*ey000-gy000*ex010);
                    (*dg)[i000] +=E010*(fx010*(-ez010)+fz010*ex010);
                    (*dg)[i010] +=E010*(fx010*(e110-e000)-fy000*(-e110+e011)+fz010*(e000-e011));
                    (*dg)[i011] +=E010*(-fy000*ez010+fz010*ey000);
                    (*dg)[i110] +=E010*(fx010*(-ey000)-fy000*(-ex010));
                    (*de)[i000] +=E010*(fx010*gz010+fz010*(-gx010));
                    (*de)[i010] +=E010*(fx010*(g000-g110)-fy000*(-g011+g110)+fz010*(g011-g000));
                    (*de)[i011] +=E010*(-fy000*(-gz010)+fz010*(-gy000));
                    (*de)[i110] +=E010*(fx010*gy000-fy000*gx010);
                }
                if (J011<eps_jacobian)
                {
                    E011=zeta_jacobian*(J011-eps_jacobian);
                    (*df)[i001] +=E011*(gx010*ez011-gz011*ex010);
                    (*df)[i010] +=E011*(-gy001*ez011+gz011*ey001);
                    (*df)[i011] +=E011*(gy001*ez011-gz011*ey001-gx010*ez011+gz011*ex010-gx010*ey001+gy001*ex010);
                    (*df)[i111] +=E011*(gx010*ey001-gy001*ex010);
                    (*dg)[i001] +=E011*(fx010*(-ez011)+fz011*ex010);
                    (*dg)[i010] +=E011*(-fy001*(-ez011)+fz011*(-ey001));
                    (*dg)[i011] +=E011*(fx010*(e111-e001)-fy001*(e111-e010)+fz011*(-e001+e010));
                    (*dg)[i111] +=E011*(fx010*(-ey001)-fy001*(-ex010));
                    (*de)[i001] +=E011*(fx010*gz011+fz011*(-gx010));
                    (*de)[i010] +=E011*(-fy001*gz011+fz011*gy001);
                    (*de)[i011] +=E011*(fx010*(g001-g111)-fy001*(g010-g111)+fz011*(-g010+g001));
                    (*de)[i111] +=E011*(fx010*gy001-fy001*gx010);
                }
                if (J100<eps_jacobian)
                {
                    E100=zeta_jacobian*(J100-eps_jacobian);
                    (*df)[i000] +=E100*(-gx100*ey100+gy100*ex100);
                    (*df)[i100] +=E100*(-gy100*ez000+gz000*ey100+gx100*ez000-gz000*ex100+gx100*ey100-gy100*ex100);
                    (*df)[i101] +=E100*(gy100*ez000-gz000*ey100);
                    (*df)[i110] +=E100*(-gx100*ez000+gz000*ex100);
                    (*dg)[i000] +=E100*(fx100*ey100-fy100*ex100);
                    (*dg)[i100] +=E100*(fx100*(e000-e110)-fy100*(e000-e101)+fz000*(-e110+e101));
                    (*dg)[i101] +=E100*(-fy100*ez000+fz000*ey100);
                    (*dg)[i110] +=E100*(fx100*ez000+fz000*(-ex100));
                    (*de)[i000] +=E100*(fx100*(-gy100)-fy100*(-gx100));
                    (*de)[i100] +=E100*(fx100*(g110-g000)-fy100*(g101-g000)+fz000*(-g101+g110));
                    (*de)[i101] +=E100*(-fy100*(-gz000)+fz000*(-gy100));
                    (*de)[i110] +=E100*(fx100*(-gz000)+fz000*gx100);
                }
                if (J101<eps_jacobian)
                {
                    E101=zeta_jacobian*(J101-eps_jacobian);
                    (*df)[i001] +=E101*(-gx100*ey101+gy101*ex100);
                    (*df)[i100] +=E101*(-gy101*ez001+gz001*ey101);
                    (*df)[i101] +=E101*(gy101*ez001-gz001*ey101+gx100*ez001-gz001*ex100+gx100*ey101-gy101*ex100);
                    (*df)[i111] +=E101*(-gx100*ez001+gz001*ex100);
                    (*dg)[i001] +=E101*(fx100*ey101-fy101*ex100);
                    (*dg)[i100] +=E101*(-fy101*(-ez001)+fz001*(-ey101));
                    (*dg)[i101] +=E101*(fx100*(e001-e111)-fy101*(-e001+e100)+fz001*(e111-e100));
                    (*dg)[i111] +=E101*(fx100*ez001+fz001*(-ex100));
                    (*de)[i001] +=E101*(fx100*(-gy101)-fy101*(-gx100));
                    (*de)[i100] +=E101*(-fy101*gz001+fz001*gy101);
                    (*de)[i101] +=E101*(fx100*(g111-g001)-fy101*(-g100+g001)+fz001*(g100-g111));
                    (*de)[i111] +=E101*(fx100*(-gz001)+fz001*gx100);
                }
                if (J110<eps_jacobian)
                {
                    E110=zeta_jacobian*(J110-eps_jacobian);
                    (*df)[i010] +=E110*(-gx110*ey100+gy100*ex110);
                    (*df)[i100] +=E110*(gx110*ez010-gz010*ex110);
                    (*df)[i110] +=E110*(-gy100*ez010+gz010*ey100-gx110*ez010+gz010*ex110+gx110*ey100-gy100*ex110);
                    (*df)[i111] +=E110*(gy100*ez010-gz010*ey100);
                    (*dg)[i010] +=E110*(fx110*ey100-fy100*ex110);
                    (*dg)[i100] +=E110*(fx110*(-ez010)+fz010*ex110);
                    (*dg)[i110] +=E110*(fx110*(-e010+e100)-fy100*(e010-e111)+fz010*(e100-e111));
                    (*dg)[i111] +=E110*(-fy100*ez010+fz010*ey100);
                    (*de)[i010] +=E110*(fx110*(-gy100)-fy100*(-gx110));
                    (*de)[i100] +=E110*(fx110*gz010+fz010*(-gx110));
                    (*de)[i110] +=E110*(fx110*(-g100+g010)-fy100*(g111-g010)+fz010*(g111-g100));
                    (*de)[i111] +=E110*(-fy100*(-gz010)+fz010*(-gy100));
                }
                if (J111<eps_jacobian)
                {
                    E111=zeta_jacobian*(J111-eps_jacobian);
                    (*df)[i011] +=E111*(-gx110*ey101+gy101*ex110);
                    (*df)[i101] +=E111*(gx110*ez011-gz011*ex110);
                    (*df)[i110] +=E111*(-gy101*ez011+gz011*ey101);
                    (*df)[i111] +=E111*(gy101*ez011-gz011*ey101-gx110*ez011+gz011*ex110+gx110*ey101-gy101*ex110);
                    (*dg)[i011] +=E111*(fx110*ey101-fy101*ex110);
                    (*dg)[i101] +=E111*(fx110*(-ez011)+fz011*ex110);
                    (*dg)[i110] +=E111*(-fy101*(-ez011)+fz011*(-ey101));
                    (*dg)[i111] +=E111*(fx110*(-e011+e101)-fy101*(-e011+e110)+fz011*(-e101+e110));
                    (*de)[i011] +=E111*(fx110*(-gy101)-fy101*(-gx110));
                    (*de)[i101] +=E111*(fx110*gz011+fz011*(-gx110));
                    (*de)[i110] +=E111*(-fy101*gz011+fz011*gy101);
                    (*de)[i111] +=E111*(fx110*(-g101+g011)-fy101*(-g110+g011)+fz011*(-g110+g101));
                }
                
                Vind++;
            }
    return 1;
}

/********************************************************************************************************************************/	

int VolumeMatching3D(float **f, float **g, float **e, float *TV, int r, int c, int s, int maxitnum, int norestriction, char *RestrictionMap, const FLOATING epsilon)
{
    int i, j, k, ind, rc=r*c, rcs=r*c*s, itnum;
    float *nf = NULL, *ng = NULL, *ne = NULL;
    float *df = NULL, *dg = NULL, *de = NULL;
    float *H = NULL, Hmark, nH, alpha=1.0;

    (void) j;
    (void) k;
    (void) ind;
    (void) rc;
    
    H=(float *)malloc(sizeof(float)*(maxitnum+1));
    nf=(float *)malloc(sizeof(float)*rcs);
    ng=(float *)malloc(sizeof(float)*rcs);
    ne=(float *)malloc(sizeof(float)*rcs);
    df=(float *)malloc(sizeof(float)*rcs);
    dg=(float *)malloc(sizeof(float)*rcs);
    de=(float *)malloc(sizeof(float)*rcs);
    
    H[0]=VolumeMatchingCost3D(*f,*g,*e,TV,r,c,s);
    
    itnum=0;
    while ((itnum<maxitnum)&&(alpha>0.000000001)&&(max_voxel_error>epsilon))
    {
        Hmark=H[itnum];
        VolumeMatchingGrads3D((*f),(*g),(*e),TV,r,c,s,&df,&dg,&de);
        for (i=0;i<rcs;i++)
        {
            nf[i]=(*f)[i]-alpha*df[i];
            ng[i]=(*g)[i]-alpha*dg[i];
            ne[i]=(*e)[i]-alpha*de[i];
        }
        nH=VolumeMatchingCost3D(nf,ng,ne,TV,r,c,s);

	/** printf("%d, %f\n", itnum, H[itnum]); */
         if (nH<Hmark)
        {
            if (itnum==0)
            {
                while (nH<Hmark)
                {
                    if (VERBOSITY) {printf("/"); fflush(stdout);}
                    Hmark=nH;
                    alpha*=2.0;
                    for (i=0;i<rcs;i++)
                    {
                        nf[i]=(*f)[i]-alpha*df[i];
                        ng[i]=(*g)[i]-alpha*dg[i];
                        ne[i]=(*e)[i]-alpha*de[i];
                    }
                    nH=VolumeMatchingCost3D(nf,ng,ne,TV,r,c,s);
                }
                if (VERBOSITY) {printf("\\"); fflush(stdout);}
                alpha*=0.5;
                for (i=0;i<rcs;i++)
                {
                    nf[i]=(*f)[i]-alpha*df[i];
                    ng[i]=(*g)[i]-alpha*dg[i];
                    ne[i]=(*e)[i]-alpha*de[i];
                }
                H[itnum+1]=Hmark;
            } else
            {
                H[itnum+1]=nH;
                alpha*=1.05;
            }
        } else
        {
            while ((nH>Hmark)&&(alpha>.000000001))
            {
                if (VERBOSITY) {printf("\\"); fflush(stdout);}
                alpha*=0.5;
                for (i=0;i<rcs;i++)
                {
                    nf[i]=(*f)[i]-alpha*df[i];
                    ng[i]=(*g)[i]-alpha*dg[i];
                    ne[i]=(*e)[i]-alpha*de[i];
                }
                nH=VolumeMatchingCost3D(nf,ng,ne,TV,r,c,s);
            }
            H[itnum+1]=nH;
        }
        if (alpha>.000000001)
        {
            if (norestriction)
                for (i=0;i<rcs;i++)
                {
                    (*f)[i]=nf[i];
                    (*g)[i]=ng[i];
                    (*e)[i]=ne[i];
                    df[i]=0.0;
                    dg[i]=0.0;
                    de[i]=0.0;
                }
            else
                for (i=0;i<rcs;i++)
                {
                    if (RestrictionMap[i])
                    {
                        (*f)[i]=nf[i];
                        (*g)[i]=ng[i];
                        (*e)[i]=ne[i];
                    }
                    df[i]=0.0;
                    dg[i]=0.0;
                    de[i]=0.0;
                }
            itnum++;
        verbose_printf(true,
                       "Iteration %5d:  "
                       "total error %6e  "
                       "max voxel error %6e  "
                       "eta %6e\n",
                       itnum, H[itnum], max_voxel_error, alpha);
        }
        if (VERBOSITY) {printf("."); fflush(stdout);}
    }
    
    if (VERBOSITY)
    {
        printf("\n");
        for (i=0;i<maxitnum+1;i++)
        {
            printf("%12.4f\n",H[i]);
            fflush(stdout);
        }
    }
    /* iterations have been exhaused, freeing memory and returning*/
    free(nf);
    free(ng);
    free(ne);
    free(df);
    free(dg);
    free(de);
    free(H);

    return 1;
}

/********************************************************************************************************************************/	

void Usage(char *programName)
{
    int i,n=0;
    char *emptyspace = NULL;
    
    while (programName[n]!='\0')
        n++;
    emptyspace=(char *)malloc(sizeof(char)*n);
    for (i=0;i<n;i++)
        emptyspace[i]=' ';
    emptyspace[n]='\0';
	fprintf(stderr,"%s usage:\n",programName);
	fprintf(stderr,"%s finds the deformation fields that agree with the prescribed values for cubic volumes to the extent possible. ",programName);
    fprintf(stderr,"Specifically, it reads an initial deformation field and the desired set of cubic volumes, and employs a ");
    fprintf(stderr,"gradient descent algorithm to minimize the squared error between the prescribed set of cubic columes to the average of corner ");
    fprintf(stderr,"Jacobians of the deformation field. Note that it reads one desired volume per cubic patch defined by 8 voxels ");
    fprintf(stderr,"at its vertices; thus, for a r x c x s deformation field, it expects to read (r-1)(c-1)(s-1) cubic volume values. ");
    fprintf(stderr,"It assumes the deformation field is defined on an isotropic grid with equal spacing in x, y, and z directions. %s takes the following options:\n",programName);
    fprintf(stderr,"%s -d <displacement field> \n",emptyspace);
    fprintf(stderr,"%s -v <prescribed volumes> \n",emptyspace);
    fprintf(stderr,"%s -o <output file name [defaultoutput.dat]> \n",emptyspace);
    fprintf(stderr,"%s -r <number of rows (y dim) [256]> \n",emptyspace);
    fprintf(stderr,"%s -c <number of columns (x dim) [256]> \n",emptyspace);
    fprintf(stderr,"%s -s <number of slices (z dim) [198]> \n",emptyspace);
    fprintf(stderr,"%s -I <number of iterations [1000]> \n",emptyspace);
    fprintf(stderr,"%s -R <restrictions file (char)> \n",emptyspace);
    fprintf(stderr,"%s -V <verbosity [1]> \n",emptyspace);
    
    free(emptyspace);
}

/********************************************************************************************************************************/	

/* returns the index of the first argument that is not an option; i.e.
   does not start with a dash or a slash
*/
int HandleOptions(int argc,char *argv[])
{
	int i,firstnonoption=0;
	float tau,popt;

    (void) firstnonoption;
    (void) tau;
    (void) popt;

	for (i=1; i< argc;i++) {
		if (argv[i][0] == '/' || argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* An argument -? means help is requested */
				case '?':
					Usage(argv[0]);
					break;
				/** case 'h': */
				/** case 'H': */
				/**     if (!stricmp(argv[i]+1,"help")) { */
				/**         Usage(argv[0]); */
				/**         break; */
				/**     } */
					/* If the option -h means anything else
					 * in your application add code here
					 * Note: this falls through to the default
					 * to print an "unknow option" message
					*/
				/* add your option switches here */
				case 'd':
					strcpy((char *)&Deformationfieldfilename,argv[i+1]);
                    NODEFORMATION=0;
					break;
				case 'v':
					strcpy((char *)&Volumefilename,argv[i+1]);
					break;
				case 'R':
					strcpy((char *)&Restrictionfilename,argv[i+1]);
                    NORESTRICTION=0;
					break;
				case 'o':
					strcpy((char *)&Outputfilename,argv[i+1]);
					NOOUTFILE=0;
					break;
				case 'r':
					r=atoi(argv[i+1]);
                    if (r<2)
                    {
                        fprintf(stderr,"number of rows too small (%d).\n",r);
                        fprintf(stderr,"bailing out ..\n");
                        return 0;
                    }
					break;
				case 'c':
					c=atoi(argv[i+1]);
                    if (c<2)
                    {
                        fprintf(stderr,"number of columns too small (%d).\n",c);
                        fprintf(stderr,"bailing out ..\n");
                        return 0;
                    }
					break;
				case 's':
					s=atoi(argv[i+1]);
                    if (s<2)
                    {
                        fprintf(stderr,"number of slices too small (%d).\n",s);
                        fprintf(stderr,"bailing out ..\n");
                        return 0;
                    }
					break;
				case 'I':
					MAXITNUM=atoi(argv[i+1]);
					break;
				case 'V':
					VERBOSITY=atoi(argv[i+1]);
					break;
				default:
					fprintf(stderr,"unknown option %s\n",argv[i]);
					break;
			}
		}
	}
	return 1;
}



/********************************************************************************************************************************/	

/*!
 * \brief Wrapper function.
 */
void volume_matching_3d(
        const size_t nx,               /*!< Width of the image  */
        const size_t ny,               /*!< Length of the image */
        const size_t nz,               /*!< Depth of the image  */
        const FLOATING dx,             /*!< x spacing */
        const FLOATING dy,             /*!< y spacing */
        const FLOATING dz,             /*!< z spacing */
        const FLOATING J[nz][ny][nx],  /*!< Target Jacobian */
        const bool mask[nz][ny][nx],   /*!< Body mask */
        const FLOATING epsilon,        /*!< Tolerance on the Jacobian per voxel */
        const FLOATING tolerance,      /*!< Jacobian tolerance on background */
        FLOATING eta,                  /*!< Initial step length for the optimisation */
        const FLOATING alpha,          /*!< Step length increase coefficient */
        const FLOATING beta,           /*!< Step length decrease coefficient */
        const FLOATING gamma,          /*!< Armijo-Goldstein parameter */
        const FLOATING delta,          /*!< Jacobian regularisation threshold */
        const FLOATING zeta,           /*!< Jacobian regularisation weight */
        const bool strict,             /*!< Always improve maximum voxel error */
        const size_t it_max,           /*!< Maximum number of iterations */
        FLOATING field[3][nz][ny][nx]  /*!< Resulting displacement field */
        )
{
    verbose_printf(DISPTOOLS_DEBUG,
                   "%s\n"
                   "nx:        %lu\n"
                   "ny:        %lu\n"
                   "nz:        %lu\n"
                   "dx:        %f\n"
                   "dy:        %f\n"
                   "dz:        %f\n"
                   "alpha:     %f\n"
                   "beta:      %f\n"
                   "gamma:     %f\n"
                   "delta:     %f\n"
                   "epsilon:   %f\n"
                   "zeta:      %f\n"
                   "eta:       %f\n"
                   "tolerance: %f\n"
                   "strict:    %d\n"
                   "it_max:    %lu\n",
                   __func__,
                   nx, ny, nz,
                   dx, dy, dz,
                   alpha, beta, gamma, delta,
                   epsilon, zeta, eta,
                   tolerance,
                   strict,
                   it_max);

    float *f = NULL, *g = NULL, *e = NULL, *TV = NULL;
    char *RestrictionMap = NULL;
	int i,j,k,ind,L;

    (void) mask;
    (void) tolerance;
    (void) eta;
    (void) alpha;
    (void) beta;
    (void) gamma;
    (void) strict;

    r = nx;
    c = ny;
    s = nz;
    NORESTRICTION = 1;
    MAXITNUM = it_max;
    VERBOSITY = 0;
    eps_jacobian = delta; 
    zeta_jacobian = zeta;

    L=(r-1)*(c-1)*(s-1);
    f=(float *)calloc(r*c*s,sizeof(float));
    g=(float *)calloc(r*c*s,sizeof(float));
    e=(float *)calloc(r*c*s,sizeof(float));
    ind=0;
    for (k=1;k<s+1;k++)
        for (j=1;j<c+1;j++)
            for (i=1;i<r+1;i++)
            {
                // Add displacement component (initial guess)
                f[ind]=(float)j + field[Y][k-1][j-1][i-1] / dy;
                g[ind]=(float)i + field[X][k-1][j-1][i-1] / dx;
                e[ind]=(float)k + field[Z][k-1][j-1][i-1] / dz;
                ind++;
            }
	
    TV=(float *)malloc(sizeof(float)*L);
    for (k=0;k<s-1;k++)
        for (j=0;j<c-1;j++)
            for (i=0;i<r-1;i++)
                TV[k*(c-1)*(r-1) + j*(r-1) + i] = J[k][j][i];
    
    VolumeMatching3D(&f,&g,&e,TV,r,c,s,MAXITNUM,NORESTRICTION,RestrictionMap,epsilon);

    // Convert from the tool's coordinate system to ITK's.
    // Also convert from deformation field to displacement field, i.e.
    // subtract the identity transformation, and scale back to the input
    // grid spacing.
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
                field[X][z][y][x] = (g[z*nx*ny + y*nx + x] - (x+1)) * dx;
                field[Y][z][y][x] = (f[z*nx*ny + y*nx + x] - (y+1)) * dy;
                field[Z][z][y][x] = (e[z*nx*ny + y*nx + x] - (z+1)) * dz;
            }
        }
    }
            
    free(f);
    free(g);
    free(e);
    free(TV);
    free(RestrictionMap);
}


