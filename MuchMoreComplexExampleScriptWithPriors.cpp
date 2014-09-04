#include "preamble.h" 
#ifdef  Global015prior  
#include <stdlib.h> 
#include <math.h> 
#include <string.h>
#include "filzbach.h"
#include <iostream> 
#include <boost/array.hpp> 
#include <boost/numeric/odeint.hpp>
using namespace std;
using namespace boost::numeric::odeint;
void pause(){PAUSE}  
void fake_data(); 
void read_data(); 
void setup_parameters();
void fit_model();
void final_output();
double aTSRot;
 double aTLRot;
 double aTOCop;
 double aTDapR;
 double aTAsp;
 double aTMes;
 double aTDapP;
 double aTHol;
 double aTBos;
 double aJSRot;
 double aJLRot;
 double aJOCop;
 double aJDapR;
 double aJAsp;
 double aJMes;
 double aJDapP;
 double aJHol;
 double aJBos;
 double aRBac;
 double aRDia;
 double aRCr;
 double aRSGr;
 double aRLGr;
 double aRSGo;
 double aRLGo;
 double aRInC;

 double	start1Bac	=	0.0050	;
double	start1Dia	=	0.0181	;
double	start1Cr	=	0.0070	;
double	start1SGr	=	0.0010	;
double	start1LGr	=	0.0205	;
double	start1SGo	=	0.0156	;
double	start1LGo	=	0.0620	;
double	start1InC	=	0.0173	;
double	start1SRot	=	0.9854	;
double	start1LRot	=	0.1002	;
double	start1OCop	=	22.6953	;
double	start1DapR	=	0.1705	;
double	start1Asp	=	8.5935	;
double	start1Mes	=	0.7572	;
double	start1DapP	=	1.8839	;
double	start1Hol	=	9.6428	;
double	start1Bos	=	0.1002	;
double	start2Bac	=	0.0218	;
double	start2Dia	=	0.0010	;
double	start2Cr	=	0.0702	;
double	start2SGr	=	0.0010	;
double	start2LGr	=	0.0069	;
double	start2SGo	=	0.0249	;
double	start2LGo	=	0.0450	;
double	start2InC	=	0.0035	;
double	start2SRot	=	3.2703	;
double	start2LRot	=	0.1026	;
double	start2OCop	=	54.6490	;
double	start2DapR	=	1.0775	;
double	start2Asp	=	0.6661	;
double	start2Mes	=	0.1005	;
double	start2DapP	=	1.9573	;
double	start2Hol	=	0.3979	;
double	start2Bos	=	0.0996	;
double	start3Bac	=	0.0291	;
double	start3Dia	=	0.0119	;
double	start3Cr	=	0.0079	;
double	start3SGr	=	0.0010	;
double	start3LGr	=	0.0010	;
double	start3SGo	=	0.0171	;
double	start3LGo	=	0.0120	;
double	start3InC	=	0.01090	;
double	start3SRot	=	2.5080	;
double	start3LRot	=	1.4960	;
double	start3OCop	=	49.7455	;
double	start3DapR	=	0.1023	;
double	start3Asp	=	0.1015	;
double	start3Mes	=	0.8759	;
double	start3DapP	=	4.5128	;
double	start3Hol	=	1.7653	;
double	start3Bos	=	0.10	;
double	start4Bac	=	0.0330	;
double	start4Dia	=	0.0147	;
double	start4Cr	=	0.0428	;
double	start4SGr	=	0.0039	;
double	start4LGr	=	0.0010	;
double	start4SGo	=	0.0010	;
double	start4LGo	=	0.0180	;
double	start4InC	=	0.0010	;
double	start4SRot	=	1.3421	;
double	start4LRot	=	0.0431	;
double	start4OCop	=	0.1119	;
double	start4DapR	=	0.3765	;
double	start4Asp	=	14.7573	;
double	start4Mes	=	71.8387	;
double	start4DapP	=	0.1542	;
double	start4Hol	=	92.4456	;
double	start4Bos	=	8.3233	;
double	start5Bac	=	0.0090	;
double	start5Dia	=	0.0030	;
double	start5Cr	=	0.0282	;
double	start5SGr	=	0.0010	;
double	start5LGr	=	0.0116	;
double	start5SGo	=	0.0107	;
double	start5LGo	=	0.0803	;
double	start5InC	=	0.0010	;
double	start5SRot	=	4.8803	;
double	start5LRot	=	0.0986	;
double	start5OCop	=	70.6117	;
double	start5DapR	=	3.1821	;
double	start5Asp	=	0.0999	;
double	start5Mes	=	1.3518	;
double	start5DapP	=	4.4872	;
double	start5Hol	=	18.7869	;
double	start5Bos	=	0.3588	;

double fT =  0.9 ;
double fJ =  0.9 ;
      double fR =  0.9 ;
      double d =  0.005 ;
      double K =  0.91 ;  
      double sdZ =  2.5 ; 
      double sdP =  1.4 ;
      double W =  1 ;
      double q =  1.2 ;
      double PhytoTotal; 
 double Q0; double W0;
double Q1; double W1;
double Q2; double W2;
double Q3; double W3;
double Q4; double W4;
double Q5; double W5;
double Q6; double W6;
double Q7; double W7;
double Q8; double W8;
double Q9; double W9;
double Q10; double W10;
double Q11; double W11;
double Q12; double W12;
double Q13; double W13;
double Q14; double W14;
double Q15; double W15;
double Q16; double W16;
 
// Note that these are mass ^ -0.15 
const double imBac= 8.814	;
const double imDia= 3.127	;
const double imCr= 	5.068	;
const double imSGr= 4.568	;
const double imLGr= 3.127	;
const double imSGo=	5.068	;
const double imLGo= 3.588	;
const double imInC= 3.588	;
const double imSRot= 2.115	;
const double imLRot= 0.784	;
const double imOCop=  0.864	;
const double imDapR=  0.682	;
const double imAsp=  0.673	;
const double imMes=  0.515	;
const double imDapP= 0.550	;
const double imHol= 0.538	;
const double imBos= 1.031	;

 
                      
                      typedef boost::array< double , 17 > state_type; // Create internal array the size of what needs to be tracked
                      typedef runge_kutta_cash_karp54< state_type > error_stepper_type; // Type of stepper used see http://headmyshoulder.github.io/odeint-v2/index.html docs for details
                      typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
                      double abs_err = 1e-12 , rel_err = 1e-06 , a_x = 1.0 , a_dxdt = 1.0;
                      controlled_stepper_type controlled_stepper(
                      default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
                      int timerow = 0; 

 void ODE( const state_type &x , state_type &dxdt , double t )
{W0=x[0]/W;Q0=pow(W0,q);W1=x[1]/W;Q1=pow(W1,q);W2=x[2]/W;Q2=pow(W2,q);W3=x[3]/W;Q3=pow(W3,q);W4=x[4]/W;Q4=pow(W4,q);W5=x[5]/W;Q5=pow(W5,q);W6=x[6]/W;Q6=pow(W6,q);W7=x[7]/W;Q7=pow(W7,q);W8=x[8]/W;Q8=pow(W8,q);W9=x[9]/W;Q9=pow(W9,q);W10=x[10]/W;Q10=pow(W10,q);W11=x[11]/W;Q11=pow(W11,q);W12=x[12]/W;Q12=pow(W12,q);W13=x[13]/W;Q13=pow(W13,q);W14=x[14]/W;Q14=pow(W14,q);W15=x[15]/W;Q15=pow(W15,q);W16=x[16]/W;Q16=pow(W16,q);PhytoTotal = 0+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]; 
 /*  Bac     */ dxdt[0] = fR*aRBac*imBac*x[0]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJSRot*imSRot*2*x[8] * Q0/(1+d*x[8]+(Q0+Q2+Q3+Q5+0))) // Impact of SRot
                           -(fJ*aJLRot*imLRot*2*x[9] * Q0/(1+d*x[9]+(Q0+Q2+Q3+Q5+0))) // Impact of LRot
                           -(fJ*aJDapR*imDapR*2*x[11] * Q0/(1+d*x[11]+(Q0+Q1+Q2+Q4+Q8+0))) // Impact of DapR
                           -(fJ*aJDapP*imDapP*2*x[14] * Q0/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q0/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           -(fJ*aJBos*imBos*2*x[16] * Q0/(1+d*x[16]+(Q0+Q2+Q3+0))) // Impact of Bos
                           ;
 /*  Dia     */ dxdt[1] = fR*aRDia*imDia*x[1]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJOCop*imOCop*2*x[10] * Q1/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJDapR*imDapR*2*x[11] * Q1/(1+d*x[11]+(Q0+Q1+Q2+Q4+Q8+0))) // Impact of DapR
                           -(fJ*aJAsp*imAsp*2*x[12] * Q1/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJDapP*imDapP*2*x[14] * Q1/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q1/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           ;
 /*  Cr      */ dxdt[2] = fR*aRCr*imCr*x[2]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJSRot*imSRot*2*x[8] * Q2/(1+d*x[8]+(Q0+Q2+Q3+Q5+0))) // Impact of SRot
                           -(fJ*aJLRot*imLRot*2*x[9] * Q2/(1+d*x[9]+(Q0+Q2+Q3+Q5+0))) // Impact of LRot
                           -(fJ*aJOCop*imOCop*2*x[10] * Q2/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJDapR*imDapR*2*x[11] * Q2/(1+d*x[11]+(Q0+Q1+Q2+Q4+Q8+0))) // Impact of DapR
                           -(fJ*aJAsp*imAsp*2*x[12] * Q2/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJDapP*imDapP*2*x[14] * Q2/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q2/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           -(fJ*aJBos*imBos*2*x[16] * Q2/(1+d*x[16]+(Q0+Q2+Q3+0))) // Impact of Bos
                           ;
 /*  SGr     */ dxdt[3] = fR*aRSGr*imSGr*x[3]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJSRot*imSRot*2*x[8] * Q3/(1+d*x[8]+(Q0+Q2+Q3+Q5+0))) // Impact of SRot
                           -(fJ*aJLRot*imLRot*2*x[9] * Q3/(1+d*x[9]+(Q0+Q2+Q3+Q5+0))) // Impact of LRot
                           -(fJ*aJOCop*imOCop*2*x[10] * Q3/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJAsp*imAsp*2*x[12] * Q3/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJDapP*imDapP*2*x[14] * Q3/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q3/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           -(fJ*aJBos*imBos*2*x[16] * Q3/(1+d*x[16]+(Q0+Q2+Q3+0))) // Impact of Bos
                           ;
 /*  LGr     */ dxdt[4] = fR*aRLGr*imLGr*x[4]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJOCop*imOCop*2*x[10] * Q4/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJDapR*imDapR*2*x[11] * Q4/(1+d*x[11]+(Q0+Q1+Q2+Q4+Q8+0))) // Impact of DapR
                           -(fJ*aJAsp*imAsp*2*x[12] * Q4/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJDapP*imDapP*2*x[14] * Q4/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q4/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           ;
 /*  SGo     */ dxdt[5] = fR*aRSGo*imSGo*x[5]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJSRot*imSRot*2*x[8] * Q5/(1+d*x[8]+(Q0+Q2+Q3+Q5+0))) // Impact of SRot
                           -(fJ*aJLRot*imLRot*2*x[9] * Q5/(1+d*x[9]+(Q0+Q2+Q3+Q5+0))) // Impact of LRot
                           -(fJ*aJOCop*imOCop*2*x[10] * Q5/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJAsp*imAsp*2*x[12] * Q5/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJDapP*imDapP*2*x[14] * Q5/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           -(fJ*aJHol*imHol*2*x[15] * Q5/(1+d*x[15]+(Q0+Q1+Q2+Q3+Q4+Q5+0))) // Impact of Hol
                           ;
 /*  LGo     */ dxdt[6] = fR*aRLGo*imLGo*x[6]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJDapP*imDapP*2*x[14] * Q6/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           ;
 /*  InC     */ dxdt[7] = fR*aRInC*imInC*x[7]*(1 -(PhytoTotal/K)) 
                              -(fJ*aJDapP*imDapP*2*x[14] * Q7/(1+d*x[14]+(Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7+0))) // Impact of DapP
                           ;
 /*  SRot    */ dxdt[8] =  fJ*aJSRot*imSRot*x[8] *(0+Q0+Q2+Q3+Q5)/(1 + d*x[8] +(0+Q0+Q2+Q3+Q5))//  Gain from predation 
                           -(aTSRot*fT*imSRot*x[8]) // Loss from Respiration 
                            -(fJ*aJOCop*imOCop*2*x[10] * Q8/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJDapR*imDapR*2*x[11] * Q8/(1+d*x[11]+(Q0+Q1+Q2+Q4+Q8+0))) // Impact of DapR
                           -(fJ*aJAsp*imAsp*2*x[12] * Q8/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJMes*imMes*2*x[13] * Q8/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;
 /*  LRot    */ dxdt[9] =  fJ*aJLRot*imLRot*x[9] *(0+Q0+Q2+Q3+Q5)/(1 + d*x[9] +(0+Q0+Q2+Q3+Q5))//  Gain from predation 
                           -(aTLRot*fT*imLRot*x[9]) // Loss from Respiration 
                            -(fJ*aJOCop*imOCop*2*x[10] * Q9/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJAsp*imAsp*2*x[12] * Q9/(1+d*x[12]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+0))) // Impact of Asp
                           -(fJ*aJMes*imMes*2*x[13] * Q9/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;
 /*  OCop    */ dxdt[10] =  fJ*aJOCop*imOCop*x[10] *(0+Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10)/(1 + d*x[10] +(0+Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10))//  Gain from predation 
                           -(aTOCop*fT*imOCop*x[10]) // Loss from Respiration 
                            -(fJ*aJOCop*imOCop*2*x[10] * Q10/(1+d*x[10]+(Q1+Q2+Q3+Q4+Q5+Q8+Q9+Q10+0))) // Impact of OCop
                           -(fJ*aJMes*imMes*2*x[13] * Q10/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;
 /*  DapR    */ dxdt[11] =  fJ*aJDapR*imDapR*x[11] *(0+Q0+Q1+Q2+Q4+Q8)/(1 + d*x[11] +(0+Q0+Q1+Q2+Q4+Q8))//  Gain from predation 
                           -(aTDapR*fT*imDapR*x[11]) // Loss from Respiration 
                            -(fJ*aJMes*imMes*2*x[13] * Q11/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;
 /*  Asp     */ dxdt[12] =  fJ*aJAsp*imAsp*x[12] *(0+Q1+Q2+Q3+Q4+Q5+Q8+Q9)/(1 + d*x[12] +(0+Q1+Q2+Q3+Q4+Q5+Q8+Q9))//  Gain from predation 
                           -(aTAsp*fT*imAsp*x[12]) // Loss from Respiration 
                            -(fJ*aJMes*imMes*2*x[13] * Q12/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;
 /*  Mes     */ dxdt[13] =  fJ*aJMes*imMes*x[13] *(0+Q8+Q9+Q10+Q11+Q12+Q16)/(1 + d*x[13] +(0+Q8+Q9+Q10+Q11+Q12+Q16))//  Gain from predation 
                           -(aTMes*fT*imMes*x[13]) // Loss from Respiration 
                            ;
 /*  DapP    */ dxdt[14] =  fJ*aJDapP*imDapP*x[14] *(0+Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7)/(1 + d*x[14] +(0+Q0+Q1+Q2+Q3+Q4+Q5+Q6+Q7))//  Gain from predation 
                           -(aTDapP*fT*imDapP*x[14]) // Loss from Respiration 
                            ;
 /*  Hol     */ dxdt[15] =  fJ*aJHol*imHol*x[15] *(0+Q0+Q1+Q2+Q3+Q4+Q5)/(1 + d*x[15] +(0+Q0+Q1+Q2+Q3+Q4+Q5))//  Gain from predation 
                           -(aTHol*fT*imHol*x[15]) // Loss from Respiration 
                            ;
 /*  Bos     */ dxdt[16] =  fJ*aJBos*imBos*x[16] *(0+Q0+Q2+Q3)/(1 + d*x[16] +(0+Q0+Q2+Q3))//  Gain from predation 
                           -(aTBos*fT*imBos*x[16]) // Loss from Respiration 
                            -(fJ*aJMes*imMes*2*x[13] * Q16/(1+d*x[13]+(Q8+Q9+Q10+Q11+Q12+Q16+0))) // Impact of Mes
                           ;

 } 
 void read_data()
{table_read("mydata1","./workspace/Pauldata91.txt",17);
table_addcolumn("mydata1", "foundBac"); 
table_addcolumn("mydata1", "foundDia"); 
table_addcolumn("mydata1", "foundCr"); 
table_addcolumn("mydata1", "foundSGr"); 
table_addcolumn("mydata1", "foundLGr"); 
table_addcolumn("mydata1", "foundSGo"); 
table_addcolumn("mydata1", "foundLGo"); 
table_addcolumn("mydata1", "foundInC"); 
table_addcolumn("mydata1", "foundSRot"); 
table_addcolumn("mydata1", "foundLRot"); 
table_addcolumn("mydata1", "foundOCop"); 
table_addcolumn("mydata1", "foundDapR"); 
table_addcolumn("mydata1", "foundAsp"); 
table_addcolumn("mydata1", "foundMes"); 
table_addcolumn("mydata1", "foundDapP"); 
table_addcolumn("mydata1", "foundHol"); 
table_addcolumn("mydata1", "foundBos"); 
table_read("mydata2","./workspace/Pauldata92.txt",17);
table_addcolumn("mydata2", "foundBac"); 
table_addcolumn("mydata2", "foundDia"); 
table_addcolumn("mydata2", "foundCr"); 
table_addcolumn("mydata2", "foundSGr"); 
table_addcolumn("mydata2", "foundLGr"); 
table_addcolumn("mydata2", "foundSGo"); 
table_addcolumn("mydata2", "foundLGo"); 
table_addcolumn("mydata2", "foundInC"); 
table_addcolumn("mydata2", "foundSRot"); 
table_addcolumn("mydata2", "foundLRot"); 
table_addcolumn("mydata2", "foundOCop"); 
table_addcolumn("mydata2", "foundDapR"); 
table_addcolumn("mydata2", "foundAsp"); 
table_addcolumn("mydata2", "foundMes"); 
table_addcolumn("mydata2", "foundDapP"); 
table_addcolumn("mydata2", "foundHol"); 
table_addcolumn("mydata2", "foundBos"); 
table_read("mydata3","./workspace/Pauldata93alt.txt",17);
table_addcolumn("mydata3", "foundBac"); 
table_addcolumn("mydata3", "foundDia"); 
table_addcolumn("mydata3", "foundCr"); 
table_addcolumn("mydata3", "foundSGr"); 
table_addcolumn("mydata3", "foundLGr"); 
table_addcolumn("mydata3", "foundSGo"); 
table_addcolumn("mydata3", "foundLGo"); 
table_addcolumn("mydata3", "foundInC"); 
table_addcolumn("mydata3", "foundSRot"); 
table_addcolumn("mydata3", "foundLRot"); 
table_addcolumn("mydata3", "foundOCop"); 
table_addcolumn("mydata3", "foundDapR"); 
table_addcolumn("mydata3", "foundAsp"); 
table_addcolumn("mydata3", "foundMes"); 
table_addcolumn("mydata3", "foundDapP"); 
table_addcolumn("mydata3", "foundHol"); 
table_addcolumn("mydata3", "foundBos"); 
table_read("mydata4","./workspace/Pauldata94.txt",17);
table_addcolumn("mydata4", "foundBac"); 
table_addcolumn("mydata4", "foundDia"); 
table_addcolumn("mydata4", "foundCr"); 
table_addcolumn("mydata4", "foundSGr"); 
table_addcolumn("mydata4", "foundLGr"); 
table_addcolumn("mydata4", "foundSGo"); 
table_addcolumn("mydata4", "foundLGo"); 
table_addcolumn("mydata4", "foundInC"); 
table_addcolumn("mydata4", "foundSRot"); 
table_addcolumn("mydata4", "foundLRot"); 
table_addcolumn("mydata4", "foundOCop"); 
table_addcolumn("mydata4", "foundDapR"); 
table_addcolumn("mydata4", "foundAsp"); 
table_addcolumn("mydata4", "foundMes"); 
table_addcolumn("mydata4", "foundDapP"); 
table_addcolumn("mydata4", "foundHol"); 
table_addcolumn("mydata4", "foundBos"); 
table_read("mydata5","./workspace/Pauldata95.txt",17);
table_addcolumn("mydata5", "foundBac"); 
table_addcolumn("mydata5", "foundDia"); 
table_addcolumn("mydata5", "foundCr"); 
table_addcolumn("mydata5", "foundSGr"); 
table_addcolumn("mydata5", "foundLGr"); 
table_addcolumn("mydata5", "foundSGo"); 
table_addcolumn("mydata5", "foundLGo"); 
table_addcolumn("mydata5", "foundInC"); 
table_addcolumn("mydata5", "foundSRot"); 
table_addcolumn("mydata5", "foundLRot"); 
table_addcolumn("mydata5", "foundOCop"); 
table_addcolumn("mydata5", "foundDapR"); 
table_addcolumn("mydata5", "foundAsp"); 
table_addcolumn("mydata5", "foundMes"); 
table_addcolumn("mydata5", "foundDapP"); 
table_addcolumn("mydata5", "foundHol"); 
table_addcolumn("mydata5", "foundBos"); 
return;} 
 
void setup_parameters()
{   /* min - max - initial - norm/lnorm -  fixed? - display?  */ 
parameter_create("aRBac",0.01,20,8.146,1,0,1);
parameter_create("aRDia",0.01,20,15.614,1,0,1); 
parameter_create("aRCr",0.01,20,19.259,1,0,1);
parameter_create("aRSGr",0.01,20,18.927,1,0,1);
parameter_create("aRLGr",0.01,20,2.225,1,0,1);
parameter_create("aRSGo",0.01,20,15.875,1,0,1);
parameter_create("aRLGo",0.01,20,6.98,1,0,1);
parameter_create("aRInC",0.01,20,6.181,1,0,1);
parameter_create("aJSRot",0.01,40,16.79,1,0,1);
parameter_create("aJLRot",0.01,40,14.107,1,0,1);
parameter_create("aJOCop",0.01,40,2.844,1,0,1);
parameter_create("aJDapR",0.01,40,4.103,1,0,1);
parameter_create("aJAsp",0.01,40,5.631,1,0,1);
parameter_create("aJMes",0.01,40,2.647,1,0,1);
parameter_create("aJDapP",0.01,40,32.828,1,0,1);
parameter_create("aJHol",0.01,40,36.797,1,0,1);
parameter_create("aJBos",0.01,40,6.746,1,0,1);
parameter_create("aTSRot",0.00001,2,0.22,1,0,1);
parameter_create("aTLRot",0.00001,2,0.351,1,0,1);
parameter_create("aTOCop",0.00001,2,0.065,1,0,1);
parameter_create("aTDapR",0.00001,2,1.331,1,0,1);
parameter_create("aTAsp",0.00001,2,1.192,1,0,1);
parameter_create("aTMes",0.00001,2,0.308,1,0,1);
parameter_create("aTDapP",0.00001,2,0.842,1,0,1);
parameter_create("aTHol",0.000001,2,1.35,1,0,1);
parameter_create("aTBos",0.000001,2,1.575,1,0,1);

parameter_create("fT",0.05,0.07,0.06,0,1,1);
parameter_create("fR",0.12,0.14,0.13,0,1,1);
parameter_create("K",0.1,10,0.3,0,1,1);
parameter_create("sdZ",1.5,8,1.537505,0,-1,1);
parameter_create("sdP",1,6,1.467582,0,-1,1);
parameter_create("fJ",0.14,0.16,0.14,0,1,1);
parameter_create("d",1e-05,100,10,0,1,1);
parameter_create("W",0.1,5,0.15,0,1,1);
parameter_create("q",1,2.5,1.4,0,1,1);

parameter_delay("sdZ", 8000);
parameter_delay("sdP", 8000);

parameter_addprior("aRBac",9,0.145);
parameter_addprior("aRDia",9,0.145);
parameter_addprior("aRCr",9,0.145);
parameter_addprior("aRSGr",9,0.145);
parameter_addprior("aRLGr",9,0.145);
parameter_addprior("aRSGo",9,0.145);
parameter_addprior("aRLGo",9,0.145);
parameter_addprior("aRInC",9,0.145);

parameter_addprior("aJSRot",10,0.125);
parameter_addprior("aJLRot",10,0.125);
parameter_addprior("aJOCop",10,0.125);
parameter_addprior("aJDapR",10,0.125);
parameter_addprior("aJAsp",10,0.125);
parameter_addprior("aJMes",10,0.125);
parameter_addprior("aJDapP",10,0.125);
parameter_addprior("aJHol",10,0.125);
parameter_addprior("aJBos",10,0.125);

parameter_addprior("aTSRot",0.110,1.6);
parameter_addprior("aTLRot",0.110,1.6);
parameter_addprior("aTOCop",0.110,1.6);
parameter_addprior("aTDapR",0.110,1.6);
parameter_addprior("aTAsp",0.110,1.6);
parameter_addprior("aTMes",0.110,1.6);
parameter_addprior("aTDapP",0.110,1.6);
parameter_addprior("aTHol",0.110,1.6);
parameter_addprior("aTBos",0.110,1.6);


parameter_showall();
                PAUSE
                return;
}  
 void write_Hy1( const state_type &x , const double t )
{
table_writevalue_multichain("mydata1", "foundBac", timerow,x[0]); 
table_writevalue_multichain("mydata1", "foundDia", timerow,x[1]); 
table_writevalue_multichain("mydata1", "foundCr", timerow,x[2]); 
table_writevalue_multichain("mydata1", "foundSGr", timerow,x[3]); 
table_writevalue_multichain("mydata1", "foundLGr", timerow,x[4]); 
table_writevalue_multichain("mydata1", "foundSGo", timerow,x[5]); 
table_writevalue_multichain("mydata1", "foundLGo", timerow,x[6]); 
table_writevalue_multichain("mydata1", "foundInC", timerow,x[7]); 
table_writevalue_multichain("mydata1", "foundSRot", timerow,x[8]); 
table_writevalue_multichain("mydata1", "foundLRot", timerow,x[9]); 
table_writevalue_multichain("mydata1", "foundOCop", timerow,x[10]); 
table_writevalue_multichain("mydata1", "foundDapR", timerow,x[11]); 
table_writevalue_multichain("mydata1", "foundAsp", timerow,x[12]); 
table_writevalue_multichain("mydata1", "foundMes", timerow,x[13]); 
table_writevalue_multichain("mydata1", "foundDapP", timerow,x[14]); 
table_writevalue_multichain("mydata1", "foundHol", timerow,x[15]); 
table_writevalue_multichain("mydata1", "foundBos", timerow,x[16]); 
timerow++; 
} 
void write_Hy2( const state_type &x , const double t )
{
table_writevalue_multichain("mydata2", "foundBac", timerow,x[0]); 
table_writevalue_multichain("mydata2", "foundDia", timerow,x[1]); 
table_writevalue_multichain("mydata2", "foundCr", timerow,x[2]); 
table_writevalue_multichain("mydata2", "foundSGr", timerow,x[3]); 
table_writevalue_multichain("mydata2", "foundLGr", timerow,x[4]); 
table_writevalue_multichain("mydata2", "foundSGo", timerow,x[5]); 
table_writevalue_multichain("mydata2", "foundLGo", timerow,x[6]); 
table_writevalue_multichain("mydata2", "foundInC", timerow,x[7]); 
table_writevalue_multichain("mydata2", "foundSRot", timerow,x[8]); 
table_writevalue_multichain("mydata2", "foundLRot", timerow,x[9]); 
table_writevalue_multichain("mydata2", "foundOCop", timerow,x[10]); 
table_writevalue_multichain("mydata2", "foundDapR", timerow,x[11]); 
table_writevalue_multichain("mydata2", "foundAsp", timerow,x[12]); 
table_writevalue_multichain("mydata2", "foundMes", timerow,x[13]); 
table_writevalue_multichain("mydata2", "foundDapP", timerow,x[14]); 
table_writevalue_multichain("mydata2", "foundHol", timerow,x[15]); 
table_writevalue_multichain("mydata2", "foundBos", timerow,x[16]); 
timerow++; 
} 
void write_Hy3( const state_type &x , const double t )
{
table_writevalue_multichain("mydata3", "foundBac", timerow,x[0]); 
table_writevalue_multichain("mydata3", "foundDia", timerow,x[1]); 
table_writevalue_multichain("mydata3", "foundCr", timerow,x[2]); 
table_writevalue_multichain("mydata3", "foundSGr", timerow,x[3]); 
table_writevalue_multichain("mydata3", "foundLGr", timerow,x[4]); 
table_writevalue_multichain("mydata3", "foundSGo", timerow,x[5]); 
table_writevalue_multichain("mydata3", "foundLGo", timerow,x[6]); 
table_writevalue_multichain("mydata3", "foundInC", timerow,x[7]); 
table_writevalue_multichain("mydata3", "foundSRot", timerow,x[8]); 
table_writevalue_multichain("mydata3", "foundLRot", timerow,x[9]); 
table_writevalue_multichain("mydata3", "foundOCop", timerow,x[10]); 
table_writevalue_multichain("mydata3", "foundDapR", timerow,x[11]); 
table_writevalue_multichain("mydata3", "foundAsp", timerow,x[12]); 
table_writevalue_multichain("mydata3", "foundMes", timerow,x[13]); 
table_writevalue_multichain("mydata3", "foundDapP", timerow,x[14]); 
table_writevalue_multichain("mydata3", "foundHol", timerow,x[15]); 
table_writevalue_multichain("mydata3", "foundBos", timerow,x[16]); 
timerow++; 
} 
void write_Hy4( const state_type &x , const double t )
{
table_writevalue_multichain("mydata4", "foundBac", timerow,x[0]); 
table_writevalue_multichain("mydata4", "foundDia", timerow,x[1]); 
table_writevalue_multichain("mydata4", "foundCr", timerow,x[2]); 
table_writevalue_multichain("mydata4", "foundSGr", timerow,x[3]); 
table_writevalue_multichain("mydata4", "foundLGr", timerow,x[4]); 
table_writevalue_multichain("mydata4", "foundSGo", timerow,x[5]); 
table_writevalue_multichain("mydata4", "foundLGo", timerow,x[6]); 
table_writevalue_multichain("mydata4", "foundInC", timerow,x[7]); 
table_writevalue_multichain("mydata4", "foundSRot", timerow,x[8]); 
table_writevalue_multichain("mydata4", "foundLRot", timerow,x[9]); 
table_writevalue_multichain("mydata4", "foundOCop", timerow,x[10]); 
table_writevalue_multichain("mydata4", "foundDapR", timerow,x[11]); 
table_writevalue_multichain("mydata4", "foundAsp", timerow,x[12]); 
table_writevalue_multichain("mydata4", "foundMes", timerow,x[13]); 
table_writevalue_multichain("mydata4", "foundDapP", timerow,x[14]); 
table_writevalue_multichain("mydata4", "foundHol", timerow,x[15]); 
table_writevalue_multichain("mydata4", "foundBos", timerow,x[16]); 
timerow++; 
} 
void write_Hy5( const state_type &x , const double t )
{
table_writevalue_multichain("mydata5", "foundBac", timerow,x[0]); 
table_writevalue_multichain("mydata5", "foundDia", timerow,x[1]); 
table_writevalue_multichain("mydata5", "foundCr", timerow,x[2]); 
table_writevalue_multichain("mydata5", "foundSGr", timerow,x[3]); 
table_writevalue_multichain("mydata5", "foundLGr", timerow,x[4]); 
table_writevalue_multichain("mydata5", "foundSGo", timerow,x[5]); 
table_writevalue_multichain("mydata5", "foundLGo", timerow,x[6]); 
table_writevalue_multichain("mydata5", "foundInC", timerow,x[7]); 
table_writevalue_multichain("mydata5", "foundSRot", timerow,x[8]); 
table_writevalue_multichain("mydata5", "foundLRot", timerow,x[9]); 
table_writevalue_multichain("mydata5", "foundOCop", timerow,x[10]); 
table_writevalue_multichain("mydata5", "foundDapR", timerow,x[11]); 
table_writevalue_multichain("mydata5", "foundAsp", timerow,x[12]); 
table_writevalue_multichain("mydata5", "foundMes", timerow,x[13]); 
table_writevalue_multichain("mydata5", "foundDapP", timerow,x[14]); 
table_writevalue_multichain("mydata5", "foundHol", timerow,x[15]); 
table_writevalue_multichain("mydata5", "foundBos", timerow,x[16]); 
timerow++; 
} 
double probBac; double realBac; double foundBac; 
double probDia; double realDia; double foundDia; 
double probCr; double realCr; double foundCr; 
double probSGr; double realSGr; double foundSGr; 
double probLGr; double realLGr; double foundLGr; 
double probSGo; double realSGo; double foundSGo; 
double probLGo; double realLGo; double foundLGo; 
double probInC; double realInC; double foundInC; 
double probSRot; double realSRot; double foundSRot; 
double probLRot; double realLRot; double foundLRot; 
double probOCop; double realOCop; double foundOCop; 
double probDapR; double realDapR; double foundDapR; 
double probAsp; double realAsp; double foundAsp; 
double probMes; double realMes; double foundMes; 
double probDapP; double realDapP; double foundDapP; 
double probHol; double realHol; double foundHol; 
double probBos; double realBos; double foundBos; 
void likelihood()
{
set_metr_ltotnew(0.0); /* set sum over log-Likelihood to zero */
set_metr_number_ok(0);
int numdata;
timerow = 0;
aRBac=cv("aRBac");
aRDia=cv("aRDia");
aRCr=cv("aRCr");
aRSGr=cv("aRSGr");
aRLGr=cv("aRLGr");
aRSGo=cv("aRSGo");
aRLGo=cv("aRLGo");
aRInC=cv("aRInC");
aTSRot=cv("aTSRot");
aTLRot=cv("aTLRot");
aTOCop=cv("aTOCop");
aTDapR=cv("aTDapR");
aTAsp=cv("aTAsp");
aTMes=cv("aTMes");
aTDapP=cv("aTDapP");
aTHol=cv("aTHol");
aTBos=cv("aTBos");
aJSRot=cv("aJSRot");
aJLRot=cv("aJLRot");
aJOCop=cv("aJOCop");
aJDapR=cv("aJDapR");
aJAsp=cv("aJAsp");
aJMes=cv("aJMes");
aJDapP=cv("aJDapP");
aJHol=cv("aJHol");
aJBos=cv("aJBos");

fR=cv("fR");
fT=cv("fT");
K=cv("K");
sdZ=cv("sdZ");
sdP=cv("sdP");
fJ=cv("fJ");
d=cv("d");
W=cv("W");
q=cv("q");
 
  timerow = 0;
 
                state_type x1 = {start1Bac,start1Dia,start1Cr,start1SGr,start1LGr,start1SGo,start1LGo,start1InC,start1SRot,start1LRot,start1OCop,start1DapR,start1Asp,start1Mes,start1DapP,start1Hol,start1Bos}; 
integrate_const(controlled_stepper,ODE,x1,0.0,1.4,0.1,write_Hy1 ); 
numdata = 14;
                 for(int ii = 0; ii < numdata; ii++)
{ realBac = table_getvalue("mydata1", "Bac",ii);
realDia = table_getvalue("mydata1", "Dia",ii);
realCr = table_getvalue("mydata1", "Cr",ii);
realSGr = table_getvalue("mydata1", "SGr",ii);
realLGr = table_getvalue("mydata1", "LGr",ii);
realSGo = table_getvalue("mydata1", "SGo",ii);
realLGo = table_getvalue("mydata1", "LGo",ii);
realInC = table_getvalue("mydata1", "InC",ii);
realSRot = table_getvalue("mydata1", "SRot",ii);
realLRot = table_getvalue("mydata1", "LRot",ii);
realOCop = table_getvalue("mydata1", "OCop",ii);
realDapR = table_getvalue("mydata1", "DapR",ii);
realAsp = table_getvalue("mydata1", "Asp",ii);
realMes = table_getvalue("mydata1", "Mes",ii);
realDapP = table_getvalue("mydata1", "DapP",ii);
realHol = table_getvalue("mydata1", "Hol",ii);
realBos = table_getvalue("mydata1", "Bos",ii);
foundBac = table_getvalue("mydata1", "foundBac",ii);
foundDia = table_getvalue("mydata1", "foundDia",ii);
foundCr = table_getvalue("mydata1", "foundCr",ii);
foundSGr = table_getvalue("mydata1", "foundSGr",ii);
foundLGr = table_getvalue("mydata1", "foundLGr",ii);
foundSGo = table_getvalue("mydata1", "foundSGo",ii);
foundLGo = table_getvalue("mydata1", "foundLGo",ii);
foundInC = table_getvalue("mydata1", "foundInC",ii);
foundSRot = table_getvalue("mydata1", "foundSRot",ii);
foundLRot = table_getvalue("mydata1", "foundLRot",ii);
foundOCop = table_getvalue("mydata1", "foundOCop",ii);
foundDapR = table_getvalue("mydata1", "foundDapR",ii);
foundAsp = table_getvalue("mydata1", "foundAsp",ii);
foundMes = table_getvalue("mydata1", "foundMes",ii);
foundDapP = table_getvalue("mydata1", "foundDapP",ii);
foundHol = table_getvalue("mydata1", "foundHol",ii);
foundBos = table_getvalue("mydata1", "foundBos",ii);
inc_metr_ltotnew(log(probBac=lognormal_density(realBac, foundBac,sdP)));
inc_metr_ltotnew(log(probDia=lognormal_density(realDia, foundDia,sdP)));
inc_metr_ltotnew(log(probCr=lognormal_density(realCr, foundCr,sdP)));
inc_metr_ltotnew(log(probSGr=lognormal_density(realSGr, foundSGr,sdP)));
inc_metr_ltotnew(log(probLGr=lognormal_density(realLGr, foundLGr,sdP)));
inc_metr_ltotnew(log(probSGo=lognormal_density(realSGo, foundSGo,sdP)));
inc_metr_ltotnew(log(probLGo=lognormal_density(realLGo, foundLGo,sdP)));
inc_metr_ltotnew(log(probInC=lognormal_density(realInC, foundInC,sdP)));
inc_metr_ltotnew(log(probSRot=lognormal_density(realSRot, foundSRot,sdZ)));
inc_metr_ltotnew(log(probLRot=lognormal_density(realLRot, foundLRot,sdZ)));
inc_metr_ltotnew(log(probOCop=lognormal_density(realOCop, foundOCop,sdZ)));
inc_metr_ltotnew(log(probDapR=lognormal_density(realDapR, foundDapR,sdZ)));
inc_metr_ltotnew(log(probAsp=lognormal_density(realAsp, foundAsp,sdZ)));
inc_metr_ltotnew(log(probMes=lognormal_density(realMes, foundMes,sdZ)));
inc_metr_ltotnew(log(probDapP=lognormal_density(realDapP, foundDapP,sdZ)));
inc_metr_ltotnew(log(probHol=lognormal_density(realHol, foundHol,sdZ)));
inc_metr_ltotnew(log(probBos=lognormal_density(realBos, foundBos,sdZ)));
inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on 
 } 
  timerow = 0;
 
                state_type x2 = {start2Bac,start2Dia,start2Cr,start2SGr,start2LGr,start2SGo,start2LGo,start2InC,start2SRot,start2LRot,start2OCop,start2DapR,start2Asp,start2Mes,start2DapP,start2Hol,start2Bos}; 
integrate_const(controlled_stepper,ODE,x2,0.0,1.4,0.1,write_Hy2 ); 
numdata = 14;
                 for(int ii = 0; ii < numdata; ii++)
{ realBac = table_getvalue("mydata2", "Bac",ii);
realDia = table_getvalue("mydata2", "Dia",ii);
realCr = table_getvalue("mydata2", "Cr",ii);
realSGr = table_getvalue("mydata2", "SGr",ii);
realLGr = table_getvalue("mydata2", "LGr",ii);
realSGo = table_getvalue("mydata2", "SGo",ii);
realLGo = table_getvalue("mydata2", "LGo",ii);
realInC = table_getvalue("mydata2", "InC",ii);
realSRot = table_getvalue("mydata2", "SRot",ii);
realLRot = table_getvalue("mydata2", "LRot",ii);
realOCop = table_getvalue("mydata2", "OCop",ii);
realDapR = table_getvalue("mydata2", "DapR",ii);
realAsp = table_getvalue("mydata2", "Asp",ii);
realMes = table_getvalue("mydata2", "Mes",ii);
realDapP = table_getvalue("mydata2", "DapP",ii);
realHol = table_getvalue("mydata2", "Hol",ii);
realBos = table_getvalue("mydata2", "Bos",ii);
foundBac = table_getvalue("mydata2", "foundBac",ii);
foundDia = table_getvalue("mydata2", "foundDia",ii);
foundCr = table_getvalue("mydata2", "foundCr",ii);
foundSGr = table_getvalue("mydata2", "foundSGr",ii);
foundLGr = table_getvalue("mydata2", "foundLGr",ii);
foundSGo = table_getvalue("mydata2", "foundSGo",ii);
foundLGo = table_getvalue("mydata2", "foundLGo",ii);
foundInC = table_getvalue("mydata2", "foundInC",ii);
foundSRot = table_getvalue("mydata2", "foundSRot",ii);
foundLRot = table_getvalue("mydata2", "foundLRot",ii);
foundOCop = table_getvalue("mydata2", "foundOCop",ii);
foundDapR = table_getvalue("mydata2", "foundDapR",ii);
foundAsp = table_getvalue("mydata2", "foundAsp",ii);
foundMes = table_getvalue("mydata2", "foundMes",ii);
foundDapP = table_getvalue("mydata2", "foundDapP",ii);
foundHol = table_getvalue("mydata2", "foundHol",ii);
foundBos = table_getvalue("mydata2", "foundBos",ii);
inc_metr_ltotnew(log(probBac=lognormal_density(realBac, foundBac,sdP)));
inc_metr_ltotnew(log(probDia=lognormal_density(realDia, foundDia,sdP)));
inc_metr_ltotnew(log(probCr=lognormal_density(realCr, foundCr,sdP)));
inc_metr_ltotnew(log(probSGr=lognormal_density(realSGr, foundSGr,sdP)));
inc_metr_ltotnew(log(probLGr=lognormal_density(realLGr, foundLGr,sdP)));
inc_metr_ltotnew(log(probSGo=lognormal_density(realSGo, foundSGo,sdP)));
inc_metr_ltotnew(log(probLGo=lognormal_density(realLGo, foundLGo,sdP)));
inc_metr_ltotnew(log(probInC=lognormal_density(realInC, foundInC,sdP)));
inc_metr_ltotnew(log(probSRot=lognormal_density(realSRot, foundSRot,sdZ)));
inc_metr_ltotnew(log(probLRot=lognormal_density(realLRot, foundLRot,sdZ)));
inc_metr_ltotnew(log(probOCop=lognormal_density(realOCop, foundOCop,sdZ)));
inc_metr_ltotnew(log(probDapR=lognormal_density(realDapR, foundDapR,sdZ)));
inc_metr_ltotnew(log(probAsp=lognormal_density(realAsp, foundAsp,sdZ)));
inc_metr_ltotnew(log(probMes=lognormal_density(realMes, foundMes,sdZ)));
inc_metr_ltotnew(log(probDapP=lognormal_density(realDapP, foundDapP,sdZ)));
inc_metr_ltotnew(log(probHol=lognormal_density(realHol, foundHol,sdZ)));
inc_metr_ltotnew(log(probBos=lognormal_density(realBos, foundBos,sdZ)));
inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on 
 } 
  timerow = 0;
 
                state_type x3 = {start3Bac,start3Dia,start3Cr,start3SGr,start3LGr,start3SGo,start3LGo,start3InC,start3SRot,start3LRot,start3OCop,start3DapR,start3Asp,start3Mes,start3DapP,start3Hol,start3Bos}; 
integrate_const(controlled_stepper,ODE,x3,0.0,1.4,0.1,write_Hy3 ); 
numdata = 14;
                 for(int ii = 0; ii < numdata; ii++)
{ realBac = table_getvalue("mydata3", "Bac",ii);
realDia = table_getvalue("mydata3", "Dia",ii);
realCr = table_getvalue("mydata3", "Cr",ii);
realSGr = table_getvalue("mydata3", "SGr",ii);
realLGr = table_getvalue("mydata3", "LGr",ii);
realSGo = table_getvalue("mydata3", "SGo",ii);
realLGo = table_getvalue("mydata3", "LGo",ii);
realInC = table_getvalue("mydata3", "InC",ii);
realSRot = table_getvalue("mydata3", "SRot",ii);
realLRot = table_getvalue("mydata3", "LRot",ii);
realOCop = table_getvalue("mydata3", "OCop",ii);
realDapR = table_getvalue("mydata3", "DapR",ii);
realAsp = table_getvalue("mydata3", "Asp",ii);
realMes = table_getvalue("mydata3", "Mes",ii);
realDapP = table_getvalue("mydata3", "DapP",ii);
realHol = table_getvalue("mydata3", "Hol",ii);
realBos = table_getvalue("mydata3", "Bos",ii);
foundBac = table_getvalue("mydata3", "foundBac",ii);
foundDia = table_getvalue("mydata3", "foundDia",ii);
foundCr = table_getvalue("mydata3", "foundCr",ii);
foundSGr = table_getvalue("mydata3", "foundSGr",ii);
foundLGr = table_getvalue("mydata3", "foundLGr",ii);
foundSGo = table_getvalue("mydata3", "foundSGo",ii);
foundLGo = table_getvalue("mydata3", "foundLGo",ii);
foundInC = table_getvalue("mydata3", "foundInC",ii);
foundSRot = table_getvalue("mydata3", "foundSRot",ii);
foundLRot = table_getvalue("mydata3", "foundLRot",ii);
foundOCop = table_getvalue("mydata3", "foundOCop",ii);
foundDapR = table_getvalue("mydata3", "foundDapR",ii);
foundAsp = table_getvalue("mydata3", "foundAsp",ii);
foundMes = table_getvalue("mydata3", "foundMes",ii);
foundDapP = table_getvalue("mydata3", "foundDapP",ii);
foundHol = table_getvalue("mydata3", "foundHol",ii);
foundBos = table_getvalue("mydata3", "foundBos",ii);
inc_metr_ltotnew(log(probBac=lognormal_density(realBac, foundBac,sdP)));
inc_metr_ltotnew(log(probDia=lognormal_density(realDia, foundDia,sdP)));
inc_metr_ltotnew(log(probCr=lognormal_density(realCr, foundCr,sdP)));
inc_metr_ltotnew(log(probSGr=lognormal_density(realSGr, foundSGr,sdP)));
inc_metr_ltotnew(log(probLGr=lognormal_density(realLGr, foundLGr,sdP)));
inc_metr_ltotnew(log(probSGo=lognormal_density(realSGo, foundSGo,sdP)));
inc_metr_ltotnew(log(probLGo=lognormal_density(realLGo, foundLGo,sdP)));
inc_metr_ltotnew(log(probInC=lognormal_density(realInC, foundInC,sdP)));
inc_metr_ltotnew(log(probSRot=lognormal_density(realSRot, foundSRot,sdZ)));
inc_metr_ltotnew(log(probLRot=lognormal_density(realLRot, foundLRot,sdZ)));
inc_metr_ltotnew(log(probOCop=lognormal_density(realOCop, foundOCop,sdZ)));
inc_metr_ltotnew(log(probDapR=lognormal_density(realDapR, foundDapR,sdZ)));
inc_metr_ltotnew(log(probAsp=lognormal_density(realAsp, foundAsp,sdZ)));
inc_metr_ltotnew(log(probMes=lognormal_density(realMes, foundMes,sdZ)));
inc_metr_ltotnew(log(probDapP=lognormal_density(realDapP, foundDapP,sdZ)));
inc_metr_ltotnew(log(probHol=lognormal_density(realHol, foundHol,sdZ)));
inc_metr_ltotnew(log(probBos=lognormal_density(realBos, foundBos,sdZ)));
inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on 
 } 
  timerow = 0;
 
                state_type x4 = {start4Bac,start4Dia,start4Cr,start4SGr,start4LGr,start4SGo,start4LGo,start4InC,start4SRot,start4LRot,start4OCop,start4DapR,start4Asp,start4Mes,start4DapP,start4Hol,start4Bos}; 
integrate_const(controlled_stepper,ODE,x4,0.0,1.4,0.1,write_Hy4 ); 
numdata = 14;
                 for(int ii = 0; ii < numdata; ii++)
{ realBac = table_getvalue("mydata4", "Bac",ii);
realDia = table_getvalue("mydata4", "Dia",ii);
realCr = table_getvalue("mydata4", "Cr",ii);
realSGr = table_getvalue("mydata4", "SGr",ii);
realLGr = table_getvalue("mydata4", "LGr",ii);
realSGo = table_getvalue("mydata4", "SGo",ii);
realLGo = table_getvalue("mydata4", "LGo",ii);
realInC = table_getvalue("mydata4", "InC",ii);
realSRot = table_getvalue("mydata4", "SRot",ii);
realLRot = table_getvalue("mydata4", "LRot",ii);
realOCop = table_getvalue("mydata4", "OCop",ii);
realDapR = table_getvalue("mydata4", "DapR",ii);
realAsp = table_getvalue("mydata4", "Asp",ii);
realMes = table_getvalue("mydata4", "Mes",ii);
realDapP = table_getvalue("mydata4", "DapP",ii);
realHol = table_getvalue("mydata4", "Hol",ii);
realBos = table_getvalue("mydata4", "Bos",ii);
foundBac = table_getvalue("mydata4", "foundBac",ii);
foundDia = table_getvalue("mydata4", "foundDia",ii);
foundCr = table_getvalue("mydata4", "foundCr",ii);
foundSGr = table_getvalue("mydata4", "foundSGr",ii);
foundLGr = table_getvalue("mydata4", "foundLGr",ii);
foundSGo = table_getvalue("mydata4", "foundSGo",ii);
foundLGo = table_getvalue("mydata4", "foundLGo",ii);
foundInC = table_getvalue("mydata4", "foundInC",ii);
foundSRot = table_getvalue("mydata4", "foundSRot",ii);
foundLRot = table_getvalue("mydata4", "foundLRot",ii);
foundOCop = table_getvalue("mydata4", "foundOCop",ii);
foundDapR = table_getvalue("mydata4", "foundDapR",ii);
foundAsp = table_getvalue("mydata4", "foundAsp",ii);
foundMes = table_getvalue("mydata4", "foundMes",ii);
foundDapP = table_getvalue("mydata4", "foundDapP",ii);
foundHol = table_getvalue("mydata4", "foundHol",ii);
foundBos = table_getvalue("mydata4", "foundBos",ii);
inc_metr_ltotnew(log(probBac=lognormal_density(realBac, foundBac,sdP)));
inc_metr_ltotnew(log(probDia=lognormal_density(realDia, foundDia,sdP)));
inc_metr_ltotnew(log(probCr=lognormal_density(realCr, foundCr,sdP)));
inc_metr_ltotnew(log(probSGr=lognormal_density(realSGr, foundSGr,sdP)));
inc_metr_ltotnew(log(probLGr=lognormal_density(realLGr, foundLGr,sdP)));
inc_metr_ltotnew(log(probSGo=lognormal_density(realSGo, foundSGo,sdP)));
inc_metr_ltotnew(log(probLGo=lognormal_density(realLGo, foundLGo,sdP)));
inc_metr_ltotnew(log(probInC=lognormal_density(realInC, foundInC,sdP)));
inc_metr_ltotnew(log(probSRot=lognormal_density(realSRot, foundSRot,sdZ)));
inc_metr_ltotnew(log(probLRot=lognormal_density(realLRot, foundLRot,sdZ)));
inc_metr_ltotnew(log(probOCop=lognormal_density(realOCop, foundOCop,sdZ)));
inc_metr_ltotnew(log(probDapR=lognormal_density(realDapR, foundDapR,sdZ)));
inc_metr_ltotnew(log(probAsp=lognormal_density(realAsp, foundAsp,sdZ)));
inc_metr_ltotnew(log(probMes=lognormal_density(realMes, foundMes,sdZ)));
inc_metr_ltotnew(log(probDapP=lognormal_density(realDapP, foundDapP,sdZ)));
inc_metr_ltotnew(log(probHol=lognormal_density(realHol, foundHol,sdZ)));
inc_metr_ltotnew(log(probBos=lognormal_density(realBos, foundBos,sdZ)));
inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on 
 } 
  timerow = 0;
 
                state_type x5 = {start5Bac,start5Dia,start5Cr,start5SGr,start5LGr,start5SGo,start5LGo,start5InC,start5SRot,start5LRot,start5OCop,start5DapR,start5Asp,start5Mes,start5DapP,start5Hol,start5Bos}; 
integrate_const(controlled_stepper,ODE,x5,0.0,1.4,0.1,write_Hy5 ); 
numdata = 14;
                 for(int ii = 0; ii < numdata; ii++)
{ realBac = table_getvalue("mydata5", "Bac",ii);
realDia = table_getvalue("mydata5", "Dia",ii);
realCr = table_getvalue("mydata5", "Cr",ii);
realSGr = table_getvalue("mydata5", "SGr",ii);
realLGr = table_getvalue("mydata5", "LGr",ii);
realSGo = table_getvalue("mydata5", "SGo",ii);
realLGo = table_getvalue("mydata5", "LGo",ii);
realInC = table_getvalue("mydata5", "InC",ii);
realSRot = table_getvalue("mydata5", "SRot",ii);
realLRot = table_getvalue("mydata5", "LRot",ii);
realOCop = table_getvalue("mydata5", "OCop",ii);
realDapR = table_getvalue("mydata5", "DapR",ii);
realAsp = table_getvalue("mydata5", "Asp",ii);
realMes = table_getvalue("mydata5", "Mes",ii);
realDapP = table_getvalue("mydata5", "DapP",ii);
realHol = table_getvalue("mydata5", "Hol",ii);
realBos = table_getvalue("mydata5", "Bos",ii);
foundBac = table_getvalue("mydata5", "foundBac",ii);
foundDia = table_getvalue("mydata5", "foundDia",ii);
foundCr = table_getvalue("mydata5", "foundCr",ii);
foundSGr = table_getvalue("mydata5", "foundSGr",ii);
foundLGr = table_getvalue("mydata5", "foundLGr",ii);
foundSGo = table_getvalue("mydata5", "foundSGo",ii);
foundLGo = table_getvalue("mydata5", "foundLGo",ii);
foundInC = table_getvalue("mydata5", "foundInC",ii);
foundSRot = table_getvalue("mydata5", "foundSRot",ii);
foundLRot = table_getvalue("mydata5", "foundLRot",ii);
foundOCop = table_getvalue("mydata5", "foundOCop",ii);
foundDapR = table_getvalue("mydata5", "foundDapR",ii);
foundAsp = table_getvalue("mydata5", "foundAsp",ii);
foundMes = table_getvalue("mydata5", "foundMes",ii);
foundDapP = table_getvalue("mydata5", "foundDapP",ii);
foundHol = table_getvalue("mydata5", "foundHol",ii);
foundBos = table_getvalue("mydata5", "foundBos",ii);
inc_metr_ltotnew(log(probBac=lognormal_density(realBac, foundBac,sdP)));
inc_metr_ltotnew(log(probDia=lognormal_density(realDia, foundDia,sdP)));
inc_metr_ltotnew(log(probCr=lognormal_density(realCr, foundCr,sdP)));
inc_metr_ltotnew(log(probSGr=lognormal_density(realSGr, foundSGr,sdP)));
inc_metr_ltotnew(log(probLGr=lognormal_density(realLGr, foundLGr,sdP)));
inc_metr_ltotnew(log(probSGo=lognormal_density(realSGo, foundSGo,sdP)));
inc_metr_ltotnew(log(probLGo=lognormal_density(realLGo, foundLGo,sdP)));
inc_metr_ltotnew(log(probInC=lognormal_density(realInC, foundInC,sdP)));
inc_metr_ltotnew(log(probSRot=lognormal_density(realSRot, foundSRot,sdZ)));
inc_metr_ltotnew(log(probLRot=lognormal_density(realLRot, foundLRot,sdZ)));
inc_metr_ltotnew(log(probOCop=lognormal_density(realOCop, foundOCop,sdZ)));
inc_metr_ltotnew(log(probDapR=lognormal_density(realDapR, foundDapR,sdZ)));
inc_metr_ltotnew(log(probAsp=lognormal_density(realAsp, foundAsp,sdZ)));
inc_metr_ltotnew(log(probMes=lognormal_density(realMes, foundMes,sdZ)));
inc_metr_ltotnew(log(probDapP=lognormal_density(realDapP, foundDapP,sdZ)));
inc_metr_ltotnew(log(probHol=lognormal_density(realHol, foundHol,sdZ)));
inc_metr_ltotnew(log(probBos=lognormal_density(realBos, foundBos,sdZ)));
inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on 
 }
                   return;
} 
void final_output()
{ /* run the most likely set of parameters */
               params_set_to_posterior_mean();
               likelihood();char fname1[100];
                                 get_filzbach_path(fname1, 100);
                                 strcat_s(fname1,"_my_outputYEAR1.txt");
                                 table_output("mydata1",fname1);char fname2[100];
                                 get_filzbach_path(fname2, 100);
                                 strcat_s(fname2,"_my_outputYEAR2.txt");
                                 table_output("mydata2",fname2);char fname3[100];
                                 get_filzbach_path(fname3, 100);
                                 strcat_s(fname3,"_my_outputYEAR3.txt");
                                 table_output("mydata3",fname3);char fname4[100];
                                 get_filzbach_path(fname4, 100);
                                 strcat_s(fname4,"_my_outputYEAR4.txt");
                                 table_output("mydata4",fname4);char fname5[100];
                                 get_filzbach_path(fname5, 100);
                                 strcat_s(fname5,"_my_outputYEAR5.txt");
                                 table_output("mydata5",fname5); 
return;
}
               
               /* ************************************************* */
               /* The control function that calls everything else.  */
               /* ************************************************* */
               
               int main()
{
               atexit(pause);
               // set the likelihood function pointer
               pfn_likelihood = &likelihood;
               initialize_filzbach();
               name_analysis("Priors015"); 
            read_data();
             setup_parameters();
             set_chains(2);
             runmcmc(10000, 20000, 0, 0); // Burn-in, Bayes, MLE adaption, MLE search
             final_output();
}
             #endif
             