// Example C++ function "CppModTemplate" dynamically loaded into "load.edp"
// ------------------------------------------------------------------------
#include <ff++.hpp>
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function
#include <math.h>
#include <iostream>
#include <algorithm>
//#include "MeshPoint.hpp"
using namespace Fem2D;

  
 Matrice_Creuse<double> *Dmx;
 Matrice_Creuse<double> *Dmy;
 Matrice_Creuse<double> *Dmz;
 Matrice_Creuse<double> *Dpx;
 Matrice_Creuse<double> *Dpy;
 Matrice_Creuse<double> *Dpz;
 
 Matrice_Creuse<double> *Dmx0;
 Matrice_Creuse<double> *Dmy0;
 Matrice_Creuse<double> *Dmz0;
 Matrice_Creuse<double> *Dpx0;
 Matrice_Creuse<double> *Dpy0;
 Matrice_Creuse<double> *Dpz0;
 
 Matrice_Creuse<double> *Dmxx;
 Matrice_Creuse<double> *Dmyy;
 Matrice_Creuse<double> *Dmzz;
 Matrice_Creuse<double> *Dpxx;
 Matrice_Creuse<double> *Dpyy;
 Matrice_Creuse<double> *Dpzz;
 Matrice_Creuse<double> *D0xx;
 Matrice_Creuse<double> *D0yy;
 Matrice_Creuse<double> *D0zz;

double dx,dy,dz;
int Nreg, Nupwind;
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double G(const double &u1,const double &u2,const double &u3,const double &u4,const double &u5,const double &u6)
{
  double  p1=max(0.,u1),  p2=min(0.,u2),  p3=max(0.,u3),  p4=min(0.,u4), p5=max(0.,u5),  p6=min(0.,u6);
  return sqrt(p1*p1+p2*p2+p3*p3+p4*p4+p5*p5+p6*p6);
}

inline double Plus(const double &u) 
{
  double  p1=max(0.,u);
  return p1;
}

inline double Moins(const double &u)
{
  double  p1=min(0.,u);
  return p1;
}

inline double Minmod(const double &a, const double &b)
{
  double sa= a > 0. ? 1. : -1.;
  double sb= b > 0. ? 1. : -1.;
  return    max(0.,sa*sb)*sa*min(abs(a),abs(b));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Derivatives1(Matrice_Creuse<double> *const &dpx,Matrice_Creuse<double> *const &dpy,Matrice_Creuse<double> *const &dpz,Matrice_Creuse<double> *const &dmx,Matrice_Creuse<double> *const &dmy,Matrice_Creuse<double> *const &dmz)
{
 Dmx = dmx;
 Dmy = dmy;
 Dmz = dmz;
 Dpx = dpx;
 Dpy = dpy;
 Dpz = dpz;
return 0;
}

double Derivatives21(Matrice_Creuse<double> *const &dpxx,Matrice_Creuse<double> *const &dpyy,Matrice_Creuse<double> *const &dpzz,Matrice_Creuse<double> *const &dmxx,Matrice_Creuse<double> *const &dmyy,Matrice_Creuse<double> *const &dmzz)
{
 Dmxx = dmxx;
 Dmyy = dmyy;
 Dmzz = dmzz;
 Dpxx = dpxx;
 Dpyy = dpyy;
 Dpzz = dpzz;
return 0;
}

double Derivatives22(Matrice_Creuse<double> *const &d0xx,Matrice_Creuse<double> *const &d0yy,Matrice_Creuse<double> *const &d0zz)
{
 D0xx = d0xx;
 D0yy = d0yy;
 D0zz = d0zz;
return 0;
}

double Derivatives3(Matrice_Creuse<double> *const &dpx0,Matrice_Creuse<double> *const &dpy0,Matrice_Creuse<double> *const &dpz0,Matrice_Creuse<double> *const &dmx0,Matrice_Creuse<double> *const &dmy0,Matrice_Creuse<double> *const &dmz0)
{
 Dmx0 = dmx0;
 Dmy0 = dmy0;
 Dmz0 = dmz0;
 Dpx0 = dpx0;
 Dpy0 = dpy0;
 Dpz0 = dpz0;
return 0;
}

double Parameters(const double & dx_, const double & dy_, const double & dz_,const long & Nupwind_,const long & Nreg_)
{
dx = dx_;
dy = dy_;
dz = dz_;
Nreg = Nreg_;
Nupwind = Nupwind_;
return 0;
}

double HJUpwind(KN<double> *const & psi0,KN<double> *const & psi,KN<double> *const & v,const double & dtReini,const double & dtUpwind)
{
  
   KN<double> &Psi0(*psi0);
   KN<double> &Psi(*psi);
   KN<double> &V(*v);
   int nm=Psi.N();
   KN<double> DmxPsi0(nm), DpxPsi0(nm), DmyPsi0(nm), DpyPsi0(nm), DmzPsi0(nm), DpzPsi0(nm),DmxPsi(nm), DpxPsi(nm), DmyPsi(nm), DpyPsi(nm), DmzPsi(nm), DpzPsi(nm),DmxxPsi(nm), DpxxPsi(nm), DmyyPsi(nm), DpyyPsi(nm), DmzzPsi(nm), DpzzPsi(nm) , D0xxPsi(nm), D0yyPsi(nm),D0zzPsi(nm),A(nm),B(nm),C(nm),D(nm),E(nm),F(nm),GammaP(nm),GammaM(nm),S0(nm),Sp(nm),Sm(nm),Vp(nm),Vm(nm);
   
   
   //Psi  = Psi0;
   ffassert(Dmx->A);
   ffassert(Dmy->A);
   ffassert(Dmz->A);
   ffassert(Dpx->A);
   ffassert(Dpy->A);
   ffassert(Dpz->A);
   
   ffassert(Dmx0->A);
   ffassert(Dmy0->A);
   ffassert(Dmz0->A);
   ffassert(Dpx0->A);
   ffassert(Dpy0->A);
   ffassert(Dpz0->A);
   
   ffassert(Dmxx->A);
   ffassert(Dmyy->A);
   ffassert(Dmzz->A);
   ffassert(Dpxx->A);
   ffassert(Dpyy->A);
   ffassert(Dpzz->A);
   ffassert(D0xx->A);
   ffassert(D0yy->A);
   ffassert(D0zz->A);
			
   for(int i=0;i<Nupwind;i++){	
		Psi0 = Psi;
		if (i % 2 == 0) //Frequency of the reinitialization
	{
	DmxPsi0=*(Dmx0->A)*Psi0;
    DpxPsi0=*(Dpx0->A)*Psi0;
	DmyPsi0=*(Dmy0->A)*Psi0;
	DpyPsi0=*(Dpy0->A)*Psi0;
	DmzPsi0=*(Dmz0->A)*Psi0;
	DpzPsi0=*(Dpz0->A)*Psi0;
			
		for(int j=0;j<Nreg;j++)
		{
		// Computation of the differentials of Psi
			
			DmxPsi=*(Dmx->A)*Psi;
			DpxPsi=*(Dpx->A)*Psi;
			DmyPsi=*(Dmy->A)*Psi;
			DpyPsi=*(Dpy->A)*Psi;
			DmzPsi=*(Dmz->A)*Psi;
			DpzPsi=*(Dpz->A)*Psi;
			
			DmxxPsi=*(Dmxx->A)*Psi;
			DpxxPsi=*(Dpxx->A)*Psi;
			D0xxPsi=*(D0xx->A)*Psi;
			DmyyPsi=*(Dmyy->A)*Psi;
			DpyyPsi=*(Dpyy->A)*Psi;
			D0yyPsi=*(D0yy->A)*Psi;
			DmzzPsi=*(Dmzz->A)*Psi;
			DpzzPsi=*(Dpzz->A)*Psi;
			D0zzPsi=*(D0zz->A)*Psi;
			
			for(int ii=0;ii<nm;ii++)
			{  
			A(ii)=DmxPsi(ii)+dx/2.*Minmod(DmxxPsi(ii),D0xxPsi(ii));
			B(ii)=DpxPsi(ii)-dx/2.*Minmod(DpxxPsi(ii),D0xxPsi(ii));
			C(ii)=DmyPsi(ii)+dy/2.*Minmod(DmyyPsi(ii),D0yyPsi(ii));
			D(ii)=DpyPsi(ii)-dy/2.*Minmod(DpyyPsi(ii),D0yyPsi(ii));
			E(ii)=DmzPsi(ii)+dz/2.*Minmod(DmzzPsi(ii),D0zzPsi(ii));
			F(ii)=DpzPsi(ii)-dz/2.*Minmod(DpzzPsi(ii),D0zzPsi(ii));
			
		// Gamma functions	
			GammaP(ii)=G(A(ii),B(ii),C(ii),D(ii),E(ii),F(ii));
			GammaM(ii)=G(B(ii),A(ii),D(ii),C(ii),F(ii),E(ii));
		// S functions
			S0(ii)=Psi0(ii)/sqrt(Psi0(ii)*Psi0(ii)+max(dx,max(dy,dz))*(DmxPsi0(ii)*DmxPsi0(ii)+DpxPsi0(ii)*DpxPsi0(ii)+DmyPsi0(ii)*DmyPsi0(ii)+DpyPsi0(ii)*DpyPsi0(ii)+DmzPsi0(ii)*DmzPsi0(ii)+DpzPsi0(ii)*DpzPsi0(ii))/20.+1.e-10);
			Sp(ii)=Plus(S0(ii));
			Sm(ii)=Moins(S0(ii));
			
			Psi(ii)=Psi(ii)-dtReini*(Sm(ii)*GammaM(ii)+Sp(ii)*GammaP(ii)-S0(ii));
			}			
		}
		
	}
			DmxPsi=*(Dmx->A)*Psi;
			DpxPsi=*(Dpx->A)*Psi;
			DmyPsi=*(Dmy->A)*Psi;
			DpyPsi=*(Dpy->A)*Psi;
			DmzPsi=*(Dmz->A)*Psi;
			DpzPsi=*(Dpz->A)*Psi;
			
			DmxxPsi=*(Dmxx->A)*Psi;
			DpxxPsi=*(Dpxx->A)*Psi;
			D0xxPsi=*(D0xx->A)*Psi;
			DmyyPsi=*(Dmyy->A)*Psi;
			DpyyPsi=*(Dpyy->A)*Psi;
			D0yyPsi=*(D0yy->A)*Psi;
			DmzzPsi=*(Dmzz->A)*Psi;
			DpzzPsi=*(Dpzz->A)*Psi;
			D0zzPsi=*(D0zz->A)*Psi;
			
			for(int ii=0;ii<nm;ii++)
			{  
			A(ii)=DmxPsi(ii)+dx/2.*Minmod(DmxxPsi(ii),D0xxPsi(ii));
			B(ii)=DpxPsi(ii)-dx/2.*Minmod(DpxxPsi(ii),D0xxPsi(ii));
			C(ii)=DmyPsi(ii)+dy/2.*Minmod(DmyyPsi(ii),D0yyPsi(ii));
			D(ii)=DpyPsi(ii)-dy/2.*Minmod(DpyyPsi(ii),D0yyPsi(ii));
			E(ii)=DmzPsi(ii)+dz/2.*Minmod(DmzzPsi(ii),D0zzPsi(ii));
			F(ii)=DpzPsi(ii)-dz/2.*Minmod(DpzzPsi(ii),D0zzPsi(ii));
			
		// Gamma functions	
			GammaP(ii)=G(A(ii),B(ii),C(ii),D(ii),E(ii),F(ii));
			GammaM(ii)=G(B(ii),A(ii),D(ii),C(ii),F(ii),E(ii));
		
			Vm(ii)=Moins(V(ii));
			Vp(ii)=Plus(V(ii));
	
			Psi(ii)=Psi(ii)-dtUpwind*(Vm(ii)*GammaM(ii)+Vp(ii)*GammaP(ii));	
			}
  }
   
  return 0.;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EuclideanDistance(KN<double> *const & P, KN<double> *const & Q)
{
   KN<double> &PP(*P);
   KN<double> &QQ(*Q);
   double d = sqrt((PP(0)-QQ(0))*(PP(0)-QQ(0))+(PP(1)-QQ(1))*(PP(1)-QQ(1))+(PP(2)-QQ(2))*(PP(2)-QQ(2)));
   return d;
}

double VolumeTetra(KN<double> *const &  P1, KN<double> *const & P2,KN<double> *const & P3, KN<double> *const &  P4)
{  
  double a4 = EuclideanDistance(P1,P2);//a4 d12
  double a6 = EuclideanDistance(P1,P3);//a6 d13
  double a1 = EuclideanDistance(P1,P4);//a1 d14
  double a5 = EuclideanDistance(P2,P3);//a5 d23
  double a2 = EuclideanDistance(P2,P4);//a2 d24
  double a3 = EuclideanDistance(P3,P4);//a3 d34
  double E1 = (pow(a1,2)*pow(a5,2)*(pow(a2,2)+pow(a3,2)+pow(a4,2)+pow(a6,2)-pow(a1,2)-pow(a5,2)));
  double E2 = (pow(a2,2)*pow(a6,2)*(pow(a1,2)+pow(a3,2)+pow(a4,2)+pow(a5,2)-pow(a2,2)-pow(a6,2)));
  double E3 = (pow(a3,2)*pow(a4,2)*(pow(a1,2)+pow(a2,2)+pow(a5,2)+pow(a6,2)-pow(a3,2)-pow(a4,2)));
  double E4 = pow(a1,2)*pow(a2,2)*pow(a4,2)+pow(a2,2)*pow(a3,2)*pow(a5,2)+pow(a1,2)*pow(a3,2)*pow(a6,2)+pow(a4,2)*pow(a5,2)*pow(a6,2);
  return sqrt((1./144)*(E1+E2+E3-E4));
}

double AreaQuad(KN<double> *const &  P1, KN<double> *const & P2,KN<double> *const & P3, KN<double> *const &  P4)
{
  double a = EuclideanDistance(P1,P2);
  double b = EuclideanDistance(P2,P3);
  double c = EuclideanDistance(P3,P4);
  double d = EuclideanDistance(P4,P1);
  double s = (a+b+c+d)/2;
  return sqrt((s-a)*(s-b)*(s-c)*(s-d));
}

double Height(KN<double> *const &  P,KN<double> *const &  P1, KN<double> *const & P2,KN<double> *const & P3)
{
  KN<double> PP(*P);
  KN<double> PP1(*P1);
  KN<double> PP2(*P2);
  KN<double> PP3(*P3);
  
  KN<double> P10,P12,P13,N(*P);
  P12 = PP2-PP1;
  P13 = PP3-PP1;
  P10 = PP-PP1;
  N(0) = (P12(1)*P13(2)-P12(2)*P13(1));
  N(1) = -(P12(0)*P13(2)-P12(2)*P13(0));
  N(2) = (P12(0)*P13(1)-P12(1)*P13(0));
  
  double NormN;
  NormN = sqrt(pow(N(0),2)+pow(N(1),2)+pow(N(2),2));
  return abs(P10(0)*N(0)+P10(1)*N(1)+P10(2)*N(2))/NormN;
}

double VolumeTetraCut(KN<double> *const & P12,KN<double> *const & P24, KN<double> *const & P34,KN<double> *const & P13,KN<double> *const & P1, KN<double> *const & P4)
{
  KN<double> PP12(*P12);
  KN<double> PP24(*P24);
  KN<double> PP34(*P34);
  KN<double> PP13(*P13);
  
  KN<double> PP124, PP134;
  PP124 = (PP12+PP24);
  PP134 = (PP13+PP34);
  
  KN<double> P124(*P12),P134(*P13); 
  P124(0)=PP124(0)/2;P124(1)=PP124(1)/2;P124(2)=PP124(2)/2;
  P134(0)=PP134(0)/2;P134(1)=PP134(1)/2;P134(2)=PP134(2)/2;
  
  KN<double> * const Q124 = &P124;
  KN<double> * const Q134 = &P134;
  
  double R1,R2,T1;
  R1  = (1./3)*AreaQuad(P12,Q124,Q134,P13)*Height(P1,P12,Q124,Q134);
  R2  = (1./3)*AreaQuad(Q124,P24,P34,Q134)*Height(P4,Q124,P24,P34);
  T1 = VolumeTetra(Q124,Q134,P4,P1);
  return (R1+R2+T1);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
KN<double>  ElementV1;
KN<double>  ElementV2;
KN<double>  ElementV3;
KN<double>  ElementV4;
KN<double>  ElementArea;
KN<double>  Vx;
KN<double>  Vy;
KN<double>  Vz;

double SetElementInfo(KN<double> *const &EV1,KN<double> *const &EV2,KN<double> *const &EV3,KN<double> *const &EV4,KN<double> *const &Area)
{
 KN<double> ElementV1c(*EV1);
 KN<double> ElementV2c(*EV2);
 KN<double> ElementV3c(*EV3);
 KN<double> ElementV4c(*EV4);
 KN<double> ElementAreac(*Area);
 ElementV1=ElementV1c;
 ElementV2=ElementV2c;
 ElementV3=ElementV3c;
 ElementV4=ElementV4c;
 ElementArea=ElementAreac;
return 0;
}

double SetNodeInfo(KN<double> *const &vx,KN<double> *const &vy,KN<double> *const &vz)
{
KN<double> Vxc(*vx);
KN<double> Vyc(*vy);
KN<double> Vzc(*vz);
Vx = Vxc;
Vy = Vyc;
Vz = Vzc;
return 0;
}

KN<double> Case1Split(const int & i,const double & SumSign,KN<double>* const & signpsi,KN<double>* const & psi)
{
int N = 4;
KN<double> SignPsi(*signpsi);
KN<double> Psi(*psi);
KN<double> NodesTetra(12);

int pos = static_cast<int>(abs(SumSign+1));
int neg = static_cast<int>(abs(SumSign-1));
int Positive[pos];
int Negative[neg];
int ipos=0,ineg=0;

int index[4];
index[0] = static_cast<int>(ElementV1(i));
index[1] = static_cast<int>(ElementV2(i));
index[2] = static_cast<int>(ElementV3(i));
index[3] = static_cast<int>(ElementV4(i));

for(int k=0; k<N; k++)
	{
	if(SignPsi(k)==1){Positive[ipos]=k;ipos+=1;}
	if(SignPsi(k)==-1){Negative[ineg]=k;ineg+=1;}
	}
double t = 0;
if(SumSign>0)//Only one negative node //Psi1*t+Psi2*(1-t)=0 //Psi2/(Psi2-Psi1)=t
	{
	NodesTetra(0) = Vx(index[Negative[0]]);
	NodesTetra(1) = Vy(index[Negative[0]]);
	NodesTetra(2) = Vz(index[Negative[0]]);
	for(int k=1; k<N; k++)
		{
		t = Psi(Positive[k-1])/(Psi(Positive[k-1])-Psi(Negative[0]));
		NodesTetra(3*k) = t*Vx(index[Negative[0]]) +(1-t)*Vx(index[Positive[k-1]]); 
		NodesTetra(3*k+1) = t*Vy(index[Negative[0]]) +(1-t)*Vy(index[Positive[k-1]]); 
		NodesTetra(3*k+2) = t*Vz(index[Negative[0]]) +(1-t)*Vz(index[Positive[k-1]]); 
		}
	}
else //Only one positive node
{
NodesTetra(0) = Vx(index[Positive[0]]);
NodesTetra(1) = Vy(index[Positive[0]]);
NodesTetra(2) = Vz(index[Positive[0]]);
for(int k=1; k<N; k++)
		{
		t = Psi(Negative[k-1])/(Psi(Negative[k-1])-Psi(Positive[0]));
		NodesTetra(3*k) = t*Vx(index[Positive[0]]) +(1-t)*Vx(index[Negative[k-1]]); 
		NodesTetra(3*k+1) = t*Vy(index[Positive[0]]) +(1-t)*Vy(index[Negative[k-1]]); 
		NodesTetra(3*k+2) = t*Vz(index[Positive[0]]) +(1-t)*Vz(index[Negative[k-1]]); 
		}
}
return NodesTetra;
}

KN<double> Case2Split(const int & i,const double & SumSign,KN<double>* const & signpsi,KN<double>* const & psi)
{
int N = 4;
KN<double> SignPsi(*signpsi);
KN<double> Psi(*psi);
KN<double> NodesPrism(24);//18

int Positive[2];
int Negative[2];
int ipos=0,ineg=0;

int index[4];
index[0] = static_cast<int>(ElementV1(i));
index[1] = static_cast<int>(ElementV2(i));
index[2] = static_cast<int>(ElementV3(i));
index[3] = static_cast<int>(ElementV4(i));

for(int k=0; k<N; k++)
	{
	if(SignPsi(k)==1){Positive[ipos]=k;ipos+=1;}
	if(SignPsi(k)==-1){Negative[ineg]=k;ineg+=1;}
	}
double t = 0;
NodesPrism(0) = Vx(index[Negative[0]]);
NodesPrism(1) = Vy(index[Negative[0]]);
NodesPrism(2) = Vz(index[Negative[0]]);
for(int k=0; k<2; k++)
{
	t = Psi(Positive[k])/(Psi(Positive[k])-Psi(Negative[0]));
	//if(t>1){cout<<"HELP"<<endl;}
	NodesPrism(3*(k+1)) = t*Vx(index[Negative[0]]) +(1-t)*Vx(index[Positive[k]]); 
	NodesPrism(3*(k+1)+1) = t*Vy(index[Negative[0]]) +(1-t)*Vy(index[Positive[k]]); 
	NodesPrism(3*(k+1)+2) = t*Vz(index[Negative[0]]) +(1-t)*Vz(index[Positive[k]]); 
}
NodesPrism(9) = Vx(index[Negative[1]]);
NodesPrism(10) = Vy(index[Negative[1]]);
NodesPrism(11) = Vz(index[Negative[1]]);
for(int k=0; k<2; k++)
{
	t = Psi(Positive[k])/(Psi(Positive[k])-Psi(Negative[1]));
	//if(t>1){cout<<"HELP"<<endl;}
	NodesPrism(3*k+12) = t*Vx(index[Negative[1]]) +(1-t)*Vx(index[Positive[k]]); 
	NodesPrism(3*k+12+1) = t*Vy(index[Negative[1]]) +(1-t)*Vy(index[Positive[k]]); 
	NodesPrism(3*k+12+2) = t*Vz(index[Negative[1]]) +(1-t)*Vz(index[Positive[k]]); 
}
NodesPrism(18) = Vx(index[Positive[0]]);
NodesPrism(19) = Vy(index[Positive[0]]);
NodesPrism(20) = Vz(index[Positive[0]]);
NodesPrism(21) = Vx(index[Positive[1]]);
NodesPrism(22) = Vy(index[Positive[1]]);
NodesPrism(23) = Vz(index[Positive[1]]);
return NodesPrism;
}

double DensityCpp(const double & eps,KN<double> *const &density,KN<double> *const &psi,KN<double> *const &signpsi)
{
   KN<double> &Density(*density);
   KN<double> &Psi(*psi);
   KN<double> &SignPsi(*signpsi);
   int nt=Density.N();
   
   KN<double> PsiVal(4);
   KN<double> SPsiVal(4);
   KN<double> NodesTetra(12);
   KN<double> NodesPrism(24);//18
   
   for (int i=0;i<nt;i++)
		{
		int i1 = static_cast<int>(ElementV1(i));
		int i2 = static_cast<int>(ElementV2(i));
		int i3 = static_cast<int>(ElementV3(i));
		int i4 = static_cast<int>(ElementV4(i));
		
		PsiVal(0) = Psi(i1);PsiVal(1) = Psi(i2);PsiVal(2) = Psi(i3);PsiVal(3) = Psi(i4);
		SPsiVal(0) = SignPsi(i1);SPsiVal(1) = SignPsi(i2);SPsiVal(2) = SignPsi(i3);SPsiVal(3) = SignPsi(i4);
		//SPsiVal=(SignPsi(i1),SignPsi(i2),SignPsi(i3),SignPsi(i4)); NE PAS INITIALISER COMME CA!
		double SumSign = (SPsiVal(0)+SPsiVal(1)+SPsiVal(2)+SPsiVal(3));
		double ProdSign = SPsiVal(0)*SPsiVal(1)*SPsiVal(2)*SPsiVal(3);
		// double sPsiVal[]={Psi(i1),Psi(i2),Psi(i3),Psi(i4)};
		// double sSPsiVal[]={SignPsi(i1),SignPsi(i2),SignPsi(i3),SignPsi(i4)};
		// double sSumSign = (sSPsiVal[0]+sSPsiVal[1]+sSPsiVal[2]+sSPsiVal[3]);
		// double sProdSign = sSPsiVal[0]*sSPsiVal[1]*sSPsiVal[2]*sSPsiVal[3];
		
			if(ProdSign ==0.){cout<<"Warning! Iso-0 matching with one or more mesh nodes"<<endl;}
			else
			{
				if(SumSign==4){Density(i) = eps;} //Empty case
				else if (SumSign==-4){Density(i)=1.;} // Full case
				else if(abs(SumSign)==2)// simpler cut case
				{
					NodesTetra = Case1Split(i,SumSign,&SPsiVal,&PsiVal);
					KN<double> P1(3),P2(3),P3(3),P4(3);
					for (int j=0;j<3;j++){P1(j)=NodesTetra(j);}
					for (int j=0;j<3;j++){P2(j)=NodesTetra(3*1+j);}
					for (int j=0;j<3;j++){P3(j)=NodesTetra(3*2+j);}
					for (int j=0;j<3;j++){P4(j)=NodesTetra(3*3+j);}
					double voltetra = VolumeTetra(&P1,&P2,&P3,&P4);
					//if(SumSign>0){Density(i)=0*eps*(1.-voltetra/ElementArea(i))+(voltetra/ElementArea(i));}
					//else{Density(i)=(1.-voltetra/ElementArea(i))+0*eps*(voltetra/ElementArea(i));}
					if(SumSign>0){Density(i)=max(eps,voltetra/ElementArea(i));}
					else{Density(i)=max(eps,1.-voltetra/ElementArea(i));}
				} 
				else // harder cut case (SumSign==0)
				{
					if(SumSign!=0){cout<<"ALERT! Case not considered"<<endl;}
					NodesPrism = Case2Split(i,SumSign,&SPsiVal,&PsiVal);
					KN<double> P1(3),P12(3),P13(3),P4(3),P24(3),P34(3);
					for (int j=0;j<3;j++){P1(j)=NodesPrism(j);}
					for (int j=0;j<3;j++){P12(j)=NodesPrism(3*1+j);}
					for (int j=0;j<3;j++){P13(j)=NodesPrism(3*2+j);}
					for (int j=0;j<3;j++){P4(j)=NodesPrism(3*3+j);}
					for (int j=0;j<3;j++){P24(j)=NodesPrism(3*4+j);}
					for (int j=0;j<3;j++){P34(j)=NodesPrism(3*5+j);}
					
					KN<double> P2(3),P3(3);
					for (int j=0;j<3;j++){P2(j)=NodesPrism(3*6+j);}
					for (int j=0;j<3;j++){P3(j)=NodesPrism(3*7+j);}
					
					//double volprism = VolumeTetraCut(&P12,&P24,&P34,&P13,&P1,&P4);//CES FORMULES NE MARCHENT PAS! 
					//double volprism2 = VolumeTetraCut(&P12,&P13,&P34,&P24,&P2,&P3);
					
					double tt1=VolumeTetra(&P24,&P13,&P4,&P34);
					double tt2=VolumeTetra(&P12,&P13,&P1,&P24);
					double tt3=VolumeTetra(&P13,&P24,&P1,&P4);
					double volprism3 = (tt1+tt2+tt3);
					
					// double tt4=VolumeTetra(&P13,&P24,&P3,&P34);
					// double tt5=VolumeTetra(&P12,&P24,&P2,&P13);
					// double tt6=VolumeTetra(&P24,&P13,&P2,&P3);
					// double volprism4 = (tt4+tt5+tt6);
					
					Density(i) = max(eps,volprism3/ElementArea(i));
					//Density(i)=0*eps*(1.-volprism/ElementArea(i))+(volprism/ElementArea(i));//CORRIGER! 
					//if((volprism+volprism2)!=ElementArea(i))
					//{
					//cout<<"Alert! "<<(volprism+volprism2)<<" vs "<<(volprism3+volprism4)<<" vs "<<ElementArea(i)<<endl;
					//cout<<"Alert! "<<volprism<<" vs "<<volprism2<<" vs "<<volprism3<<endl;
					//}
					//if(Density(i)>1){
					//cout<<" P12 "<<P12<<" P24 "<<P24<<cout<<" P34 "<<P34<<" P13 "<<P13<<cout<<" P1 "<<P1<<" P4 "<<P4<<" P3 "<<P3<<" P2 "<<P2<<endl;
					//cout<<"Volume 1: "<<volprism<<"Volume 2: "<<volprism2<<" Add "<<(volprism+volprism2)<<" real "<<voltetra<<" Freefem "<<ElementArea(i)<<endl;
					//cout<<"Psi: "<<PsiVal<<endl;
					//}
				} 
			}
		}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   add the function name to the freefem++ table 
class Init { public:
  Init();
};
Init init;
Init::Init(){
  // Add function with 3 arguments
  Global.Add("G","(",new OneOperator6_<double, double, double , double, double, double>(G)); 
  Global.Add("Minmod","(",new OneOperator2_<double, double >(Minmod)); 
  Global.Add("Plus","(",new OneOperator1_<double >(Plus)); 
  Global.Add("Moins","(",new OneOperator1_<double >(Moins)); 
  Global.Add("Schema","(",new OneOperator5_<double,KN<double>*, KN<double>*, KN<double>*,double,double>(HJUpwind)); 
  Global.Add("Parameters","(",new OneOperator5_<double,double,double,double,long,long>(Parameters)); 
  Global.Add("Derivatives1","(",new OneOperator6_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives1)); 
  Global.Add("Derivatives21","(",new OneOperator6_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives21)); 
  Global.Add("Derivatives22","(",new OneOperator3_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives22)); 
  Global.Add("Derivatives3","(",new OneOperator6_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives3)); 
  Global.Add("EuclideanDistance","(",new OneOperator2_<double, KN<double>*, KN<double>*>(EuclideanDistance)); 
  Global.Add("VolumeTetra","(",new OneOperator4_<double, KN<double>*, KN<double>*, KN<double>*,KN<double>*>(VolumeTetra)); 
  Global.Add("AreaQuad","(",new OneOperator4_<double, KN<double>*, KN<double>*, KN<double>*,KN<double>*>(AreaQuad)); 
  Global.Add("Height","(",new OneOperator4_<double, KN<double>*, KN<double>*, KN<double>*,KN<double>*>(Height)); 
  Global.Add("VolumeTetraCut","(",new OneOperator6_<double, KN<double>*, KN<double>*, KN<double>*,KN<double>*,KN<double>*,KN<double>*>(VolumeTetraCut)); 
  Global.Add("SetElementInfo","(",new OneOperator5_<double, KN<double>*, KN<double>*, KN<double>*,KN<double>*,KN<double>*>(SetElementInfo)); 
  Global.Add("SetNodeInfo","(",new OneOperator3_<double, KN<double>*, KN<double>*,KN<double>*>(SetNodeInfo)); 
  Global.Add("DensityCpp","(",new OneOperator4_<double,double, KN<double>*, KN<double>*,KN<double>*>(DensityCpp)); 
  //Global.Add("Case1Split","(",new OneOperator5_<double,int,double,KN<double>*,KN<double>*,KN<double>*>(Case1Split)); 
}
