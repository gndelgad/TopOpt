// Example C++ function "CppModTemplate" dynamically loaded into "load.edp"
// ------------------------------------------------------------------------
#include <ff++.hpp>
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function
using namespace Fem2D;

  
 Matrice_Creuse<double> *Dmx;
 Matrice_Creuse<double> *Dmy;
 Matrice_Creuse<double> *Dpx;
 Matrice_Creuse<double> *Dpy;
 
 Matrice_Creuse<double> *Dmx0;
 Matrice_Creuse<double> *Dmy0;
 Matrice_Creuse<double> *Dpx0;
 Matrice_Creuse<double> *Dpy0;
 
 Matrice_Creuse<double> *Dmxx;
 Matrice_Creuse<double> *Dmyy;
 Matrice_Creuse<double> *Dpxx;
 Matrice_Creuse<double> *Dpyy;
 Matrice_Creuse<double> *D0xx;
 Matrice_Creuse<double> *D0yy;

double dx,dy;
int Nreg, Nupwind;
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline double G(const double &u1,const double &u2,const double &u3,const double &u4)
{
  double  p1=max(0.,u1),  p2=min(0.,u2),  p3=max(0.,u3),  p4=min(0.,u4);
  return sqrt(p1*p1+p2*p2+p3*p3+p4*p4);
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

//////////////////////
inline double AdaptMeshLevelSet(const double &b)
{
  //const char * command= "/home/delgado/TOP/T4/MMG/mmg-5.2.1-Linux-4.2.0-1-amd64/bin/mmg3d_O3 -in cube -out adpcube.mesh -sol LevelSet.sol -hausd 0.001 -hmax 0.05 -hmin 0.03 -hgrad 1.3 -nr -ls";
  //-ls -nr
  const char * command= "mmg2d_O3 -in cube.mesh -out adpcube.mesh -sol LevelSet.sol -hausd 0.001 -ls -hmax 0.05 -hmin 0.02 -hgrad 1.3 -nr";//   
  system(command);
  return b;
}
/////////////////////////////////////////////////////////////////////////////
inline double ReinitCharac(const double &u)
{
  
  //const char * command="C:/Users/gabriel.delgado-keef/Dropbox/Projet_TOP/T4/MeshDist/sources_Windows7/Dist.exe Sh.LevelSet.mesh";
  const char * command="Dist_W7 Th.LevelSet.mesh";
  system(command);
  
   std::ifstream    inFile("Th.LevelSet.sol");
   std::ofstream    outFile("Th.LevelSet.chi.sol");
   outFile << inFile.rdbuf();
  
  return 1.;
}

inline double AdvectionCharac(const double &u)
{
 
  string s;
   {
        std::ostringstream oss;
        oss << u;
        s = oss.str();
    }
  //const char * command="C:/Users/Gabriel/Desktop/advection/sourcesAdvect/Advect.exe -t "+buffer+" Sh.LevelSet.mesh";//-t 0.001
  //string commands = "C:/Users/gabriel.delgado-keef/Dropbox/Projet_TOP/T4/Advect/sourcesAdvect_Windows7/Advect.exe -t "+s+" Sh.LevelSet.mesh";//-r -rk4 
  string commands = "Advect_W7 -t "+s+" Th.LevelSet.mesh";//-r -rk4
  
  size_t size = commands.size() + 1;
  char * buffer = new char[ size ];
  strncpy( buffer, commands.c_str(), size );
  cout << buffer << '\n'; // "une chaîne de caractères" 
  // libérer la mémoire
  //delete [] buffer;
  system(buffer);
  
  std::ifstream    inFile("Th.LevelSet.new.chi.sol");
  std::ofstream    outFile("Th.LevelSet.chi.sol");
  outFile << inFile.rdbuf();
  
  return 1.;
}

double ReadSolution(KN<double> *const & psi)
{
  ifstream fp ("Th.LevelSet.chi.sol");
  //ifstream fp ("Sh.LevelSet.sol");
  string line;
  setprecision(16);
  KN<double> &Psi(*psi);
  
  for (int j=0;j<Psi.N()+8;++j)
  {
   getline(fp,line);
   if(j>=8){Psi(j-8)=atof(line.c_str());
   }
  }
  fp.close();
  
  return 1.;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Derivatives1(Matrice_Creuse<double> *const &dpx,Matrice_Creuse<double> *const &dpy,Matrice_Creuse<double> *const &dmx,Matrice_Creuse<double> *const &dmy)
{
 Dmx = dmx;
 Dmy = dmy;
 Dpx = dpx;
 Dpy = dpy;
return 0;
}

double Derivatives2(Matrice_Creuse<double> *const &dpxx,Matrice_Creuse<double> *const &dpyy,Matrice_Creuse<double> *const &dmxx,Matrice_Creuse<double> *const &dmyy,Matrice_Creuse<double> *const &d0xx,Matrice_Creuse<double> *const &d0yy)
{
 Dmxx = dmxx;
 Dmyy = dmyy;
 Dpxx = dpxx;
 Dpyy = dpyy;
 D0xx = d0xx;
 D0yy = d0yy;
return 0;
}

double Derivatives3(Matrice_Creuse<double> *const &dpx0,Matrice_Creuse<double> *const &dpy0,Matrice_Creuse<double> *const &dmx0,Matrice_Creuse<double> *const &dmy0)
{
 Dmx0 = dmx0;
 Dmy0 = dmy0;
 Dpx0 = dpx0;
 Dpy0 = dpy0;
return 0;
}

double Parameters(const double & dx_, const double & dy_,const long & Nupwind_,const long & Nreg_)
{
dx = dx_;
dy = dy_;
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
   KN<double> DmxPsi0(nm), DpxPsi0(nm), DmyPsi0(nm), DpyPsi0(nm),DmxPsi(nm), DpxPsi(nm), DmyPsi(nm), DpyPsi(nm),DmxxPsi(nm), DpxxPsi(nm), DmyyPsi(nm), DpyyPsi(nm), D0xxPsi(nm), D0yyPsi(nm),A(nm),B(nm),C(nm),D(nm),GammaP(nm),GammaM(nm),S0(nm),Sp(nm),Sm(nm),Vp(nm),Vm(nm);
   
   
   //Psi  = Psi0;
   ffassert(Dmx->A);
   ffassert(Dmy->A);
   ffassert(Dpx->A);
   ffassert(Dpy->A);
   
   ffassert(Dmx0->A);
   ffassert(Dmy0->A);
   ffassert(Dpx0->A);
   ffassert(Dpy0->A);
   
   ffassert(Dmxx->A);
   ffassert(Dmyy->A);
   ffassert(Dpxx->A);
   ffassert(Dpyy->A);
   ffassert(D0xx->A);
   ffassert(D0yy->A);
   
			
   for(int i=0;i<Nupwind;i++){	
		Psi0 = Psi;
		if (i % 2 == 0) //Frequency of the reinitialization
	{
	DmxPsi0=*(Dmx0->A)*Psi0;
    DpxPsi0=*(Dpx0->A)*Psi0;
	DmyPsi0=*(Dmy0->A)*Psi0;
	DpyPsi0=*(Dpy0->A)*Psi0;
			
		for(int j=0;j<Nreg;j++)
		{
		// Computation of the differentials of Psi
			
			DmxPsi=*(Dmx->A)*Psi;
			DpxPsi=*(Dpx->A)*Psi;
			DmyPsi=*(Dmy->A)*Psi;
			DpyPsi=*(Dpy->A)*Psi;
			
			DmxxPsi=*(Dmxx->A)*Psi;
			DpxxPsi=*(Dpxx->A)*Psi;
			D0xxPsi=*(D0xx->A)*Psi;
			DmyyPsi=*(Dmyy->A)*Psi;
			DpyyPsi=*(Dpyy->A)*Psi;
			D0yyPsi=*(D0yy->A)*Psi;
			
			for(int ii=0;ii<nm;ii++)
			{  
			A(ii)=DmxPsi(ii)+dx/2.*Minmod(DmxxPsi(ii),D0xxPsi(ii));
			B(ii)=DpxPsi(ii)-dx/2.*Minmod(DpxxPsi(ii),D0xxPsi(ii));
			C(ii)=DmyPsi(ii)+dy/2.*Minmod(DmyyPsi(ii),D0yyPsi(ii));
			D(ii)=DpyPsi(ii)-dy/2.*Minmod(DpyyPsi(ii),D0yyPsi(ii));
			
		// Gamma functions	
			GammaP(ii)=G(A(ii),B(ii),C(ii),D(ii));
			GammaM(ii)=G(B(ii),A(ii),D(ii),C(ii));
		// S functions
			S0(ii)=Psi0(ii)/sqrt(Psi0(ii)*Psi0(ii)+max(dx,dy)*(DmxPsi0(ii)*DmxPsi0(ii)+DpxPsi0(ii)*DpxPsi0(ii)+DmyPsi0(ii)*DmyPsi0(ii)+DpyPsi0(ii)*DpyPsi0(ii))/20.+1.e-10);
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
			
			DmxxPsi=*(Dmxx->A)*Psi;
			DpxxPsi=*(Dpxx->A)*Psi;
			D0xxPsi=*(D0xx->A)*Psi;
			DmyyPsi=*(Dmyy->A)*Psi;
			DpyyPsi=*(Dpyy->A)*Psi;
			D0yyPsi=*(D0yy->A)*Psi;
			
			for(int ii=0;ii<nm;ii++)
			{  
			A(ii)=DmxPsi(ii)+dx/2.*Minmod(DmxxPsi(ii),D0xxPsi(ii));
			B(ii)=DpxPsi(ii)-dx/2.*Minmod(DpxxPsi(ii),D0xxPsi(ii));
			C(ii)=DmyPsi(ii)+dy/2.*Minmod(DmyyPsi(ii),D0yyPsi(ii));
			D(ii)=DpyPsi(ii)-dy/2.*Minmod(DpyyPsi(ii),D0yyPsi(ii));
			
		// Gamma functions	
			GammaP(ii)=G(A(ii),B(ii),C(ii),D(ii));
			GammaM(ii)=G(B(ii),A(ii),D(ii),C(ii));
		
			Vm(ii)=Moins(V(ii));
			Vp(ii)=Plus(V(ii));
	
			Psi(ii)=Psi(ii)-dtUpwind*(Vm(ii)*GammaM(ii)+Vp(ii)*GammaP(ii));	
			}
  }
   
   // Psi0=Psi;
   ////Second reinitialization
   	// {
	// DmxPsi0=*(Dmx0->A)*Psi0;
    // DpxPsi0=*(Dpx0->A)*Psi0;
	// DmyPsi0=*(Dmy0->A)*Psi0;
	// DpyPsi0=*(Dpy0->A)*Psi0;
			
		// for(int j=0;j<Nreg;j++)
		// {
		////Computation of the differentials of Psi
			
			// DmxPsi=*(Dmx->A)*Psi;
			// DpxPsi=*(Dpx->A)*Psi;
			// DmyPsi=*(Dmy->A)*Psi;
			// DpyPsi=*(Dpy->A)*Psi;
			
			// DmxxPsi=*(Dmxx->A)*Psi;
			// DpxxPsi=*(Dpxx->A)*Psi;
			// D0xxPsi=*(D0xx->A)*Psi;
			// DmyyPsi=*(Dmyy->A)*Psi;
			// DpyyPsi=*(Dpyy->A)*Psi;
			// D0yyPsi=*(D0yy->A)*Psi;
			
			// for(int ii=0;ii<nm;ii++)
			// {  
			// A(ii)=DmxPsi(ii)+dx/2.*Minmod(DmxxPsi(ii),D0xxPsi(ii));
			// B(ii)=DpxPsi(ii)-dx/2.*Minmod(DpxxPsi(ii),D0xxPsi(ii));
			// C(ii)=DmyPsi(ii)+dy/2.*Minmod(DmyyPsi(ii),D0yyPsi(ii));
			// D(ii)=DpyPsi(ii)-dy/2.*Minmod(DpyyPsi(ii),D0yyPsi(ii));
			
		////Gamma functions	
			// GammaP(ii)=G(A(ii),B(ii),C(ii),D(ii));
			// GammaM(ii)=G(B(ii),A(ii),D(ii),C(ii));
		////S functions
			// S0(ii)=Psi0(ii)/sqrt(Psi0(ii)*Psi0(ii)+max(dx,dy)*(DmxPsi0(ii)*DmxPsi0(ii)+DpxPsi0(ii)*DpxPsi0(ii)+DmyPsi0(ii)*DmyPsi0(ii)+DpyPsi0(ii)*DpyPsi0(ii))/20.+1.e-10);
			// Sp(ii)=Plus(S0(ii));
			// Sm(ii)=Moins(S0(ii));
			
			// Psi(ii)=Psi(ii)-dtReini*(Sm(ii)*GammaM(ii)+Sp(ii)*GammaP(ii)-S0(ii));
			// }			
		// }
		
	// }
  return 0.;
}

//   add the function name to the freefem++ table 
class Init { public:
  Init();
};
Init init;
Init::Init(){
  // Add function with 3 arguments
  Global.Add("AdaptMeshLevelSet","(",new OneOperator1_<double>(AdaptMeshLevelSet));
  Global.Add("ReadSolution","(",new OneOperator1_<double,KN<double>*>(ReadSolution)); 
  Global.Add("AdvectionCharac","(",new OneOperator1_<double>(AdvectionCharac)); 
  Global.Add("ReinitCharac","(",new OneOperator1_<double>(ReinitCharac)); 
  Global.Add("G","(",new OneOperator4_<double, double, double , double >(G)); 
  Global.Add("Minmod","(",new OneOperator2_<double, double >(Minmod)); 
  Global.Add("Plus","(",new OneOperator1_<double >(Plus)); 
  Global.Add("Moins","(",new OneOperator1_<double >(Moins)); 
  Global.Add("Schema","(",new OneOperator5_<double,KN<double>*, KN<double>*, KN<double>*,double,double>(HJUpwind)); 
  Global.Add("Parameters","(",new OneOperator4_<double,double,double,long,long>(Parameters)); 
  Global.Add("Derivatives1","(",new OneOperator4_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives1)); 
  Global.Add("Derivatives2","(",new OneOperator6_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives2)); 
  Global.Add("Derivatives3","(",new OneOperator4_<double,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *,Matrice_Creuse<double> *>(Derivatives3)); 
}



