#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
void step(cmplx* psi0, cmplx* psi1, const double k, const int Nx, const double dx, const double dt);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = ((xmax-xmin)/(Nx-1));
	const double dt =  dx/10 ; //hab ich jetzt einfach so wie in lab 12 gewählt, ka wieso genau 
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5); //damit rundet er richtig, da int abschneidet

	const double lambda = 10;
	const double omega = 0.2;
	const double alpha = sqrt(omega); //aus Bedingungen aus Aufgabe
	const double k = omega*omega;
	
	stringstream strm;

	cmplx* psi0 = new cmplx[Nx]; //dynamisches array
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];
	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  
		  step(psi0,psi1,k,Nx,dx,dt);
		  
		  t+=dt;
		  
		  h = psi0;
		  psi0 = psi1;
		  psi1 = h;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
	delete[] psi0;
	delete[] psi1;
	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  ); //analytische Lösung aus Aufgabenstellung
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl; 
	 // Ausgabe: x, Betragsquadrat von phi, Re und Im von phi, und Btrag-,Im und Re der analytischen Lösung 
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 ); // aus Aufgabenstellung
	}
}

//--------------------------------------
void step(cmplx* psi0, cmplx* psi1, const double k, const int Nx, const double dx, const double dt){
  
    cmplx d[Nx];
    cmplx a = cmplx(0,-dt/(4*dx*dx)); // aus Aufgabenstellung. Muss kein Vektor sein, da a überall gleich 
    cmplx aStern = cmplx(0,dt/(4*dx*dx));
    //u tilde mache ich nicht, da u tilde immer = a ist 
    cmplx psiS[Nx]; //psi tilde 
    cmplx dStern[Nx]; // d* 
    for(int i = 0 ; i<Nx ; i++){
     double x = i*dx;
      double Imd = dt/(2*dx*dx) + dt*0.25*k*x*x; //Imaginärteil von d 
      d[i] = cmplx(1,Imd);
      dStern[i] = cmplx(1,-Imd);
    }
    
    cmplx dS[Nx]; // d schlange in upper triangular matrix  
   
    //--------forward substitution
    dS[0] = d[0];
    psiS[0] = psi0[0];
    for (int i = 1; i < Nx; i++){
      dS[i] = d[i]-a/d[i-1]*a ;                 //aus Berechungen siehe Aufzeichnungen
      psiS[i] = psi0[i] - a/d[i-1]*psiS[i-1];
    }
    
    //--------backward substitution
  psi1[Nx-1] = 1.0/(dS[Nx-1]) * (aStern*psiS[Nx-2] + dStern[Nx-1]*psiS[Nx-1]);
  
  for (int i = Nx-2; i>0; i--){ // soll von unten nach oben gehen
    psi1[i] = 1.0/(dS[i]) * (aStern*psiS[i-1] + dStern[i]*psiS[i] + aStern*psiS[i+1] - a*psi1[i+1]);
  }
  
  psi1[0] = 1.0/(dS[0]) * (dStern[0]*psiS[0] + aStern*psiS[1] - a*psi1[1]);
    
  
}
