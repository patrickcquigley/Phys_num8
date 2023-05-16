#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double> > vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void
triangular_solve(vector<T> const& diag,
                 vector<T> const& lower,
                 vector<T> const& upper,
                 vector<T> const& rhs,
                 vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// Potentiel V(x) :
double
V(double const& x, double const& m, double const& delta, double const& w1, double const& w2)

{
// TODO: Initialiser le potentiel  
double vl;
double vr;

vl=0.5*m*pow(w1,2)*pow((x+delta),2);
vr=0.5*m*pow(w1,2)*pow((x-delta),2);

return min(vl,vr);
       


}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule entre les points de maillage nL et nR,
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.
double
prob(vec_cmplx const& psi, int nL, int nR, double dx);
double
E(vec_cmplx const& psi,
  vec_cmplx const& diagH,
  vec_cmplx const& lowerH,
  vec_cmplx const& upperH,
  double const& dx);
double
xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double
x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double
pmoy(vec_cmplx const& psi, double const& dx, double const& hbar);
double
p2moy(vec_cmplx const& psi, double const& dx, double const& hbar);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx
normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.

int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double w1 = configFile.get<double>("w1");
    double w2 = configFile.get<double>("w2");
    double delta = configFile.get<double>("delta");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double

    // Parametres numeriques :

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");
    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        
        x[i] =xL+i*(xR-xL)/Npoints;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);
  
    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
        double x0 = configFile.get<double>("x0");
        double k0 = .0; //à modifier par l'élève
        double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);
        for (int i(0); i < Npoints; ++i)
            psi[i] = exp(i*k0*x[i])*exp(-((x[i]-x0)*(x[i]-x0))/(2*pow(sigma0,2)));
        // Modifications des valeurs aux bords :
        psi[0] = complex<double>(0., 0.);
        psi[Npoints - 1] = complex<double>(0., 0.);
        // Normalisation :
        psi = normalize(psi, dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

  complex<double> a = .0; // Coefficient complexe a de l'equation (4.100), à modifier par l'élève
    
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        complex<double> b = .0; // Coefficient complexe b de l'equation (4.100), à modifier par l'élève
	// TODO: éléments de matrices diagonales
        dH[i] = .0;
        dA[i] = .0;
        dB[i] = .0;
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
	// TOTO: éléments de matrices sous- et sur- diagonales
        aH[i] = .0;
        cH[i] = .0;
        aA[i] = .0;
        cA[i] = .0;
        aB[i] = .0;
        cB[i] = .0;
    }

    // Conditions aux limites: psi nulle aux deux bords
    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
    //dA[0] = ....
    //... 

    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(x[i], m, delta, w1, w2) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(15);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // Boucle temporelle :
    double t = 0;
    
    //TODO: intersection des deux paraboles (x=x_c)
    double x_local_max = 0; 
    
    // TODO: Nombre d'intervalles entre xL et x_local_max (x=x_c)
    unsigned int Nx0 = 1;
    
    for (int nt = 0.; nt <= Nsteps; nt += 1) {
        // Ecriture de |psi|^2 and phase information
        for (int i(0); i < Npoints; ++i){
            fichier_psi << pow(abs(psi[i]), 2) 
            		 << " " 
            		 << real(psi[i]) 
            		 << " " 
            		 << imag(psi[i]) << " ";
            		 }
            fichier_psi << endl;

        // Ecriture des observables :
        fichier_observables << t << " " 
        		     << prob(psi, 0, Nx0, dx)
                            << " " // probabilite que la particule soit en x < x0
                            << prob(psi, Nx0, Npoints, dx)
                            << " " // probabilite que la particule soit en x >= x0
                            << E(psi, dH, aH, cH, dx) << " " // Energie
                            << xmoy (psi, x,  dx) << " "       // Position moyenne
                            << x2moy(psi, x,  dx) << " "      // Position^2 moyenne
                            << pmoy (psi, dx, hbar) << " "    // Quantite de mouvement moyenne
                            << p2moy(psi, dx, hbar) << endl; // (Quantite de mouvement)^2 moyenne
        // Only do writing of data on the last step
       if (nt < Nsteps) {
            // Calcul du membre de droite :
            vec_cmplx psi_tmp(Npoints, 0.);
            // TODO: Programmer l'algorithme semi-implicite
            psi = psi_tmp; // à modifier!
            
            t += dt;
        }
    } // Fin de la boucle temporelle
    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}

double
prob(vec_cmplx const& psi, int nL, int nR, double dx)
{
    // TODO: calculer la probabilite de trouver la particule entre les points nL et nR
    double resultat(0.);
    return resultat;
}

double
E(vec_cmplx const& psi,
  vec_cmplx const& diagH,
  vec_cmplx const& lowerH,
  vec_cmplx const& upperH,
  double const& dx)
{   
    // TODO: calcul de l’énergie de la particule, moyenne de l’Hamiltonien
    vec_cmplx psi_tmp(psi.size(), 0.);
    double resultat(0.);

    // TODO: calculer H(psi)*psi
    
    
    // TODO: calculer Integrale de conj(psi)* H(psi)*psi dx
    
    return resultat;
}

double
xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
    // TODO: calcul de la position moyenne de la particule <x>
    double resultat(0.);

    for (int i(0); i < psi.size(); ++i){
    resultat+=real(dx*(conj(psi[i])*x[i]*psi[i]+conj(psi[i+1])*x[i+1]*psi[i+1])*0.5);
    }

    return resultat;
}

double
x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
    // TODO: calcul de la x^2 moyenne de la particule <x^2>
    double resultat(0.);
    return resultat;
}

double
pmoy(vec_cmplx const& psi, double const& dx, double const& hbar)
{
    // TODO: calcul de la quantité de mouvement moyenne de la particule <p>
    double resultat(0.);
    return resultat;
}

double
p2moy(vec_cmplx const& psi, double const& dx, double const& hbar)
{
    // TODO: calcul de la p^2 moyenne de la particule <p^2>
    double resultat(0.);
    return resultat;
}

vec_cmplx
normalize(vec_cmplx const& psi, double const& dx)
{
    // TODO: Normalisation de la fonction d'onde initiale psi
    vec_cmplx psi_norm(psi.size(), 0.);
    return psi_norm;
}
