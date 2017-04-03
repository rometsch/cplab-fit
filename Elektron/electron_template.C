// Elektronenstreuung in Materie
// Zu ersetzende oder fehlende Programmteile sind gekennzeichnet

// Durch Kompilieren des Programms kann die Ausf�hrungszeit gewaltig verk�rzt werden (s. ROOT-Manual)

#include<iostream>
#include<fstream>
#include<math.h>
#include <TF1.h>
#include <TList.h>
#include <TView.h>
#include <TPolyline3D.h>
#include <TCanvas.h>
#include <TVector3.h>

using namespace std;

// constants

const double pi = 3.14159265358;      // pi
const double me = 0.511;              // electron mass [MeV]

const int N_MAX = 10000;
const int N_ELEC = 50;           // number of electrons
const double X_MAX = 100;        // x-size of the detector
const double Y_MAX = 100;        // y-size of the detector
const double E_0 = 20;           // energy of electrons
const double E_MIN = 0.02;       // minimum energy of electrons
TF1 *f_length;
TF1 *f_theta;
TF1 *f_phi;

void calc_new_direction( double otheta, double ophi,
                      double beta, double gamma,
                      double* angles, double* vec);
void electron_template(void);

class electron: public TObject {
public :
    electron(double x1, double y1, double z1, double E1 , double theta1, double phi1);
    ~electron();
    double x;
    double y;
    double z;
    double E; // kin. energy
    double theta;
    double phi;
};

electron::electron(double x1, double y1, double z1, double E1 , double theta1, double phi1) {
    x=x1; y=y1; z=z1; E=E1; theta=theta1; phi=phi1;
}

electron::~electron(){};

int nextcoord(electron *e1, electron *e2)
{

// Keep starting position of struck electron e2

    e2->x = e1->x;
    e2->y = e1->y;
    e2->z = e1->z;
    double otheta = e1->theta;  // old theta
    double ophi   = e1->phi;    // old phi

// calculate new coordinates

// Moeller scattering
    double n_theta = f_theta->GetRandom();
    double n_phi = f_phi->GetRandom();
    double len= f_length->GetRandom()*e1->E/10; // scattering length (scales with E)
    double p = sqrt(pow((e1->E+me),2)-pow(me,2));      // Momentum of e1 in target system

// transform angles into target system (done)
    double gamma = ( e1->E + 2*me )/sqrt( 2*pow(me,2) + 2*(e1->E + me)*me );
    double st = sin(n_theta);
    double ct = cos(n_theta);
    double theta1 = atan( st/gamma/( 1 + ct ) );
    double theta2 = atan( st/gamma/( 1 - ct ) );

    double phi1 = n_phi;
    double phi2 = n_phi;

// calculate new momenta and energies (stimmt nicht!)
    double p1 = 2*me*( e1->E + 2*me )*p*cos(theta1)/( pow(e1->E + 2*me,2) - pow(p*cos(theta1),2) );
    double p2 = 2*me*( e1->E + 2*me )*p*cos(theta2)/( pow(e1->E + 2*me,2) - pow(p*cos(theta2),2) );
    // double p1 = 0.9*p;
    // double p2 = 0.1*p;

// stimmt!
    e1->E = (sqrt(pow(p1,2)+pow(me,2))-me)-0.05*len/e1->E;
    e2->E = (sqrt(pow(p2,2)+pow(me,2))-me)-0.05*len/e1->E;

// stop if no energy left
    if (e1->E <= E_MIN) {return 0;}

  // calculate direction of scattered e1 in target system
  // z-axis aligns with direction of e1
    TVector3 r1(
      sin(theta1)*cos(phi1),
      sin(theta1)*sin(phi1),
      cos(theta1)
    );

// Rotate direction vector into reference system
// in which the z-axis aligns with the former direction of e1
    r1.RotateY(otheta);
    r1.RotateZ(ophi);

// New position and direction of e1
    e1->x += len*r1.X();
    e1->y += len*r1.Y();
    e1->z += len*r1.Z();
    e1->theta = r1.Theta();
    e1->phi = r1.Phi();

// calculate direction of scattered e2 in target system
// z-axis aligns with direction of e2
  TVector3 r2(
    sin(theta2)*cos(phi2),
    sin(theta2)*sin(phi2),
    cos(theta2)
  );

// Rotate direction vector into reference system
// in which the z-axis aligns with the former direction of e2
  r2.RotateY(otheta);
  r2.RotateZ(ophi);

// Direction of e2
    e2->theta = r2.Theta();
    e2->phi   = r2.Phi();


// stop if electron leaves detector
    if ((e1->x < 0) || (e1->x>X_MAX)) {cout << "detector left" << endl; return 0;}
    if ((e1->y < 0) || (e1->y>Y_MAX)) {cout << "detector left" << endl; return 0;}

  // return value 1 means success
    return 1;
}

void electron_template(void)
{
    ofstream outfile;
    ofstream logfile;
    TList electrons;    // List of electron objects (see root User's guide)
    electron* e1;       // first electron (whose trajectory is followed)
    electron* e2;       // second electron (which is just produced and traced afterwards)

// define graph
    TView* view;
    TPolyLine3D* traj = NULL;
    TCanvas* c1 = new TCanvas("SimCanvas", "SimCanvas", 1);
    int k = 0;

    // define functions for random number generator (kann so bleiben)
    f_length = new TF1("f_length", "exp(-x)" ,0 ,100);   // mean free path
    f_phi = new TF1("f_phi", "1" ,0 , 2*pi);             // azimuthal angle
    f_theta = new TF1("f_theta","pow(sin(x/2),-4)",0.2 ,pi/2); // Rutherford,
    // The cross section cannot be integrated to 0 degrees
    c1->cd();
    view = TView::CreateView(1);
    view->SetRange(0,0,0,100,100,100);
    view->ShowAxis();

    for (int n = 0; n < N_ELEC; n++) {
	// loop for electrons
	cout << "electron #" << n << endl;
	electrons.Add(new electron(50, 50 ,0 ,E_0, 0, 0));  // push new electron on stack
	for (int i=0; i<N_MAX; i++ ) {
	    e1 = (electron*)electrons.First();              // pop next starting point from stack
	    if (e1==NULL) break;
	    electrons.Remove(e1);                           // remove it from the stack
	    traj = new TPolyLine3D(1024);
	    k = 0;
	    traj->SetPoint(k, e1->x, e1->y, e1->z);
	    traj->Draw("ogl");
	    k++;
	    for (int j=0; j<N_MAX; j++) {
		e2 = new electron(0, 0 ,0 ,0 , 0, 0);
		int ret=nextcoord(e1, e2);
		traj->SetPoint(k, e1->x, e1->y, e1->z);
		traj->Draw("ogl");
		k++;
		if (ret == 0) break;
		if (e2->E > E_MIN) electrons.Add(e2);     // put secondary electron on the stack
	    }  // end while
      }
      c1->Modified();
      c1->Update();
  }  // end for
}

void calc_new_direction( double otheta, double ophi,
                      double beta, double gamma,
                      double* angles, double* vec) {
    /*  Rotate the unit vector specified by otheta and ophi
    *   by the euler angles alpha=0, beta, gamma.
    *   Store a unit vector of the new direction in vec
    *   and the azimuthal and lognitudonal in angles. */
    // Calculate unit vector.
    double st = sin(otheta);
    double sp = sin(ophi);
    double ct = cos(otheta);
    double cp = cos(ophi);
    double x = st*sp;
    double y = st*cp;
    double z = ct;

    // Calculate rotation matrix.
    double cb = cos(beta);
    double sb = sin(beta);
    double cg = cos(gamma);
    double sg = sin(gamma);
    double R11 = cb*cg;
    double R12 = sg;
    double R13 = -sb*cg;
    double R21 = -cg*sb;
    double R22 = cg;
    double R23 = sb*sg;
    double R31 = sb;
    double R32 = 0;
    double R33 = cg;

    // Rotate unit vector.
    vec[0] = R11*x + R12*y + R13*z;
    vec[1] = R21*x + R22*y + R23*z;
    vec[2] = R31*x + R32*y + R33*z;

    // Calculate angles from rotated unit vector.
    // theta
    angles[0] = acos( vec[2] );
    // phi
    if ( x == 0 && y == 0 ) {
      angles[1] = 0;
    } else {
      angles[1] = atan2(vec[1],vec[0]);
      if ( angles[1] < 0 ) angles[1] += 2*pi;
    }

}
