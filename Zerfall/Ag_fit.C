{
gROOT->Reset();

Int_t ik;
Float_t v;

// Histogramm definieren
h1f = new TH1F("h1f","Decay of silver",120,0,600);
h1f->GetXaxis()->SetTitle("t/s");

// Daten einlesen
ifstream in;
in.open("Agdecay.dat");    
string s;
getline(in,s);  // 1.- 3. Zeile ueberlesen
getline(in,s);
getline(in,s);
for (int i=0; i<120; i++) {
    in >> ik >> v;
    h1f->SetBinContent(i,v);
}
in.close();

// Diagramm definieren    
c1 = new TCanvas("c1","Silver decay",200,10,700,900);
c1->SetFillColor(18);
gStyle->SetOptFit(1); 
gStyle->SetOptStat(0); 
h1f->SetFillColor(45);
h1f->Draw();

// Fitfunktion definieren
TF1 *agfit = new TF1("agfit","[0]+[2]*[1]*exp(-x*[2]) + [4]*[3]*exp(-x*[4])", 0,600);
agfit->SetParameter(0, 5);
agfit->SetParameter(1,1096);
agfit->SetParameter(2,0.03);
agfit->SetParameter(3,200);
agfit->SetParameter(4,5e-3);

h1f->Fit("agfit");
c1->Update();

}
