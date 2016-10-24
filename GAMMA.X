{

gStyle->SetOptStat(111111111);

G->Draw("gamma.gen>>(500,500,1000,15000,0,15000)");

G->Print("*");

G->Scan("*");

}