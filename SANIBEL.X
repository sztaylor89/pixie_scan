{
gStyle->SetOptStat(00000000);
V->Draw("vandle.qdc:vandle.tof-983>>(600,50,350,50000,0,50000)","vandle.vtype==1 && vandle.bsnrl>16 && vandle.bsnrr>16 && vandle.bqdcl>50 && vandle.bqdcr>50 && tape.move==0","COLZ");
}
