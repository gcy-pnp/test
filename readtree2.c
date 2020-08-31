//读取root文件中的tree，并处理分析

TH1D *hTOF;
TH1D *tdiff;
TH1D *dtd;
TH1D *qdiff;
TH1D *dqd;

void readtree2()
{
    const Double_t D = 500.;      //cm, distance between target and the scin.(Center)
    const Double_t L = 100.;      //cm, half length of the scin.
    const Double_t dD = 5.;       //cm, thickness of the scin.
    const Double_t TRes = 1.;     //ns, time resolution(FWHM) of the scintillator.
    const Double_t Lambda = 380.; //cm, attenuation lenght of the scin.
    const Double_t QRes = 0.1;    //relative energy resolution(FWHM) of the scin.
    const Double_t Vsc = 7.5;     //ns/cm, speed of light in the scin.
    const Double_t En0 = 100;     //MeV, average neutron energy
    const Double_t EnRes = 50.;   //MeV, energy spread of neutron(FWHM)
    const Double_t Eg0 = 1;       //MeV, gamma energy
    const Double_t Rg = 0.3;      //ratio of gamma,ratio of neutron 1-Rg

    // 1.打开文件，得到TTree指针
    TFile *ipf = new TFile("tree.root");     //打开ROOT文件
    TTree *tree = (TTree *)ipf->Get("tree"); //得到名字为“tree”的TTree指针

    //2. 声明tree的Branch变量
    Double_t x;
    Double_t e;
    int pid;
    Double_t tof, ctof;
    Double_t tu, td;
    Double_t qu, qd;
    Double_t tx, d, ntof, Tof0, qx, ce;

    //3. 将变量指向对应Branch的地址
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("tof", &tof);
    tree->SetBranchAddress("ctof", &ctof);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("tu", &tu);
    tree->SetBranchAddress("td", &td);
    tree->SetBranchAddress("qu", &qu);
    tree->SetBranchAddress("qd", &qd);

    TFile *opf = new TFile("test.root", "recreate");
    //Histogram
    hTOF = new TH1D("hTOF", "Time of flight", 1000, 0, 100);
    tdiff = new TH1D("tdiff", "td-tu", 140, -20, 50);
    qdiff = new TH1D("qdiff", "log(qd/qu)", 50, -1, 1);
    //4. 逐事件读取tree的branch数据
    Long64_t nentries = tree->GetEntries();
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        tree->GetEntry(jentry);
        hTOF->Fill(ctof);
        tdiff->Fill(td - tu);
        qdiff->Fill(TMath::Log(qu / qd));
    }
    dtd = new TH1D("dtd", "dt/dx", 141, -20.25, 50.25);
    for (int i = 1; i < tdiff->GetNbinsX(); i++)
    {
        Double_t df = tdiff->GetBinContent(i + 1) - tdiff->GetBinContent(i);
        dtd->Fill(tdiff->GetBinLowEdge(i + 1), df);
    }

    hTOF->Fit("gaus");
    dtd->Fit("gaus");
    dtd->Sumw2(0);

    dqd = new TH1D("dqd", "dq/dx", 51, -1.02, 1.02);
    for (int i = 1; i < qdiff->GetNbinsX(); i++)
    {
        Double_t df2 = qdiff->GetBinContent(i + 1) - qdiff->GetBinContent(i);
        dqd->Fill(qdiff->GetBinLowEdge(i + 1), df2);
    }
    TF1 *f3 = new TF1("f3", "gaus", -0.7, -0.3);
    TF1 *f4 = new TF1("f4", "gaus", 0.3, 0.7);
    f3->SetParameter(0, 600);
    f3->SetParameter(1, -0.5);
    f3->SetParameter(2, 0.05);
    f4->SetParameter(0, -600);
    f4->SetParameter(1, 0.5);
    f4->SetParameter(2, 0.05);
    dqd->Fit(f3, "R");
    dqd->Fit(f4, "R+");
    dqd->Sumw2(0);

    ipf->Close();
    hTOF->Write();
    dtd->Write();
    dqd->Write();
    opf->Close();
}