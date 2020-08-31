//读取root文件中的tree，并处理分析

TCanvas *c1;
TCanvas *c2;
TCanvas *c3;
TH1D *hTOF;
TH1D *tdiff;
TH1D *dtd;
TH1D *qdiff;
TH1D *dqd;
TH1D *htxx;
TH2D *hgtofx;
TH1D *hgctof;

void readtree()
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
    c1 = new TCanvas("c1", "c1");
    c2 = new TCanvas("c2", "c2");
    c3 = new TCanvas("c3", "c3");

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
    TF1 *f1 = new TF1("f1", "gaus", -14, -9);
    TF1 *f2 = new TF1("f2", "gaus", 39.5, 43);
    f2->SetParameter(0, -300);
    f2->SetParameter(1, 42);
    f2->SetParameter(2, 0.2);
    dtd->Fit(f1, "R");
    dtd->Fit(f2, "R+");
    dtd->Sumw2(0);
    Double_t p1, p2, p3, p4;
    p1 = f1->GetParameter(1);
    p2 = f2->GetParameter(1);

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
    p3 = f3->GetParameter(1);
    p4 = f4->GetParameter(1);

    c1->cd();
    dtd->Draw();
    c1->Draw();
    c2->cd();
    dqd->Draw();
    c2->Draw();

    Double_t tof0 = 3.333 * ((D + dD / 2) * 0.01);
    hgtofx = new TH2D("hgtofx", "hgtofx", 100, -120, 120, 100, 39, 48);
    htxx = new TH1D("htxx", "htx-x", 100, -20, 20);
    hgctof = new TH1D("hgctof", "hgctof", 100, 39, 48);
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    { //对每个事件进行遍历
        tree->GetEntry(jentry);
        tx = 2 * L / (p2 - p1) * (td - tu - (p1 + p2) / 2);
        htxx->Fill(tx - x);
        if (ctof > 42 && ctof < 44.5)
        {
            hgtofx->Fill(tx, ctof);
            if (abs(tx) < 5)
                hgctof->Fill(ctof); //gamma hits the center of the det.
        }
    }
    c3->cd();
    hgctof->Draw();
    TF1 *f5 = new TF1("f5", "gaus");
    hgctof->Fit(f5);
    Tof0 = f5->GetParameter(1);
    c3->Draw();

    cout << "p1"
         << "\t\t"
         << "p2"
         << "\t"
         << "p3"
         << "\t\t"
         << "p4"
         << "\t\t"
         << "Tof0" << endl;
    cout << p1 << "\t" << p2 << "\t" << p3 << "\t" << p4 << "\t" << Tof0 << endl;


}