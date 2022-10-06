#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <assert.h>

#include "Rtypes.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TThread.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TButton.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"

#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGTab.h"
#include "TGLayout.h"

#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"
#include "TString.h"
#include "TGraph.h"
#include "TTimeStamp.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TEnv.h"
#include "TTimeStamp.h"


using namespace std;

static TApplication *app = new TApplication("App", NULL, NULL);

enum tabname {
  charges, posXY, posScan, temp
};

//_____________________________________________________________________________
class TabWidget : public TGMainFrame {
private:
  TGTab *fTab;
  TRootEmbeddedCanvas *fCanvas[5];
  TGLayoutHints *fLcan;
  TGHorizontalFrame *fHframe0;
public:
  TabWidget();
  virtual ~TabWidget();
  inline TCanvas *GetCanvas(Int_t tabnum) { return fCanvas[tabnum]->GetCanvas(); }
  inline Int_t GetCurrent() const { return fTab->GetCurrent(); }

};

//_____________________________________________________________________________
void init_plots();
void reader_thread(void *arg = NULL);
void updater_thread(void *arg = NULL);
void analyze(int i);
void SetCharges(TH1F *h, int i);
void SetPosXY(TH2F *hEven,TH2F *hOdd,int i);
void SetTemp(TGraph *g);
bool Stop(unsigned int i);

//_____________________________________________________________________________
int main(int argc, char *argv[])
{
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetLabelSize(0.03,"xyz");
  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  init_plots();

  TThread reader("Reader", &reader_thread, NULL);
  reader.Run();
  TThread updater("Updater", &updater_thread, NULL);
  updater.Run();
  app->Run();

  return 0;
}

bool Stop(unsigned int i)
{
  sleep(i);
  return true;
}


//_____________________________________________________________________________
TabWidget::TabWidget() : TGMainFrame(gClient->GetRoot(), 1200, 800)
{
  int cout_check = 0;
  cout << cout_check++ << endl;
  SetCleanup(kDeepCleanup);
  fLcan = new TGLayoutHints(kLHintsCenterX|kLHintsCenterY|kLHintsExpandX|kLHintsExpandY,5,5,5,5);
  fHframe0 = new TGHorizontalFrame(this, 0, 0, 0);
  fHframe0->Resize(1200, 700);
  AddFrame(fHframe0, fLcan);
  fTab = new TGTab(fHframe0, 0, 0);
  fTab->AddTab("charges");
  fTab->AddTab("position scan even & odd");
  fTab->AddTab("position scan");
  fTab->AddTab("temp");
  fHframe0->AddFrame(fTab, fLcan);
  for(Int_t i=0;i<4;i++) {
    TGCompositeFrame *frame = fTab->GetTabContainer(i);
    fCanvas[i] = new TRootEmbeddedCanvas(Form("Canvas%d",i), frame, 1200, 700);
    frame->AddFrame(fCanvas[i]);
    frame->AddFrame(fCanvas[i], fLcan);
    TCanvas *canvas = fCanvas[i]->GetCanvas();
    canvas->SetGrid();
    canvas->SetFillColor(0);
  }
  cout << cout_check++ << endl;

  SetWindowName("DAQ Online Display");
  cout << cout_check++ << endl;
  MapSubwindows();
  cout << cout_check++ << endl;
  MapWindow();
  cout << cout_check++ << endl;
  Resize(1200, 800);
  cout << cout_check++ << endl;
}


//_____________________________________________________________________________
TabWidget::~TabWidget()
{
  Cleanup();
}

//_____________________________________________________________________________
const Int_t nhfadc_charges = 1;
const Int_t nhfadc_posXY = 2;
const Int_t nhfadc_posScan = 32;
static TH1F *hfadc_charges;
static TH2F *hfadc_posXYeven, *hfadc_posXYodd;
static TGraph *g_temp;
static TRandom *rnd;

Int_t runno;

static TabWidget *tab;


//_____________________________________________________________________________
void SetCharges(TH1F *h, int i)
{
  h->Reset();
  vector<double> *charge = 0;
  TFile *inputFile;
  TTree *tree1;
  TBranch *branch1 = 0;
  inputFile = TFile::Open(Form("./20201216_1924_cosmicRayTst.root"));
  tree1 = (TTree*)inputFile->Get("T_Event");
  tree1->SetBranchAddress("charges",&charge,&branch1);
  tree1->GetEntry(i);
  for(int j=0;j<charge->size();j++){
    h->Fill(charge->at(j));
  }
  cout << "Event number : " << i << endl;
  cout << "# of elements : " << charge->size() << endl;
  inputFile->Close();
  return;
}

void SetPosXY(TH2F *hEven,TH2F *hOdd,int i)
{
  hEven->Reset();
  hOdd->Reset();
  vector<double> *X=0, *Y=0;
  vector<int> *cellIDs=0;
  TFile *inputFile;
  TTree *treeX, *treeY;
  TBranch *branchX, *branchY, *branchcellIDs;
  inputFile = TFile::Open(Form("./20201216_1924_cosmicRayTst.root"));
  treeX = (TTree*)inputFile->Get("T_Event");
  treeY = (TTree*)inputFile->Get("T_Event");
  treeX->SetBranchAddress("posX",&X,&branchX);
  treeY->SetBranchAddress("posY",&Y,&branchY);
  treeX->SetBranchAddress("cellIDs",&cellIDs,&branchcellIDs);
  assert(treeX->GetEntries()==treeY->GetEntries());
  assert(X->size()==Y->size());
  treeX->GetEntry(i);
  treeY->GetEntry(i);

  for(int j=0;j<X->size();j++){
    int layer = cellIDs->at(j)%1000000;
    if(layer%2 == 0){
      hEven->Fill(X->at(j),Y->at(j));
    }else if(layer%2 != 0){
      hOdd->Fill(X->at(j),Y->at(j));
    }
  }
  return;
}

void SetTemp(TGraph *g, int cellID)
{
  vector<double> *temp = nullptr;
  vector<int> *cellIDs = nullptr;
  TFile *inputFile;
  TTree *tree;
  inputFile = TFile::Open(Form("./20201216_1924_cosmicRayTst.root"));
  tree = (TTree*)inputFile->Get("T_Event");
  tree->SetBranchAddress("temp",&temp);
  tree->SetBranchAddress("cellIDs",&cellIDs);

  int count = 0;

  for(int ientry=0;ientry<tree->GetEntries();ientry++)
  {
    tree->GetEntry(ientry);

    for(int i=0;i<temp->size();i++)
    {
      // if(cellIDs->at(i) == cellID){
      if(cellIDs->at(i) == 11050020){
        g->SetPoint(count,ientry,temp->at(i));
        // cout << ientry << endl;
        count++;
      }
    }
  }

  return;
}

void SetPedestal(TH1F *g, int i)
{
  vector<double> *temp = nullptr;
  vector<int> *cellIDs = nullptr;
  TFile *inputFile;
  TTree *tree;
  inputFile = TFile::Open(Form("./20201216_1924_cosmicRayTst.root"));
  tree = (TTree*)inputFile->Get("T_Event");
  tree->SetBranchAddress("temp",&temp);
  tree->SetBranchAddress("cellIDs",&cellIDs);

  int count = 0;

  for(int ientry=0;ientry<tree->GetEntries();ientry++)
  {
    tree->GetEntry(ientry);

    for(int i=0;i<temp->size();i++)
    {
      // if(cellIDs->at(i) == cellID){
      if(cellIDs->at(i) == 11050020){
        g->SetPoint(count,ientry,temp->at(i));
        // cout << ientry << endl;
        count++;
      }
    }
  }

  return;
}

//_____________________________________________________________________________
void updater_thread(void *arg)
{
  while (true) {

    TCanvas *charges_can = tab->GetCanvas(charges);
    TCanvas *posXY_can = tab->GetCanvas(posXY);
    TCanvas *temp_can = tab->GetCanvas(temp);

    charges_can->cd()->Modified();
    if(tab->GetCurrent() == charges) charges_can->Update();

    for(int i=0;i<nhfadc_posXY;i++){
      posXY_can->cd(i+1)->Modified();
    }
    if(tab->GetCurrent() == posXY) posXY_can->Update();

    temp_can->cd()->Modified();
    if(tab->GetCurrent() == temp) temp_can->Update();

  }
}


//_____________________________________________________________________________
void analyze(int i) {

      SetCharges(hfadc_charges,i);
      SetPosXY(hfadc_posXYeven,hfadc_posXYodd,i);
      SetTemp(g_temp,i);
}


//_____________________________________________________________________________
void reader_thread(void *arg)
{
  int i = 0;
  istream::int_type ch;
  while ((ch = cin.get()) != EOF) {
    cout << "Event # ? : ";
    cin >> i;
    analyze(i);
  }
}


//_____________________________________________________________________________
void init_plots() {

  TFile *inputFile;
  TTree *tree;
  inputFile = TFile::Open(Form("./20201216_1924_cosmicRayTst.root"));
  tree = (TTree*)inputFile->Get("T_Event");
  const int tree_entry = tree->GetEntries();


  rnd = new TRandom();
  rnd->SetSeed();

    hfadc_charges = new TH1F("charge","charge" , 100, 200, 600);
    hfadc_posXYeven = new TH2F("posXYeven","X vs Y", 5, -112.5, 112.5,45,-112.5,112.5);
    hfadc_posXYodd = new TH2F("posXYodd","X vs Y", 45, -112.5, 112.5,5,-112.5,112.5);
    g_temp = new TGraph();
    g_temp->SetTitle("temperature monitor");
    g_temp->GetXaxis()->SetTitle("event number");
    g_temp->GetYaxis()->SetTitle("temperature [C]");
  cout << "aa" << endl;
  tab = new TabWidget();
  cout << "aa" << endl;
  TCanvas *charges_can = tab->GetCanvas(charges);
  TCanvas *posXY_can = tab->GetCanvas(posXY);

    int cout_check = 0;
    cout << cout_check++ << endl;

  TCanvas *temp_can = tab->GetCanvas(temp);
  cout << cout_check++ << endl;

//  charges_can->Divide(2);
  posXY_can->Divide(2,1,1.e-7,1.e-7);
  cout << cout_check++ << endl;

  charges_can->cd(0)->SetGrid();
  hfadc_charges->Draw();
  charges_can->Update();
  charges_can->Flush();

  posXY_can->cd(1)->SetGrid();
  hfadc_posXYeven->Draw("colz");
  posXY_can->cd(2)->SetGrid();
  hfadc_posXYodd->Draw("colz");
  posXY_can->cd(0);
  posXY_can->Update();
  posXY_can->Flush();
  cout << cout_check++ << endl;

  temp_can->cd(0);
  cout << cout_check++ << endl;
  g_temp->Draw();
  cout << cout_check++ << endl;
  temp_can->Update();
  cout << cout_check++ << endl;
  temp_can->Flush();
    cout << cout_check++ << endl;

}
