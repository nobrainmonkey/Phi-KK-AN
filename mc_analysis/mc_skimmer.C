// Run from mc_analysis/ after ./build_upcdst.sh:  root -l -b -q mc_skimmer.C
#include "Rtypes.h" // defines R__ADD_INCLUDE_PATH / R__LOAD_LIBRARY
// clang-format off
R__ADD_INCLUDE_PATH(../mc/upcdst/include)
R__LOAD_LIBRARY(lib/libUPCDst.so)
// clang-format on

#include "stdio.h"
#include "iostream"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"

#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCVertex.h"
#include "StUPCMCEvent.h"
#include "StUPCMCTrack.h"
#include "StUPCMCVertex.h"

StUPCTrack *bestRecoTrack(const StUPCMCTrack *mc)
{
    StUPCTrack *best = nullptr;
    for (Int_t k = 0; k < mc->getNumberOfRecoTracks(); ++k)
    {
        StUPCTrack *rc = mc->getRecoTrack(k);
        if (!rc) continue;
        if (!best || rc->getQaTruth() > best->getQaTruth()
            || (rc->getQaTruth() == best->getQaTruth() && rc->getNhitsFit() > best->getNhitsFit()))
            best = rc;
    }
    return best;   // nullptr if the MC track was not reconstructed
}

void mc_skimmer(const char *input = "../mc/phi_kk_flat_rap_flat_pt_4.root")
{
    TFile *inputFile = TFile::Open(input);
    TTree *upcTree = (TTree *)inputFile->Get("mUPCTree");
    StUPCEvent *evt = nullptr;
    upcTree->SetBranchAddress("mUPCEvent", &evt);
    Int_t nEntries = upcTree->GetEntries();

    TH1F* evt_zdc_coin = new TH1F("zdc_coin","zdc_coin",120,0,6000);

    for (Int_t e = 1; e < nEntries; e++)
    {
        upcTree->GetEntry(e);
        evt_zdc_coin->Fill(evt->getZDCCoincRate());
        Int_t nRCTracks = evt->getNumberOfTracks();
        // sample RC loop
        for (Int_t trc = 0; trc < nRCTracks; trc++)
        {
            StUPCTrack *rc_trk = evt->getTrack(trc);
            // RC to MC association rcmc_trk mean reco's monte carlo track
            StUPCMCTrack *rcmc_trk = rc_trk->getMCTrack();
            if (!rcmc_trk) continue;
        }
        Int_t nMCTracks = evt->getNumberOfMCTracks();
        // sample MC loop
        for (Int_t tmc = 0; tmc < nMCTracks; tmc++)
        {
            StUPCMCTrack *mc_trk = evt->getMCTrack(tmc);

            // MC to RC association, mcrc_trk mean mc's reco track
            // with no argument, it returns the first mathced reco track. 
            StUPCTrack *mcrc_trk=mc_trk->getRecoTrack();
        }
    } // event loop end

    evt_zdc_coin->Draw();
    evt_zdc_coin->SaveAs("zdc_coin.pdf");
}
