#include "../../../events/PlotStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
#include "cuts.h"
#include <fstream>
#include <vector>

// Run after pid_plots.C: root -l -b -q tof_pool.C
// Everything is in this function: read TOF tracks, generate phi decays,
// sample PID, then repeat at the measured daughter momenta.
void tof_pool() {
  const int nPBins = 100; // 0-1 GeV/c in 10 MeV/c bins
  const double pMax = 1.0;
  const double dp = pMax / nPBins;
  const double tofM2Min = .22, tofM2Max = .28;
  const int minDonors = 30, drawsPerPair = 4000;
  const int generatedPairs = 2000000;
  // Data comparison: 0.3-0.5, 0.5-0.7, 0.7-0.9 GeV/c.
  const double ptSplit1 = 0.5;
  const double ptSplit2 = 0.7;
  double ptEdges[4] = {pairPtMin, ptSplit1, ptSplit2, pairPtMax};
  TRandom3 rng(7182026);
  // Separate streams: changing the generated sample does not change the
  // random donor draws in the measured-daughter comparison.
  TRandom3 toyRng(7182027);
  TH1::AddDirectory(false);
  gSystem->mkdir("figures", true);

  TFile calibration(correctionFile);
  TH1D *shift = (TH1D *)calibration.Get("hMeanNSigmaKaon");
  TFile input(kaonFile);
  TTree *tree = (TTree *)input.Get("KaonTree");
  TFile dataFile("pid_plots.root");
  TTree *pairs = (TTree *)dataFile.Get("pairs");
  if (!shift || !tree || !pairs || !pairs->GetBranch("y")) {
    printf("Missing calibration, KaonTree, or pairs with rapidity. Run "
           "pid_plots.C "
           "first.\n");
    return;
  }
  double correctionMin = shift->GetBinCenter(1);
  double correctionMax = shift->GetBinCenter(shift->GetNbinsX());
  std::vector<double> *trackPt = 0, *trackEta = 0, *trackK = 0;
  std::vector<double> *trackE = 0, *trackPi = 0, *trackM2 = 0;
  std::vector<Short_t> *trackCharge = 0;
  std::vector<UShort_t> *trackHits = 0;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("track_pt", 1);
  tree->SetBranchStatus("track_eta", 1);
  tree->SetBranchStatus("track_nSigmaKaon", 1);
  tree->SetBranchStatus("track_nSigmaElectron", 1);
  tree->SetBranchStatus("track_nSigmaPion", 1);
  tree->SetBranchStatus("track_m2", 1);
  tree->SetBranchStatus("track_charge", 1);
  tree->SetBranchStatus("track_nHits", 1);
  tree->SetBranchAddress("track_pt", &trackPt);
  tree->SetBranchAddress("track_eta", &trackEta);
  tree->SetBranchAddress("track_nSigmaKaon", &trackK);
  tree->SetBranchAddress("track_nSigmaElectron", &trackE);
  tree->SetBranchAddress("track_nSigmaPion", &trackPi);
  tree->SetBranchAddress("track_m2", &trackM2);
  tree->SetBranchAddress("track_charge", &trackCharge);
  tree->SetBranchAddress("track_nHits", &trackHits);

  // Same indices in K, E, Pi and rawK always refer to the SAME TOF track.
  // First index: 0 = negative charge, 1 = positive charge.
  std::vector<float> poolK[2][nPBins], poolE[2][nPBins];
  std::vector<float> poolPi[2][nPBins], poolRawK[2][nPBins];
  // Number of donors that could pass the pair's kaon cut.
  int kaonDonors[2][nPBins] = {};
  long donorCount = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    for (unsigned int i = 0; i < trackPt->size(); ++i) {
      double p = (*trackPt)[i] * std::cosh((*trackEta)[i]);
      if (!std::isfinite(p))
        continue;
      if (p < 0.1 || p >= pMax)
        continue;
      if ((*trackHits)[i] < hitsMin)
        continue;
      if ((*trackM2)[i] <= tofM2Min || (*trackM2)[i] >= tofM2Max)
        continue;
      if (!std::isfinite((*trackK)[i]) || !std::isfinite((*trackE)[i]) ||
          !std::isfinite((*trackPi)[i]))
        continue;
      double correctionP = p;
      if (correctionP < correctionMin)
        correctionP = correctionMin;
      if (correctionP > correctionMax)
        correctionP = correctionMax;
      double k = (*trackK)[i] - shift->Interpolate(correctionP);
      int q = 0; // Negative track.
      if ((*trackCharge)[i] > 0)
        q = 1; // Positive track.
      // For example, p = 0.437 GeV/c goes in bin 43 (0.43-0.44).
      int b = int(p / dp);
      poolK[q][b].push_back(k);
      poolE[q][b].push_back((*trackE)[i]);
      poolPi[q][b].push_back((*trackPi)[i]);
      poolRawK[q][b].push_back((*trackK)[i]);
      if (k * k <= kkMax)
        ++kaonDonors[q][b];
      ++donorCount;
    }
  }
  printf("SL23e TOF pool: %ld tracks; %.3f GeV/c momentum bins\n", donorCount,
         dp);
  TFile output("tof_pool.root", "RECREATE");
  TH1D donorsMinus("donors_minus", ";p [GeV/c];Tracks", nPBins, 0, pMax);
  TH1D donorsPlus("donors_plus", ";p [GeV/c];Tracks", nPBins, 0, pMax);
  for (int b = 0; b < nPBins; ++b) {
    donorsMinus.SetBinContent(b + 1, poolK[0][b].size());
    donorsPlus.SetBinContent(b + 1, poolK[1][b].size());
  }

  // Read the nominal mass-window data before electron rejection.
  Double_t mass, pt, y, p1, p2, kk, ee, pipi;
  Short_t charge1, charge2;
  pairs->SetBranchAddress("mass", &mass);
  pairs->SetBranchAddress("pt", &pt);
  pairs->SetBranchAddress("y", &y);
  pairs->SetBranchAddress("p1", &p1);
  pairs->SetBranchAddress("p2", &p2);
  pairs->SetBranchAddress("charge1", &charge1);
  pairs->SetBranchAddress("charge2", &charge2);
  pairs->SetBranchAddress("kk", &kk);
  pairs->SetBranchAddress("ee", &ee);
  pairs->SetBranchAddress("pipi", &pipi);

  std::vector<double> parentPt;
  std::vector<double> parentY;
  TTree parentData("parent_data",
                   "Nominal mass-window pairs used for parent pT/y");
  parentData.Branch("pt", &pt);
  parentData.Branch("y", &y);
  TH2D dataParentPtY("data_parent_pt_y", ";Pair p_{T} [GeV/c];Pair y", 30,
                     pairPtMin, pairPtMax, 20, -pairYMax, pairYMax);
  TH2D generatedParentPtY("generated_parent_pt_y", ";Pair p_{T} [GeV/c];Pair y",
                          30, pairPtMin, pairPtMax, 20, -pairYMax, pairYMax);
  dataParentPtY.Sumw2();
  generatedParentPtY.Sumw2();
  for (Long64_t entry = 0; entry < pairs->GetEntries(); ++entry) {
    pairs->GetEntry(entry);
    if (mass <= massMin || mass >= massMax)
      continue;
    if (kk > kkMax)
      continue;
    if (pipi < pipiMin)
      continue;
    parentPt.push_back(pt);
    parentY.push_back(y);
    parentData.Fill();
    dataParentPtY.Fill(pt, y);
  }
  int nParents = parentPt.size();
  if (nParents == 0) {
    printf("No nominal mass-window pairs to sample.\n");
    return;
  }

  // Generated phi->KK: take pT and y from the SAME randomly chosen data pair.
  // Keep generating the mass and decay angles.
  TH1D generated("generated", "", 3, ptEdges);
  TH1D covered("covered", "", 3, ptEdges);
  TH1D toyBefore("toy_without_ee", "", 3, ptEdges);
  TH1D toyAfter("toy_with_ee", "", 3, ptEdges);
  generated.Sumw2();
  covered.Sumw2();
  toyBefore.Sumw2();
  toyAfter.Sumw2();
  for (int trial = 0; trial < generatedPairs; ++trial) {
    double mass;
    do {
      mass = toyRng.BreitWigner(1.019461, .004249);
    } while (mass <= massMin || mass >= massMax);
    int parentIndex = toyRng.Integer(nParents);
    double pt = parentPt[parentIndex];
    double y = parentY[parentIndex];
    generatedParentPtY.Fill(pt, y);
    double phi = toyRng.Uniform(-TMath::Pi(), TMath::Pi());
    double mt = std::hypot(mass, pt);
    TLorentzVector parent(pt * std::cos(phi), pt * std::sin(phi),
                          mt * std::sinh(y), mt * std::cosh(y));
    double costh;
    do {
      costh = toyRng.Uniform(-1, 1);
    } while (toyRng.Uniform() >
             1 - costh * costh); // sin^2(theta), helicity frame
    double angle = toyRng.Uniform(-TMath::Pi(), TMath::Pi());
    double pstar = .5 * std::sqrt(mass * mass - 4 * .493677 * .493677);
    TVector3 direction(pstar * std::sqrt(1 - costh * costh) * std::cos(angle),
                       pstar * std::sqrt(1 - costh * costh) * std::sin(angle),
                       pstar * costh);
    direction.RotateUz(parent.Vect().Unit());
    TLorentzVector plus(direction, mass / 2), minus(-direction, mass / 2);
    plus.Boost(parent.BoostVector());
    minus.Boost(parent.BoostVector());
    generated.Fill(pt);
    // These two restrictions come from the original pair-tree skim.
    if (plus.Pt() >= 1 || minus.Pt() >= 1)
      continue;
    int b1 = int(plus.P() / dp), b2 = int(minus.P() / dp);
    if (b1 < 0 || b1 >= nPBins || b2 < 0 || b2 >= nPBins)
      continue;
    if (kaonDonors[1][b1] < minDonors || kaonDonors[0][b2] < minDonors)
      continue;
    covered.Fill(pt);
    int i = toyRng.Integer(poolK[1][b1].size());
    int j = toyRng.Integer(poolK[0][b2].size());
    double k1 = poolK[1][b1][i], k2 = poolK[0][b2][j];
    double e1 = poolE[1][b1][i], e2 = poolE[0][b2][j];
    double pi1 = poolPi[1][b1][i], pi2 = poolPi[0][b2][j];
    double raw1 = poolRawK[1][b1][i], raw2 = poolRawK[0][b2][j];
    if (raw1 * raw1 + raw2 * raw2 >= 100)
      continue;
    double kk = k1 * k1 + k2 * k2;
    double ee = e1 * e1 + e2 * e2;
    double pipi = pi1 * pi1 + pi2 * pi2;
    if (kk > kkMax)
      continue;
    if (pipi < pipiMin)
      continue;
    toyBefore.Fill(pt);
    if (ee >= eeMin && ee >= eeFactor * kk)
      toyAfter.Fill(pt);
  }
  TH1D toyLoss("toy_electron_loss", "", 3, ptEdges);
  for (int b = 1; b <= 3; ++b) {
    double n = toyBefore.GetBinContent(b);
    if (n == 0)
      continue;
    double f = 1 - toyAfter.GetBinContent(b) / n;
    toyLoss.SetBinContent(b, 100 * f);
    toyLoss.SetBinError(b, 100 * std::sqrt(f * (1 - f) / n));
  }

  // Direct comparison: use each measured pair's daughter momenta and charges.
  // Count only the mass window. No fits, sidebands or background subtraction.
  // Bin 0 holds the total; bins 1-3 hold the pair-pT intervals.
  int total[4] = {};
  int matched[4] = {};
  int rejected[4] = {};
  double sumLoss[4] = {};
  double sumLoss2[4] = {};
  double sumMCVariance[4] = {};
  TH1D dataE("data_ee", "", 60, 0, 120);
  TH1D poolEHistogram("predicted_ee", "", 60, 0, 120);
  dataE.Sumw2();
  poolEHistogram.Sumw2();
  double predictedLoss = -1;
  bool hasPool = false, rejectedElectron = false;
  int acceptedDraws = 0;
  TTree match("matched_pairs", "Raw mass-window pairs and TOF predictions");
  match.Branch("pt", &pt);
  match.Branch("p1", &p1);
  match.Branch("p2", &p2);
  match.Branch("covered", &hasPool);
  match.Branch("accepted_draws", &acceptedDraws);
  match.Branch("predicted_loss", &predictedLoss);
  match.Branch("rejected_electron", &rejectedElectron);
  for (Long64_t entry = 0; entry < pairs->GetEntries(); ++entry) {
    pairs->GetEntry(entry);
    if (mass <= massMin || mass >= massMax)
      continue;
    if (kk > kkMax)
      continue;
    if (pipi < pipiMin)
      continue;
    int ptBin;
    if (pt < ptSplit1)
      ptBin = 1; // 0.3-0.5 GeV/c.
    else if (pt < ptSplit2)
      ptBin = 2; // 0.5-0.7 GeV/c.
    else
      ptBin = 3; // 0.7-0.9 GeV/c.

    ++total[0];
    ++total[ptBin];
    rejectedElectron = false;
    if (ee < eeMin || ee < eeFactor * kk)
      rejectedElectron = true;

    int b1 = int(p1 / dp);
    int b2 = int(p2 / dp);
    int q1 = 0;
    int q2 = 0;
    if (charge1 > 0)
      q1 = 1;
    if (charge2 > 0)
      q2 = 1;

    hasPool = true;
    if (b1 < 0 || b1 >= nPBins)
      hasPool = false;
    if (b2 < 0 || b2 >= nPBins)
      hasPool = false;
    if (hasPool) {
      if (kaonDonors[q1][b1] < minDonors)
        hasPool = false;
      if (kaonDonors[q2][b2] < minDonors)
        hasPool = false;
    }
    acceptedDraws = 0;
    int failedDraws = 0;
    double eeBins[62] = {}; // 60 histogram bins plus underflow and overflow.
    if (hasPool) {
      for (int trial = 0; trial < drawsPerPair; ++trial) {
        int i = rng.Integer(poolK[q1][b1].size());
        int j = rng.Integer(poolK[q2][b2].size());
        double k1 = poolK[q1][b1][i], k2 = poolK[q2][b2][j];
        double e1 = poolE[q1][b1][i], e2 = poolE[q2][b2][j];
        double pi1 = poolPi[q1][b1][i], pi2 = poolPi[q2][b2][j];
        double raw1 = poolRawK[q1][b1][i], raw2 = poolRawK[q2][b2][j];
        if (raw1 * raw1 + raw2 * raw2 >= 100)
          continue;
        double sampledKK = k1 * k1 + k2 * k2;
        double sampledEE = e1 * e1 + e2 * e2;
        double sampledPiPi = pi1 * pi1 + pi2 * pi2;
        if (sampledKK > kkMax)
          continue;
        if (sampledPiPi < pipiMin)
          continue;
        ++acceptedDraws;
        if (sampledEE < eeMin || sampledEE < eeFactor * sampledKK)
          ++failedDraws;
        ++eeBins[poolEHistogram.FindBin(sampledEE)];
      }
    }
    if (acceptedDraws < minDonors)
      hasPool = false;
    predictedLoss = -1;
    if (hasPool)
      predictedLoss = double(failedDraws) / acceptedDraws;
    match.Fill();
    if (!hasPool)
      continue;
    // Count this pair in the total (index 0) and in its pT bin.
    ++matched[0];
    ++matched[ptBin];
    if (rejectedElectron) {
      ++rejected[0];
      ++rejected[ptBin];
    }
    sumLoss[0] += predictedLoss;
    sumLoss[ptBin] += predictedLoss;
    sumLoss2[0] += predictedLoss * predictedLoss;
    sumLoss2[ptBin] += predictedLoss * predictedLoss;
    double samplingVariance =
        predictedLoss * (1 - predictedLoss) / acceptedDraws;
    sumMCVariance[0] += samplingVariance;
    sumMCVariance[ptBin] += samplingVariance;
    dataE.Fill(ee);
    for (int b = 0; b < 62; ++b) {
      double weight = eeBins[b] / acceptedDraws;
      poolEHistogram.AddBinContent(b, weight);
      double oldError = poolEHistogram.GetBinError(b);
      double newError = std::sqrt(oldError * oldError + weight * weight);
      poolEHistogram.SetBinError(b, newError);
    }
  }
  TH1D dataLoss("data_electron_loss", "", 3, ptEdges);
  TH1D predicted("predicted_electron_loss", "", 3, ptEdges);
  std::ofstream report("tof_pool.txt");
  report << "SL23e data and correction; same-track K/e/pi PID tuples.\n"
         << donorCount << " donors; momentum bins " << dp
         << " GeV/c; charge matched.\n"
         << "TOF m2: " << tofM2Min << " to " << tofM2Max
         << "; nHits >= " << hitsMin << ".\n"
         << "No TPC cut on donors. Sparse momentum bins are left uncovered.\n"
         << "KK <= " << kkMax << ", pipi >= " << pipiMin << "; reject ee < "
         << eeMin << " OR ee < " << eeFactor << "*KK.\n"
         << "Raw counts in " << massMin << " < M_KK < " << massMax
         << "; no background subtraction.\n\n";
  for (int b = 0; b <= 3; ++b) {
    double n = matched[b];
    if (n == 0)
      continue;
    double observed = rejected[b] / n;
    double prediction = sumLoss[b] / n;
    double observedError = std::sqrt(observed * (1 - observed) / n);
    // Spread among the measured pairs, plus the finite donor-sampling error.
    double momentumVariance = sumLoss2[b] - n * prediction * prediction;
    if (momentumVariance < 0)
      momentumVariance = 0;
    double predictionError = std::sqrt(momentumVariance + sumMCVariance[b]) / n;
    if (b > 0) {
      dataLoss.SetBinContent(b, 100 * observed);
      dataLoss.SetBinError(b, 100 * observedError);
      predicted.SetBinContent(b, 100 * prediction);
      predicted.SetBinError(b, 100 * predictionError);
    }
    if (b == 0)
      report << "All pT";
    else
      report << "pT bin " << b;
    report << ": covered " << matched[b] << "/" << total[b]
           << "; observed rejected " << rejected[b] << "; predicted rejected "
           << sumLoss[b] << "; observed loss " << 100 * observed << " +/- "
           << 100 * observedError << "%; predicted loss " << 100 * prediction
           << " +/- " << 100 * predictionError << "%\n";
  }
  report << "\npT bins: " << pairPtMin << "-" << ptSplit1 << ", " << ptSplit1
         << "-" << ptSplit2 << ", " << ptSplit2 << "-" << pairPtMax
         << " GeV/c.\n"
         << "Prediction errors describe the data momentum mixture and finite "
            "draws.\n"
         << "Data and prediction use the same momentum sample; their errors "
            "are correlated.\n"
         << "This is a conditional electron-cut test, not an absolute PID "
            "efficiency.\n"
         << "Generated-phi check: " << generatedPairs
         << " decays; joint parent pT/y sampled from " << nParents
         << " nominal mass-window data pairs before electron rejection,\n"
         << "sin^2(theta) helicity distribution, natural width 4.249 MeV.\n";
  for (int b = 1; b <= 3; ++b) {
    report << "Generated phi, pT " << ptEdges[b - 1] << "-" << ptEdges[b]
           << ": before electron cuts " << toyBefore.GetBinContent(b)
           << "; after electron cuts " << toyAfter.GetBinContent(b) << "; loss "
           << toyLoss.GetBinContent(b) << " +/- " << toyLoss.GetBinError(b)
           << "%\n";
  }
  printf("TOF coverage %d/%d; observed loss %d/%d; predicted loss %.3f/%d\n",
         matched[0], total[0], rejected[0], matched[0], sumLoss[0], matched[0]);

  PlotStyle::Options opt;
  opt.title = "Run-19 Au+Au";
  opt.subtitle = "#sqrt{s_{NN}} = 200 GeV";
  opt.note = "ZDCMB; incoherent region";
  opt.xlabel = "Pair p_{T} [GeV/c]";
  opt.ylabel = "Pairs rejected by electron PID [%]";
  opt.showEntries = false;
  opt.legendX1 = .54;
  opt.legendTextSize = .026;
  opt.YRange(-1, 45);
  TCanvas *c = PlotStyle::Draw(
      {{&dataLoss, "ULS data", "E0"},
       {&predicted, "TOF kaons, matched p and charge", "E0"},
       {&toyLoss, "Generated #phi #rightarrow K^{+}K^{-}, TOF PID", "E0"}},
      opt);
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(.030);
  text.DrawLatex(.18, .75, "1.00 < M_{KK} < 1.04 GeV/c^{2}; |y_{KK}| < 0.5");
  text.DrawLatex(.18, .69,
                 Form("#chi^{2}_{KK} #leq %.0f; #chi^{2}_{#pi#pi} #geq %.0f; "
                      "no extra track-p_{T} cut",
                      kkMax, pipiMin));
  text.DrawLatex(.18, .63,
                 Form("Electron PID: #chi^{2}_{ee} #geq %.0f AND "
                      "#chi^{2}_{ee} #geq %.0f#chi^{2}_{KK}",
                      eeMin, eeFactor));
  text.DrawLatex(.18, .57,
                 Form("TOF coverage: %d/%d ULS pairs", matched[0], total[0]));
  text.DrawLatex(
      .18, .51,
      "Generated #phi: data parent (p_{T}, y); sin^{2}(#theta) helicity decay");
  c->Print("figures/tof_electron_loss.png");
  c->Print("figures/tof_electron_loss.pdf");
  delete c;
  opt.note = "TOF-tagged kaons; 70-80% centrality";
  opt.YRange(-1, std::max(35., 1.8 * toyLoss.GetMaximum()));
  c = PlotStyle::Draw(
      {{&toyLoss, "#phi #rightarrow K^{+}K^{-} with TOF-kaon PID", "E0"}}, opt);
  text.DrawLatex(.18, .75, "1.00 < M_{KK} < 1.04 GeV/c^{2}; |y_{KK}| < 0.5");
  text.DrawLatex(.18, .69,
                 Form("#chi^{2}_{KK} #leq %.0f; #chi^{2}_{#pi#pi} #geq %.0f",
                      kkMax, pipiMin));
  text.DrawLatex(.18, .63,
                 Form("Electron PID: #chi^{2}_{ee} #geq %.0f AND "
                      "#chi^{2}_{ee} #geq %.0f#chi^{2}_{KK}",
                      eeMin, eeFactor));
  text.DrawLatex(.18, .57,
                 "Data parent (p_{T}, y); sin^{2}(#theta) helicity decay");
  c->Print("figures/tof_generated_phi.png");
  c->Print("figures/tof_generated_phi.pdf");
  delete c;
  opt.note = "ZDCMB; incoherent region";
  opt.xlabel = "#chi^{2}_{ee}";
  opt.ylabel = "Pairs / 2";
  opt.hasYRange = false;
  opt.XRange(0, 40);
  opt.headroom = 2.4;
  c = PlotStyle::Draw(
      {{&dataE, "ULS data", "E1"},
       {&poolEHistogram, "TOF kaons, matched p and charge", "E1"}},
      opt);
  text.DrawLatex(.18, .75, "1.00 < M_{KK} < 1.04 GeV/c^{2}; |y_{KK}| < 0.5");
  text.DrawLatex(.18, .69,
                 Form("#chi^{2}_{KK} #leq %.0f; #chi^{2}_{#pi#pi} #geq %.0f; "
                      "no electron rejection",
                      kkMax, pipiMin));
  text.DrawLatex(
      .18, .63,
      Form("0.3 < p_{T}^{KK} < 0.9 GeV/c; %d matched ULS pairs", matched[0]));
  c->Print("figures/tof_ee_response.png");
  c->Print("figures/tof_ee_response.pdf");
  delete c;
  output.cd();
  donorsMinus.Write();
  donorsPlus.Write();
  parentData.Write();
  dataParentPtY.Write();
  generatedParentPtY.Write();
  generated.Write();
  covered.Write();
  toyBefore.Write();
  toyAfter.Write();
  toyLoss.Write();
  match.Write();
  dataLoss.Write();
  predicted.Write();
  dataE.Write();
  poolEHistogram.Write();
  output.Close();
}
