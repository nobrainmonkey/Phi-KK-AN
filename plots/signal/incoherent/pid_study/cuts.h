#ifndef INCOHERENT_PID_CUTS_H
#define INCOHERENT_PID_CUTS_H

// Shared by both macros. Change the cuts here, then rerun both macros.
const char *pairFile = "../../../../data/phi-kk-mudst.root"; // SL23e
const char *kaonFile = "../../../../data/kaon_mudst.root";   // SL23e
const char *correctionFile = "../../../track/PID/nsigmakaon_mean_vs_p.root";

const double pairPtMin = 0.3, pairPtMax = 0.9;
const double pairYMax = 0.5;
const double vzMax = 100, vrMax = 5;
const int hitsMin = 20;
// BBC east = west = 0; exactly two primary tracks; no extra track-pT cut.

const double massMin = 1.00, massMax = 1.04; // GeV/c^2
const double kkMax = 10;                     // Corrected kaon chi2
const double pipiMin = 15;
const double eeMin = 10, eeFactor = 2; // ee >= eeMin AND ee >= eeFactor * kk

#endif
