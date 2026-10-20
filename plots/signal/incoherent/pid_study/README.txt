Incoherent PID study

Edit cuts.h for the input files, event cuts, pair cuts, mass window and PID cuts.
Both inputs and the kaon correction are SL23e.

Run in this directory, in this order:
  root -l -b -q pid_plots.C
  root -l -b -q tof_pool.C

pid_plots.C
  Reads uls_data directly from the original pair ROOT file.
  Makes KK versus ee, pipi and pp plots, the four 1-D chi2 plots, and the
  kaon/electron/pion/proton nSigma plots. Kaon nSigma and KK are corrected.
  All PID distributions use 1.00 < M_KK < 1.04 GeV/c^2.
  In 1-D, blue is before PID and orange has KK and pion cuts applied,
  before electron cuts. The 2-D plots use ROOT COLZ, with cut lines overlaid.
  The nSigma plots include both tracks. Errors account for the two tracks
  belonging to the same counted pair.
  The mass comparison holds KK and pion cuts fixed, then turns BOTH electron
  requirements on/off. Yields are raw counts inside the mass window.
  Binning and display ranges are written next to the histogram definitions.
  Also draws corrected nSigmaK(track 1) versus nSigmaK(track 2), with the
  nominal KK <= 10 circle, and compares pion rejection inside/above the window.
  The pion mass comparison and electron-rejected mass spectrum use raw counts.
  Nominal PID: corrected KK <= 10, pipi >= 15, no electron or proton rejection.

tof_pool.C
  Reads the SL23e TOF-kaon file directly. No intermediate pool file is needed.
  Pool settings, momentum bins, number of decays and draws are at the top.
  Keeps each donor's K/e/pi responses together and matches momentum and charge.
  First generates phi -> KK; then samples at the actual data daughter momenta.
  Sparse momentum bins stay uncovered.
  The data comparison uses the SAME covered raw mass-window pairs for both
  curves. These raw counts include any background in the window.
  The generated sample takes parent pT and y together from one randomly
  chosen nominal mass-window data pair, before electron rejection.
  It samples all 298 selected pairs, preserving their pT/y correlation.
  Mass and decay angles are still generated; daughters receive TOF donors.
  Independent random streams keep this check separate from the direct
  measured-daughter prediction. Both use the same three pair-pT bins.

Outputs
  figures/             PNG and PDF plots, using the existing PlotStyle.h.
  pid_plots.root       All PID histograms, scatter points, mass histograms,
                       window yields and the selected-pair tree.
  window_yields.txt    Raw counts with and without the electron requirements.
  tof_pool.root        Pool occupancy, parent_data tree, data/generated parent
                       pT/y histograms, generated-phi histograms, matched-pair
                       predictions and the data/predicted electron losses.
  tof_pool.txt         Pool coverage and predicted/observed raw-count losses.

No fitting or sideband subtraction. The kaon-cut table compares 8, 10 and 15.
After changing cuts.h, rerun both macros.
