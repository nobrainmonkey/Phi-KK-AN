= Dataset and Trigger

This analysis uses the Run-19 Au+Au dataset collected at $sqrt(s_(N N)) = 200  "GeV"$.
Run-19 is the first Au+Au dataset collected after the full installation of the inner Time Projection Chamber (iTPC) upgrade.
The improved inner-sector acceptance allows charged tracks to be reconstructed down to transverse momenta of approximately $p_T tilde.op 60 "MeV"\/c$.
This is particularly important for this analysis because the decay kaons are soft.
The $phi.alt$ carries very little momentum, especially in coherent production, and in the $phi.alt$ rest frame each kaon has a momentum
$ p^* = sqrt((m_(phi.alt) / 2)^2 - m_K^2) approx 127  "MeV"\/c, $ <eq:kaon_momentum>
so in the laboratory frame the kaon transverse momenta are of order $100  "MeV"\/c$.

The dataset used in this analysis has the following configuration:
- Production library: `SL23e`;
- Data format: `MuDst`.

No dedicated UPC trigger was available for this dataset.
Instead, events are selected from the Zero Degree Calorimeter (ZDC) minimum-bias (zdcmb) trigger sample.

In UPCs of heavy nuclei, virtual photon exchanges can electromagnetically excite one or both nuclei, most commonly through the Giant Dipole Resonance (GDR).
The excited nucleus subsequently de-excites through the emission of one or more forward neutrons, which remain close to beam rapidity and are detected by the ZDCs.
The zdcmb trigger therefore provides a baseline sample of UPC events accompanied by mutual electromagnetic dissociation.

The trigger configuration is summarized below:
- Offline trigger ID: `700001`;
- Trigger requirements: ZDC-East, ZDC-West, and ZDC TAC coincidence, with no minimum Time-of-Flight (TOF) requirement;
- Number of recorded events: $138.38$ million.

The absence of a TOF requirement is essential for this analysis: as discussed above, the decay kaons are largely too soft to reach the TOF, so a trigger requiring TOF hits would remove most of the coherent signal.

The analyzed dataset was recorded with ZDC coincidence rates ranging from $1.5 "kHz"$ to $5.0 "kHz"$, as shown in @fig:run19_zdcmb_rate. This corresponds to a clean, low-luminosity running condition. At this collision rate, the mean time separation between collisions is $200$–$700 thin upright(mu)"s"$, which is much longer than the TPC drift time of $~40 thin upright(mu)"s"$.

#figure(image("../../plots/events/all_bbc_zdc_coin.pdf", width: 75%), caption: [ZDC and BBC coincidence rates for the Run-19 Au+Au zdcmb trigger.])<fig:run19_zdcmb_rate>

In the Run-19 trigger configuration, there is no cut on the ZDC ADC values, so no ZDC selection efficiency is applied in this analysis. The unattenuated ZDC ADC values that saturate the 12-bit ADC range are placed in the overflow bin, which appears as a spike at the right edge of the ZDC ADC distribution in @fig:run19_zdcua.
#figure(image("../../plots/events/all_zdc_ua.pdf", width: 75%), caption: [ZDC unattenuated ADC values for the Run-19 Au+Au zdcmb trigger.])<fig:run19_zdcua>



== Luminosity

The integrated luminosity for the zdcmb trigger is obtained from the #link("https://www.star.bnl.gov/protected/common/common2019/trigger2019/lumiAu200GeV/")[STAR luminosity database for Run-19], where it is determined from the ZDC coincidence rate and the ZDC coincidence cross section. The corresponding net cross section, extracted as the slope of the live ZDC coincidences in each run, $chevron.l R_"ZDC" chevron.r times t_"run" times "livetime"$, versus the database luminosity, is $12.2 "b"$ (@fig:run19_lumi_vs_zdc_rate); here $t_"run"$ is the run duration. The number of recorded events versus the integrated luminosity for each run is shown in @fig:run19_lumi_vs_events and exhibits a strong linear correlation, indicating that the luminosity determination is consistent for all analyzed runs. The integrated luminosity for the zdcmb trigger before bad-run removal is $L_"int" = 13.2666 thin upright(mu)"b"^(-1)$, with a prescale factor of 1.

#figure(
  image("../../plots/events/zdc_rate_vs_lumi.pdf", width: 100%),
  caption: [Run-by-run integrated luminosity compared to the ZDC coincidence rate for the Run-19 Au+Au zdcmb trigger. Left: live ZDC coincidences per run versus integrated luminosity, with a proportional fit. Right: their ratio, the implied ZDC coincidence cross section, versus the mean ZDC coincidence rate.],
)<fig:run19_lumi_vs_zdc_rate>

#figure(
  image("../../plots/events/lumi_vs_events.pdf", width: 100%),
  caption: [Run-by-run integrated luminosity compared to the number of events recorded for the Run-19 Au+Au zdcmb trigger. Left: events per run versus integrated luminosity, with a proportional fit. Right: events per unit luminosity for each run.],
)<fig:run19_lumi_vs_events>

== Bad Runs

Prior to applying event-level selections, runs with known detector-quality issues are removed from the analysis.
#link("https://drupal.star.bnl.gov/STAR/pwg/common/bes-ii-run-qa/Run19-AuAu-200-GeV-COL")[The standard STAR BES-II Run QA bad-run list for Run-19 Au+Au collisions at $sqrt(s_(N N)) = 200  "GeV"$] is used to exclude periods with degraded TPC performance.
@tab:bad_run_lumi lists the run-by-run integrated luminosities for the bad runs.

#figure(
  table(
    columns: 2,
    align: center,
    stroke: none,
    table.hline(),
    table.header([Run number], [Integrated luminosity \[$upright(mu)"b"^(-1)$\]]),
    table.hline(),
    [20191005], [0.145248],
    [20191015], [0.000152],
    [20192001], [0.040369],
    [20193001], [0.404860],
    [20193019], [0.000067],
    table.hline(),
    [Total], [0.590695],
    table.hline(),
  ),
  caption: [Integrated luminosity associated with the BES-II QA bad runs removed from the Run-19 Au+Au dataset.],
) <tab:bad_run_lumi>

The nominal integrated trigger luminosity after bad-run removal for this analysis is
$ L_"nominal" = 12.6759 thin upright(mu)"b"^(-1) $ <eq:lumi_nominal>
