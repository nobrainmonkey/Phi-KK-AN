= Candidate Selection and Raw Yields <sec:candidates>
This section combines the event selection (@sec:event_selection) with the track requirements of @sec:tracks to select coherent and incoherent $phi.alt$ candidates.
In both cases the raw yield is the number of ULS pairs in the mass window $1.00 < M_(K K) < 1.04 thin "GeV"\/c^2$, before background subtraction, acceptance and efficiency corrections, and luminosity normalization.
The cutflow tables start from the preselected ULS pairs of @sec:tracks and count pairs over the full mass range, $M_(K K) < 1.2 thin "GeV"\/c^2$, and in the mass window.
An $N-1$ distribution shows one variable with all other nominal requirements applied; no mass-window requirement is applied to these distributions.
In the figures, dashed lines mark the cut values of the plotted variable, or the mass window in mass spectra.

== Coherent and Incoherent Regions <sec:pt_regions>

In coherent production the photon couples to the nucleus as a whole, so the pair $p_T$ is limited by the nuclear size, $p_T^(K K) lt.tilde planck c \/ R_"Au"$; incoherent production on individual nucleons extends to much larger $p_T^(K K)$.
We define the two regions by a cut on the pair $p_T$ alone:
$ "coherent:" p_T^(K K) <= 0.1 thin "GeV"\/c, wide "incoherent:" 0.3 < p_T^(K K) < 0.9 thin "GeV"\/c. $ <eq:pt_regions>
The upper limit of the incoherent region follows from the electron-rejection study of @sec:incoh_electron.

@fig:pt_regions shows the pair $p_T$ before any kaon identification.
The ULS spectrum has a pronounced coherent peak below $0.1 thin "GeV"\/c$ that is absent in the LS pairs, and the two regions are well separated by the $p_T$ cut alone.
This separation depends only on the nuclear form factor and not on the produced vector meson, so it holds even though the sample is still dominated by pions at this stage.

#figure(image("../../plots/track/qa/pair/pair_pt_regions.pdf", width: 60%), caption: [Pair $p_T$ of ULS and LS pairs, with the coherent and incoherent regions of @eq:pt_regions.]) <fig:pt_regions>

== Coherent Candidates <sec:coh_selection>

=== Selection

@tab:coh_selection lists the nominal coherent selection, and @tab:coh_cutflow gives the cumulative ULS pair counts.
The pair, kaon-identification, and track-momentum requirements are motivated in @sec:coh_pair, @sec:coh_pid, and @sec:coh_trackpt.
@fig:coh_nm1_vz and @fig:coh_nm1_vr show the $N-1$ distributions of the vertex coordinates.

#figure(
  kind: table,
  table(columns: (1fr, 3fr), align: (left, left), stroke: none,
    table.hline(), table.header([Level], [Requirement]), table.hline(),
    [Event], [BBC veto; $N_"prim" = 2$; $abs(V_z) <= 100 thin "cm"$; $V_r <= 5 thin "cm"$],
    [Pair], [ULS; $p_T^(K K) <= 0.1 thin "GeV"\/c$; $abs(y_(K K)) < 0.5$],
    [Track quality], [Both tracks: $N_"hits" >= 20$],
    [PID], [$chi^2_(K K) <= 15$],
    [Track momentum], [Both tracks: $p_T <= 0.20 thin "GeV"\/c$],
    [Mass window], [$1.00 < M_(K K) < 1.04 thin "GeV"\/c^2$],
    table.hline(),
  ),
  caption: [Nominal coherent selection.],
) <tab:coh_selection>

#figure(
  kind: table,
  [#include "../../plots/signal/coherent/cutflow.typ"],
  caption: [Cumulative ULS pair counts for the coherent selection.],
) <tab:coh_cutflow>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/coherent/nm1_vz.pdf", width: 100%), caption: [$N-1$ distribution of $V_z$ for coherent candidates.]) <fig:coh_nm1_vz>],
  [#figure(image("../../plots/signal/coherent/nm1_vr.pdf", width: 100%), caption: [$N-1$ distribution of $V_r$ for coherent candidates.]) <fig:coh_nm1_vr>],
)

=== Pair Transverse Momentum and Rapidity <sec:coh_pair>

The requirement $p_T^(K K) <= 0.1 thin "GeV"\/c$ selects the coherent region of @eq:pt_regions (@fig:coh_nm1_pair_pt).

The pair-rapidity window is chosen from the mass spectra in four signed-rapidity intervals (@fig:coh_rap_0 to @fig:coh_rap_3), with all other nominal requirements applied except the track-momentum requirement.
Mass-window pairs occur only in the three inner intervals, $-0.5 <= y_(K K) <= 0.5$; the interval $0.5 < y_(K K) <= 1$ contains none.
The full rapidity distribution (@fig:coh_rapidity_full) shows that no mass-window pair lies at $abs(y_(K K)) > 0.5$ on either side.
We therefore use $abs(y_(K K)) < 0.5$.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/coherent/nm1_pair_pt.pdf", width: 100%), caption: [$N-1$ distribution of $p_T^(K K)$ for coherent candidates.]) <fig:coh_nm1_pair_pt>],
  [#figure(image("../../plots/signal/coherent_rap_range/pair_rapidity.pdf", width: 100%), caption: [Pair rapidity inside and outside the mass window, before the track-momentum requirement.]) <fig:coh_rapidity_full>],
  [#figure(image("../../plots/signal/coherent_rap_range/mass_rapidity_bin0.pdf", width: 100%), caption: [ULS mass spectrum for $-0.5 <= y_(K K) < -0.25$, without the rapidity and track-momentum requirements.]) <fig:coh_rap_0>],
  [#figure(image("../../plots/signal/coherent_rap_range/mass_rapidity_bin1.pdf", width: 100%), caption: [Same as @fig:coh_rap_0 for $abs(y_(K K)) <= 0.25$.]) <fig:coh_rap_1>],
  [#figure(image("../../plots/signal/coherent_rap_range/mass_rapidity_bin2.pdf", width: 100%), caption: [Same as @fig:coh_rap_0 for $0.25 < y_(K K) <= 0.5$.]) <fig:coh_rap_2>],
  [#figure(image("../../plots/signal/coherent_rap_range/mass_rapidity_bin3.pdf", width: 100%), caption: [Same as @fig:coh_rap_0 for $0.5 < y_(K K) <= 1$.]) <fig:coh_rap_3>],
)

=== Kaon Identification <sec:coh_pid>

In the coherent region the decay kaons have momenta below about $0.3 thin "GeV"\/c$ (@fig:coh_track_p), where the kaon and pion $dif E slash dif x$ bands are well separated (@fig:dedx_p).
A loose requirement, $chi^2_(K K) <= 15$, is therefore applied.
@fig:coh_nm1_pid shows its $N-1$ distribution.

#figure(image("../../plots/signal/coherent/nm1_pid.pdf", width: 60%), caption: [$N-1$ distribution of $chi^2_(K K)$ for coherent candidates.]) <fig:coh_nm1_pid>

=== Track Transverse Momentum <sec:coh_trackpt>

The pair $p_T^(K K)$ is the magnitude of the vector sum of the two track transverse momenta, so a low pair $p_T$ does not imply that each track has low $p_T$.
@fig:coh_nm1_trackpt shows the larger of the two track momenta, and @fig:coh_mass_trackpt the mass spectrum before and after requiring $p_T <= 0.20 thin "GeV"\/c$ for both tracks.
For nearly back-to-back tracks, the track-momentum requirement is correlated with $M_(K K)$ and suppresses the higher-mass population, while it keeps all mass-window pairs.
We adopt $p_T <= 0.20 thin "GeV"\/c$ for both tracks.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/coherent/nm1_track_pt_max.pdf", width: 100%), caption: [$N-1$ distribution of the larger track $p_T$ for coherent candidates.]) <fig:coh_nm1_trackpt>],
  [#figure(image("../../plots/signal/coherent/mass_track_pt_selection.pdf", width: 100%), caption: [Coherent mass spectrum before and after the track-momentum requirement.]) <fig:coh_mass_trackpt>],
)

=== Mass Spectrum and Raw Yield

@fig:coh_mass shows the mass spectrum after the nominal coherent selection.
The number of ULS pairs in the mass window is the coherent raw yield,
$ N_"raw"^"coh" = 61. $ <eq:nraw_coh>
@fig:coh_pair_pt, @fig:coh_pair_y, and @fig:coh_pair_phi show the pair distributions of the mass-window pairs; their track distributions are collected in @sec:app_coh.

#figure(image("../../plots/signal/coherent/mass.pdf", width: 75%), caption: [Invariant-mass spectrum after the nominal coherent selection.]) <fig:coh_mass>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/coherent/pair_pt.pdf", width: 100%), caption: [$p_T^(K K)$ of the coherent mass-window pairs.]) <fig:coh_pair_pt>],
  [#figure(image("../../plots/signal/coherent/pair_rapidity.pdf", width: 100%), caption: [$y_(K K)$ of the coherent mass-window pairs.]) <fig:coh_pair_y>],
  [#figure(image("../../plots/signal/coherent/pair_phi.pdf", width: 100%), caption: [Pair azimuth of the coherent mass-window pairs.]) <fig:coh_pair_phi>],
)


== Incoherent Candidates <sec:incoh_selection>

=== Selection

@tab:incoh_selection lists the nominal incoherent selection, and @tab:incoh_cutflow gives the cumulative ULS pair counts.
The event, vertex, and track-quality requirements and the rapidity interval are those of the coherent selection, and there is no track-momentum requirement beyond the preselection.
The higher parent $p_T$ shifts the daughter momenta upward, where the TPC response bands approach each other, so the PID is chosen in a separate study (@sec:incoh_pid, @sec:incoh_electron).
@fig:incoh_nm1_pt and @fig:incoh_nm1_y show the $N-1$ distributions of the pair requirements.

#figure(
  kind: table,
  table(columns: (1fr, 3fr), align: (left, left), stroke: none,
    table.hline(), table.header([Level], [Requirement]), table.hline(),
    [Event], [BBC veto; $N_"prim" = 2$; $abs(V_z) <= 100 thin "cm"$; $V_r <= 5 thin "cm"$],
    [Pair], [ULS; $0.3 < p_T^(K K) < 0.9 thin "GeV"\/c$; $abs(y_(K K)) < 0.5$],
    [Track quality], [Both tracks: $N_"hits" >= 20$],
    [PID], [$chi^2_(K K) <= 10$; $chi^2_(pi pi) >= 15$; no electron or proton rejection],
    [Track momentum], [No requirement beyond the preselection],
    [Mass window], [$1.00 < M_(K K) < 1.04 thin "GeV"\/c^2$],
    table.hline(),
  ),
  caption: [Nominal incoherent selection.],
) <tab:incoh_selection>

#figure(
  kind: table,
  [#include "../../plots/signal/incoherent/cutflow.typ"],
  caption: [Cumulative ULS pair counts for the incoherent selection.],
) <tab:incoh_cutflow>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/incoherent/nm1_pair_pt.pdf", width: 100%), caption: [$N-1$ distribution of $p_T^(K K)$ for incoherent candidates.]) <fig:incoh_nm1_pt>],
  [#figure(image("../../plots/signal/incoherent/nm1_pair_y.pdf", width: 100%), caption: [$N-1$ distribution of $y_(K K)$ for incoherent candidates.]) <fig:incoh_nm1_y>],
)

=== Kaon Identification <sec:incoh_pid>

The PID requirements are chosen from the ULS mass-window pairs that pass all nominal requirements except the PID.
Their single-track $n sigma$ and pair $chi^2$ distributions are shown in @fig:incoh_nsig_k to @fig:incoh_chi2_pp, before and after the nominal PID; the error bars of the single-track distributions account for the correlation of the two tracks of a pair.
The corrected kaon response contains a concentration near zero, together with a broad population away from the kaon expectation; the pion, electron, and proton hypotheses show how these candidates populate the other response bands.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/nsigma_kaon.pdf", width: 100%), caption: [$n sigma_K^"corr"$ of tracks in incoherent mass-window pairs, before and after the nominal PID.]) <fig:incoh_nsig_k>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/nsigma_electron.pdf", width: 100%), caption: [Same as @fig:incoh_nsig_k for $n sigma_e$.]) <fig:incoh_nsig_e>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/nsigma_pion.pdf", width: 100%), caption: [Same as @fig:incoh_nsig_k for $n sigma_pi$.]) <fig:incoh_nsig_pi>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/nsigma_proton.pdf", width: 100%), caption: [Same as @fig:incoh_nsig_k for $n sigma_p$.]) <fig:incoh_nsig_p>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/chi2_kk_full.pdf", width: 100%), caption: [$chi^2_(K K)$ of incoherent mass-window pairs, before and after the nominal PID.]) <fig:incoh_chi2_kk>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/chi2_ee_full.pdf", width: 100%), caption: [Same as @fig:incoh_chi2_kk for $chi^2_(e e)$.]) <fig:incoh_chi2_ee>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/chi2_pipi_full.pdf", width: 100%), caption: [Same as @fig:incoh_chi2_kk for $chi^2_(pi pi)$.]) <fig:incoh_chi2_pipi>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/chi2_pp_full.pdf", width: 100%), caption: [Same as @fig:incoh_chi2_kk for $chi^2_(p p)$.]) <fig:incoh_chi2_pp>],
)

@fig:incoh_kaon_plane shows the corrected kaon response of track 1 against that of track 2.
The central population lies near the origin, separated from the broad population at negative $n sigma_K^"corr"$.
The requirement
$ chi^2_(K K) = (n sigma_(K,1)^"corr")^2 + (n sigma_(K,2)^"corr")^2 <= 10 $
selects a circle of radius $sqrt(10)$ around the kaon expectation, which retains the central population without extending far into the surrounding tails.
@tab:incoh_kaon_choice compares nearby boundaries: tightening to 8 removes 11 mass-window pairs, while loosening to 15 adds 9 mass-window pairs and 70 pairs outside the window.
We therefore use $chi^2_(K K) <= 10$.

#figure(
  image("../../plots/signal/incoherent/pid_study/figures/nsigma_kaon_pair.pdf", width: 100%),
  caption: [$n sigma_(K,1)^"corr"$ versus $n sigma_(K,2)^"corr"$ for incoherent mass-window pairs before PID, with circles at $chi^2_(K K) = 10$ (solid) and 15 (dashed).],
) <fig:incoh_kaon_plane>

#figure(
  kind: table,
  table(columns: (1.2fr, 1fr, 1fr, 1fr), align: (left, right, right, right), stroke: none,
    table.hline(),
    table.header([Kaon requirement], [All masses], [Mass window], [Outside window]),
    table.hline(),
    [$chi^2_(K K) <= 8$], [505], [287], [218],
    [$chi^2_(K K) <= 10$], [540], [298], [242],
    [$chi^2_(K K) <= 15$], [619], [307], [312],
    table.hline(),
  ),
  caption: [ULS pair counts for three kaon requirements, without pion or electron rejection.],
) <tab:incoh_kaon_choice>

All mass-window pairs already satisfy $chi^2_(pi pi) >= 15$; this mild pion requirement is kept for background control and removes 6 pairs above the mass window.

=== Electron Rejection <sec:incoh_electron>

At the larger daughter momenta of the incoherent sample, a kaon-compatible TPC response can also overlap the electron hypothesis.
@fig:incoh_electron_plane shows $chi^2_(K K)$ versus $chi^2_(e e)$ in the mass window; some candidates that pass the kaon requirement would fail one of the electron requirements
$ chi^2_(e e) >= 10, wide chi^2_(e e) >= 2 chi^2_(K K). $ <eq:electron_cuts>
The mass spectrum retains a clear enhancement near the $phi.alt$ mass without electron rejection (@fig:incoh_mass_ecut).
Applying both electron requirements removes 13 of the 298 mass-window pairs, $(4.4 plus.minus 1.2)%$, including pairs in the $phi.alt$ peak (@fig:incoh_mass_erej).
Since the $phi.alt$ is close to the $K K$ threshold, electron pairs reconstructed under the kaon mass hypothesis can also populate this region.
We therefore test with TOF-identified kaons whether this loss is expected from the kaon $dif E slash dif x$ response itself.

#figure(
  image("../../plots/signal/incoherent/pid_study/figures/kk_vs_ee.pdf", width: 100%),
  caption: [$chi^2_(K K)$ versus $chi^2_(e e)$ in the incoherent mass window before PID, with the kaon requirement (vertical line) and the electron requirements of @eq:electron_cuts.],
) <fig:incoh_electron_plane>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/mass_electron_cuts.pdf", width: 100%), caption: [Incoherent mass spectrum with and without the electron requirements.]) <fig:incoh_mass_ecut>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/mass_rejected_electron.pdf", width: 100%), caption: [Mass spectrum of the pairs removed by the electron requirements.]) <fig:incoh_mass_erej>],
)

==== Prediction from TOF-Identified Kaons <sec:incoh_tof_pool>

The prediction asks how often true kaons with the measured daughter momenta would fail the electron requirements.
It uses kaon candidates from the peripheral sample (@sec:nsigk_tof) with $0.22 < m^2 < 0.28 thin "GeV"^2\/c^4$, $N_"hits" >= 20$, and $0.1 <= p < 1.0 thin "GeV"\/c$, without any TPC $n sigma$ requirement: 5.2 million tracks, binned by charge and in $0.01 thin "GeV"\/c$ of momentum.
Their $n sigma_K$ is corrected in the same way as for the UPC candidates.

For each mass-window pair, one donor track is drawn for each daughter, matched in charge and momentum bin.
A donor supplies its measured $(n sigma_K^"corr", n sigma_e, n sigma_pi)$ together, which preserves their correlations and the measured response tails; the two daughters are drawn independently.
For each pair, 4000 donor pairs are drawn.
Among the draws that pass the preselection $chi^2_(K K) < 100$ and the nominal kaon and pion requirements, the fraction that fails either electron requirement gives the rejection probability of pair $i$ and the expected number of rejected pairs,
$ P_i^"loss" = N_i^"rejected" / N_i^"accepted", wide N_"predicted"^"rejected" = sum_i P_i^"loss". $
Pairs with fewer than 30 kaon donors in either daughter's momentum bin, or fewer than 30 accepted draws, are excluded from both the data and the prediction.

The pool covers 268 of the 298 mass-window pairs.
In this common sample the electron requirements reject 13 pairs, against 14.6 predicted: a loss of $(4.9 plus.minus 1.3)%$ in data and $(5.4 plus.minus 1.2)%$ predicted.
The prediction uncertainty includes the spread of the measured momenta and the finite number of draws.
The dependence on pair $p_T$ also agrees: nearly all of the rejection occurs at $0.7 < p_T^(K K) < 0.9 thin "GeV"\/c$ (@tab:incoh_tof_loss, @fig:incoh_tof_loss).

#figure(
  kind: table,
  table(columns: (1.3fr, 1fr, 1fr, 1fr), align: (left, right, right, right), stroke: none,
    table.hline(),
    table.header([$p_T^(K K)$ [GeV/$c$]], [Covered pairs], [Observed loss], [Predicted loss]),
    table.hline(),
    [$0.3$--$0.5$], [51 / 73], [0], [0.00],
    [$0.5$--$0.7$], [118 / 126], [1], [1.27],
    [$0.7$--$0.9$], [99 / 99], [12], [13.30],
    table.hline(stroke: 0.5pt),
    [Full interval], [268 / 298], [13], [14.57],
    table.hline(),
  ),
  caption: [Observed and predicted electron-requirement losses, in pairs, for the mass-window pairs covered by the TOF-identified kaons.],
) <tab:incoh_tof_loss>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/tof_electron_loss.pdf", width: 100%), caption: [Electron-requirement loss versus $p_T^(K K)$ for data, the TOF-kaon prediction, and generated $phi.alt$ decays.]) <fig:incoh_tof_loss>],
  [#figure(image("../../plots/signal/incoherent/pid_study/figures/tof_ee_response.pdf", width: 100%), caption: [Measured and predicted $chi^2_(e e)$ distributions of the covered mass-window pairs.]) <fig:incoh_tof_ee>],
)

As a kinematic cross-check, 2 million $phi.alt -> K^+ K^-$ decays are generated with the parent $(p_T, y)$ drawn jointly from the 298 mass-window pairs before electron rejection (including any background in the window), a Breit-Wigner mass restricted to the mass window, a uniform parent azimuth, and a $sin^2 theta$ decay distribution in the helicity frame.
The daughters receive donor responses at their generated momenta and charges.
Unlike the prediction above, which uses the measured daughter momenta, this check generates new daughter momenta from the parent kinematics; its loss, evaluated in the part of the generated sample covered by the TOF-identified kaons, is shown with the data and the prediction in @fig:incoh_tof_loss.

The transfer assumes that charge- and momentum-matched TOF-identified kaons represent the TPC response of the decay kaons.
The agreement with the observed loss, including its concentration at high pair $p_T$, supports kaon-response overlap as the explanation for the removed candidates in the $phi.alt$ peak.
We therefore apply no electron rejection.

The loss grows with $p_T^(K K)$ as the daughters reach the momentum range where the kaon and electron responses overlap: in $0.7 < p_T^(K K) < 0.9 thin "GeV"\/c$ the prediction from the measured daughters gives $(13.4 plus.minus 2.9)%$ and the generated decays $19.5%$.
This rising loss sets the upper limit $p_T^(K K) < 0.9 thin "GeV"\/c$ of the incoherent region (@eq:pt_regions).

=== Mass Spectrum and Raw Yield

@fig:incoh_mass shows the mass spectrum after the nominal incoherent selection.
The number of ULS pairs in the mass window is the incoherent raw yield,
$ N_"raw"^"incoh" = 298. $ <eq:nraw_incoh>
@fig:incoh_pair_pt, @fig:incoh_pair_y, and @fig:incoh_pair_phi show the pair distributions of the mass-window pairs; their track distributions are collected in @sec:app_incoh.

#figure(image("../../plots/signal/incoherent/mass.pdf", width: 75%), caption: [Invariant-mass spectrum after the nominal incoherent selection.]) <fig:incoh_mass>

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/signal/incoherent/pair_pt.pdf", width: 100%), caption: [$p_T^(K K)$ of the incoherent mass-window pairs.]) <fig:incoh_pair_pt>],
  [#figure(image("../../plots/signal/incoherent/pair_rapidity.pdf", width: 100%), caption: [$y_(K K)$ of the incoherent mass-window pairs.]) <fig:incoh_pair_y>],
  [#figure(image("../../plots/signal/incoherent/pair_phi.pdf", width: 100%), caption: [Pair azimuth of the incoherent mass-window pairs.]) <fig:incoh_pair_phi>],
)
