= Track Selection and Particle Identification <sec:tracks>
After the event selection, each candidate event contains a single primary vertex with exactly two primary tracks.
This section characterizes these tracks through their tracking quality, kinematics, and specific energy loss ($dif E slash dif x$) in the TPC, and establishes the track-quality requirement and the kaon identification used by the candidate selections of @sec:candidates.
Since the final state is reconstructed from only two tracks in an otherwise empty detector, the track sample is already relatively clean.
We therefore apply a minimal set of track-quality requirements, chosen to retain as many good tracks as possible while rejecting poorly reconstructed ones.

Particle identification is based on the TPC $dif E slash dif x$.
For a track $i$ and a particle hypothesis $X in {e, pi, K, p}$, the normalized deviation of the measured energy loss from the expected one is
$ n sigma_X^((i)) = 1 / sigma_(dif E slash dif x) ln (dif E slash dif x)_"meas"^((i)) / (dif E slash dif x)_X^((i)), $ <eq:nsigma>
where $(dif E slash dif x)_X^((i))$ is the Bichsel expectation for hypothesis $X$ at the momentum of track $i$ and $sigma_(dif E slash dif x)$ is the $dif E slash dif x$ resolution.
For the two tracks of a pair, the compatibility with the hypothesis that both tracks are of type $X$ is quantified by
$ chi^2_(X X) = (n sigma_X^((1)))^2 + (n sigma_X^((2)))^2, $ <eq:chi2_pair>
which is evaluated for the kaon ($chi^2_(K K)$), pion ($chi^2_(pi pi)$), electron ($chi^2_(e e)$), and proton ($chi^2_(p p)$) hypotheses.

The analysis tuple contains the pairs that pass the following preselection:
- zdcmb trigger (700001), bad runs removed;
- exactly two primary tracks associated with the pair's vertex;
- pair $dif E slash dif x$ compatibility with the kaon--kaon hypothesis, $chi^2_(K K) < 100$ (@eq:chi2_pair), computed from the uncorrected $n sigma_K$ (@sec:nsigk_tof);
- track $p_T < 1 thin "GeV"\/c$ and pair $p_T^(K K) < 1 thin "GeV"\/c$;
- pair invariant mass, computed with the kaon mass hypothesis, $M_(K K) < 1.2 thin "GeV"\/c^2$;
- no TOF requirement.

Unless stated otherwise, the distributions in this section also require the BBC veto, $N_"prim" = 2$, and $abs(V_z) <= 100 thin "cm"$ (@sec:event_selection).
The two tracks form an unlike-sign (ULS) pair, the signal sample, or a like-sign (LS) pair, which is kept as a background reference.


== Track Quality              // nHitsFit, nHitsDedx, DCA, |η| — common

=== Primary-Vertex Reconstruction <sec:minuit_vertex>

Primary vertices in this production (SL23e) are reconstructed with the Minuit vertex finder (`StMinuitVertexFinder`), which fits the vertex position from the distances of closest approach of TPC global tracks.
Only global tracks passing `StMinuitVertexFinder::accept()` enter the fit:
```cpp
bool StMinuitVertexFinder::accept(StTrack* track) const
{
  return (track &&
          track->flag() >= 0 &&
          track->fitTraits().numberOfFitPoints() >= mMinNumberOfFitPointsOnTrack &&
          !track->topologyMap().trackFtpc() &&
          finite(track->length()) &&
          track->geometry()->helix().valid());
}
```
The cut values come from the `Calibrations/tracker/PrimaryVertexCuts` database table (@tab:minuit_cuts).

#figure(
  table(
    columns: 3,
    align: (left, center, left),
    stroke: none,
    table.hline(),
    table.header([Parameter], [Value], [Requirement]),
    table.hline(stroke: 0.5pt),
    [`MinNumberOfFitPointsOnTrack`], [20], [global-track fit points $N_"fit"^"glob" >= 20$],
    [`RImpactMax`], [6 cm], [transverse distance of the global helix to the beamline $< 6 thin "cm"$],
    [`DcaZMax`], [6 cm], [track $z$ at closest approach within $6 thin "cm"$ of the vertex seed],
    [`MinTrack`], [2], [at least two accepted tracks per vertex],
    [`ZMin`, `ZMax`], [$-200$, $205$ cm], [range of the vertex search along $z$],
    [`MinTrackPt`], [0], [not applied],
    [`MinFracOfPossFitPointsOnTrack`], [0], [not applied],
    table.hline(),
  ),
  caption: [Track and vertex requirements of the Minuit vertex finder for the Run-19 Au+Au production.],
) <tab:minuit_cuts>

Since each selected vertex has exactly two tracks, both must pass $N_"fit"^"glob" >= 20$.
The primary refit adds the vertex as a fit point, giving the sharp edge at $N_"fit" = 21$ (@fig:nhitsfit_lin, @fig:nhitsfit_log); the small tail below it consists of tracks that lost hits in this refit.
We therefore adopt
$ N_"hits" >= 20 $ <eq:nhits_cut>
as the nominal track-quality selection, with $N_"hits" = N_"fit" - 1$ the number of TPC fit points.
Its efficiency is evaluated with embedding.

The vertex finder applies *no minimum $p_T$ requirement*: `MinTrackPt` is zero and is not read by the Minuit finder.
The only implicit low-$p_T$ limit is that strongly curling tracks may not reach 20 fit points; this is included in the embedding efficiency.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/track/qa/track/track_nhitsfit.pdf", width: 100%), caption: [Number of fit points of primary tracks in ULS and LS pairs.]) <fig:nhitsfit_lin>],
  [#figure(image("../../plots/track/qa/track/track_nhitsfit_log.pdf", width: 100%), caption: [Same as @fig:nhitsfit_lin on a logarithmic scale.]) <fig:nhitsfit_log>],
)

=== Inner TPC Upgrade and the Hit Requirement <sec:itpc>

The full iTPC (@sec:dataset) increased the inner-sector pad rows from 13 to 40 and the maximum number of points per track from 45 to 72.
The iTPC is fully operational in this sample: $N_"fit"$ extends well beyond 45 (@fig:nhitsfit_lin), $N_"poss"$ peaks at the full inner-sector and full-TPC lengths (@fig:nhitsposs), and nearly all tracks have iTPC hits (@fig:itpc_rows).

This makes $N_"hits" >= 20$ a mild requirement here, although it is tighter than the 15 hits common in STAR UPC analyses.
Many tracks, mostly at large $abs(eta)$, leave through the endcap without reaching the outer sectors; with 40 inner rows they still collect 20 hits, which produces the edge at 20 in @fig:itpc_rows, whereas with 13 inner rows they could not.

The fit-point threshold is not new: `MinNumberOfFitPointsOnTrack` is 20 also in 2012, 2014, and 2018.
What changed is `MinTrack`, which was 5 in those years and is 2 here, so the Minuit finder can now reconstruct exclusive two-track events.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/track/qa/track/track_nhitsposs.pdf", width: 100%), caption: [Number of possible points of primary tracks in ULS and LS pairs.]) <fig:nhitsposs>],
  [#figure(image("../../plots/track/qa/track/track_itpc_rows.pdf", width: 100%), caption: [Number of iTPC pad rows hit per track.]) <fig:itpc_rows>],
)

== Particle Identification

@fig:dedx_p shows $dif E slash dif x$ versus momentum for all tracks (ULS and LS).
The curves give the expected $dif E slash dif x$ for each hypothesis, taken as the mean $dif E slash dif x$ of tracks with $abs(n sigma_X) < 0.3$, i.e. the reference used in @eq:nsigma.
At low momentum the kaon band is well separated from the pion band, while the two merge above $p approx 1 thin "GeV"\/c$.

The sample is nevertheless dominated by pions.
The preselection $chi^2_(K K) < 100$ only limits each track to $abs(n sigma_K) < 10$ (@fig:nsigk_p), which at these momenta still accepts most of the pion band.
This is also visible in the one-dimensional distributions: $n sigma_pi$ peaks sharply at zero, while $n sigma_K$ is dominated by a broad pion contribution at negative values (@fig:nsig_1d_k, @fig:nsig_1d_pi).
A dedicated kaon selection is therefore required (@sec:coh_pid, @sec:incoh_pid).
@fig:nsigk_p also shows that the kaon band is centered at positive $n sigma_K$ rather than at zero; this shift is measured with TOF-identified kaons in @sec:nsigk_tof.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/track/qa/track2d/track_p_dedx_expected_all.pdf", width: 100%), caption: [$dif E slash dif x$ versus momentum for all tracks, with the expected $dif E slash dif x$ for each hypothesis.]) <fig:dedx_p>],
  [#figure(image("../../plots/track/qa/track2d/track_p_nsigma_kaon_all.pdf", width: 100%), caption: [$n sigma_K$ versus momentum for the same tracks.]) <fig:nsigk_p>],
  [#figure(image("../../plots/track/qa/track/track_nsigma_kaon.pdf", width: 100%), caption: [$n sigma_K$ of tracks in ULS and LS pairs.]) <fig:nsig_1d_k>],
  [#figure(image("../../plots/track/qa/track/track_nsigma_pion.pdf", width: 100%), caption: [$n sigma_pi$ of the same tracks.]) <fig:nsig_1d_pi>],
)
=== $n sigma_K$ Calibration with TOF-Identified Kaons <sec:nsigk_tof>

The $n sigma_K$ shift is measured with kaons identified independently of $dif E slash dif x$ by the TOF.
The $phi.alt -> K^+ K^-$ UPC sample itself yields few TOF-identified kaons: the decay kaons are soft and mostly at large $abs(eta)$, outside the TOF acceptance, and the sample is small.
We therefore use kaons from peripheral hadronic Au+Au collisions in the same Run-19 data (zdcmb trigger), referred to below as the peripheral sample, which cover the full momentum range with high statistics.
The TOF-identified kaons are selected as follows:
- centrality 70--80% from `StRefMultCorr` (`getCentralityBin9() == 0`), $abs(V_z) < 50 thin "cm"$, $V_r < 2 thin "cm"$, bad runs removed;
- primary tracks of the highest-ranked vertex with $N_"hits" >= 20$ and $N_"hits" \/ N_"poss" > 0.52$;
- TOF-matched with $beta < 1$ and $0.22 < m^2 < 0.27 thin "GeV"^2\/c^4$, where $m^2 = p^2 (1\/beta^2 - 1)$ (the kaon tuple extends to $0.28 thin "GeV"^2\/c^4$);
- $abs(n sigma_K) < 5$.

The TOF mass spectrum (@fig:tof_m2_all) shows well-separated pion, kaon, and proton peaks, and the $m^2$ window selects the kaon peak (@fig:tof_m2_kaon).
The selected tracks form a single band in $n sigma_K$ over the full momentum range (@fig:tof_nsigk_p), with only small pion and proton contamination at low momentum, confirming that the sample consists of genuine kaons.

The tracks are divided into momentum bins of $0.1 thin "GeV"\/c$ between 0.1 and $2.5 thin "GeV"\/c$, and the $n sigma_K$ distribution in each bin is fitted with a Gaussian in $[-5, 5]$.
The fitted mean $chevron.l n sigma_K chevron.r (p)$ (@fig:tof_nsigk_mean) is positive at all momenta and largest at low momentum, consistent with the shift seen in @fig:nsigk_p.
We correct for it with
$ n sigma_K^"corr" = n sigma_K - chevron.l n sigma_K chevron.r (p). $ <eq:nsigk_corr>
All kaon requirements from here on use $n sigma_K^"corr"$, and $chi^2_(K K)$ denotes the value computed from it.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/track/PID/track_m2_all.pdf", width: 100%), caption: [TOF $m^2$ of all TOF-matched tracks in the peripheral sample.]) <fig:tof_m2_all>],
  [#figure(image("../../plots/track/PID/track_m2.pdf", width: 100%), caption: [$m^2$ of the kaon candidates in the tuple ($0.22 < m^2 < 0.28 thin "GeV"^2\/c^4$).]) <fig:tof_m2_kaon>],
  [#figure(image("../../plots/track/PID/track_p_nsigma_kaon.pdf", width: 100%), caption: [$n sigma_K$ versus momentum for the same kaon candidates, before the $abs(n sigma_K) < 5$ requirement.]) <fig:tof_nsigk_p>],
  [#figure(image("../../plots/track/PID/nsigmakaon_mean_vs_p.pdf", width: 100%), caption: [Mean $n sigma_K$ of TOF-identified kaons versus momentum, from Gaussian fits.]) <fig:tof_nsigk_mean>],
)

The Gaussian fit in each momentum bin up to $p = 1.6 thin "GeV"\/c$ is shown in @fig:tof_nsigk_fits.
The correction is applied in this range, which contains essentially all tracks of the $phi.alt$ sample; at higher momentum the kaon and pion bands merge and the fits become unreliable.

#figure(
  kind: image,
  table(
    columns: (1fr, 1fr, 1fr),
    stroke: none,
    inset: 2pt,
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_00.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_01.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_02.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_03.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_04.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_05.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_06.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_07.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_08.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_09.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_10.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_11.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_12.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_13.pdf", width: 100%),
    image("../../plots/track/PID/fits/nsigmakaon_fit_pbin_14.pdf", width: 100%),
  ),
  caption: [$n sigma_K$ of TOF-identified kaons in each momentum bin with the Gaussian fit (red) used for @fig:tof_nsigk_mean.],
) <fig:tof_nsigk_fits>

=== Cross-Check with Kaons from $phi.alt -> K^+ K^-$ <sec:nsigk_phi>

The $n sigma_K$ shift is cross-checked in situ with kaons from $phi.alt$ decays in the UPC sample, selected as follows:
- BBC veto, $N_"prim" = 2$, $abs(V_z) <= 100 thin "cm"$, $V_r <= 5 thin "cm"$, pair rapidity $abs(y_(K K)) < 0.5$, and $abs(n sigma_K) <= 5$ for both tracks;
- tag and probe: a track (probe) is used when the other track (tag) satisfies the pion rejection $abs(n sigma_pi) >= 3$, so that no pion rejection is applied to the track whose $n sigma_K$ is measured.

The non-$phi.alt$ background in the mass window $1.00 < M_(K K) < 1.04 thin "GeV"\/c^2$ is described by probes from the sideband $1.06 < M_(K K) < 1.10 thin "GeV"\/c^2$.
The sideband is scaled to the background in the window obtained from a fit of the probe mass spectrum with a Voigtian (PDG width) plus a linear background (@fig:phi_probe_mass).
In each momentum bin, the $n sigma_K$ distribution of the window probes is fitted with a binned Poisson likelihood using a Gaussian for the kaons plus the scaled sideband distribution as a fixed background template.

The resulting $chevron.l n sigma_K chevron.r (p)$ agrees with the TOF-identified kaons within the statistical uncertainties (@fig:phi_nsigk_mean).
Since the TOF sample is far more precise and covers the full momentum range, the TOF result (@sec:nsigk_tof) is used as the nominal correction.

#table(
  columns: (1fr, 1fr),
  stroke: none,
  align: top,
  [#figure(image("../../plots/track/PID/phi/mass.pdf", width: 100%), caption: [$M_(K K)$ of tag-and-probe pairs with the fit; dashed lines mark the mass window and the sideband.]) <fig:phi_probe_mass>],
  [#figure(image("../../plots/track/PID/phi/mean_vs_p.pdf", width: 100%), caption: [Mean $n sigma_K$ of kaons from $phi.alt$ decays compared with the TOF-identified kaons.]) <fig:phi_nsigk_mean>],
)
