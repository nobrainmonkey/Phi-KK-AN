#import "../helpers.typ": figure-pair

= Efficiencies
== Vertex Selection Efficiency
The vertex selection efficiency is the fraction of well-reconstructed vertices from the luminous region that lie inside the fiducial region $|V_z| < 100 "cm"$. It is measured with two samples: all UPC-like events (@sec:vz_upclike), whose $V_z$ distribution has a broad background plateau, and a clean sample of coherent di-pion events (@sec:vz_pion). Both give an efficiency consistent with unity, and $epsilon_(V_z) = 1.00$ is used as the nominal value (@sec:vz_nominal).

Both samples start from the following events, referred to as UPC-like events:
- #raw("zdcmb") triggered, excluding bad runs.
- $#raw("bbcEast") = 0  "and"  #raw("bbcWest") = 0$.
- Exactly two primary tracks matched to the same vertex.

=== Using All UPC-like Events <sec:vz_upclike>
For all UPC-like events, there is no kinematic or particle-identification restriction on the reconstructed vertices, so secondary interactions are not rejected. @fig:Vr_all and @fig:Vx_Vy_all show signatures of interactions at $V_r approx 3.8 "cm"$, which correspond to the beryllium beam pipe radius.

A second, sharp structure appears at $V_r = 10 "cm"$ and extends up to $V_r = 10 sqrt(2) approx 14.1 "cm"$. In @fig:Vx_Vy_all these vertices lie on the edges of the square $|V_x| = 10 "cm"$ or $|V_y| = 10 "cm"$, which is consistent with vertices pushed to the boundary of the transverse region allowed in the vertex reconstruction rather than interactions at a physical radius. These vertices carry no information on the luminous region and are not considered further.

#figure-pair(
  [#figure(image("../../plots/events/vz_vr.pdf", width: 100%), caption: [$V_r$ distribution for all UPC-like events.]) <fig:Vr_all>],
  [#figure(image("../../plots/events/vz_vxvy.pdf", width: 100%), caption: [$V_x$ versus $V_y$ distribution for all UPC-like events.]) <fig:Vx_Vy_all>],
)

To suppress secondary interactions for the efficiency study only, the following additional vertex selection criteria are implemented:
- $V_r < 3 "cm"$.
- Unlike-sign (ULS) pairs only.

@fig:Vx_Vy_vr3_uls shows that the additional vertex quality cuts remove this background. @fig:Vz_all and @fig:Vz_vr3 show the $V_z$ distribution before and after the vertex quality cuts.

#figure(image("../../plots/events/vz_vxvy_vr3_uls.pdf", width: 50%), placement: auto, caption: [$V_x$ versus $V_y$ distribution for unlike-sign pairs in UPC-like events with $V_r < 3 "cm"$.]) <fig:Vx_Vy_vr3_uls>
#figure-pair(
  [#figure(image("../../plots/events/vz_vz.pdf", width: 100%), caption: [$V_z$ distribution for all UPC-like events.]) <fig:Vz_all>],
  [#figure(image("../../plots/events/vz_vz_vr3.pdf", width: 100%), caption: [$V_z$ distribution for UPC-like events with $V_r < 3 "cm"$.]) <fig:Vz_vr3>],
)

The $V_z$ distribution after the vertex quality cuts is fitted with the sum of two Gaussians,
$ f(V_z) = A_1 exp(-(V_z - mu_1)^2 / (2 sigma_1^2)) + A_2 exp(-(V_z - mu_2)^2 / (2 sigma_2^2)), $ <eq:vz_double_gaus>
in 5 cm bins over $|V_z| < 190 "cm"$, which stays inside the $plus.minus 200 "cm"$ edge of the vertex reconstruction. The narrow Gaussian ($sigma_1$) describes the luminous region, while the wide Gaussian ($sigma_2$) absorbs the broad plateau seen in @fig:Vz_vr3. The fraction of a Gaussian inside the fiducial region is
$ F(mu, sigma) = 1/2 [ "erf"((100 - mu) / (sqrt(2) sigma)) - "erf"((-100 - mu) / (sqrt(2) sigma)) ], $ <eq:gaus_fraction>
from which three estimates of the fraction inside $|V_z| < 100 "cm"$ are formed:
- $f_"narrow" = F(mu_1, sigma_1)$, from the narrow Gaussian only;
- $f_"both" = (1 - w_2) F(mu_1, sigma_1) + w_2 F(mu_2, sigma_2)$, from both Gaussians, where $w_2 = A_2 sigma_2 \/ (A_1 sigma_1 + A_2 sigma_2)$ is the area fraction of the wide Gaussian;
- $f_"count"$, the fraction of entries with $|V_z| < 100 "cm"$, counted directly from the histogram.

The fits to the unlike-sign pairs (@fig:Vz_fit_vr3_uls) and, as a cross-check, to the unlike-sign minus like-sign (ULS $-$ LS) distribution (@fig:Vz_fit_vr3_sub) are summarized in @tab:vz_fit.

#figure(
  table(
    columns: 7,
    align: (left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    table.header(
      [Sample], [$sigma_1$ [cm]], [$sigma_2$ [cm]], [$w_2$], [$f_"narrow"$], [$f_"both"$], [$f_"count"$],
    ),
    table.hline(),
    [ULS],       [19.9], [124.3], [0.372], [1.0000], [0.843], [0.879],
    [ULS $-$ LS], [19.0], [64.1],  [0.242], [1.0000], [0.971], [0.975],
    table.hline(),
  ),
  caption: [Double-Gaussian fit results for UPC-like events with $V_r < 3 "cm"$. All fractions refer to $|V_z| < 100 "cm"$.],
) <tab:vz_fit>

The narrow Gaussian is stable between the two samples ($sigma_1 approx 19$–$20 "cm"$, $mu_1 approx -5.4 "cm"$) and lies entirely inside the fiducial region, $f_"narrow" = 1.0000$. The wide Gaussian, in contrast, changes substantially when the like-sign pairs are subtracted ($sigma_2$ from 124 to 64 cm, $w_2$ from 0.37 to 0.24), showing that it is driven by the residual background plateau rather than by the luminous region. The lower values of $f_"both"$ and $f_"count"$ therefore reflect background that extends beyond $|V_z| = 100 "cm"$ and are not used as the efficiency. The vertex selection efficiency from this method is taken from the narrow Gaussian, $epsilon_(V_z) approx 1.000$. Both fits have large $chi^2\/"ndf"$ (964/70 and 340/70), as two Gaussians do not fully describe the peak-plus-plateau shape; this does not affect $f_"narrow"$, which is set by $sigma_1 approx 20 "cm"$ relative to the 100 cm cut.

#figure-pair(
  [#figure(image("../../plots/events/vz_fit_vr3_uls.pdf", width: 100%), caption: [Double-Gaussian fit to the $V_z$ distribution of unlike-sign pairs in UPC-like events with $V_r < 3 "cm"$. Dashed vertical lines mark $|V_z| = 100 "cm"$.]) <fig:Vz_fit_vr3_uls>],
  [#figure(image("../../plots/events/vz_fit_vr3_sub.pdf", width: 100%), caption: [Double-Gaussian fit to the unlike-sign minus like-sign $V_z$ distribution in UPC-like events with $V_r < 3 "cm"$. Dashed vertical lines mark $|V_z| = 100 "cm"$.]) <fig:Vz_fit_vr3_sub>],
)

=== Using Di-pion Events <sec:vz_pion>
To avoid the background plateau of the UPC-like events, the luminous region is also measured with a clean sample of coherently photoproduced $pi^+ pi^-$ pairs. In addition to the UPC-like event definition, the following selection is applied:
- $#raw("TofMult") <= 6$ and at least one track matched to the TOF.
- Both tracks satisfy $|n sigma_pi| < 3$.
- Pair transverse momentum $p_T^(pi pi) < 0.1 "GeV"\/c$ (coherent region).
- $0.60 < M_(pi pi) < 0.95 "GeV"\/c^2$ around the $rho^0$ pole, for the vertex distributions only.
The vertex coordinates are those of the vertex the two pions are matched to.

==== Di-pion Mass Spectrum
The mass spectrum in this section is not corrected for acceptance or efficiency and is not used to extract a cross section. It only serves to check that the selected di-pion sample is not background dominated and that its shape agrees with the published coherent $rho^0$ result.

@fig:pion_mass compares the unlike-sign and like-sign invariant-mass spectra of the selected pairs. The like-sign pairs are about two orders of magnitude below the unlike-sign pairs over the whole mass range. The unlike-sign minus like-sign spectrum is fitted over $0.6 < M_(pi pi) < 1.2 "GeV"\/c^2$ with the modified Söding parametrization used by STAR for coherent $rho^0$ photoproduction (Phys. Rev. C 96 (2017) 054904),
$ (d N)/(d M_(pi pi)) prop lr(|A_rho (sqrt(M_(pi pi) M_rho Gamma_rho)) / (M_(pi pi)^2 - M_rho^2 + i M_rho Gamma_rho) + B_(pi pi) + C_omega e^(i phi_omega) (sqrt(M_(pi pi) M_omega Gamma_(omega -> pi pi))) / (M_(pi pi)^2 - M_omega^2 + i M_omega Gamma_omega)|)^2 + f_p, $ <eq:soding>
with mass-dependent widths as in that reference and a linear remnant background $f_p$. The spectrum uses 2.5 MeV$\/c^2$ bins in $0.72 < M_(pi pi) < 0.84 "GeV"\/c^2$, where the $rho^0$–$omega$ interference is located, and 10 MeV$\/c^2$ bins elsewhere. The fit and its components are shown in @fig:pion_mass_fit and the parameters are compared with the published values in @tab:pion_mass_fit.

#figure(image("../../plots/events/vz_pion_mass.pdf", width: 75%), placement: auto, caption: [Invariant-mass spectra of unlike-sign and like-sign pion pairs passing the di-pion selection, without vertex or mass-window requirements.]) <fig:pion_mass>
#figure(image("../../plots/events/vz_pion_mass_fit.pdf", width: 85%), placement: auto, caption: [Modified Söding fit to the unlike-sign minus like-sign $pi^+ pi^-$ spectrum, with the $rho^0$, $omega$, interference, non-resonant and remnant-background components.]) <fig:pion_mass_fit>

#figure(
  table(
    columns: 3,
    align: (left, right, right),
    stroke: none,
    table.hline(),
    table.header([Parameter], [This analysis], [STAR, PRC 96, 054904]),
    table.hline(),
    [$M_rho$ [GeV$\/c^2$]],      [0.7759], [$0.7762 plus.minus 0.0006$],
    [$Gamma_rho$ [GeV$\/c^2$]],  [0.151],  [$0.156 plus.minus 0.001$],
    [$B_(pi pi) \/ A_rho$ [$("GeV"\/c^2)^(-1\/2)$]], [$-0.742$], [$-0.787$],
    [$C_omega \/ A_rho$],        [0.295],  [0.358],
    [$M_omega$ [GeV$\/c^2$]],    [0.7829], [$0.7824 plus.minus 0.0008$],
    [$Gamma_omega$ [MeV$\/c^2$]], [15.8],  [$17 plus.minus 2$],
    [$phi_omega$ [rad]],         [1.33],   [$1.46 plus.minus 0.11$],
    [$chi^2\/"ndf"$],            [85.5/86], [–],
    table.hline(),
  ),
  caption: [Modified Söding fit to the coherent unlike-sign minus like-sign $pi^+ pi^-$ spectrum, compared with the published fit. The published ratios are computed from the quoted amplitudes.],
) <tab:pion_mass_fit>

The $rho^0$ mass and width and the $omega$ parameters are close to the published values. As the spectrum is uncorrected and no systematic uncertainties are evaluated, this comparison is qualitative: it shows that the lineshape is reproduced. As in the published fit, no mass resolution is included in the model, so the fitted $Gamma_omega$ absorbs the detector resolution and is about twice the natural width. The ratio $B_(pi pi) \/ A_rho = -0.742$ is also close to the published $-0.787$, showing that the non-resonant $pi^+ pi^-$ component is described by $B_(pi pi)$ rather than absorbed into the remnant background $f_p$. Inside the mass window $0.60 < M_(pi pi) < 0.95 "GeV"\/c^2$ there are 56788 unlike-sign and 301 like-sign pairs ($N_"LS" \/ N_"ULS" = 0.0053$), and the remnant background $f_p$ amounts to 2415 pairs, giving a $pi^+ pi^-$ purity of $1 - (N_"LS" + f_p) \/ N_"ULS" = 0.952$. The selected sample is therefore dominated by coherent $pi^+ pi^-$ production.

==== Vertex Distributions and Efficiency
@fig:pion_vr and @fig:pion_vxvy show the transverse position of the di-pion vertices. Unlike the UPC-like events, there is no ring at $V_r approx 3.8 "cm"$: the vertices are concentrated around the beam spot, with only a small population at the $V_r = 10 "cm"$ reconstruction boundary. The $V_z$ distribution (@fig:pion_vz) shows no plateau, and the like-sign contribution is negligible.

#figure-pair(
  [#figure(image("../../plots/events/vz_pion_vr.pdf", width: 100%), caption: [$V_r$ distribution of coherent $pi^+ pi^-$ pairs.]) <fig:pion_vr>],
  [#figure(image("../../plots/events/vz_pion_vxvy.pdf", width: 100%), caption: [$V_x$ versus $V_y$ distribution of coherent $pi^+ pi^-$ pairs.]) <fig:pion_vxvy>],
)
#figure(image("../../plots/events/vz_pion_vz.pdf", width: 65%), placement: auto, caption: [$V_z$ distribution of coherent $pi^+ pi^-$ pairs.]) <fig:pion_vz>

The $V_z$ distributions are fitted with the same double-Gaussian model as for the UPC-like events, for unlike-sign pairs and for unlike-sign minus like-sign pairs, with and without the $V_r < 3 "cm"$ requirement (@tab:pion_vz_fit). The fits for unlike-sign pairs without a $V_r$ requirement and for unlike-sign minus like-sign pairs with $V_r < 3 "cm"$ are shown in @fig:pion_vz_fit_uls and @fig:pion_vz_fit_vr3_sub.

#figure(
  table(
    columns: 8,
    align: (left, left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    table.header(
      [Sample], [$V_r$ cut], [$sigma_1$ [cm]], [$sigma_2$ [cm]], [$w_2$], [$chi^2\/"ndf"$], [$f_"both"$], [$f_"count"$],
    ),
    table.hline(),
    [ULS],        [none],          [18.6], [33.3], [0.155], [81/57], [0.9996], [0.9990],
    [ULS $-$ LS], [none],          [18.5], [32.4], [0.169], [79/60], [0.9996], [0.9993],
    [ULS],        [$< 3 "cm"$],    [18.5], [32.2], [0.168], [56/48], [0.9997], [0.9995],
    [ULS $-$ LS], [$< 3 "cm"$],    [18.5], [32.0], [0.172], [56/48], [0.9997], [0.9995],
    table.hline(),
  ),
  caption: [Double-Gaussian fit results for coherent $pi^+ pi^-$ pairs. In all cases $f_"narrow" = 1.0000$. All fractions refer to $|V_z| < 100 "cm"$.],
) <tab:pion_vz_fit>

In contrast to the UPC-like events, the wide Gaussian is the same in all four cases ($sigma_2 approx 32$–$33 "cm"$, $w_2 approx 0.16$–$0.17$): it does not depend on the like-sign subtraction or on the $V_r$ requirement and is therefore part of the luminous region rather than background. Both Gaussians are thus included, and the fits describe the data well ($chi^2\/"ndf" approx 1.2$–$1.4$). The wide Gaussian lies almost entirely inside the fiducial region, so $f_"both"$, $f_"narrow"$ and $f_"count"$ agree within $0.001$. The vertex selection efficiency from the di-pion sample is $epsilon_(V_z) = 0.9996$, with a spread of $0.0006$ among the variations in @tab:pion_vz_fit, consistent with the narrow-Gaussian result from the UPC-like events.

#figure-pair(
  [#figure(image("../../plots/events/vz_pion_fit_uls.pdf", width: 100%), caption: [Double-Gaussian fit to the $V_z$ distribution of coherent unlike-sign $pi^+ pi^-$ pairs without a $V_r$ requirement. Dashed vertical lines mark $|V_z| = 100 "cm"$.]) <fig:pion_vz_fit_uls>],
  [#figure(image("../../plots/events/vz_pion_fit_vr3_sub.pdf", width: 100%), caption: [Double-Gaussian fit to the unlike-sign minus like-sign $V_z$ distribution of coherent $pi^+ pi^-$ pairs with $V_r < 3 "cm"$. Dashed vertical lines mark $|V_z| = 100 "cm"$.]) <fig:pion_vz_fit_vr3_sub>],
)

=== Nominal Vertex Selection Efficiency <sec:vz_nominal>
Both methods give a vertex selection efficiency consistent with unity: the narrow Gaussian of the UPC-like events gives $f_"narrow" = 1.000$, and the coherent $pi^+ pi^-$ sample gives $epsilon_(V_z) = 0.9996$, with all variations within $0.001$ of unity. We therefore use $epsilon_(V_z) = 1.00$ as the nominal value.

== BBC Veto Selection Efficiency
The BBC veto (@eq:bbc_dsm_veto) rejects events with activity in the small BBC tiles. Besides hadronic background, it also removes genuine UPC events that are accompanied by additional BBC activity, either from the same collision (e.g. electromagnetic dissociation of the nuclei) or from pile-up. The corresponding efficiency is measured with coherently photoproduced $pi^+ pi^-$ pairs, as the ratio of the $rho^0$ yields with and without the veto.

The di-pion selection of the vertex study is used, with two changes: the $#raw("TofMult") <= 6$ requirement is removed, since TOF activity is correlated with BBC activity and such a cut would bias the veto efficiency, and no $V_z$ requirement is applied. The selection is therefore exactly two primary tracks on one vertex, at least one track matched to the TOF, $|n sigma_pi| < 3$ for both tracks and $p_T^(pi pi) < 0.1 "GeV"\/c$ as a proxy for coherent selection.

The efficiency can only be measured on a sample that remains signal dominated without the BBC veto. This is ensured by the selection itself. Exactly two primary tracks identified as pions are required, and $pi^+ pi^-$ photoproduction (dominantly $rho^0$) is the most abundant exclusive final state in UPC. In addition, $p_T^(pi pi) < 0.1 "GeV"\/c$ restricts the sample to the coherent peak, a narrow window in pair transverse momentum where hadronic and incoherent contributions are small. As a result, even without the BBC veto the like-sign pairs amount to only $N_"LS" \/ N_"ULS" = 1.0%$ in $0.6 < M_(pi pi) < 1.2 "GeV"\/c^2$, compared with $0.64%$ with the veto.

The unlike-sign minus like-sign mass spectra of all selected events, of the events passing the BBC veto, and of the events rejected by it are fitted over $0.6 < M_(pi pi) < 1.2 "GeV"\/c^2$ with the modified Söding parametrization (@eq:soding), using the same binning and fit procedure as in the vertex study (@fig:bbc_pre, @fig:bbc_post, @fig:bbc_rejected). All three fits describe the data well and their shape parameters agree within uncertainties (@tab:bbc_fit), so the veto changes the yield but not the shape of the spectrum.

None of the samples is background dominated. The remnant background $f_p$ amounts to 5.8% of the unlike-sign minus like-sign pairs for all events and 4.8% for events passing the veto. The slightly larger background without the veto is expected, since the veto preferentially removes hadronic background: only about half of the like-sign pairs (750 to 379) survive the veto, compared with 82% of the signal. Two yields are extracted from each fit: the integral of the fit without the remnant background $f_p$ over the fit range, and the integral of the $rho^0$ Breit-Wigner from $2 m_pi$ to $M_rho + 5 Gamma_rho$.

#figure-pair(
  [#figure(image("../../plots/events/BBC_veto_pre.pdf", width: 100%), caption: [Modified Söding fit to the unlike-sign minus like-sign $pi^+ pi^-$ spectrum without the BBC veto, with all fit components.]) <fig:bbc_pre>],
  [#figure(image("../../plots/events/BBC_veto_post.pdf", width: 100%), caption: [Modified Söding fit to the unlike-sign minus like-sign $pi^+ pi^-$ spectrum with the BBC veto, with all fit components.]) <fig:bbc_post>],
)
#figure(image("../../plots/events/BBC_veto_rejected.pdf", width: 60%), placement: auto, caption: [Modified Söding fit to the unlike-sign minus like-sign $pi^+ pi^-$ spectrum of the events rejected by the BBC veto, with all fit components.]) <fig:bbc_rejected>

#figure(
  table(
    columns: 4,
    align: (left, right, right, right),
    stroke: none,
    table.hline(),
    table.header([], [All events], [Pass BBC veto], [Rejected by BBC veto]),
    table.hline(),
    [$M_rho$ [GeV$\/c^2$]],       [$0.7756 plus.minus 0.0012$], [$0.7759 plus.minus 0.0013$], [$0.7749 plus.minus 0.0025$],
    [$Gamma_rho$ [GeV$\/c^2$]],   [$0.150 plus.minus 0.003$],   [$0.151 plus.minus 0.004$],   [$0.143 plus.minus 0.006$],
    [$B_(pi pi) \/ A_rho$ [$("GeV"\/c^2)^(-1\/2)$]], [$-0.740 plus.minus 0.017$], [$-0.742 plus.minus 0.018$], [$-0.744 plus.minus 0.040$],
    [$C_omega \/ A_rho$],         [$0.283 plus.minus 0.049$],   [$0.295 plus.minus 0.056$],   [$0.249 plus.minus 0.103$],
    [$Gamma_omega$ [MeV$\/c^2$]], [$16.2 plus.minus 4.1$],      [$16.0 plus.minus 4.5$],      [$17.0 plus.minus 8.2$],
    [$chi^2\/"ndf"$],             [81/86],  [85/86],  [87/86],
    table.hline(),
    [Yield (all $-$ $f_p$)],      [$67586 plus.minus 793$],  [$55730 plus.minus 722$], [$11722 plus.minus 315$],
    [$rho^0$ BW yield],           [$72574 plus.minus 1110$], [$59688 plus.minus 990$], [$12818 plus.minus 479$],
    table.hline(),
  ),
  placement: auto,
  caption: [Modified Söding fits to the coherent $pi^+ pi^-$ spectrum for all selected events, for events passing the BBC veto, and for events rejected by it. Uncertainties are statistical, from the fits; the uncertainties of $B_(pi pi) \/ A_rho$ and $C_omega \/ A_rho$ include the correlation with $A_rho$.],
) <tab:bbc_fit>

The fit to the rejected events (@fig:bbc_rejected) shows that these are predominantly genuine UPC-like di-pions rather than hadronic background. Their spectrum shows the same $rho^0$ lineshape, with $M_rho$, $Gamma_rho$ and $B_(pi pi) \/ A_rho$ consistent with those of the accepted events and the $rho^0$–$omega$ interference visible at $M_(pi pi) approx 0.78 "GeV"\/c^2$. The background in this sample is larger but still small: $N_"LS" \/ N_"ULS" = 2.7%$ and $f_p$ amounts to 10% of the unlike-sign minus like-sign pairs. The BBC veto therefore removes genuine UPC signal.

The efficiency is computed from the two disjoint samples as $epsilon_"BBC" = Y_"pass" \/ (Y_"pass" + Y_"rejected")$. Since the samples do not overlap, their fit uncertainties are independent and are propagated directly; a binomial uncertainty on $Y_"pass" \/ Y_"all"$ would neglect the uncertainty of the fit decomposition, which dominates the yield uncertainties. The efficiency is
- $epsilon_"BBC" = 0.826 plus.minus 0.004$ from the yield without $f_p$,
- $epsilon_"BBC" = 0.823 plus.minus 0.006$ from the $rho^0$ Breit-Wigner yield,
- $epsilon_"BBC" = 0.817 plus.minus 0.002$ from a direct unlike-sign minus like-sign count over the fit range.
The sum of the passing and rejected yields reproduces the yield of all events to 0.2%, and the ratio $Y_"pass" \/ Y_"all"$ gives consistent values (0.825 and 0.822 for the two fit-based yields). The direct count is lower because it includes the remnant background, which is preferentially removed by the veto. The spectra with and without the veto are compared in @fig:bbc_compare.

#figure(image("../../plots/events/BBC_veto_compare.pdf", width: 65%), placement: auto, caption: [Unlike-sign minus like-sign $pi^+ pi^-$ spectra with and without the BBC veto, with the resulting efficiencies.]) <fig:bbc_compare>

To check whether the efficiency depends on the kinematics, it is computed in bins of $M_(pi pi)$ and of the pair rapidity $y_(pi pi)$ within the $rho^0$ window $0.60 < M_(pi pi) < 0.95 "GeV"\/c^2$, from unlike-sign minus like-sign counts with binomial uncertainties (@fig:bbc_eff_mass, @fig:bbc_eff_rapidity). The efficiency does not vary with the kinematics: compared with the efficiency of the whole window, the bins give $chi^2 = 8.4$ for 7 bins in $M_(pi pi)$ (0.809–0.823) and $chi^2 = 2.4$ for 8 bins in $y_(pi pi)$ (0.814–0.821). The binned counts lie about 0.008 below the nominal value because, like the integrated count, they include the remnant background that the veto removes preferentially. The pair transverse momentum is not used for this check. The selection is restricted to the coherent region, where the sample is clean; towards larger $p_T^(pi pi)$ the hadronic background grows, so a dependence on $p_T^(pi pi)$ would reflect the changing background fraction rather than the efficiency for the signal.

#figure-pair(
  [#figure(image("../../plots/events/BBC_veto_eff_mass.pdf", width: 100%), caption: [BBC veto efficiency versus $M_(pi pi)$ from unlike-sign minus like-sign counts. The solid line and band show the nominal value and its uncertainty, the dashed line the efficiency from the integrated count over the window.]) <fig:bbc_eff_mass>],
  [#figure(image("../../plots/events/BBC_veto_eff_rapidity.pdf", width: 100%), caption: [BBC veto efficiency versus the pair rapidity $y_(pi pi)$, for $0.60 < M_(pi pi) < 0.95 "GeV"\/c^2$. Lines as in @fig:bbc_eff_mass.]) <fig:bbc_eff_rapidity>],
)


=== Nominal BBC Veto Efficiency
We use $epsilon_"BBC" = 0.826 plus.minus 0.004$ as the nominal value, from the yield without $f_p$, where the uncertainty is statistical, from the fits. The $rho^0$ Breit-Wigner yield gives a consistent result.
