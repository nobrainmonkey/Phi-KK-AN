= Event Selection

UPCs are characterized by large impact parameters, which strongly suppress direct hadronic interactions between the two nuclei.
In photon-induced vector-meson production, the exchanged photon is colorless and the interaction with the target proceeds through a color-singlet exchange.
Consequently, no color connection is established between the two beam nuclei, suppressing particle production between the centrally produced vector meson and the outgoing nuclear remnants.
This gives rise to large rapidity gaps on either side of the produced $phi.alt$ meson.

We exploit this topology to suppress hadronic Au+Au backgrounds.
The rapidity-gap requirement is implemented through a veto on forward charged-particle activity in the Beam-Beam Counters (BBCs), while the exclusive nature of the final state is further selected using low central-event multiplicity requirements.

== Vertex Selection

The reconstructed primary vertex is used to identify the collision location and to ensure that selected events originate within the nominal STAR interaction region.
Since the geometrical acceptance and tracking performance of the TPC depend on the vertex position, particularly along the beam direction, events reconstructed far from the detector center may exhibit reduced or non-uniform tracking efficiency.
Vertex-quality requirements are therefore applied to define a well-controlled collision region and to suppress poorly reconstructed or non-collision vertices before the remaining event-level selections are imposed.

For this analysis, we require the reconstructed event vertex to satisfy
$ abs(V_z) < 100 thin "cm". $ <eq:vz_cut>
The integrated luminosity is determined using the ZDC-based luminosity measurement, whose longitudinal acceptance is not restricted in the same way as the TPC tracking acceptance.
Because the effective fiducial region of the TPC is smaller than the longitudinal interaction region sampled by the ZDC luminosity measurement, the $V_z$ requirement introduces an additional fiducial acceptance on the analyzed event sample.
The corresponding efficiency associated with this restricted vertex region, including the $V_z$-dependent vertex reconstruction and tracking efficiency, is therefore evaluated separately using embedding and incorporated into the analysis correction.

== BBC Veto

The BBCs are forward scintillator detectors located on the east and west sides of STAR, approximately $3.7 thin "m"$ from the interaction point.
The small (inner) BBC tiles are available in this dataset and cover the pseudorapidity range $3.3 < abs(eta) < 5.0$.
The BBCs are therefore sensitive to charged-particle production in the forward region.
In this analysis, BBC activity is used as a veto to select the large rapidity gaps characteristic of UPC events and to suppress backgrounds from hadronic Au+Au collisions.

In past STAR UPC triggers, the BBC-veto online trigger was utilized for rapidity gap enforcement with the small BBC tiles.
For each of the 18 small-tile PMTs on a given side, a BBC hit was defined by:
$ "ADC"_"PMT" > 20 $ <eq:bbc_adc>
together with the timing requirement
$ 100 <= "TDC"_"PMT" <= 2300. $ <eq:bbc_tdc>
The qualifying PMT ADC signals (those with a BBC hit) were summed separately in the East and West,
$ S_"East" equiv sum_(e in "hits") "ADC"_e, quad S_"West" equiv sum_(w in "hits") "ADC"_w, $ <eq:bbc_adc_sum>
and the UPC-Main online veto threshold was 50 ADC channels on each side,
$ S_"East" <= 50 quad and quad S_"West" <= 50. $ <eq:bbc_adc_sum_cut>

The zdcmb trigger used in this analysis does not require this veto, but the BBC conditions are still evaluated at Level 0 and recorded with every event.
The online BBC veto decision is obtained directly from the Level-0 trigger information stored in `StTriggerData`.
For the Run-19 UPC zdcmb configuration, the East and West BBC conditions correspond to global DSM bits 17 and 18, respectively.
Since `lastDSM(1)` contains global trigger bits 16--31, these correspond to bits 1 and 2 of the returned 16-bit word.
The quantities are extracted as
$ #raw("bbcEast") = (#raw("lastDSM(1)") ">>" 1) thick \& thick 1,
  wide
  #raw("bbcWest") = (#raw("lastDSM(1)") ">>" 2) thick \& thick 1. $ <eq:bbc_dsm_bits>
A value of one indicates that the corresponding BBC activity condition fired, i.e. that the small-tile BBC signal exceeded the online veto threshold on that side.
An event therefore passes the BBC veto when neither the East nor the West bit is set.
This analysis applies exactly this online requirement offline,
$ #raw("bbcEast") = 0 wide "and" wide #raw("bbcWest") = 0. $ <eq:bbc_dsm_veto>

== Event Multiplicity Selection

To select events consistent with the exclusive two-track final state of $phi.alt -> K^+ K^-$, an event-multiplicity requirement is applied in addition to the BBC veto.
The corresponding selections are:

- Each event must contain exactly two primary tracks.
- The two primary tracks must be associated with the same reconstructed primary vertex.

These requirements select events with a single two-track primary vertex and suppress events containing additional charged-particle activity, as expected from hadronic interactions or pileup.
