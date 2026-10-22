# UPCDst for Run-19 φ→K⁺K⁻ embedding

Source of `StUPCFilterMaker`, which converts STAR embedding MuDst into the
UPCDst format, together with a minimal standalone reader for its output.

## Layout

| Path | Content |
|---|---|
| `StRoot/StUPCFilterMaker/` | filter maker and the UPCDst data classes (`.h` + `.cxx`) |
| `include/` | headers of the data classes and `star-upcLinkDef.h` |
| `RunFilterMaker.C`, `config_run_19_embedding.in` | STAR-side production macro and its embedding config |
| `tests/` | MC-truth tests (STAR environment) |
| `../../mc_analysis/build_upcdst.sh` | builds the data classes for a local ROOT |
| `../../mc_analysis/mc_skimmer.C` | example reader: event, RC and MC loops |

## Producing UPCDst (STAR environment)

Compile `StRoot/` with `cons`, then

```sh
root -l -b -q 'RunFilterMaker.C("files.list", 1000, "out.UPCDst.root", "config_run_19_embedding.in")'
```

`is_mc 2` in the config stores the full MC event (tracks, vertices) and the
MC↔RC associations, and bypasses the trigger selection.

## Reading UPCDst locally (no STAR needed)

The data classes depend only on ROOT. From `mc_analysis/`:

```sh
./build_upcdst.sh           # → lib/libUPCDst.so (+ .rootmap, .pcm)
root -l -b -q mc_skimmer.C  # default input: ../mc/phi_kk_flat_rap_flat_pt_4.root
```

The skimmer loads the library and headers with

```cpp
// clang-format off
R__ADD_INCLUDE_PATH(../mc/upcdst/include)
R__LOAD_LIBRARY(lib/libUPCDst.so)
// clang-format on
```

The `clang-format` guards keep the formatter from turning the paths into
`lib / libUPCDst.so`. These macros only act in interpreted mode; for ACLiC
(`mc_skimmer.C+`) use `gSystem->AddIncludePath` and `gSystem->Load` instead.

The tree is `mUPCTree` with a single branch `mUPCEvent` (`StUPCEvent`):

```cpp
StUPCEvent *evt = nullptr;
tree->SetBranchAddress("mUPCEvent", &evt);
for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);                       // fills evt
    for (Int_t j = 0; j < evt->getNumberOfTracks(); ++j)   StUPCTrack   *rc = evt->getTrack(j);
    for (Int_t j = 0; j < evt->getNumberOfMCTracks(); ++j) StUPCMCTrack *mc = evt->getMCTrack(j);
}
```

## Classes

```
StUPCEvent              reconstructed event: trigger, ZDC/BBC/VPD, multiplicities
 ├── StUPCTrack[]       primary tracks: pT/η/φ, nHits, DCA, dE/dx nσ, TOF/BEMC flags
 ├── StUPCVertex[]      primary vertices: position, errors, nTracks
 └── StUPCMCEvent       generator truth
      ├── StUPCMCTrack[]   true particles (inherits TParticle): PDG, kinematics, decay chain
      └── StUPCMCVertex[]  true production/decay vertices
```

MC↔RC links, stored in both directions:

| From | To | Call |
|---|---|---|
| `StUPCTrack` | `StUPCMCTrack` | `getMCTrack()` — `nullptr` if unmatched |
| `StUPCMCTrack` | `StUPCTrack` | `getNumberOfRecoTracks()`, `getRecoTrack(k)` |
| `StUPCVertex` | `StUPCMCVertex` | `getMCVertex()` |
| `StUPCMCVertex` | `StUPCVertex` | `getNumberOfRecoVertices()`, `getRecoVertex(k)` |

## Things to know

- **Tracks are primary only.** Every stored track has `kPrimary`. pT/η/φ are
  vertex-constrained; DCA and helix come from the global track.
- **Vertices with > 4 primary tracks are not written.** A kaon at such a
  vertex has no RC match; this is part of the efficiency.
- **Multiple vertices per event are common.** Require pair daughters to share
  `getVertexId()`.
- **Unmatched RC tracks** carry `idTruth` ≥ 10000 (embedding background) and
  have `getMCTrack() == nullptr`.
- **MC tracks 0 and 1 are the φ daughters.** Secondaries store `getPt() = -999`;
  a primary kaon along the beam (η ≈ 6) can too. Use the `TParticle` methods
  (`Pt()`, `Eta()`, `Px/Py/Pz/Energy()`), which are always filled.
- **One MC track can match several RC tracks** (split tracks, ~0.07 % of
  reconstructed kaons). `getRecoTrack()` returns the first in write order,
  not the best. `bestRecoTrack()` in `mc_skimmer.C` picks highest `qaTruth`,
  then most `nHitsFit`.
- **Efficiency:** fill numerator and denominator with truth kinematics and
  count each MC particle at most once in the numerator.
