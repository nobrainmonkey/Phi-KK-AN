# Coherent phi → K+K− analysis

All rapidity-slice code, plots, ROOT outputs, reports and analysis logs are in [rapidity_slice/](rapidity_slice/README.md).

From this workspace, run the simple C++/ROOT extraction with:

```sh
root -l -b -q 'rapidity_slice/coherent_four_bins.C()'
root -l -b -q 'rapidity_slice/coherent_four_bins.C(0.2)'
```

The macro also runs from inside `rapidity_slice/`; default input/output paths are resolved relative to its own location. Outputs go to `rapidity_slice/results/four_bins_root/`.

The original ROOT input files remain in `data/`, and the shared Python environment remains in `.venv/`.
