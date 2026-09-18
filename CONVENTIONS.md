# Data and analysis conventions

Reference for anyone reading the code or re-using the Tambopata (PE-TNR) flux
products. It documents the variable-naming contract inherited from the upstream
EcoFlux/Biomet.net pipeline, the sign conventions, and the analytical rules the
manuscript follows.

The flux variables were produced by an upstream MATLAB/Biomet.net pipeline that
reads a raw instrument database not distributed with this paper. That pipeline
is not part of this repository, so the three definitions a reader most needs are
reproduced verbatim in §0 below. The public naming reference is the
[EcoFlux/Biomet.net pipeline documentation](https://ecoflux-lab.github.io/PipelineDocumentation/5_5_Full_Doc_Third_Stage_Cleaning_And_Ameriflux_Output.html#div_id_RunAmerifluxStage).

---

## 0. Upstream definitions, quoted verbatim

These are the three trace definitions from the second-stage configuration that
determine how PAR, the storage term and NEE were computed.

### PAR is derived from shortwave, not measured

```
[Trace]
	variableName    = 'PPFD_IN_1_1_1'
	Evaluate        = 'PPFD_IN_1_1_1 = SW_IN_1_1_1*2.3-1.22;' % correcting the old PAR data from old sensor, based on a linear correlation with SW_IN_1_1_1
	title           = 'incoming PAR'
	units           = '\mu mol/m^2/s'
[End]
```

The original PAR sensor had a calibration fault, so PAR across the whole record
is `2.30 × SW_IN − 1.22`. Values are **not** truncated at zero. Any regression of
the archived PPFD against SW_IN is therefore circular and returns R² ≈ 1; the
R² ≈ 0.97 reported in the paper comes from the original sensor fit, which cannot
be recovered from the archived product.

### The storage term switches source in 2024

```
[Trace]
	variableName    = 'SC'
	Evaluate        = 'yearIn=year(clean_tv(10));
                       if (yearIn < 2024),
                           SC=SC_model;
                       else,
                          SC = SC_profiler;
                       end'
	title           = 'Estimate of storage CO2 flux'
	units           = 'umol/m^2/s'
[End]
```

So `SC` is the **neural-network reconstruction** (`SC_model`, labelled
"Model of storage CO2 flux v31") for every year before 2024, and the
**profile measurement** (`SC_profiler`) from 2024 onward. This is a year switch,
not a "measured where available" rule — roughly 77% of the record uses the
modelled term, and 2017–2020 has no profile measurement at all. The
consequences for the diel cycle are quantified by
`code/20_storage_diel_comparison.R` (Fig S10).

### NEE

```
[Trace]
	variableName    = 'NEE'
	Evaluate        = 'NEE = sum([FC,SC],2,''omitnan'');
                       NEE(all(isnan(FC)&isnan(SC),2)) = NaN;'
	title           = 'net ecosystem exchange'
	units           = '\mu mol m^{-2} s^{-1}'
[End]
```

`NEE = FC + SC`, with NaN only where both terms are missing. Note that `omitnan`
means a half-hour with a valid `FC` and a missing `SC` yields `NEE = FC`.

One further detail from the export step: the manuscript series comes from the
stem `*_PI_SC_JSZ_MAD_RP_uStar_orig` → `*_PI_SC_JSZ_MAD`, with **no fallback to
raw names**. The first-stage configuration applies a `[-100, 100]` min/max bound
to `SC_model` that is *not* re-applied on export.

---

## 1. Flux-variable naming

### Sequential third-stage suffixes

Third-stage suffixes are **cumulative** and describe transformations applied
sequentially from left to right:

| Suffix | Meaning |
|---|---|
| *(none)* | standard third-stage cleaning, including configured wind-sector and precipitation filtering |
| `_PI_SC` | storage correction added, so NEE = FC + SC |
| `_JSZ` | sliding-window z-score outlier filtering added |
| `_MAD` | median-absolute-deviation outlier filtering added |
| `_RP` | REddyProc processing added, including the configured u\* filtering and flux-partitioning branch |

Not every generated variable carries every optional suffix. The stem actually
used by the manuscript is `*_PI_SC_JSZ_MAD_RP_uStar_orig` → `*_PI_SC_JSZ_MAD`
(see §0).

### REddyProc gap-filling suffixes

Within one exact processing stem and one exact u\* branch:

| Suffix | Meaning |
|---|---|
| `*_orig` | non-gap-filled values supplied to the gap-filling procedure |
| `*_f` | combined series: original values plus filled gaps |
| `*_fqc` | gap-fill quality flag — `0` original, `1` most reliable, `2` medium, `3` least reliable |
| `*_fmeth` | gap-filling method |
| `*_fwin` | full window length used for gap filling |
| `*_fnum`, `*_fsd` | contributing observations and uncertainty |

u\* branches appear as `_uStar`, `_U05`, `_U50`, `_U95`.

**Never pair `_orig`, `_f`, `_fqc`, `_fmeth` or `_fwin` columns from different
processing stems or different u\* branches.** The meaning of an absent `_uStar`
suffix depends on the installed Biomet.net/REddyProc version — check rather than
infer.

### "Observed" is not a synonym for "non-gap-filled"

A non-gap-filled NEE value can still contain a **model-derived storage term**.
Temporal gap filling and modelled storage are separate processing dimensions and
must never be collapsed into a single flag. Carry two independent row-level
indicators:

- `nee_gapfilled` — whether the value came from temporal gap filling
- `storage_source` — one of `measured_profile`, `modeled_profile`,
  `single_point`, `none`, `unknown`

At this site storage is `modeled_profile` before 2024 and `measured_profile`
from 2024 onward — see the switch quoted in §0.

### Consistency checks worth repeating

1. `*_f` equals `*_orig` wherever `*_fqc == 0` and both are finite. Report the
   mismatch count, the maximum absolute mismatch, and the offending timestamps.
2. Wherever `*_fqc > 0` and `*_f` is finite, the matching `*_orig` should be
   missing. Investigate exceptions.
3. Keep the gap-fill quality class; do not reduce it to a boolean too early.
4. Derive gap-run lengths from consecutive missing `*_orig` with finite `*_f`.
5. The manuscript retains gap-filled runs only where the gap does not exceed
   **24 consecutive half-hours**. Apply that as a *separate* variable rather than
   overwriting the original `*_f` series.

---

## 2. Sign conventions

- **NEE** — negative daytime NEE is net ecosystem CO2 **uptake**; positive
  nighttime NEE is net ecosystem CO2 **release**.
- **Partitioning** — `Reco > 0` denotes respiratory release and `GEP < 0`
  photosynthetic uptake, so that `NEE = Reco + GEP` holds at every half-hour.
- **Carbon balance sums** — reported in g C m⁻² d⁻¹, positive = source.
- **SWC** — reported as integrated 0–100 cm soil-water **storage in centimetres**,
  not volumetric water content.
- **GEP ≡ GPP** in these outputs.
- `hour` is the half-hour-of-day bin (48 per day) in local time (America/Lima).

Positive GEP values are retained rather than clipped: GEP is the residual
`NEE − Reco`, so its error is two-sided and clipping at a physical bound would
bias the retained estimator.

---

## 3. Analytical rules the manuscript follows

**Provenance.** Distinguish measured, gap-filled, modelled, calculated and
derived quantities. Preserve original row identifiers and timestamps. Never
silently drop unusual values or negative nighttime NEE; apply only documented
quality-control rules.

**Interpretation.** Raw boxplots and LOESS curves are descriptive *marginal*
relationships. Do not describe a predictor as an independent control without a
conditional analysis. Assess common environmental support before comparing
seasons, and do not extrapolate conditional curves into unsupported combinations
of drivers.

**Language.** Avoid causal language unless the design supports it. Prefer
"conditional predictor" to "dominant regulator", and "resistance" or
"persistence of uptake" to "resilience" unless recovery was actually analysed.

**Statistics.** Control explicitly for PAR and local time in daytime Ta–VPD
analyses. Evaluate temporal autocorrelation. Report predictor covariance, common
support and model identifiability. Inspect concurvity or collinearity in any
model containing both Ta and VPD. Do not select models on p-values alone, and
distinguish exploratory from confirmatory results — including failed, singular
or non-identifiable fits.

**Resampling.** Use whole dates or multi-day blocks for validation and
bootstrapping. Never randomly split half-hourly observations from the same date
between training and test sets.

---

## 4. Code conventions

- Every figure and table is reproducible from a script, and every reported value
  is written to a machine-readable output.
- Scripts fail clearly when a required variable is absent rather than silently
  producing a shorter result.
- Random seeds are set and reported; session information is recorded.
- Paths resolve through `code/paths.R`; no analysis script hard-codes an
  absolute path.
- Code comments begin with a lower-case letter.
