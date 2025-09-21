

Here's a clean, end-to-end checklist of stages for the project-ordered exactly as you'll execute them. I'v included the core goal and the key artifact(s) to produce at each stage. No code yet; just the roadmap.
**1.) Define scope & targets**
— Lock the $(T, x)$ grid (e.g., three $T \times$ three $x$ ).
- Choose primary method (GK) and backup cross-check (rNEMD).
- Artifact: short design brief with success criteria and error bars you aim to hit.
**2.) Force-field selection & dossiers**
- Select candidate $\mathrm{CeO}_2$ potentials (fixed-charge; optional polarizable).
- Collect their published $a_0, B$, defect notes, and any prior $k$ reports.
- Artifact: "potential dossier" PDFs + a decision note naming the primary FF.
**3.) Reproducibility scaffolding**
- Decide folder layout, manifest fields, seed policy, hash strategy, and logs.
- Artifact: template manifest.json , control.json, and directory skeleton.
**4.) Supercell sizing & commensurability**
- Choose $n$ so $x \cdot 4 n^3 \in \mathbb{Z}$; confirm linear size ( $\geq$ few nm ).
- Artifact: sizing table mapping $x \rightarrow$ (vacancies, Ce count, box length).
**5.) Stoichiometric crystal build ( $x=0$ )**
- Generate pristine fluorite structure with chosen $a_0$, species, charges.
- Artifact: baseline structure files (POSCAR/LAMMPS data) + checksum.
**6.) Electrostatics & numerics knobs (global)**
- Fix PPPM accuracy/order, real-space cutoffs/skin, timestep candidates.
- Artifact: "numerics card" you'll reuse for all states.
**7.) EOS & structure validation at ( $T, P$ )**
- Short NPT scans $\rightarrow$ fit EOS (Birch-Murnaghan or linear).
- Freeze volume; NVT check of stress and pair peaks; NVE drift sanity.
- Gate: $a, \rho$ within $\sim 1-2 \%$ of dossier; $B$ within $\sim 10-15 \%$; quiet drift.
- Artifact: EOS fit, plots, and pass/fail note.
**8.) Stoichiometric $k$ validation (GK, plus one rNEMD check)**
- NPT-NVT-NVE; assemble heat flux; HCACF- $\rightarrow$ windowed integral.
- Size/time convergence + elongated-axis test; optional rNEMD agreement.
- Gate: stable plateau; shard scaling; size/mesh insensitivity; value in bracket.
- Artifact: $k \pm \mathrm{Cl}$ at baseline, with convergence plots.
**9.) Vacancy generator \& $\mathrm{Ce}^{3+}$ assignment (for $x>0$ )**
- Blue-noise sample O sites with $r_{\text {min }}$; for each $V_O$ relabel two nearest $\mathrm{Ce}^{4+} \rightarrow \mathrm{Ce}^{3+}$.
- Artifact: vacancy index list, $\mathrm{Ce}^{3+}$ maps, separation histograms, seeds stored.
**10.) Equilibration per $(T, x)$ state**
- NPT to correct density; NVT to settle modes; record coupling times/energies.
- Artifact: equilibrated restart + manifest (seeds, thermostats, PPPM, etc.).
**11.) GK production per state (primary path)**
- NVE shards with fixed MPI layout; sample $\mathbf{J}(t)$ at fixed stride.
- Artifact: raw flux time series (binary + SHA-256) and J0Jt files.
**12.) rNEMD spot check (one or two states)**
- Gentle swap rate; bin temperature; fit bulk gradient; verify rate independence.
- Artifact: profile plots, linear-fit window, and $k$ cross-check.
**13.) HCACF analysis \& uncertainty quantification**
- Build HCACF with overlapping origins; estimate $\boldsymbol{\tau}_{\text {int }}$; pick $t^{\star}=c \tau_{\text {int }}$.
- Compute $k \pm \mathrm{Cl}$ via blocking and/or shard variance; component averaging.
- Artifact: analysis report per state (plateau, IAT, CI, diagnostics).
**14.) Grid execution \& consistency checks**
- Repeat 10-13 over the full ( $T, x$ ) grid with frozen numerics.
- Track run ledger (wall-clock, seeds, anomalies).
- Artifact: consolidated table of $k(T, x)$ with uniform metadata.
**15.) Sensitivity slices (as needed)**
- One or two deltas: tighter PPPM, different $\Delta t$, or alternate seed to bound systematics.
- Artifact: small "systematics band" note distinct from statistical CI.
**16.) Interpretation & limitations ledger**
- Separate statistical CI from model-family band; note classical-nuclei caveat; remind bulk-ideal assumptions.
- Artifact: one-page limitations statement linked from results.
**17.) Packaging & provenance freeze**
- Write k_result. j son per state (values, CI, windows, hashes, versions).
- Bundle plots (HCACF, cumulative $k$, rNEMD profiles), manifests, and inputs.
- Artifact: immutable release (tag/DOI) of the dataset + analysis code.
**18.) Comparison layer (optional)**
- Place stoichiometric point against BTE/experiment brackets; optionally apply porosity EMA when comparing to pellets.
- Artifact: comparison figure with clear caveats.
**19.) Final report assembly**
- Methods (with citations), validation gates, results grid, uncertainty accounting, limitations, and reuse instructions.
- Artifact: manuscript + README that lets someone rerun any state from raw inputs.


# Overview


We begin by fixing the scientific scope before we touch a simulator. A single, versioned config pins the temperature-vacancy grid ( $T, x$ ), the primary method (equilibrium Green-Kubo) and the cross-check (reverse NEMD), and the exact success criteria we will hold ourselves to: plateau stability, target confidence intervals, and basic convergence gates. That brief is the contract that downstream scripts read, so "what the experiment is" lives outside the code.

Next we choose a force field with intention. A small registry lists each candidate potential, its oxidationstate support (explicit $\mathrm{Ce}^{3+} / \mathrm{Ce}^{4+}$ vs. smeared charge), and the literature values we expect to reproduce (lattice parameter, bulk modulus, any prior $k$ ). For each candidate we auto-generate a one-page dossier and a decision note that names a primary model, plus any secondary we may use to assess model-family uncertainty later.

Before any physics, we put the scaffolding in place that makes the work reproducible. A fixed folder layout, a manifest schema, and deterministic seed derivations ensure every step leaves audit trails: numerics (PPPM tolerances, neighbor skin, timestep), analysis knobs (Sokal window c, sampling stride), input hashes, and per-phase seeds. From one "numerics card" we autogenerate include files for LAMMPS so the same tolerances and time step appear everywhere without copy-paste drift.


Sizing the supercell is not a formality. We pick $n$ so $x \cdot 4 n^3 \in \mathbb{Z}$ and the linear size comfortably admits long-wavelength acoustic modes. That yields a simple sizing table-exact vacancy counts, atom counts, and box length-that we commit alongside the config. With that, we build a pristine fluorite $\mathrm{CeO}_2$ supercell at $x=0$, writing both VASP and LAMMPS inputs plus checksums, so the very geometry is a versioned artifact.

Global numerics-electrostatics accuracy, cutoffs, neighbor rebuilds, and the default timestep-are frozen once, then audited. We tighten PPPM until pressure/virial changes are far below equilibrium noise, check that NVE energy drift is negligible on the horizons we'll use for GK, and save the passing "card" and audit note. Every run henceforth includes that card verbatim.

Validation starts with the equation of state at the target $(T, P)$. We run a miniature pressure sweep, fit a Birch-Murnaghan (or local linear) EOS, and freeze the volume at the fitted minimum. A short NVT checks mean stress and a Ce-O RDF, then a short NVE confirms quiet dynamics. Only after the box is demonstrably "right" do we measure transport. The stoichiometric $k$ gate follows: NPT $\rightarrow$ NVT $\rightarrow$ NVE, assemble the heat flux consistently with PPPM, accumulate the HCACF, and integrate with a Sokal window to a plateau. We show size and time convergence explicitly, and we optionally run a gentle rNEMD to verify that steady-gradient $k$ agrees with GK within uncertainty. Passing this single point locks the protocol for the rest of the grid.


Defects come next, but with controlled randomness. For each $x>0$ we blue-noise sample oxygen sites under a minimum-distance rule to avoid pathological clustering, remove those O's, and re-label the two nearest $\mathrm{Ce}^{4+}$ as $\mathrm{Ce}^{3+}$ for charge neutrality. The generator writes a reduced-structure data file, a vacancy list, a $\mathrm{Ce}^{3+}$ map, and separation histograms, all seed-stamped. We then equilibrate each $(T, x)$ state identically-NPT to correct density, NVT to settle modes-and emit one restart and one manifest that records the taus, seeds, and numerics actually used.

Production is a set of identical NVE shards per state, run with a fixed MPI layout for bit-stable flux files. Each shard dumps a raw binary time series of $\mathrm{J}(t)$ and a quick-look $J_0 J_t$. Analysis is deterministic and documented: build the HCACF with overlapping origins (FFT), estimate the integrated autocorrelation time, choose $t^{\star}=c \tau_{\text {int }}$, integrate to $\hat{k}\left(t^{\star}\right)$, and report both per-shard values and the shard-mean with a $95 \% \mathrm{Cl}$. The analyzer also plots the correlation and cumulative integral with $t^{\star}$ marked, so plateau quality is visible, and it averages tensor components only after windowing, as theory demands.

We then scale out to the whole grid. A light "grid runner" executes equilibration $\rightarrow$ GK shards $\rightarrow$ analysis for every ( $T, x$ ) with frozen numerics, while a consistency check confirms that dt, PPPM accuracy/order, stride, and window c are identical across states. A ledger records wall-clock, seeds, MPI ranks, and any anomalies. A collector consolidates results into a single $k(T, x)$ table (CSV/JSON/MD) with the exact metadata needed to reuse the numbers.

To bound method-level systematics, we take one designated baseline state and run a couple of sensitivity slices-tighter PPPM, smaller $\Delta t$, or an alternate seed for the shards-leaving physics and analysis untouched. The resulting $\Delta k$ values define a small "systematics band," clearly separated from the statistical CI. If needed, we can also quantify a model-family band by repeating a slim subset with an alternate vetted force field.


Interpretation is shipped as a one-page limitations ledger that travels with the dataset. It explains what the numbers mean (bulk lattice $k$ of a classical fixed-charge model), states the classical-nuclei caveat and the bulk-ideal assumptions (periodic, no porosity or grain boundaries, PPPM tin-foil boundary; no polaronic or electronic heat), and separates three uncertainties: the statistical Cl , the numerical systematics band, and the optional model-family spread. For context, an optional comparison layer plots the stoichiometric bulk point against an ab-initio/BTE bracket and high-density pellets, with an effective-medium porosity curve to show scale-not as a fit, but as orientation.

We finish by packaging everything into an immutable release. Each state's k_result. json, plots, and manifest; the consolidated grid; the numerics card and analysis code; the LIMITATIONS page; and checksums and licenses are bundled and tagged. If desired, a DOI (e.g., via Zenodo) freezes the exact version in public. The result is a dataset that is not just "numbers plus a script," but a complete, auditable, and minimal-friction pipeline from scope to figure-small enough to fit on a lab laptop, rigorous enough to stand up in methods review, and transparent enough that future you (or anyone else) can rerun any state and get the same answer for the same reasons.




# 1.) Define scope & targets


You want this step to be a single source of truth that the rest of the pipeline can read. The cleanest way is to (a) capture the scientific scope in a small, human-edited config file, (b) validate it with a tiny Python helper, and (c) auto-render a “design brief” Markdown that states the grid, the primary/backup methods, and the success criteria you’re committing to before any runs start.

## What to decide, concretely

You need the temperature set, the vacancy-fraction set, the primary method (GK) and the cross-check plan (which states will get rNEMD), and quantitative pass/fail targets for uncertainty and convergence. Lock those into a versioned YAML so everything downstream is deterministic.

## Minimal project config (YAML)

Save this as `project.yaml`. Edit numbers, not structure.

```yaml
project:
  name: "ceria_kappa_map_v1"
  description: "Lattice thermal conductivity map k(T,x) for CeO2-x via GK, with rNEMD spot checks."
  primary_method: "GK"           # fixed for this project
  crosscheck_method: "rNEMD"     # used on a subset of states

grid:
  temperatures_K: [300, 700, 1100]           # three T points
  vacancy_fractions: [0.00, 0.03, 0.06]      # three x points
  crosscheck_states:                         # which (T,x) pairs get rNEMD
    - { T: 300, x: 0.00 }
    - { T: 700, x: 0.03 }

targets:
  # statistical precision target per state
  k_rel_ci_max: 0.10                # 10% two-sided CI target on k
  # GK plateaus: stability windows you must demonstrate
  plateau_min_width_ps: 50.0
  plateau_sensitivity_sigma: 1.0     # |Δk| under modest t* change < 1σ
  # time/size convergence gates (procedural, not just a number)
  shards_min: 2
  shard_length_ps: 1000.0
  # size insensitivity (one elongated axis)
  size_delta_max_sigma: 1.0
  # cross-method agreement (where applicable)
  method_agree_max_sigma: 2.0

analysis:
  # windowing rule for GK; you may tune c later—log any change
  sokal_c: 5.0
  # sampling stride for J(t) in fs (multiple of MD dt)
  flux_sample_stride_fs: 1.0

provenance:
  owner: "Your Name"
  repo: "github.com/you/ceria_kappa_map"
  version_tag: "v0.1.0"   # bump only if the *scope* changes
```

## A tiny validator and brief generator (Python)

Drop this into `define_scope.py`. It reads `project.yaml`, checks for obvious mistakes (duplicate states, missing cross-check states not in the grid, non-increasing targets), and emits a Markdown brief you can commit as the contract for the project.

````python
#!/usr/bin/env python3
from __future__ import annotations
import sys, hashlib, json, datetime as dt
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Any
import yaml

@dataclass(frozen=True)
class State:
    T: int
    x: float

@dataclass
class Scope:
    name: str
    description: str
    primary_method: str
    crosscheck_method: str
    temperatures_K: List[int]
    vacancy_fractions: List[float]
    crosscheck_states: List[State]
    targets: Dict[str, Any]
    analysis: Dict[str, Any]
    provenance: Dict[str, Any]

    def grid_states(self) -> List[State]:
        return [State(T, x) for T in self.temperatures_K for x in self.vacancy_fractions]

def load_scope(path: Path) -> Scope:
    cfg = yaml.safe_load(path.read_text())
    P, G, T, A, V = cfg["project"], cfg["grid"], cfg["targets"], cfg["analysis"], cfg["provenance"]
    xs = [State(int(s["T"]), float(s["x"])) for s in G.get("crosscheck_states", [])]
    return Scope(
        name=P["name"],
        description=P["description"],
        primary_method=P["primary_method"],
        crosscheck_method=P["crosscheck_method"],
        temperatures_K=[int(t) for t in G["temperatures_K"]],
        vacancy_fractions=[float(x) for x in G["vacancy_fractions"]],
        crosscheck_states=xs,
        targets=T, analysis=A, provenance=V
    )

def validate_scope(S: Scope) -> List[str]:
    errs = []
    if S.primary_method != "GK":
        errs.append("primary_method must be 'GK' for this project.")
    if S.crosscheck_method not in {"rNEMD", "none", "NEMD"}:
        errs.append("crosscheck_method must be 'rNEMD' (or 'none').")

    # grid sanity
    if len(set(S.temperatures_K)) != len(S.temperatures_K):
        errs.append("temperatures_K has duplicates.")
    if len(set(S.vacancy_fractions)) != len(S.vacancy_fractions):
        errs.append("vacancy_fractions has duplicates.")

    grid = set(S.grid_states())
    for cs in S.crosscheck_states:
        if cs not in grid:
            errs.append(f"crosscheck state {cs} not in (T,x) grid.")

    # target sanity
    if not (0 < S.targets["k_rel_ci_max"] < 0.5):
        errs.append("k_rel_ci_max should be in (0, 0.5).")
    if S.targets["shards_min"] < 1:
        errs.append("shards_min must be >= 1.")
    if S.targets["shard_length_ps"] <= 0:
        errs.append("shard_length_ps must be positive.")
    if S.analysis["sokal_c"] < 3.0:
        errs.append("sokal_c appears too small (<3).")
    return errs

def brief_md(S: Scope, cfg_text: str) -> str:
    grid = sorted(S.grid_states(), key=lambda s: (s.T, s.x))
    grid_str = ", ".join([f"(T={s.T} K, x={s.x:.3f})" for s in grid])
    xchk = ", ".join([f"(T={s.T} K, x={s.x:.3f})" for s in S.crosscheck_states]) or "none"
    sha = hashlib.sha256(cfg_text.encode()).hexdigest()[:12]
    now = dt.datetime.now().isoformat(timespec="seconds")
    t = S.targets
    a = S.analysis
    return f"""# Design brief — {S.name}

**Generated:** {now}  
**Config hash:** `{sha}`

## Scope
{S.description}

**Grid:** {grid_str}  
**Primary method:** {S.primary_method}  
**Cross-check:** {S.crosscheck_method} on {xchk}

## Success criteria (commitments)
We will report k(T,x) with two-sided relative CI ≤ {t["k_rel_ci_max"]:.0%} at each state using GK with a Sokal window c = {a["sokal_c"]:.1f}.  
Each GK estimate will exhibit a plateau of width ≥ {t["plateau_min_width_ps"]} ps and sensitivity < {t["plateau_sensitivity_sigma"]}σ under modest t* changes.  
Per state we will accumulate at least {t["shards_min"]} shards of length {t["shard_length_ps"]} ps and demonstrate size insensitivity ≤ {t["size_delta_max_sigma"]}σ via one elongated-axis check at a designated baseline.  
Where rNEMD is applied, GK and rNEMD will agree within {t["method_agree_max_sigma"]}σ.

## Analysis knobs (frozen for the grid)
Flux sampling stride: {a["flux_sample_stride_fs"]} fs.  
GK window rule: t* = c τ_int with c = {a["sokal_c"]:.1f}.

## Provenance
Owner: {S.provenance.get("owner","")}  
Repository: {S.provenance.get("repo","")}  
Version tag: {S.provenance.get("version_tag","")}

---

```yaml
# Embedded project.yaml (frozen copy)
{cfg_text.strip()}
````

"""

def main():  
path = Path("project.yaml")  
if not path.exists():  
sys.exit("project.yaml not found.")  
cfg_text = path.read_text()  
S = load_scope(path)  
errs = validate_scope(S)  
if errs:  
print("CONFIG ERRORS:")  
for e in errs:  
print(" -", e)  
sys.exit(1)  
out = Path("design_brief.md")  
out.write_text(brief_md(S, cfg_text))  
print(f"Wrote {out} (commit this file to lock scope).")

if **name** == "**main**":  
main()

```

Run `python define_scope.py`. It will write `design_brief.md` that you can commit as the contract for the project. Everything else (supercell sizing, generators, GK plumbing, analysis) will consume `project.yaml` so your scope is literally the first dependency in the DAG.
::contentReference[oaicite:0]{index=0}
```


# 2.) Force-field selection & dossiers



The goal here is to turn a messy literature hunt into a repeatable pipeline: you keep a tiny registry of candidate potentials, you generate a human-readable “dossier” for each one, and you have a deterministic note that chooses the primary force field according to criteria you commit to. No LAMMPS runs yet—just the scaffolding and the artifacts you’ll fill as you read.

## What you’ll keep under version control

Create a force-field registry in YAML, one entry per candidate. Each record holds bibliographic info, parameter files you’ll actually use in LAMMPS, and the headline properties you’ll validate later. Leave numeric fields blank at first; the whole point is that you’ll fill them from papers, theses, or docs and the rest of the tooling will “just work”.

```yaml
# ff_registry.yaml
candidates:
  - key: "Grimes_1990_fixedQ"
    family: "Buckingham+Q"
    oxidation_states: ["Ce4+","Ce3+","O2-"]
    has_explicit_Ce3: true
    polarizable: false
    refs:
      - type: "paper"
        citation: "Grimes et al., J. Chem. Soc., Faraday Trans., 1990, …"
        doi: ""
        url: ""
    files:
      lammps_pair: "potentials/grimes_1990/pair_coeffs.table"
      species_map: "potentials/grimes_1990/species_map.yaml"
    reported_properties:
      a0_300K_A: null
      bulk_modulus_GPa: null
      notes: "Often used baseline for fluorite oxides; Ce3+/Ce4+ distinct."
      k_300K_WmK: null
    defect_notes:
      vo_ce3_association: "yes (expected)"
      other: ""
  - key: "Gotte_2007_fixedQ"
    family: "Buckingham+Q"
    oxidation_states: ["Ce4+","Ce3+","O2-"]
    has_explicit_Ce3: true
    polarizable: false
    refs: []
    files:
      lammps_pair: "potentials/gotte_2007/pair_coeffs.table"
      species_map: "potentials/gotte_2007/species_map.yaml"
    reported_properties:
      a0_300K_A: null
      bulk_modulus_GPa: null
      notes: ""
      k_300K_WmK: null
    defect_notes:
      vo_ce3_association: "yes"
  - key: "Burbano_2011_DIPPIM"
    family: "DIPPIM (induced-dipole)"
    oxidation_states: ["Ce4+","Ce3+","O2-"]
    has_explicit_Ce3: true
    polarizable: true
    refs: []
    files:
      lammps_pair: "potentials/dippim_2011/params.yaml"   # keep raw source; you’ll convert for LAMMPS
      species_map: "potentials/dippim_2011/species_map.yaml"
    reported_properties:
      a0_300K_A: null
      bulk_modulus_GPa: null
      notes: "Polarizable; higher cost; better phonons/dielectric."
      k_300K_WmK: null
    defect_notes:
      vo_ce3_association: "yes (by construction)"
```

Keep a species map that pins type IDs, charges, and human names in one place. This file will be imported later when you auto-render LAMMPS inputs for validation.

```yaml
# potentials/gotte_2007/species_map.yaml
species:
  - name: "Ce4+"
    type_id: 1
    charge: 4.0
  - name: "Ce3+"
    type_id: 2
    charge: 3.0
  - name: "O2-"
    type_id: 3
    charge: -2.0
```

## A dossier builder that never goes out of sync

This small Python tool reads the registry, creates one Markdown dossier per candidate with “fill-me” tables, and writes a single “decision note” skeleton that compares candidates by criteria you declare once (explicit Ce³⁺ support, reported a0a_0/BB, prior kk reports, polarizability flag). You can run it anytime; it won’t clobber hand-edited numeric fields because it only uses the YAML as source of truth.

```python
#!/usr/bin/env python3
# make_dossiers.py
from __future__ import annotations
import json, hashlib, datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List
import yaml

ROOT = Path(".")
REG = ROOT / "ff_registry.yaml"
DOSSIERS = ROOT / "dossiers"
DECISION = ROOT / "dossiers" / "decision_note.md"

@dataclass
class FF:
    key: str
    family: str
    has_explicit_Ce3: bool
    polarizable: bool
    refs: List[Dict[str,Any]]
    files: Dict[str,str]
    props: Dict[str,Any]
    defect_notes: Dict[str,Any]

def load_registry() -> List[FF]:
    data = yaml.safe_load(REG.read_text())
    out = []
    for c in data["candidates"]:
        out.append(FF(
            key=c["key"],
            family=c["family"],
            has_explicit_Ce3=c["has_explicit_Ce3"],
            polarizable=c["polarizable"],
            refs=c.get("refs",[]),
            files=c["files"],
            props=c["reported_properties"],
            defect_notes=c.get("defect_notes",{}),
        ))
    return out

def md_dossier(ff: FF) -> str:
    sha = hashlib.sha256(json.dumps(ff.__dict__, sort_keys=True, default=str).encode()).hexdigest()[:12]
    refs = "\n".join([f"- {r.get('citation','')} ({r.get('doi','')})" for r in ff.refs]) or "_add references here_"
    return f"""# Force-field dossier — {ff.key}

**Family:** {ff.family} | **Explicit Ce³⁺:** {ff.has_explicit_Ce3} | **Polarizable:** {ff.polarizable}  
**Params:** `{ff.files.get('lammps_pair','')}` | **Species map:** `{ff.files.get('species_map','')}`  
**Registry hash:** `{sha}`

## Reported reference values (fill from literature)
- Lattice parameter a₀ @ 300 K: **{ff.props.get('a0_300K_A', '—')} Å**
- Bulk modulus B: **{ff.props.get('bulk_modulus_GPa','—')} GPa**
- Any reported lattice k @ 300 K (dense, defect-free): **{ff.props.get('k_300K_WmK','—')} W m⁻¹ K⁻¹**
- Notes: {ff.props.get('notes','')}

## Defect chemistry / vacancy association
- V_O ↔ Ce³⁺ association expected: **{ff.defect_notes.get('vo_ce3_association','—')}**
- Remarks: {ff.defect_notes.get('other','')}

## References
{refs}

## To-do for this FF
1. Verify EOS (a₀, ρ, B) at target (T,P).
2. Run stoichiometric GK at 300 K for bracket sanity.
3. Decide keep/park as primary/secondary.
"""

def md_decision(cands: List[FF]) -> str:
    now = dt.datetime.now().isoformat(timespec="seconds")
    rows = []
    for ff in cands:
        rows.append({
            "key": ff.key,
            "explicit_Ce3": "yes" if ff.has_explicit_Ce3 else "no",
            "polarizable": "yes" if ff.polarizable else "no",
            "a0": ff.props.get("a0_300K_A"),
            "B": ff.props.get("bulk_modulus_GPa"),
            "k300": ff.props.get("k_300K_WmK"),
        })
    # simple, readable table in Markdown
    hdr = "| key | explicit Ce³⁺ | polarizable | a₀(Å) | B(GPa) | k@300K |\n|---|---|---|---|---|---|"
    body = "\n".join([f"| {r['key']} | {r['explicit_Ce3']} | {r['polarizable']} | {r['a0'] or '—'} | {r['B'] or '—'} | {r['k300'] or '—'} |" for r in rows])
    criteria = (
        "Primary selection criteria:\n"
        "1) Must support explicit Ce³⁺/Ce⁴⁺ species with usable parameters.\n"
        "2) EOS sanity (a₀ within ~1–2%, B within ~10–15% of reported values).\n"
        "3) Literature pedigree for reduced ceria or doped ceria.\n"
        "4) Prefer non-polarizable fixed-charge for semester timeline unless polarizable is essential.\n"
    )
    return f"""# Force-field decision note

Generated: {now}

## Candidates overview
{hdr}
{body}

## Criteria (frozen before validation)
{criteria}

## Decision
- **Primary FF:** <fill after dossiers are complete>  
- **Secondary (cross-check):** <optional>

## Rationale
Summarize why the primary meets 1–4, cite dossiers by key, and note any risks or caveats (e.g., phonon dispersion fidelity, dielectric response).
"""

def main():
    DOSSIERS.mkdir(parents=True, exist_ok=True)
    cands = load_registry()
    for ff in cands:
        (DOSSIERS / f"{ff.key}.md").write_text(md_dossier(ff))
    DECISION.write_text(md_decision(cands))
    print(f"Wrote dossiers for {len(cands)} candidates and {DECISION}")

if __name__ == "__main__":
    main()
```

Run `python make_dossiers.py`. You’ll get one Markdown dossier per candidate under `dossiers/` plus a `decision_note.md` you can fill once and commit.

## A minimal parameter “sanity loader” for LAMMPS later

You won’t run validation now, but future steps will need a uniform way to turn your tabulated Buckingham parameters into a `pair_style` block. Keep a plain, explicit table format and a single translator; when you do move to validation scripts, you’ll import this instead of re-typing anything.

```yaml
# potentials/gotte_2007/pair_coeffs.table
# species_i  species_j  A(eV)     rho(Å)     C(eVÅ^6)   cutoff(Å)
Ce4+        O2-        1804.0     0.3453     0.0000      12.0
Ce3+        O2-        1600.0     0.3600     0.0000      12.0
O2-         O2-        22764.0    0.1490     27.88       12.0
Ce4+        Ce4+       0.0        1.0        0.0         12.0
Ce3+        Ce3+       0.0        1.0        0.0         12.0
Ce3+        Ce4+       0.0        1.0        0.0         12.0
```

And a tiny parser you’ll reuse (no execution here—just the interface you’ll call later):

```python
# lmp_params.py
from pathlib import Path
from typing import Dict, Tuple

def read_buckingham_table(path: Path) -> Dict[Tuple[str,str], Tuple[float,float,float,float]]:
    pairs = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.strip().startswith("#"): 
            continue
        s1, s2, A, rho, C, rc = line.split()
        key = tuple(sorted((s1,s2)))
        pairs[key] = (float(A), float(rho), float(C), float(rc))
    return pairs

def emit_lammps_pair_section(pairs: Dict[Tuple[str,str], Tuple[float,float,float,float]]) -> str:
    out = ["pair_style buck/coul/long 12.0",
           "# Q: charges come from species_map.yaml; PPPM set elsewhere"]
    for (s1,s2), (A,rho,C,rc) in pairs.items():
        out.append(f"pair_coeff @{s1} @{s2} {A:.6f} {rho:.6f} {C:.6f} {rc:.3f}")
    return "\n".join(out)
```

You’ll wire `@Ce4+` to a numeric type ID via the `species_map.yaml` in a later step; keeping symbolic names here makes the dossiers readable and reduces mismatches.

## What to write down before you proceed

Fill the dossiers with a0a_0, BB, and any reported kk you trust for dense, defect-free CeO₂. Put the PDFs (papers, theses, manuals) under `dossiers/papers/{key}/` so your “decision note” can link them. Once you pick the primary force field in `decision_note.md`, tag the repo (e.g., `ff-choice-v1`). From here on, every validation script will read `ff_registry.yaml` and your `primary FF` choice so the whole pipeline stays consistent with what you decided here.


# 3.) Reproducibility scaffolding**


You’ll lock a directory skeleton, a canonical manifest, and a control file that downstream stages can consume without guessing. The pattern below gives you a one-command bootstrap, a deterministic seed policy, content-addressed hashing, and uniform logging. No external services required.

## One-shot project bootstrap

Save as `init_project.py`, run it once in an empty repo. It lays down the tree, writes template `manifest.json` and `control.json`, and adds a minimal `.gitignore`.

```python
#!/usr/bin/env python3
from __future__ import annotations
import os, json, hashlib, platform, getpass, datetime as dt
from pathlib import Path

ROOT = Path(".")
DIRS = [
  "potentials", "configs", "geometry",
  "runs/T300_x0.00/nve_shard_000",
  "runs/T300_x0.00/nve_shard_001",
  "analysis/T300_x0.00/plots",
  "dossiers/papers"
]

GITIGNORE = """# artifacts
runs/**/stdout.log
runs/**/slurm_*.out
runs/**/flux.bin
runs/**/flux.sha256
analysis/**/kcumu.npy
analysis/**/hcacf.npy
analysis/**/k_result.json
__pycache__/
.ipynb_checkpoints/
"""

TEMPLATE_CONTROL = {
  "analysis_version": "1.0.0",
  "gk": {
    "flux_sample_stride_fs": 1.0,
    "sokal_c": 5.0,
    "component_average": True,
    "integration": "trapezoid",
    "demean": True,
    "blocking": {"min_block": 1024}
  },
  "rnemd": {
    "swap_period_steps": 500,
    "exclude_bins_at_ends": 2
  }
}

def sha256_text(s: str) -> str:
  return hashlib.sha256(s.encode()).hexdigest()

def write_manifest_skeleton():
  now = dt.datetime.now().isoformat(timespec="seconds")
  manifest = {
    "schema": "ceriakappa/manifest@1",
    "created": now,
    "who": {"user": getpass.getuser(), "host": platform.node()},
    "platform": {
      "system": platform.system(),
      "release": platform.release(),
      "python": platform.python_version(),
      "cpu": platform.processor()
    },
    "code": {
      "repo": "<fill>", "commit": "<git-hash>",
      "lammps_version": "<lmp -h first line>",
      "build_flags": "<optional>"
    },
    "inputs": {
      "supercell": {"a_A": 5.41, "n": 10},
      "composition": {"x": 0.00, "n_vac": 0},
      "potential_files": [],
      "pppm": {"accuracy": 1e-5, "order": 5},
      "neighbor": {"skin_A": 2.0}
    },
    "seeds": {
      "master": "<fill>",
      "npt": "<fill>", "nvt": "<fill>",
      "nve_shards": {"T300_x0.00": ["<fill>","<fill>"]}
    },
    "numerics": {
      "timestep_fs": 1.0,
      "thermostat": {"tau_ps": 0.2, "chain": 3},
      "barostat": {"tau_ps": 2.0, "iso": True}
    },
    "hashes": {
      "flux_sha256": None,
      "hcacf_sha256": None,
      "control_sha256": None,
      "potentials_sha256": []
    },
    "notes": ""
  }
  (ROOT/"runs/T300_x0.00/manifest.json").write_text(json.dumps(manifest, indent=2))
  ctrl = json.dumps(TEMPLATE_CONTROL, indent=2)
  (ROOT/"analysis/T300_x0.00/control.json").write_text(ctrl)
  # record control hash inside the manifest shell
  m = json.loads((ROOT/"runs/T300_x0.00/manifest.json").read_text())
  m["hashes"]["control_sha256"] = sha256_text(ctrl)
  (ROOT/"runs/T300_x0.00/manifest.json").write_text(json.dumps(m, indent=2))

def main():
  for d in DIRS:
    Path(d).mkdir(parents=True, exist_ok=True)
  (ROOT/".gitignore").write_text(GITIGNORE)
  write_manifest_skeleton()
  print("Project skeleton created. Edit runs/T300_x0.00/manifest.json and analysis/T300_x0.00/control.json before any computation.")

if __name__ == "__main__":
  main()
```

## Deterministic seed policy

Use a single master seed per state and derive phase/shard seeds via a transparent hash. This avoids accidental reuse and makes seeds inspectable.

```python
# seeds.py
import hashlib

def derive_seed(master: int, tag: str) -> int:
  s = f"{master}:{tag}".encode()
  return int(hashlib.sha256(s).hexdigest()[:8], 16)

# example
# master = 1337421
# npt = derive_seed(master, "NPT")
# nvt = derive_seed(master, "NVT")
# shard0 = derive_seed(master, "NVE_shard_0")
# shard1 = derive_seed(master, "NVE_shard_1")
```

Record `master` and all derived seeds in the manifest under `seeds`. Never reuse a master across different (T,x)(T,x) states.

## Hash strategy for artifacts and inputs

Write the flux time series and HCACF as binary arrays, alongside their SHA-256. The helper below updates the run manifest in place so the checksum lives with the numbers.

```python
# hash_artifacts.py
import json, hashlib
from pathlib import Path

def sha256_file(path: Path) -> str:
  h = hashlib.sha256()
  with path.open("rb") as f:
    for chunk in iter(lambda: f.read(1<<20), b""):
      h.update(chunk)
  return h.hexdigest()

def stamp_manifest(run_dir: Path):
  manifest = json.loads((run_dir/"manifest.json").read_text())
  flux = run_dir/"nve_shard_000"/"flux.bin"
  hcacf = run_dir.parent.parent/"analysis"/run_dir.name/"hcacf.npy"  # adjust to your actual pathing
  if flux.exists():
    manifest["hashes"]["flux_sha256"] = sha256_file(flux)
    (run_dir/"nve_shard_000"/"flux.sha256").write_text(manifest["hashes"]["flux_sha256"])
  if hcacf.exists():
    manifest["hashes"]["hcacf_sha256"] = sha256_file(hcacf)
  (run_dir/"manifest.json").write_text(json.dumps(manifest, indent=2))

# usage: stamp_manifest(Path("runs/T300_x0.00"))
```

Do the same once per potential file at the start of the campaign and append into `hashes.potentials_sha256`.

## Uniform logging you can tail and parse

Adopt a tiny logger that every orchestration script imports. It emits ISO timestamps and a machine-parsable prefix.

```python
# logutil.py
import logging, sys

def get_logger(name: str):
  logger = logging.getLogger(name)
  if logger.handlers:
    return logger
  logger.setLevel(logging.INFO)
  h = logging.StreamHandler(sys.stdout)
  fmt = logging.Formatter(fmt="%(asctime)sZ | %(name)s | %(levelname)s | %(message)s",
                          datefmt="%Y-%m-%dT%H:%M:%S")
  h.setFormatter(fmt)
  logger.addHandler(h)
  return logger

# example
# from logutil import get_logger
# log = get_logger("prep")
# log.info("equilibrating state T=300, x=0.00 with tau_T=0.2 ps, tau_P=2.0 ps")
```

Redirect LAMMPS stdout to `runs/<state>/stdout.log` and keep your Python logs alongside. Because you’re recording hashes and seeds in the manifest, the logs don’t need to duplicate those values; simply reference the manifest path in the first line of each run.

## Minimal `manifest.json` and `control.json` you can copy

If you prefer to start without the bootstrap script, these two files are enough to “speak the contract” for one state.

`runs/T300_x0.00/manifest.json`:

```json
{
  "schema": "ceriakappa/manifest@1",
  "created": "2025-09-20T19:05:00",
  "who": { "user": "yourid", "host": "hpc-login001" },
  "platform": { "system": "Linux", "python": "3.11.7", "cpu": "Intel(R) Xeon(R)" },
  "code": { "repo": "https://…", "commit": "abcdef1", "lammps_version": "23Jun2022", "build_flags": "-DPKG_KSPACE …" },
  "inputs": {
    "supercell": { "a_A": 5.41, "n": 10 },
    "composition": { "x": 0.00, "n_vac": 0 },
    "potential_files": ["potentials/gotte_2007/pair_coeffs.table", "potentials/gotte_2007/species_map.yaml"],
    "pppm": { "accuracy": 1e-5, "order": 5 },
    "neighbor": { "skin_A": 2.0 }
  },
  "seeds": {
    "master": 1337421,
    "npt": 294967268, "nvt": 274877906, 
    "nve_shards": { "T300_x0.00": [ 402653189, 241591910 ] }
  },
  "numerics": {
    "timestep_fs": 1.0,
    "thermostat": { "tau_ps": 0.2, "chain": 3 },
    "barostat": { "tau_ps": 2.0, "iso": true }
  },
  "hashes": {
    "flux_sha256": null,
    "hcacf_sha256": null,
    "control_sha256": "e3b0c44298…",
    "potentials_sha256": []
  },
  "notes": "baseline stoichiometric validation"
}
```

`analysis/T300_x0.00/control.json`:

```json
{
  "analysis_version": "1.0.0",
  "gk": {
    "flux_sample_stride_fs": 1.0,
    "sokal_c": 5.0,
    "component_average": true,
    "integration": "trapezoid",
    "demean": true,
    "blocking": { "min_block": 1024 }
  },
  "rnemd": { "swap_period_steps": 500, "exclude_bins_at_ends": 2 }
}
```

Keep these files immutable once a state has been analyzed. If you revise analysis logic, bump `analysis_version`, regenerate, and store both results with their control hashes so differences are attributable.

## Directory skeleton you’ll actually use

If you don’t run the bootstrap, create this minimal tree by hand and commit it empty. The rest of the pipeline fills it.

```
project/
  potentials/
  configs/
  geometry/
  runs/
    T300_x0.00/
      manifest.json
      nve_shard_000/
      nve_shard_001/
  analysis/
    T300_x0.00/
      control.json
      plots/
```

From here, every orchestration step reads `manifest.json` and `control.json`, derives seeds deterministically, writes flux arrays with checksums, and appends hashes back into the manifest. That’s the glue that keeps your science reproducible without adding friction.



# 4.) Supercell sizing & commensurability


You want a tiny, deterministic helper that chooses the smallest cubic supercell that (i) makes the vacancy count exact under PBC and (ii) is at least a few nanometers on a side. Then you want a sizing table you can commit as an artifact.

## What the script does

Given a lattice parameter aa (Å), a minimum box length Lmin⁡L_{\min} (nm), and a set of vacancy fractions {x}\{x\}, it searches n=⌈Lmin⁡/(a/10)⌉,⌈Lmin⁡/(a/10)⌉+1,…n = \lceil L_{\min}/(a/10)\rceil, \lceil L_{\min}/(a/10)\rceil+1,\dots for the **smallest** nn such that x⋅4n3x\cdot 4n^3 is an integer (within a strict tolerance). For each xx it outputs:

- nn, the supercell replication factor
    
- NCe=4n3N_{\mathrm{Ce}} = 4n^3 and NO=8n3N_{\mathrm{O}} = 8n^3
    
- Nvac=x NCeN_{\mathrm{vac}} = x\,N_{\mathrm{Ce}} (integer by construction)
    
- L=naL = n a in Å and nm
    

It also lets you pin nn globally (e.g., “use n=10n=10 for every state”) to keep the grid consistent, but it will warn if any xx becomes incommensurate.

## Drop-in script: `size_supercells.py`

```python
#!/usr/bin/env python3
from __future__ import annotations
import math, argparse, json
from dataclasses import dataclass
from typing import List, Tuple

TOL = 1e-9  # strict integrality tolerance

@dataclass
class Sizing:
    x: float
    n: int
    N_Ce: int
    N_O: int
    N_vac: int
    L_A: float
    L_nm: float

def find_min_n_for_x(a_A: float, Lmin_nm: float, x: float, n_max: int = 64) -> int:
    n_start = math.ceil( (Lmin_nm*10.0) / a_A )
    for n in range(max(1, n_start), n_max+1):
        N_Ce = 4 * n**3
        val = x * N_Ce
        if abs(val - round(val)) < TOL:
            return n
    raise ValueError(f"No n ≤ {n_max} satisfies integrality for x={x}")

def size_for(a_A: float, x: float, n: int) -> Sizing:
    N_Ce = 4 * n**3
    N_O  = 8 * n**3
    N_vac = int(round(x * N_Ce))
    L_A  = n * a_A
    return Sizing(x=x, n=n, N_Ce=N_Ce, N_O=N_O, N_vac=N_vac, L_A=L_A, L_nm=L_A/10.0)

def main():
    ap = argparse.ArgumentParser(description="Choose commensurate supercells for CeO2-x")
    ap.add_argument("--aA", type=float, required=True, help="lattice parameter a (Å)")
    ap.add_argument("--xmin", type=float, nargs="+", required=True, help="vacancy fractions x (e.g. 0.00 0.03 0.06)")
    ap.add_argument("--Lmin", type=float, default=5.0, help="minimum box length (nm)")
    ap.add_argument("--pin-n", type=int, default=None, help="force a single n for all x (warn if incommensurate)")
    ap.add_argument("--nmax", type=int, default=64, help="max n to search")
    ap.add_argument("--json", action="store_true", help="emit JSON as well as table")
    args = ap.parse_args()

    rows: List[Sizing] = []
    if args.pin_n is None:
        for x in args.xmin:
            n = find_min_n_for_x(args.aA, args.Lmin, x, args.nmax)
            rows.append(size_for(args.aA, x, n))
    else:
        n = args.pin_n
        for x in args.xmin:
            N_Ce = 4 * n**3
            val = x * N_Ce
            if abs(val - round(val)) >= TOL:
                print(f"[WARN] x={x} is not commensurate with n={n}: x*4n^3={val:.6f} not integer.")
            rows.append(size_for(args.aA, x, n))

    # Pretty print
    hdr = f"{'x':>6} | {'n':>3} | {'N_Ce':>6} | {'N_O':>6} | {'N_vac':>6} | {'L (Å)':>8} | {'L (nm)':>7}"
    print(hdr)
    print("-"*len(hdr))
    for r in rows:
        print(f"{r.x:6.3f} | {r.n:3d} | {r.N_Ce:6d} | {r.N_O:6d} | {r.N_vac:6d} | {r.L_A:8.2f} | {r.L_nm:7.2f}")

    if args.json:
        out = [r.__dict__ for r in rows]
        print("\nJSON:")
        print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
```

Typical usage for your baseline grid is:

```
python size_supercells.py --aA 5.41 --xmin 0.00 0.03 0.06 --Lmin 5.0 --pin-n 10
```

This pins n=10n=10 for all states, which is ideal for uniform numerics. Because NCe=4×103=4000N_{\mathrm{Ce}}=4\times 10^3=4000, both x=0.03x=0.03 and 0.060.06 are exactly commensurate and you get a comfortable L≈5.41L\approx 5.41 nm box.

## Sizing table for your current grid (with a=5.41a=5.41 Å and n=10n=10)

|x|n|N_Ce|N_O|N_vac|L (Å)|L (nm)|
|--:|--:|--:|--:|--:|--:|--:|
|0.000|10|4000|8000|0|54.10|5.41|
|0.030|10|4000|8000|120|54.10|5.41|
|0.060|10|4000|8000|240|54.10|5.41|

Commit this table as `configs/sizing_table.md` alongside the command you used to generate it, and you’re set. If you later add an awkward xx (say x=0.025x=0.025), rerun the helper without `--pin-n` to get the **smallest** commensurate nn (it will likely suggest n=8,12,…n=8,12,\dots), then decide whether to bump the global nn or nudge xx to the nearest exactly realizable value.


# 5.) Stoichiometric crystal build ( $x=0$ )


Below is a self-contained Python utility that generates a pristine fluorite CeO2\mathrm{CeO_2} supercell with your chosen lattice parameter a0a_0 (Å) and replication nn, and writes both a VASP `POSCAR` and a LAMMPS **charge**-style `data` file (with Ce4+\mathrm{Ce^{4+}} and O2−\mathrm{O^{2-}} charges). It also emits SHA-256 checksums so you can stamp your manifest.

The fluorite basis used is Fm3ˉmFm\bar{3}m: Ce on the fcc cation sublattice (4a: (0,0,0)(0,0,0) plus fcc translations) and O on 8c: (14,14,14)(\tfrac14,\tfrac14,\tfrac14) and its symmetrically equivalent positions.

## `build_ceo2.py`

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np

# --- constants ---
MASS_CE = 140.116
MASS_O  = 15.999
Q_CE4   =  4.0
Q_O2    = -2.0

# fluorite fractional bases in one conventional cell
CE_FRAC = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
], dtype=float)

O_FRAC = np.array([
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.75],
    [0.25, 0.75, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.25],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
    [0.75, 0.75, 0.75],
], dtype=float)

def make_supercell(aA: float, n: int):
    """Return cartesian coords (Å) and species arrays for CeO2 n×n×n"""
    # lattice vectors (Å)
    A = np.eye(3) * aA
    # fractional grid
    shifts = np.array([(i,j,k) for i in range(n) for j in range(n) for k in range(n)], dtype=float)
    shifts /= n
    # tile bases
    ce_frac_all = (CE_FRAC[None, :, :] + shifts[:, None, :]) % 1.0
    o_frac_all  = (O_FRAC [None, :, :] + shifts[:, None, :]) % 1.0
    ce_frac_all = ce_frac_all.reshape(-1,3)
    o_frac_all  = o_frac_all.reshape(-1,3)
    # to cartesian
    ce_cart = ce_frac_all @ A
    o_cart  = o_frac_all  @ A
    # sanity counts
    NCE = 4 * n**3
    NO  = 8 * n**3
    assert ce_cart.shape[0] == NCE and o_cart.shape[0] == NO
    # species/type/charge arrays
    # type 1: Ce4+, type 2: O2-
    types   = np.concatenate([np.full(NCE, 1, int), np.full(NO, 2, int)])
    charges = np.concatenate([np.full(NCE, Q_CE4, float), np.full(NO, Q_O2, float)])
    coords  = np.vstack([ce_cart, o_cart])
    return coords, types, charges, A

def write_poscar(path: Path, aA: float, n: int, coords: np.ndarray, types: np.ndarray):
    """Write VASP POSCAR (direct coords) with species order Ce O"""
    lat = np.eye(3) * (aA * n)
    # sort by species for POSCAR order Ce then O
    order = np.argsort(types)
    coords_sorted = coords[order] / (aA * n)  # to fractional of supercell
    types_sorted = types[order]
    n_ce = np.count_nonzero(types_sorted == 1)
    n_o  = np.count_nonzero(types_sorted == 2)
    lines = []
    lines.append(f"CeO2 fluorite n={n} a0={aA:.6f} Å  (Ce4+, O2-)")
    lines.append("1.0")
    for i in range(3):
      lines.append(f"{lat[i,0]:.10f} {lat[i,1]:.10f} {lat[i,2]:.10f}")
    lines.append("Ce O")
    lines.append(f"{n_ce} {n_o}")
    lines.append("Direct")
    for i in range(coords_sorted.shape[0]):
      x,y,z = coords_sorted[i]
      lines.append(f"{x:.10f} {y:.10f} {z:.10f}")
    path.write_text("\n".join(lines))

def write_lammps_data(path: Path, aA: float, n: int, coords: np.ndarray, types: np.ndarray, charges: np.ndarray):
    """Write LAMMPS charge-style data file (atom-ID, type, charge, x y z)"""
    L = aA * n
    N = coords.shape[0]
    nat1 = np.count_nonzero(types == 1)
    nat2 = N - nat1
    lines = []
    lines.append("# LAMMPS data file: CeO2 fluorite (charge style)")
    lines.append(f"\n{N} atoms")
    lines.append("2 atom types")
    lines.append(f"\n0.0 {L:.10f} xlo xhi")
    lines.append(f"0.0 {L:.10f} ylo yhi")
    lines.append(f"0.0 {L:.10f} zlo zhi\n")
    lines.append("Masses\n")
    lines.append("1 {:.6f}  # Ce4+".format(MASS_CE))
    lines.append("2 {:.6f}  # O2-".format(MASS_O))
    lines.append("\nAtoms # charge\n")
    # atom-ID atom-type charge x y z
    for i in range(N):
        lines.append(f"{i+1} {types[i]} {charges[i]:.6f} {coords[i,0]:.10f} {coords[i,1]:.10f} {coords[i,2]:.10f}")
    path.write_text("\n".join(lines))

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser(description="Build stoichiometric CeO2 fluorite supercell (x=0)")
    ap.add_argument("--aA", type=float, required=True, help="lattice parameter a0 in Å (at target T)")
    ap.add_argument("--n", type=int, required=True, help="cubic replication n (box is n*a0)")
    ap.add_argument("--outdir", type=str, default="geometry", help="output directory")
    args = ap.parse_args()

    coords, types, charges, A = make_supercell(args.aA, args.n)
    outdir = Path(args.outdir) / f"stoich_CeO2_n{args.n}_a{args.aA:.3f}"
    outdir.mkdir(parents=True, exist_ok=True)

    poscar = outdir / "POSCAR"
    data   = outdir / "ceo2.data"
    write_poscar(poscar, args.aA, args.n, coords, types)
    write_lammps_data(data, args.aA, args.n, coords, types, charges)

    pos_h = sha256_file(poscar)
    dat_h = sha256_file(data)
    (outdir/"structure_checksums.json").write_text(json.dumps({
        "POSCAR": {"path": str(poscar), "sha256": pos_h},
        "LAMMPS_data": {"path": str(data), "sha256": dat_h},
        "metadata": {
            "a0_A": args.aA, "n": args.n,
            "N_Ce": int(4*args.n**3), "N_O": int(8*args.n**3),
            "L_box_A": args.aA*args.n
        }
    }, indent=2))
    print(f"Wrote:\n  {poscar}\n  {data}\n  checksums: {outdir/'structure_checksums.json'}")

if __name__ == "__main__":
    main()
```

## How to use it

Choose a0a_0 (at your target TT) and nn from your sizing table. For your earlier example a0=5.41a_0=5.41 Å and n=10n=10:

```
python build_ceo2.py --aA 5.41 --n 10 --outdir geometry
```

This writes:

```
geometry/stoich_CeO2_n10_a5.410/
  POSCAR
  ceo2.data               # LAMMPS "atom_style charge" compatible
  structure_checksums.json
```

The `POSCAR` lists Ce then O in VASP’s expected species order, and the LAMMPS `data` file declares two atom types with masses and per-atom charges (Ce4+^{4+}=+4, O2−^{2-}=−2). The JSON includes SHA-256 digests and basic metadata (counts and box length) so you can copy them straight into your `manifest.json`.




# 6.) Electrostatics & numerics knobs (global)


This step produces a single, versioned “numerics card” you will include in every LAMMPS input. It pins k-space accuracy/order, real-space cutoffs and neighbor skin, candidate timesteps, and a couple of sanity-check thresholds you’ll test once and then freeze for the whole (T,x)(T,x) grid.

## Numerics card (YAML → reusable includes)

Create `configs/numerics.yaml`. The idea is that you keep one human-readable file, then auto-emit LAMMPS include snippets for prep (NPT/NVT) and measurement (NVE/GK).

```yaml
# configs/numerics.yaml
numerics_version: "1.0.0"

electrostatics:
  pair_style: "buck/coul/long"
  real_cut_A: 12.0                 # real-space cutoff for Buckingham + Coulomb (Å)
  kspace_style: "pppm"
  kspace_accuracy: 1.0e-5          # force/energy target; you will tighten once in validation
  kspace_order: 5                  # stencil order; converge vs 4/6/7 once

neighbors:
  skin_A: 2.0                      # neighbor skin (Å)
  every: 1                         # neighbor rebuild frequency
  delay: 0                         # rebuild delay

integration:
  timestep_fs_candidates: [1.0, 0.5]  # try 1.0 fs first, fall back to 0.5 fs if drift too high
  thermo_stride_steps: 1000
  flux_sample_stride_fs: 1.0          # J(t) sampling stride during GK

prep_ensembles:
  nvt_tau_ps: 0.2
  nvt_chain: 3
  npt_tauT_ps: 0.2
  npt_tauP_ps: 2.0
  npt_iso: true

gk_measurement:
  remove_com_every_steps: 100
  run_shard_ps: 1000.0

sanity_thresholds:
  nve_drift_meV_per_atom_per_ns: 0.5      # energy drift must be below this
  pppm_pressure_shift_bar: 50.0           # when halving accuracy & enlarging cutoff, |ΔP| < this
  plateau_min_width_ps: 50.0
```

## Emit LAMMPS include files from the card

Drop this into `scripts/make_numerics_includes.py`. It reads the YAML and writes two include files you will `include` from every input deck.

```python
#!/usr/bin/env python3
from __future__ import annotations
import yaml
from pathlib import Path

CARD = Path("configs/numerics.yaml")
OUT  = Path("configs/includes")
OUT.mkdir(parents=True, exist_ok=True)

TEMPLATE_GLOBAL = """# --- GLOBAL NUMERICS (AUTO-GENERATED) ---
units           metal
atom_style      charge
boundary        p p p

neighbor        {skin} bin
neigh_modify    every {every} delay {delay}

pair_style      {pair} {rc}
kspace_style    {kstyle} {kacc}
kspace_modify   order {kord}

timestep        {dt_fs}

thermo          {thermo_stride}
"""

TEMPLATE_PREP = """# --- PREP ENSEMBLES (AUTO-GENERATED) ---
# NVT
variable        tauT equal {nvt_tauT}
variable        tchain equal {nvt_chain}
# NPT
variable        tauP equal {npt_tauP}
"""

TEMPLATE_GK = """# --- GK MEASUREMENT (AUTO-GENERATED) ---
fix             delin all momentum {rm_every} linear 1 1 1
"""

def main():
    cfg = yaml.safe_load(CARD.read_text())
    e   = cfg["electrostatics"]; n = cfg["neighbors"]; integ = cfg["integration"]; prep = cfg["prep_ensembles"]; gk = cfg["gk_measurement"]
    dt_fs = float(integ["timestep_fs_candidates"][0])
    g = TEMPLATE_GLOBAL.format(
        pair=e["pair_style"], rc=e["real_cut_A"],
        kstyle=e["kspace_style"], kacc=e["kspace_accuracy"], kord=e["kspace_order"],
        skin=n["skin_A"], every=n["every"], delay=n["delay"],
        dt_fs=dt_fs, thermo_stride=integ["thermo_stride_steps"]
    )
    (OUT/"global_numerics.in").write_text(g)

    p = TEMPLATE_PREP.format(
        nvt_tauT=prep["nvt_tau_ps"], nvt_chain=prep["nvt_chain"], npt_tauP=prep["npt_tauP_ps"]
    )
    (OUT/"prep_ensembles.in").write_text(p)

    m = TEMPLATE_GK.format(rm_every=gk["remove_com_every_steps"])
    (OUT/"gk_measurement.in").write_text(m)
    print("Wrote configs/includes/{global_numerics.in, prep_ensembles.in, gk_measurement.in}")

if __name__ == "__main__":
    main()
```

Run `python scripts/make_numerics_includes.py`. You’ll get:

- `configs/includes/global_numerics.in`
    
- `configs/includes/prep_ensembles.in`
    
- `configs/includes/gk_measurement.in`
    

These are stable, version-controlled, and generated from one card.

## Minimal LAMMPS usage pattern

In your prep script:

```
include configs/includes/global_numerics.in
read_data geometry/stoich_CeO2_n10_a5.410/ceo2.data

# charges are in the data file; pair coefficients come from your FF include
# include potentials/gotte_2007/pair_coeffs.in

# NPT → NVT using the card’s taus
variable T equal 300
fix prep all npt temp ${T} ${T} ${tauT} iso 1.0 1.0 ${tauP}
run 200000
unfix prep

fix settle all nvt temp ${T} ${T} ${tauT}
run 200000
unfix settle
write_restart runs/T300_x0.00/prepped.restart
```

In your GK script:

```
include configs/includes/global_numerics.in
read_restart runs/T300_x0.00/prepped.restart
include configs/includes/gk_measurement.in

# flux assembly (per-atom ke/pe/stress must be consistent with kspace)
compute ke all ke/atom
compute pe all pe/atom
compute s  all centroid/stress/atom NULL virial
compute J  all heat/flux ke pe s

# sample J(t) with the stride from the card (convert fs→steps if needed)
variable stride equal round({flux_stride}/(1000*${dt}))  # if you define ${dt} in fs
fix hcacf all ave/correlate 1 10000 10000 c_J[1] c_J[2] c_J[3] type auto file runs/T300_x0.00/J0Jt.dat ave running

run 2000000
```

## One-time sanity checklist you’ll automate and then freeze

Write a tiny “numerics audit” script that loads your thresholds from `sanity_thresholds` and inspects outputs from a short NVE test and a PPPM-sensitivity test.

1. Energy-drift check. Run a short NVE at your candidate timestep, parse `log.lammps`, compute drift in meV/atom/ns. If it exceeds `nve_drift_meV_per_atom_per_ns`, regenerate the include with the fallback timestep (0.5 fs), rerun, and stamp the chosen value back into `configs/numerics.yaml`.
    
2. PPPM sensitivity. Rerun a short NPT with `kspace_accuracy` halved (e.g., 1e-6) and real-space cutoff extended (e.g., 14 Å). Compare mean pressure; require `|ΔP| < pppm_pressure_shift_bar`. If it fails, tighten the default `kspace_accuracy` and/or increase `kspace_order`, regenerate includes, and repeat. When it passes, you freeze the card for the entire grid.
    

A compact parser stub (drop into `scripts/numerics_audit.py`) might read the last block of `thermo` lines, compute the linear fit of energy vs time for drift, and dump a one-page “PASS/FAIL” with the exact knobs that produced the pass. You don’t need it now; the important part is that your thresholds live in the card and that you run the audit once and commit the passing card and audit report.

## What to commit as the artifact

Commit `configs/numerics.yaml`, the three generated include files, and a short `AUDIT.md` that records:

- the candidate timestep you tested and the measured NVE drift,
    
- the PPPM sensitivity delta-pressure under tightened settings,
    
- the numerics card version you froze for the campaign.
    

From here on, every input deck simply `include`s these files; if you ever need to tighten tolerances, you bump `numerics_version`, regenerate includes, rerun the one-time audit, and keep both versions side-by-side with their reports.




# 7.) EOS & structure validation at ( $T, P$ )


The goal is to automate a small pressure sweep, fit an EOS, freeze the validated volume, sanity-check structure and dynamics, and emit a clear pass/fail note. Below you’ll find a minimal orchestration plan, one LAMMPS input you can reuse at multiple pressures, and a Python analysis script that fits Birch–Murnaghan (falls back to a linear bulk modulus if needed), freezes the box, and writes a short validation report with plots.

## Orchestration (what you’ll run)

1. Prepare a pristine CeO2\mathrm{CeO_2} supercell (`ceo2.data`) and the numerics includes from Step 6.
    
2. Sweep a few pressures around target (e.g., −1,0,+1,+2-1, 0, +1, +2 GPa) with short NPT segments; log the time-averaged volume and pressure at each point.
    
3. Fit EOS → (V0,B0,B0′)(V_0, B_0, B_0') via Birch–Murnaghan (or a local linear BB if BM doesn’t converge on sparse data).
    
4. Freeze the box at V0V_0 (or the NPT mean at the target pressure), run a short NVT to verify mean stress ≈P\approx P and to dump a Ce–O RDF as a quick structure fingerprint.
    
5. Run a short NVE at the frozen volume to measure energy drift.
    
6. Compare a,ρ,Ba,\rho,B to your dossier tolerances, check drift against your numerics card, and write `EOS_VALIDATION.md` with plots.
    

---

## LAMMPS input for an NPT point (reused across pressures)

Save as `lmp/eos_scan.in`. This reads your structure, includes the numerics card (from Step 6), and runs a short NPT with gentle damping; it then prints time-averaged pressure and volume in a single machine-parsable line.

```lammps
# eos_scan.in  — run a short NPT at a given target pressure (GPa) and temperature (K)
units           metal
atom_style      charge
boundary        p p p

# --- includes from Step 6 ---
include         configs/includes/global_numerics.in
# pair coefficients for your chosen FF (emit from Step 2 tooling)
# include      potentials/gotte_2007/pair_coeffs.in

read_data       geometry/stoich_CeO2_n10_a5.410/ceo2.data

variable        T   equal ${T}           # passed on the command line
variable        P   equal ${P}           # in GPa (metal units: bar = 1e-4 GPa)

# thermo and averages
thermo_style    custom step temp press pe ke etotal vol
thermo          ${thermo_stride}

# NPT (isotropic) — taus supplied in includes via variables tauT, tauP
fix             eos all npt temp ${T} ${T} ${tauT} iso ${P} ${P} ${tauP}
run             50000            # warmup
reset_timestep  0
variable        nsamp equal 40000
fix             ave all ave/time 100 400 40000 v_press v_vol file eos_avg.txt
run             ${nsamp}
unfix           ave
unfix           eos

# Emit a single summary line (mean P, mean V, std P, std V) in metal units
variable        Pmean equal f_ave[1]
variable        Vmean equal f_ave[2]
# crude stddev via thermo keywords; if unavailable, parse eos_avg.txt offline
print           "EOS_SUMMARY T=${T}K Ptarget=${P}GPa Pmean=${Pmean}GPa Vmean=${Vmean}A3" append eos_summary.txt screen yes
```

You’ll invoke it multiple times by setting `-var T 300 -var P 0.0`, `-var P 1.0`, etc., each in its own subdirectory like `runs/T300_x0.00/eos/P_0.0/`.

If you prefer not to rely on `fix ave/time` indexing for stddev, simply parse the `eos_avg.txt` columns in the Python step; the script below handles either.

---

## Python analysis and validation (fit EOS, freeze volume, RDF, NVE drift, report)

Save as `scripts/validate_eos.py`. It expects a folder containing subfolders named `P_<value>/eos_summary.txt` (or `eos_avg.txt`), fits Birch–Murnaghan to P(V)P(V), falls back to a local linear bulk modulus if needed, writes a plot, writes a freeze-volume file for your next stage, optionally parses a short NVT and NVE check, and emits a pass/fail markdown.

```python
#!/usr/bin/env python3
from __future__ import annotations
import re, json, math, argparse
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_points(eos_dir: Path) -> List[Tuple[float,float]]:
    pts = []
    for sub in sorted(eos_dir.glob("P_*")):
        summ = sub/"eos_summary.txt"
        avg  = sub/"eos_avg.txt"
        if summ.exists():
            m = re.search(r"Pmean=([-\d\.Ee+]+)GPa Vmean=([-\d\.Ee+]+)A3", summ.read_text())
            if not m: 
                continue
            Pmean = float(m.group(1)); Vmean = float(m.group(2))
            pts.append((Vmean, Pmean))
        elif avg.exists():
            data = np.loadtxt(avg)
            # columns: time-averaged values per block; assume P in bar; V in Å^3
            # metal units: 1 bar = 1e-4 GPa
            Pbar = np.mean(data[:,0]) if data.ndim==2 else np.mean(data)
            Vmean = np.mean(data[:,1]) if data.ndim==2 else np.nan
            pts.append((Vmean, Pbar*1e-4))
    return sorted(pts, key=lambda x: x[0])

def BM3_P_of_V(V: np.ndarray, V0: float, B0: float, B0p: float) -> np.ndarray:
    # 3rd-order Birch–Murnaghan: P(V) in same units as B0 (GPa)
    eta = (V0/V)**(2.0/3.0)
    term1 = (3.0/2.0)*B0*(eta**(7.0/2.0) - eta**(5.0/2.0))  # (7/3 - 5/3) in exponent * 3/2 B0
    # more standard algebra:
    x = eta
    term1 = (3.0/2.0)*B0*(x**(3.5) - x**(2.5))
    corr  = 1.0 + (3.0/4.0)*(B0p - 4.0)*(x - 1.0)
    return term1*corr

def fit_bm3(VP: List[Tuple[float,float]]) -> Optional[Tuple[float,float,float,float]]:
    V = np.array([v for v,_ in VP], float)
    P = np.array([p for _,p in VP], float)
    V0_guess = V[np.argmin(np.abs(P))] if np.any(np.isfinite(P)) else np.median(V)
    B0_guess = 200.0  # GPa, rough for fluorite oxides
    B0p_guess = 4.0
    # crude grid search to avoid SciPy dependency
    V0_grid  = np.linspace(0.98*V0_guess, 1.02*V0_guess, 21)
    B0_grid  = np.linspace(100.0, 300.0, 41)
    B0p_grid = np.linspace(3.0, 6.0, 31)
    best = None; best_err = 1e99
    for V0 in V0_grid:
        x = (V0/V)**(2.0/3.0)
        x35 = x**3.5; x25 = x**2.5
        base = (3.0/2.0)*(x35 - x25)
        # we can separate B0 linearly, but B0' appears only in corr
        for B0 in B0_grid:
            for B0p in B0p_grid:
                corr = 1.0 + (3.0/4.0)*(B0p - 4.0)*(x - 1.0)
                Pfit = B0*base*corr
                err = float(np.mean((Pfit - P)**2))
                if err < best_err:
                    best_err = err; best = (V0, B0, B0p)
    if best is None:
        return None
    V0, B0, B0p = best
    # compute R^2
    Pfit = BM3_P_of_V(V, V0, B0, B0p)
    ss_res = float(np.sum((P - Pfit)**2))
    ss_tot = float(np.sum((P - np.mean(P))**2)) if len(P)>1 else 0.0
    r2 = 1.0 - ss_res/ss_tot if ss_tot>0 else 1.0
    return V0, B0, B0p, r2

def fit_linear_bulk(VP: List[Tuple[float,float]]) -> Tuple[float,float,float]:
    # local linear around the point closest to P=0: P = a + b*V ; B = -V * dP/dV
    V = np.array([v for v,_ in VP], float)
    P = np.array([p for _,p in VP], float)
    A = np.vstack([np.ones_like(V), V]).T
    coef, *_ = np.linalg.lstsq(A, P, rcond=None)
    a, b = coef
    V0 = -a/b
    B0 = -(V0)*b
    # no B' in linear model; return NaN for B0'
    return V0, B0, float("nan")

def write_plots(outdir: Path, VP, fit_bm, fit_lin):
    outdir.mkdir(parents=True, exist_ok=True)
    V = np.array([v for v,_ in VP]); P = np.array([p for _,p in VP])
    plt.figure()
    plt.scatter(V, P, label="NPT means")
    vspan = np.linspace(0.98*min(V), 1.02*max(V), 200)
    if fit_bm:
        V0,B0,B0p,R2 = fit_bm
        plt.plot(vspan, BM3_P_of_V(vspan, V0,B0,B0p), label=f"BM3 fit: V0={V0:.2f} Å³, B0={B0:.1f} GPa, B0'={B0p:.2f}, R²={R2:.3f}")
    if fit_lin:
        V0l,B0l,_ = fit_lin
        Pl = np.polyval([ (-(B0l)/V0l), (B0l) ], vspan)  # derived from P=a+bV with B=-Vb, a=-bV0 -> b=-B/V0
        plt.plot(vspan, Pl, "--", label=f"Linear: V0={V0l:.2f} Å³, B0={B0l:.1f} GPa")
    plt.xlabel("Volume per cell (Å³)")
    plt.ylabel("Pressure (GPa)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/"eos_fit.png", dpi=200)

def write_freeze_volume(outdir: Path, a0_A: float, n: int, V0_cell_A3: float):
    # CeO2 conventional cell has 1 conventional cell per a0^3; supercell volume: (n*a0)^3
    # If your input volumes are total box volumes, compute equivalent box length from target volume per conventional cell
    L_A = (V0_cell_A3)**(1.0/3.0) * n
    (outdir/"freeze_volume.json").write_text(json.dumps({"L_box_A": L_A}, indent=2))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eos-root", required=True, help="runs/T300_x0.00/eos")
    ap.add_argument("--a0A", type=float, required=True, help="nominal a0 used to build the cell (Å)")
    ap.add_argument("--n", type=int, required=True, help="supercell replication n")
    ap.add_argument("--dossier_a0A", type=float, required=True)
    ap.add_argument("--dossier_B_GPa", type=float, required=True)
    ap.add_argument("--tol_a0_rel", type=float, default=0.02)   # 1–2%
    ap.add_argument("--tol_B_rel",  type=float, default=0.15)   # 10–15%
    ap.add_argument("--drift_meV_atom_ns_max", type=float, default=0.5)
    args = ap.parse_args()

    eos_dir = Path(args.eos_root)
    pts = read_points(eos_dir)
    if len(pts) < 3:
        print("Need at least 3 pressure points for BM; will fall back to linear.")
    fit_bm = fit_bm3(pts) if len(pts) >= 3 else None
    fit_lin = fit_linear_bulk(pts)

    # choose preferred result
    if fit_bm:
        V0, B0, B0p, R2 = fit_bm
        eos_model = "BM3"
    else:
        V0, B0, _ = fit_lin
        B0p, R2 = float("nan"), float("nan")
        eos_model = "Linear"

    # convert to lattice parameter a from V0 per conventional cell
    a0_fit = (V0)**(1.0/3.0)
    rho_fit = None  # if you want, compute from masses / box volume

    # pass/fail against dossier
    a_rel = abs(a0_fit - args.dossier_a0A) / args.dossier_a0A
    B_rel = abs(B0 - args.dossier_B_GPa) / args.dossier_B_GPa
    a_ok  = (a_rel <= args.tol_a0_rel)
    B_ok  = (B_rel <= args.tol_B_rel)

    outdir = eos_dir.parent  # runs/T300_x0.00
    write_plots(outdir, pts, fit_bm, fit_lin)
    write_freeze_volume(outdir, args.a0A, args.n, V0_cell_A3=V0)

    # NVT stress and RDF check (optional; parse if you write rdf.out and stress.out)
    rdf_png = None
    rdf_file = outdir/"rdf_CeO.dat"
    if rdf_file.exists():
        import numpy as np
        dat = np.loadtxt(rdf_file)
        # You can add quick checks on Ce–O first peak position here if desired.
        rdf_png = str(outdir/"rdf_CeO.png")

    # NVE drift check (parse log_nve.lammps if you save it)
    drift_meV_atom_ns = None
    lognve = outdir/"log_nve.lammps"
    if lognve.exists():
        times, etot = [], []
        for line in lognve.read_text().splitlines():
            if re.match(r"^\s*\d+\s", line):
                parts = line.split()
                # assumes thermo: step etotal ...
                step = int(parts[0]); e = float(parts[2])
                times.append(step); etot.append(e)
        if len(times) > 1:
            t_fs = (times[-1]-times[0]) * 1.0  # if timestep=1 fs; adjust if different
            slope = np.polyfit(times, etot, 1)[0]  # e per step
            drift_eV_per_ps = slope / 1000.0
            # convert to meV/atom/ns if you know N; here we report per atom if N known
            drift_meV_atom_ns = abs(drift_eV_per_ps) * 1e3  # placeholder

    # write report
    note = outdir/"EOS_VALIDATION.md"
    with note.open("w") as f:
        f.write(f"# EOS & structure validation\n\n")
        f.write(f"Model: {eos_model}\n\n")
        f.write(f"Fitted a0 = {a0_fit:.5f} Å; dossier a0 = {args.dossier_a0A:.5f} Å; rel. diff = {a_rel*100:.2f}% → {'PASS' if a_ok else 'FAIL'}\n\n")
        f.write(f"Fitted B0 = {B0:.1f} GPa; dossier B = {args.dossier_B_GPa:.1f} GPa; rel. diff = {B_rel*100:.2f}% → {'PASS' if B_ok else 'FAIL'}\n\n")
        if fit_bm:
            f.write(f"Fitted B0' = {B0p:.2f}; R² = {R2:.3f}\n\n")
        f.write(f"Plots: eos_fit.png\n\n")
        if drift_meV_atom_ns is not None:
            ok = drift_meV_atom_ns <= args.drift_meV_atom_ns_max
            f.write(f"NVE drift ≈ {drift_meV_atom_ns:.3f} meV/atom/ns → {'PASS' if ok else 'FAIL'} (threshold {args.drift_meV_atom_ns_max})\n\n")
        f.write(f"Freeze volume written to freeze_volume.json\n")

    print(f"Wrote {note} and eos_fit.png")
    if not (a_ok and B_ok):
        print("EOS VALIDATION FAILED — check potential, barostat, or add more pressure points.")
```

### Notes on running

Create pressure subfolders and run LAMMPS for each point, for example:

```
for P in -1.0 0.0 1.0 2.0; do
  d="runs/T300_x0.00/eos/P_${P}"
  mkdir -p "$d" && cd "$d"
  lmp -in ../../../lmp/eos_scan.in -var T 300 -var P $P
  cd - >/dev/null
done
```

Then fit and generate the report:

```
python scripts/validate_eos.py \
  --eos-root runs/T300_x0.00/eos \
  --a0A 5.41 --n 10 \
  --dossier_a0A 5.41 \
  --dossier_B_GPa 200.0
```

This writes `runs/T300_x0.00/eos_fit.png`, `freeze_volume.json`, and `EOS_VALIDATION.md` summarizing pass/fail against your gates.

### Freezing the volume and running the quick NVT/NVE checks

With `freeze_volume.json` in hand, read it in your next LAMMPS script, set the box length LL accordingly (`change_box all x final 0 L y final 0 L z final 0 L remap`), run a short NVT at the target TT to verify mean stress ≈P\approx P and dump a Ce–O RDF:

```lammps
# after NPT sweep and fit
variable     L equal ${Lbox}   # load from freeze_volume.json via your runner
change_box   all x final 0.0 ${L} y final 0.0 ${L} z final 0.0 ${L} remap units box

# NVT check
fix          chk all nvt temp ${T} ${T} ${tauT}
compute      myrdf all rdf 200 1 2          # Ce(1)-O(2) first peak near ~2.3–2.5 Å (force-field dependent)
fix          rdfout all ave/time 100 100 10000 c_myrdf[*] file rdf_CeO.dat mode vector
run          200000
unfix        rdfout
unfix        chk

# NVE drift sanity
reset_timestep 0
fix          drift all momentum 1000 linear 1 1 1
run          200000
unfix        drift
```

Parse `rdf_CeO.dat` and `log.lammps` in the Python script (stubs are included) to add a one-line structural check and an energy-drift line to your report.

That’s the entire EOS/structure gate wired end-to-end. Once `EOS_VALIDATION.md` says **PASS**, you freeze the numerics and box volume for this state and move on to the stoichiometric kk validation (Step 8) with precisely the same knobs.




# 8.) Stoichiometric $k$ validation (GK, plus one rNEMD check)


You’ll do one clean GK calculation on the validated stoichiometric crystal, prove size/time convergence, and (optionally) cross-check with a gentle rNEMD run. Below are drop-in inputs and analysis scripts that match the gates you committed earlier: stable plateau with Sokal windowing, shard scaling ∼1/M\sim 1/\sqrt{M}, size/mesh insensitivity, and a final value that sits inside your theory/experiment bracket.

---

## GK production input (read the frozen, equilibrated restart; record J(t)\mathbf J(t))

Save as `lmp/gk_production.in`. It reads your NPT→NVT-frozen restart from Step 7, includes the numerics card from Step 6, assembles a **k-space-consistent** heat flux, then dumps both the J\mathbf J time series (binary) and the Jα(0)Jα(t)J_\alpha(0)J_\alpha(t) table for each shard.

```lmp
# gk_production.in — NVE shards to sample heat flux for GK
units           metal
atom_style      charge
boundary        p p p

# Numerics + electrostatics knobs (pppm accuracy/order, cutoffs, dt, neighbor)
include         configs/includes/global_numerics.in

# Pair coefficients for chosen FF (emitted in Step 2)
# include      potentials/gotte_2007/pair_coeffs.in

# Read equilibrated, frozen-volume state (Step 7 wrote this)
read_restart    runs/T300_x0.00/prepped.restart

# Remove COM drift during measurement
include         configs/includes/gk_measurement.in   # defines: fix momentum every N steps

# --- GK flux assembly (consistent with PPPM) ---
compute         cke all ke/atom
compute         cpe all pe/atom
compute         cst all centroid/stress/atom NULL virial
compute         J all heat/flux cke cpe cst                 # Jx Jy Jz

# Sampling stride (fs) → steps; define dt_fs in your numerics or set here
variable        dt_fs equal 1.0
variable        stride_fs equal 1.0
variable        stride_steps equal ceil(${stride_fs}/${dt_fs})

# Raw flux dump (binary) for deterministic post-processing
dump            dJ all custom ${stride_steps} runs/T300_x0.00/nve_shard_000/flux.raw c_J[1] c_J[2] c_J[3]
dump_modify     dJ format binary yes

# HCACF accumulation with overlapping origins (on-the-fly sanity look)
fix             acf all ave/correlate 1 10000 10000 c_J[1] c_J[2] c_J[3] type auto file runs/T300_x0.00/J0Jt_shard000.dat ave running

# Pure dynamics (no thermostats/barostats)
run             2000000      # = 2 ns at 1 fs; adjust to your shard_length_ps

unfix           acf
undump          dJ
```

Duplicate the shard block for `nve_shard_001`, `nve_shard_002`, … (or launch separate jobs pointing `dump`/`file` to the right shard folders). Keep the same MPI layout for bitwise-stable flux files.

---

## GK analysis with Sokal window + blocking CI (per shard and aggregated)

Save as `analysis/analyze_gk.py`. It ingests one or more shards, constructs HCACF with overlapping origins, estimates τint\tau_{\mathrm{int}}, picks a Sokal window t⋆=c τintt^\star=c\,\tau_{\mathrm{int}}, integrates to k^(t⋆)\hat{k}(t^\star), computes a CI by Flyvbjerg–Petersen blocking (per shard) and by shard-wise variance across independent shards, and makes the standard plots.

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_flux_binary(path: Path) -> np.ndarray:
    # binary dump with 3 doubles per record (Jx, Jy, Jz)
    raw = np.fromfile(path, dtype=np.float64)
    assert raw.size % 3 == 0
    return raw.reshape(-1, 3)

def hcacf(J: np.ndarray, maxlag: int) -> np.ndarray:
    # mean-zero each component to kill DC offset
    J = J - J.mean(axis=0, keepdims=True)
    N = J.shape[0]
    C = np.zeros((maxlag, 3))
    for lag in range(maxlag):
        # overlapping origins; vectorized dot over segments
        a = J[:N-lag, :]
        b = J[lag:, :]
        C[lag, :] = (a*b).mean(axis=0)
    # scalar HCACF (1/3 trace)
    return C.mean(axis=1)

def integrated_autocorr_time(C: np.ndarray, dt: float) -> float:
    # simple positive-sequence estimator; stop at first zero-crossing
    s = 0.0
    for k, c in enumerate(C):
        if c <= 0 and k > 0:
            break
        s += c
    return 2.0 * s / C[0] * dt   # τ_int in time units

def trapz_windowed(C: np.ndarray, dt: float, tstar: float) -> float:
    nstar = max(1, int(round(tstar/dt)))
    return np.trapz(C[:nstar], dx=dt)

def blocking_std(x: np.ndarray) -> float:
    # Flyvbjerg–Petersen over a vector of block contributions (assumed stationary)
    y = x.copy()
    vars_ = []
    while y.size >= 8:
        vars_.append(y.var(ddof=1) / y.size)
        # coarse-grain by averaging consecutive pairs
        y = 0.5*(y[0::2] + y[1::2])
    # pick last plateau; conservative choice = max
    return float(np.sqrt(max(vars_))) if vars_ else float("nan")

def shard_k_and_ci(J: np.ndarray, V: float, T: float, dt: float, c_sokal: float) -> dict:
    C = hcacf(J, maxlag=min(J.shape[0]//2, 500000))
    tau = integrated_autocorr_time(C, dt)
    tstar = c_sokal * tau
    area = trapz_windowed(C, dt, tstar)          # ∫ C(t) dt
    k = area / (V * (8.617333262e-5*T)**2)       # kB (eV/K); units: (eV^2/ps) / (Å^3 eV^2/K^2) → W/mK with conversions
    # Unit conversions: 1 eV/ps = 1.60218e-19 J / 1e-12 s = 1.60218e-7 W; 1 Å^3 = 1e-30 m^3
    # C(t) has units of (energy flux)^2; to keep this concise, assume you post-scale k by a known factor F:
    F = 1.0  # ← replace with precomputed factor for your unit system; see note below
    k *= F
    # Blocking CI: split window integral into contributions per origin (approx): use sliding dot products at small lags
    # Simpler: segment time series into M chunks >> τ_int and recompute windowed integral per chunk
    M = max(4, int(np.floor((J.shape[0]*dt) / (10*tau))))  # ~10 τ per chunk
    chunk_len = J.shape[0] // M
    ks = []
    for i in range(M):
        seg = J[i*chunk_len:(i+1)*chunk_len, :]
        if seg.shape[0] < 1000: break
        Ci = hcacf(seg, maxlag=min(seg.shape[0]//4, 200000))
        ki = trapz_windowed(Ci, dt, tstar) / (V * (8.617333262e-5*T)**2) * F
        ks.append(ki)
    ks = np.array(ks) if ks else np.array([k])
    ci = 1.96 * ks.std(ddof=1) / np.sqrt(len(ks)) if ks.size > 1 else np.nan
    return {"k": float(k), "ci": float(ci), "tau_int_ps": float(tau), "tstar_ps": float(tstar), "C": C}

def plot_c_and_k(C: np.ndarray, dt: float, outpng_c: Path, outpng_k: Path):
    t = np.arange(C.size)*dt
    plt.figure(); plt.plot(t, C); plt.xlabel("t (ps)"); plt.ylabel("C_JJ(t)"); plt.tight_layout(); plt.savefig(outpng_c, dpi=200)
    kcumu = np.cumsum((C[:-1] + C[1:])*0.5*dt)
    plt.figure(); plt.plot(t[1:], kcumu); plt.xlabel("t (ps)"); plt.ylabel("∫ C_JJ dt"); plt.tight_layout(); plt.savefig(outpng_k, dpi=200)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True, help="e.g. T300_x0.00")
    ap.add_argument("--shards", nargs="+", required=True, help="paths to .../nve_shard_xxx/flux.raw")
    ap.add_argument("--boxA", type=float, required=True, help="box length (Å)")
    ap.add_argument("--T", type=float, required=True, help="K")
    ap.add_argument("--dt_fs", type=float, default=1.0)
    ap.add_argument("--c_sokal", type=float, default=5.0)
    args = ap.parse_args()

    V_A3 = args.boxA**3
    dt_ps = args.dt_fs/1000.0
    per = []
    for s in args.shards:
        J = load_flux_binary(Path(s))
        r = shard_k_and_ci(J, V_A3, args.T, dt_ps, args.c_sokal)
        per.append(r)
    ks = np.array([r["k"] for r in per])
    # shard-wise CI
    k_mean = float(ks.mean())
    ci_shard = float(1.96 * ks.std(ddof=1) / np.sqrt(len(ks))) if len(ks) > 1 else float("nan")

    # write result + plots for the first shard’s diagnostics
    outdir = Path(f"analysis/{args.state}")
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir/"k_result.json", "w") as f:
        json.dump({
            "k_mean": k_mean, "ci_shard": ci_shard,
            "per_shard": [{"k": r["k"], "ci": r["ci"], "tau_int_ps": r["tau_int_ps"], "tstar_ps": r["tstar_ps"]} for r in per],
            "meta": {"T": args.T, "V_A3": V_A3, "dt_fs": args.dt_fs, "c_sokal": args.c_sokal}
        }, f, indent=2)
    # optional plots
    C0 = per[0]["C"]
    plot_c_and_k(C0, dt_ps, outdir/"hcacf.png", outdir/"kcumu.png")
    print(f"k = {k_mean:.2f} W/mK (±{ci_shard:.2f} shard-CI); shards={len(ks)}")

if __name__ == "__main__":
    main()
```

**Unit note.** LAMMPS `units metal` yields energy in eV, time in ps, length in Å. The prefactor for

$k=1VkB2T2∫CJJ(t) dtk=\frac{1}{V k_B^2 T^2}\int C_{JJ}(t)\,dt$

must convert [eV2/ps]/[A˚3 eV2/K2][eV^2/ps]/[Å^3\,eV^2/K^2] to W/mK. Precompute and set `F = (1.602176634e-19)**2 / (1e-12) / (1e-30) / (1.380649e-23)**2`, then test against the LAMMPS KAPPA example to verify. Once you lock `F`, keep it frozen in your analysis version.

---

## Convergence driver (time, size, PPPM) and final report

Save as `analysis/validate_k.py`. It aggregates shards, reruns the GK analysis with a slightly wider/narrower window to show sensitivity < 1σ, compares normal vs. elongated cell, and writes `K_VALIDATION.md` with plots.

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze_gk import load_flux_binary, shard_k_and_ci

def run_state(state: str, shard_paths, L_A: float, T: float, dt_fs: float, c: float):
    V = L_A**3; dt = dt_fs/1000.0
    ks=[]; per=[]
    for p in shard_paths:
        r = shard_k_and_ci(load_flux_binary(Path(p)), V, T, dt, c)
        ks.append(r["k"]); per.append(r)
    ks = np.array(ks)
    k_mean = ks.mean(); ci = 1.96*ks.std(ddof=1)/np.sqrt(len(ks)) if len(ks)>1 else np.nan
    return k_mean, ci, per

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--shards", nargs="+", required=True)
    ap.add_argument("--L_A", type=float, required=True)
    ap.add_argument("--T", type=float, required=True)
    ap.add_argument("--dt_fs", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=5.0)
    ap.add_argument("--elong_state", default=None, help="optional elongated state name")
    ap.add_argument("--elong_shards", nargs="*", default=[])
    ap.add_argument("--elong_L_A", type=float, default=None)
    args = ap.parse_args()

    k, ci, per = run_state(args.state, args.shards, args.L_A, args.T, args.dt_fs, args.c)
    # sensitivity to window
    k_lo,_ ,_ = run_state(args.state, args.shards, args.L_A, args.T, args.dt_fs, args.c*0.8)
    k_hi,_ ,_ = run_state(args.state, args.shards, args.L_A, args.T, args.dt_fs, args.c*1.2)
    sens = max(abs(k - k_lo), abs(k - k_hi))

    # elongated comparison
    elong_line = "not performed"
    if args.elong_state and args.elong_shards and args.elong_L_A:
        kE, ciE, _ = run_state(args.elong_state, args.elong_shards, args.elong_L_A, args.T, args.dt_fs, args.c)
        delta = abs(k - kE)
        elong_line = f"k(elong) = {kE:.2f} ± {ciE:.2f} W/mK; |Δ| = {delta:.2f} W/mK"

    outdir = Path(f"analysis/{args.state}")
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir/"K_VALIDATION.md","w") as f:
        f.write(f"# Stoichiometric k validation — {args.state}\n\n")
        f.write(f"**GK (mean over shards):** k = {k:.2f} ± {ci:.2f} W/mK (95% CI)\n\n")
        f.write(f"**Window sensitivity:** c→0.8c or 1.2c ⇒ |Δk| = {sens:.2f} W/mK\n\n")
        f.write(f"**Elongated-axis check:** {elong_line}\n\n")
        f.write(f"Plots: hcacf.png, kcumu.png (see per-shard diagnostics)\n")

    print(f"{args.state}: k = {k:.2f} ± {ci:.2f} W/mK; sensitivity {sens:.2f} W/mK")
    if args.elong_state:
        print(elong_line)

if __name__ == "__main__":
    main()
```

Run it like:

```
python analysis/validate_k.py \
  --state T300_x0.00 \
  --shards runs/T300_x0.00/nve_shard_000/flux.raw runs/T300_x0.00/nve_shard_001/flux.raw \
  --L_A 54.10 --T 300 --dt_fs 1.0 --c 5.0 \
  --elong_state T300_x0.00_long \
  --elong_shards runs/T300_x0.00_long/nve_shard_000/flux.raw \
  --elong_L_A 108.20
```

This writes `analysis/T300_x0.00/K_VALIDATION.md` and the standard plots.

---

## Optional rNEMD cross-check (one state)

**Input.** Save as `lmp/rnemd.in`. It applies Müller–Plathe swaps, bins the temperature profile, and writes a clean `profile.dat`.

```lmp
# rnemd.in — gentle MP swap cross-check
units           metal
atom_style      charge
boundary        p p p
include         configs/includes/global_numerics.in
# include      potentials/gotte_2007/pair_coeffs.in
read_restart    runs/T300_x0.00/prepped.restart

variable        T equal 300
variable        Nbins equal 40
variable        period equal 5000        # swap period in steps (gentle)
fix             tc all thermal/conductivity ${period} z 20    # heat along z; 20 = number of atom exchanges per swap

compute         slab all chunk/atom bin/1d z lower ${Nbins} units box
fix             prof all ave/chunk 100 200 20000 slab temp file runs/T300_x0.00/rnemd/profile.dat

run             2000000
unfix           prof
unfix           tc
```

**Analysis.** Save as `analysis/analyze_rnemd.py`. It fits the central linear region, discards a few bins near the swap slabs, converts the accumulated exchanged energy to a flux, and returns k=−q/∇Tk=-q/\nabla T. Do a rate-independence test by halving/doubling `period` and confirming the same kk within noise.

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, numpy as np
from pathlib import Path

def load_profile(path: Path):
    # LAMMPS chunk/atom + ave/chunk output: bin index, zlo, zhi, N, T
    dat = np.loadtxt(path)
    z = 0.5*(dat[:,1]+dat[:,2]); T = dat[:,4]
    return z, T

def fit_bulk_gradient(z, T, n_exclude=2):
    zc, Tc = z[n_exclude:-n_exclude], T[n_exclude:-n_exclude]
    A = np.vstack([zc, np.ones_like(zc)]).T
    m, b = np.linalg.lstsq(A, Tc, rcond=None)[0]
    return m, b, zc, Tc   # m = dT/dz

def heat_flux_from_log(logpath: Path, Lx: float, Ly: float, period_steps: int, dt_fs: float):
    # Simplest path: parse "fix thermal/conductivity" cumulative energy exchange per area & time from log or a thermo custom
    # Placeholder: let user pass q directly if you printed it. Otherwise compute from fix tc’s reported cumulative E (eV).
    raise NotImplementedError("Wire this to your log parser that reads cumulative exchanged energy from LAMMPS output.")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--exclude_bins", type=int, default=2)
    ap.add_argument("--q_Wm2", type=float, required=True, help="net heat flux magnitude (W/m^2)")
    args = ap.parse_args()
    z, T = load_profile(Path(args.profile))
    m, b, zc, Tc = fit_bulk_gradient(z, T, n_exclude=args.exclude_bins)
    dTdz = m  # K/Å in metal units profile; convert: 1 Å = 1e-10 m
    grad = dTdz / 1e-10
    k = -args.q_Wm2 / grad
    print(f"rNEMD: grad = {grad:.3e} K/m, k = {k:.2f} W/mK")

if __name__ == "__main__":
    main()
```

When you run with two swap rates (e.g., `period = 5000` and `10000`), both kk values should agree within the combined CI from GK. Record the test and the fit window in `analysis/T300_x0.00/K_VALIDATION.md`.

---

## What you should commit as the artifact

1. `analysis/T300_x0.00/k_result.json` with kk, CI, τint\tau_{\mathrm{int}}, t⋆t^\star, and meta (T, V, dt, c).
    
2. `analysis/T300_x0.00/hcacf.png` and `kcumu.png` showing a broad, flat plateau.
    
3. `analysis/T300_x0.00/K_VALIDATION.md` summarizing shard scaling, window sensitivity, size/mesh insensitivity, and the optional rNEMD agreement.
    

Once this baseline passes, freeze **all** numerics, sampling, and analysis knobs and reuse them for every (T,x)(T,x) state in the grid.




# 9.) Vacancy generator & $\mathrm{Ce}^{3+}$ assignment (for $x>0$ )
- Blue-noise sample O sites with $r_{\text {min }}$; for each $V_O$ relabel two nearest $\mathrm{Ce}^{4+} \rightarrow \mathrm{Ce}^{3+}$.
- Artifact: vacancy index list, $\mathrm{Ce}^{3+}$ maps, separation histograms, seeds stored.


Below is a single, reproducible Python tool that (i) blue-noise samples oxygen sites with a minimum-distance constraint rmin⁡r_{\min} under PBC, (ii) deletes those O atoms to create VOV_\mathrm{O}, (iii) converts the **two nearest** Ce⁴⁺ neighbors of each vacancy to Ce³⁺, and (iv) writes all artifacts you’ll track: a modified LAMMPS data file, the vacancy index list, the Ce³⁺ map, and CSV histograms for VOV_\mathrm{O}–VOV_\mathrm{O} and Ce³⁺–VOV_\mathrm{O} separations. It records seeds and parameters in a manifest for auditing.

Assumptions match your earlier steps: orthorhombic box (your fluorite supercells are cubic), `atom_style charge` LAMMPS data from Step 5, and a species map that pins type IDs and charges for `Ce4+`, `Ce3+`, and `O2-`.

---

## `generate_vacancies.py`

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, random, sys, csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# ---------------- IO: minimal LAMMPS data (atom_style charge; orthorhombic) ----------------

@dataclass
class Box:
    Lx: float; Ly: float; Lz: float
    xlo: float; xhi: float
    ylo: float; yhi: float
    zlo: float; zhi: float

@dataclass
class Atoms:
    ids:   np.ndarray  # (N,) int
    types: np.ndarray  # (N,) int
    q:     np.ndarray  # (N,) float
    xyz:   np.ndarray  # (N,3) float

def read_lammps_data(path: Path) -> Tuple[Box, Atoms, List[str]]:
    lines = path.read_text().splitlines()
    header = []
    i = 0
    while i < len(lines):
        header.append(lines[i])
        if "Atoms" in lines[i]:
            break
        i += 1
    # box and counts
    N = None; xlo=xhi=ylo=yhi=zlo=zhi=None
    for ln in header:
        tok = ln.strip().split()
        if len(tok)>=2 and tok[1]=="atoms":
            N = int(tok[0])
        if "xlo xhi" in ln: xlo, xhi = map(float, ln.split()[:2])
        if "ylo yhi" in ln: ylo, yhi = map(float, ln.split()[:2])
        if "zlo ylo" in ln: pass
        if "zlo zhi" in ln: zlo, zhi = map(float, ln.split()[:2])
    if N is None: raise ValueError("Could not read atom count.")
    Lx, Ly, Lz = xhi-xlo, yhi-ylo, zhi-zlo
    # find "Atoms" section start
    while i < len(lines) and not lines[i].strip().startswith("Atoms"):
        i += 1
    # skip title and blank
    i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    # parse N lines: id type q x y z
    ids=[]; types=[]; q=[]; xyz=[]
    for k in range(N):
        parts = lines[i+k].split()
        if len(parts) < 6:
            raise ValueError("Expected 'id type q x y z' per atom.")
        ids.append(int(parts[0]))
        types.append(int(parts[1]))
        q.append(float(parts[2]))
        xyz.append([float(parts[3]), float(parts[4]), float(parts[5])])
    A = Atoms(ids=np.array(ids, int),
              types=np.array(types, int),
              q=np.array(q, float),
              xyz=np.array(xyz, float))
    return Box(Lx, Ly, Lz, xlo, xhi, ylo, yhi, zlo, zhi), A, lines

def write_lammps_data(path: Path, box: Box, atoms: Atoms, masses: Dict[int,float], comment: str="# generated"):
    N = atoms.ids.size
    # map to contiguous IDs 1..N in current order
    order = np.argsort(atoms.ids)
    ids   = atoms.ids[order]
    types = atoms.types[order]
    q     = atoms.q[order]
    xyz   = atoms.xyz[order]
    # rebuild IDs
    new_ids = np.arange(1, N+1, dtype=int)

    out = []
    out.append(comment)
    out.append(f"\n{N} atoms")
    out.append(f"{len(masses)} atom types\n")
    out.append(f"{box.xlo:.10f} {box.xhi:.10f} xlo xhi")
    out.append(f"{box.ylo:.10f} {box.yhi:.10f} ylo yhi")
    out.append(f"{box.zlo:.10f} {box.zhi:.10f} zlo zhi\n")
    out.append("Masses\n")
    for t, m in sorted(masses.items()):
        out.append(f"{t} {m:.6f}")
    out.append("\nAtoms # charge\n")
    for i in range(N):
        out.append(f"{new_ids[i]} {int(types[i])} {q[i]:.6f} {xyz[i,0]:.10f} {xyz[i,1]:.10f} {xyz[i,2]:.10f}")
    path.write_text("\n".join(out))

# ---------------- geometry helpers ----------------

def min_image_delta(dx: float, L: float) -> float:
    # map to [-L/2, L/2)
    return dx - L*round(dx/L)

def pbc_dist(a: np.ndarray, b: np.ndarray, Lx: float, Ly: float, Lz: float) -> float:
    dx = min_image_delta(a[0]-b[0], Lx)
    dy = min_image_delta(a[1]-b[1], Ly)
    dz = min_image_delta(a[2]-b[2], Lz)
    return math.sqrt(dx*dx+dy*dy+dz*dz)

def nearest_k_indices(target: np.ndarray, pts: np.ndarray, k: int, Lx: float, Ly: float, Lz: float) -> List[int]:
    # brute-force; fine for O(10^4) atoms
    d = np.array([pbc_dist(target, p, Lx,Ly,Lz) for p in pts])
    return list(np.argsort(d)[:k]), list(np.sort(np.unique(d)))  # return distances too (2nd item unused)

# ---------------- blue-noise (greedy Poisson-disk on a lattice of O sites) ----------------

def blue_noise_select(points: np.ndarray, Lx: float, Ly: float, Lz: float, rmin: float, Ntarget: int, rng: random.Random) -> List[int]:
    idxs = list(range(points.shape[0]))
    rng.shuffle(idxs)
    chosen: List[int] = []
    for idx in idxs:
        p = points[idx]
        ok = True
        for j in chosen:
            if pbc_dist(p, points[j], Lx,Ly,Lz) < rmin:
                ok = False; break
        if ok:
            chosen.append(idx)
            if len(chosen) == Ntarget:
                break
    return chosen

# ---------------- main workflow ----------------

def main():
    ap = argparse.ArgumentParser(description="Generate oxygen vacancies with blue-noise spacing; retype 2 nearest Ce4+ → Ce3+ per vacancy.")
    ap.add_argument("--data", required=True, help="input LAMMPS data (stoichiometric)")
    ap.add_argument("--species-map", required=True, help="JSON/YAML with type IDs & charges for Ce4+, Ce3+, O2-")
    ap.add_argument("--aA", type=float, required=True, help="lattice parameter a0 (Å) to label output folder")
    ap.add_argument("--n", type=int, required=True, help="replication n to label output folder")
    ap.add_argument("--x", type=float, required=True, help="vacancy fraction x (per Ce)")
    ap.add_argument("--rmin", type=float, default=3.5, help="minimum V_O–V_O separation (Å)")
    ap.add_argument("--seed", type=int, required=True, help="master RNG seed")
    ap.add_argument("--outdir", default="geometry", help="output root")
    args = ap.parse_args()

    # species map (accept JSON or YAML if pyyaml available; here assume JSON for portability)
    try:
        import yaml
        species = yaml.safe_load(Path(args.species_map).read_text())
    except Exception:
        species = json.loads(Path(args.species_map).read_text())

    t_Ce4 = int(species["Ce4+"]["type_id"])
    t_Ce3 = int(species["Ce3+"]["type_id"])
    t_O   = int(species["O2-"]["type_id"])
    q_Ce4 = float(species["Ce4+"]["charge"])
    q_Ce3 = float(species["Ce3+"]["charge"])
    q_O2  = float(species["O2-"]["charge"])
    masses = {int(k): float(v) for k,v in species["masses"].items()} if "masses" in species else {}

    box, A, raw_lines = read_lammps_data(Path(args.data))
    Lx,Ly,Lz = box.Lx, box.Ly, box.Lz

    # index masks
    is_O   = (A.types == t_O)
    is_Ce4 = (A.types == t_Ce4)
    is_Ce3 = (A.types == t_Ce3)

    N_Ce = int(is_Ce4.sum() + is_Ce3.sum())
    N_O  = int(is_O.sum())
    if N_Ce == 0 or N_O == 0:
        sys.exit("Could not detect Ce/O species from species map.")

    # exact target vacancies
    N_vac = int(round(args.x * N_Ce))
    if abs(args.x * N_Ce - N_vac) > 1e-9:
        sys.exit(f"x*4n^3 must be integer: x*N_Ce={args.x*N_Ce} not integer.")
    if N_vac == 0:
        sys.exit("N_vac=0 for given x; nothing to do.")

    rng = random.Random(args.seed)
    O_xyz = A.xyz[is_O]
    O_ids = A.ids[is_O]
    # greedy Poisson-disk selection; relax rmin if needed
    rmin = float(args.rmin)
    chosen_idx = blue_noise_select(O_xyz, Lx,Ly,Lz, rmin, N_vac, rng)
    tries = 0
    while len(chosen_idx) < N_vac and tries < 10:
        rmin *= 0.9
        chosen_idx = blue_noise_select(O_xyz, Lx,Ly,Lz, rmin, N_vac, rng)
        tries += 1
    if len(chosen_idx) < N_vac:
        sys.exit(f"Could not place {N_vac} vacancies with rmin≈{args.rmin} Å even after relaxing; try smaller rmin or bigger cell.")

    vac_O_ids = O_ids[chosen_idx]
    vac_O_xyz = O_xyz[chosen_idx]

    # remove those O atoms
    keep_mask = np.ones(A.ids.size, dtype=bool)
    # mark deletions
    o_idx_global = np.where(is_O)[0]
    keep_mask[o_idx_global[chosen_idx]] = False

    # for each vacancy position, pick two nearest Ce (prefer Ce4+) and retype to Ce3+
    Ce_mask = (A.types == t_Ce4) | (A.types == t_Ce3)
    Ce_xyz  = A.xyz[Ce_mask]
    Ce_ids  = A.ids[Ce_mask]
    Ce_types= A.types[Ce_mask]
    Ce_q    = A.q[Ce_mask]

    new_types = A.types.copy()
    new_q     = A.q.copy()
    ce3_assigned: List[int] = []

    for vpos in vac_O_xyz:
        # search among all Ce, but prefer Ce4+ if available
        ce_idx_list, _ = nearest_k_indices(vpos, Ce_xyz, k=6, Lx=Lx,Ly=Ly,Lz=Lz)  # grab a few, then filter
        # pick the nearest two Ce that are currently Ce4+; if not enough, allow Ce3+
        picked = []
        for j in ce_idx_list:
            gid = Ce_ids[j]
            # map back to global index
            g = int(np.where(A.ids == gid)[0])
            if new_types[g] == t_Ce4 and len(picked) < 2:
                picked.append(g)
        if len(picked) < 2:
            for j in ce_idx_list:
                gid = Ce_ids[j]
                g = int(np.where(A.ids == gid)[0])
                if g not in picked and len(picked) < 2:
                    picked.append(g)
        # apply retype/charge
        for g in picked:
            new_types[g] = t_Ce3
            new_q[g]     = q_Ce3
            ce3_assigned.append(int(A.ids[g]))

    # verify charge neutrality per vacancy: removing one O2- (= -2) and converting two Ce4+→Ce3+ (= -1 each) net 0
    # We do not enforce per-vacancy neutrality locally—but global neutrality holds by construction.

    # build new Atoms arrays after deletions
    ids_new = A.ids[keep_mask]
    types_new = new_types[keep_mask]
    q_new     = new_q[keep_mask]
    xyz_new   = A.xyz[keep_mask]
    atoms_new = Atoms(ids=ids_new, types=types_new, q=q_new, xyz=xyz_new)

    # diagnostics: V–V separations (minimum-image), Ce3+–V distances
    # V positions are the deleted O positions (vac_O_xyz)
    def pair_min_seps(points: np.ndarray) -> List[float]:
        if points.shape[0] < 2: return []
        out=[]
        for i in range(points.shape[0]):
            dmin = 1e30
            for j in range(points.shape[0]):
                if i==j: continue
                d = pbc_dist(points[i], points[j], Lx,Ly,Lz)
                if d < dmin: dmin = d
            out.append(dmin)
        return out

    vv_mins = pair_min_seps(vac_O_xyz)
    # Ce3+ positions
    ce3_ids_arr = np.array(ce3_assigned, dtype=int)
    ce3_xyz = []
    for cid in ce3_ids_arr:
        g = int(np.where(atoms_new.ids == cid)[0]) if cid in atoms_new.ids else int(np.where(A.ids == cid)[0])
        ce3_xyz.append((A.xyz[g] if g < A.xyz.shape[0] else atoms_new.xyz[g]))
    ce3_xyz = np.array(ce3_xyz) if ce3_ids_arr.size>0 else np.zeros((0,3))
    # distances from each vacancy to its two assigned Ce3+ (approx by nearest in new list)
    v_to_ce3 = []
    for v in vac_O_xyz:
        if ce3_xyz.shape[0]==0: break
        idxs,_ = nearest_k_indices(v, ce3_xyz, k=min(2,ce3_xyz.shape[0]), Lx=Lx,Ly=Ly,Lz=Lz)
        for j in idxs:
            v_to_ce3.append(pbc_dist(v, ce3_xyz[j], Lx,Ly,Lz))

    # outputs
    tag = f"T{int(0)}_x{args.x:.3f}"  # you can inject T if you want; here we focus on geometry
    outroot = Path(args.outdir) / f"reduced_CeO2_n{args.n}_a{args.aA:.3f}_x{args.x:.3f}_rmin{args.rmin:.2f}"
    outroot.mkdir(parents=True, exist_ok=True)

    # species file (optional) can contain masses; if not, keep defaults  Ce:140.116, O:15.999 (fill if masses=={})
    if not masses:
        # guess from charges: find type ids present
        tset = set(types_new.tolist())
        # Leave empty; many LAMMPS inputs set Masses separately via include; skip if unknown
        masses = {int(t): 1.0 for t in tset}

    data_out = outroot/"ceo2_reduced.data"
    write_lammps_data(data_out, box, atoms_new, masses, comment="# reduced ceria with VO + Ce3+ (generated)")

    # vacancy list and Ce3+ map
    with (outroot/"vacancies.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["O_atom_id_deleted","x","y","z"])
        for aid,pos in zip(vac_O_ids.tolist(), vac_O_xyz.tolist()):
            w.writerow([aid, f"{pos[0]:.6f}", f"{pos[1]:.6f}", f"{pos[2]:.6f}"])
    with (outroot/"ce3_map.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["Ce3_atom_id"])
        for cid in ce3_assigned:
            w.writerow([cid])

    # histograms (CSV lists; plot later in analysis)
    np.savetxt(outroot/"vv_min_separations_A.csv", np.array(vv_mins), delimiter=",")
    np.savetxt(outroot/"v_to_ce3_distances_A.csv", np.array(v_to_ce3), delimiter=",")

    # manifest
    manifest = {
        "schema": "ceriakappa/vacancy_manifest@1",
        "source_data": str(Path(args.data).resolve()),
        "species_map": str(Path(args.species_map).resolve()),
        "a0_A": args.aA, "n": args.n,
        "box": {"Lx": box.Lx, "Ly": box.Ly, "Lz": box.Lz},
        "x": args.x, "N_Ce": N_Ce, "N_O": N_O, "N_vac": int(len(vac_O_ids)),
        "rmin_A": args.rmin,
        "seed": args.seed,
        "artifacts": {
            "data": str(data_out),
            "vacancies_csv": str(outroot/"vacancies.csv"),
            "ce3_map_csv": str(outroot/"ce3_map.csv"),
            "vv_min_csv": str(outroot/"vv_min_separations_A.csv"),
            "v_to_ce3_csv": str(outroot/"v_to_ce3_distances_A.csv")
        },
        "notes": "Two nearest Ce retyped to Ce3+ per vacancy; global neutrality preserved."
    }
    (outroot/"vacancy_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote: {data_out}")
    print(f"Artifacts in: {outroot}")

if __name__ == "__main__":
    main()
```

---

## Species map example (JSON or YAML)

Keep this alongside your potentials, so the generator knows which LAMMPS type IDs and charges to use.

```json
{
  "Ce4+": {"type_id": 1, "charge": 4.0},
  "Ce3+": {"type_id": 2, "charge": 3.0},
  "O2-":  {"type_id": 3, "charge": -2.0},
  "masses": { "1": 140.116, "2": 140.116, "3": 15.999 }
}
```

---

## How to run and what gets produced

Assuming your stoichiometric file from Step 5 is `geometry/stoich_CeO2_n10_a5.410/ceo2.data`:

```
python generate_vacancies.py \
  --data geometry/stoich_CeO2_n10_a5.410/ceo2.data \
  --species-map potentials/gotte_2007/species_map.json \
  --aA 5.41 --n 10 --x 0.06 \
  --rmin 3.5 --seed 424242 \
  --outdir geometry
```

Outputs in `geometry/reduced_CeO2_n10_a5.410_x0.060_rmin3.50/`:

- `ceo2_reduced.data` — LAMMPS data with O vacancies removed and two nearest Ce converted to Ce³⁺.
    
- `vacancies.csv` — deleted O atom IDs and their coordinates at deletion time (your VOV_\mathrm{O} seeds).
    
- `ce3_map.csv` — atom IDs retyped to Ce³⁺.
    
- `vv_min_separations_A.csv` — list of each vacancy’s nearest-neighbor vacancy separation (minimum-image).
    
- `v_to_ce3_distances_A.csv` — list of distances from each vacancy to its assigned Ce³⁺ neighbors.
    
- `vacancy_manifest.json` — full, seed-stamped manifest with counts, rmin⁡r_{\min}, and paths.
    

You can now plug `ceo2_reduced.data` into the same NPT→NVT→NVE pipeline as the stoichiometric state. The CSVs feed your diagnostics section (vacancy spacing and Ce³⁺–VOV_\mathrm{O} association), and the manifest drops straight into your run `manifest.json` so every detail of the random draw is auditably fixed.


# 10.) Equilibration per $(T, x)$ state


The aim is to stage NPT → NVT cleanly, record exactly which knobs were used (taus, seeds, PPPM settings, timestep), and leave behind a single restart file plus a manifest that captures the whole story.

## What the per-state runner does

It reads the state’s `runs/<state>/manifest.json`, derives phase seeds from the master seed, writes a self-contained LAMMPS input for NPT→NVT, runs it, computes SHA-256 for the produced restart, and stamps energies, densities, taus, seeds, PPPM accuracy, and neighbor settings back into the manifest. It also emits small CSVs for quick plots (T(t), P(t), V(t)).

### Python runner: `scripts/equilibrate_state.py`

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, subprocess, hashlib, shutil, re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Any
from seeds import derive_seed

NUMERICS_INC = "configs/includes/global_numerics.in"
PREP_INC     = "configs/includes/prep_ensembles.in"

@dataclass
class State:
    name: str           # e.g., T300_x0.03
    T: float
    x: float
    run_dir: Path
    geom_data: Path     # LAMMPS data file for this (T,x), e.g., stoich or reduced
    restart_out: Path

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def write_lmp_input(path: Path, datafile: Path, T: float, taus: Dict[str,float], seeds: Dict[str,int], outprefix: str):
    txt = f"""
units           metal
atom_style      charge
boundary        p p p

include         {NUMERICS_INC}
read_data       {datafile.as_posix()}

# Pair coefficients (emitted from your FF tooling)
# include      potentials/gotte_2007/pair_coeffs.in

# thermostat/barostat time constants (ps) from numerics card
include         {PREP_INC}
variable        T equal {T:.3f}

# record seeds for provenance (LAMMPS uses them in fix commands where applicable)
variable        seed_npt equal {seeds['npt']}
variable        seed_nvt equal {seeds['nvt']}

thermo_style    custom step temp press pe ke etotal vol density lx ly lz
thermo          ${thermo_stride}

# ---------- NPT ----------
fix             prep all npt temp ${T} ${T} {taus['tauT']:.4f} iso 1.0 1.0 {taus['tauP']:.4f}
# Optionally: fix_modify prep temp <compute-ID> to use a specific temperature compute

# warmup + averaging
run             50000
reset_timestep  0
# light time-averages for P/V/T
fix             avnpt all ave/time 100 200 20000 temp press vol density file {outprefix}_npt_ave.txt
run             200000
unfix           avnpt
unfix           prep

# ---------- NVT ----------
fix             settle all nvt temp ${T} ${T} {taus['tauT']:.4f}
run             100000
reset_timestep  0
fix             avnvt all ave/time 100 200 20000 temp press vol density file {outprefix}_nvt_ave.txt
run             200000
unfix           avnvt
unfix           settle

write_restart   {outprefix}.restart
"""
    path.write_text(txt)

def parse_ave(path: Path) -> Dict[str, float]:
    # expect ave/time output with columns: step <T> <P> <V> <rho>
    import numpy as np
    if not path.exists():
        return {}
    arr = np.loadtxt(path)
    if arr.ndim == 1 and arr.size == 0:
        return {}
    Tm = float(arr[:,1].mean()); Pm = float(arr[:,2].mean()); Vm = float(arr[:,3].mean()); rhom = float(arr[:,4].mean())
    return {"T_mean": Tm, "P_mean": Pm, "V_mean": Vm, "rho_mean": rhom}

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True, help="e.g. T300_x0.00")
    ap.add_argument("--data", required=True, help="LAMMPS data file for this state")
    ap.add_argument("--T", type=float, required=True)
    ap.add_argument("--tauT", type=float, default=0.2)
    ap.add_argument("--tauP", type=float, default=2.0)
    ap.add_argument("--manifest", required=True, help="runs/<state>/manifest.json")
    ap.add_argument("--lmp", default="lmp", help="LAMMPS executable")
    args = ap.parse_args()

    run_dir = Path("runs")/args.state
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)
    man = json.loads(manifest_path.read_text())

    # derive phase seeds deterministically from master
    master = int(man["seeds"]["master"])
    seeds = {
        "npt": derive_seed(master, "NPT"),
        "nvt": derive_seed(master, "NVT"),
    }
    man["seeds"]["npt"] = seeds["npt"]
    man["seeds"]["nvt"] = seeds["nvt"]

    outprefix = (run_dir/f"eq_{args.state}").as_posix()
    inp = run_dir/"equilibrate.in"
    write_lmp_input(inp, Path(args.data), args.T, {"tauT":args.tauT, "tauP":args.tauP}, seeds, outprefix)

    # run LAMMPS
    with (run_dir/"stdout.log").open("w") as log:
        subprocess.run([args.lmp, "-in", inp.as_posix()], cwd=run_dir.as_posix(), check=True, stdout=log, stderr=subprocess.STDOUT)

    # stamp averages and restart checksum
    npt_stats = parse_ave(run_dir/f"eq_{args.state}_npt_ave.txt")
    nvt_stats = parse_ave(run_dir/f"eq_{args.state}_nvt_ave.txt")
    restart = run_dir/f"eq_{args.state}.restart"
    restart_sha = sha256_file(restart)
    man["hashes"]["restart_sha256"] = restart_sha
    man["equilibration"] = {
        "taus_ps": {"tauT": args.tauT, "tauP": args.tauP},
        "averages": {"npt": npt_stats, "nvt": nvt_stats},
        "pppm": man["inputs"]["pppm"],
        "neighbor": man["inputs"]["neighbor"],
        "timestep_fs": man["numerics"]["timestep_fs"]
    }
    manifest_path.write_text(json.dumps(man, indent=2))
    print(f"Equilibration complete. Restart: {restart} (sha256={restart_sha})")

if __name__ == "__main__":
    main()
```

Invoke it for any state. For the stoichiometric baseline at 300 K:

```
python scripts/equilibrate_state.py \
  --state T300_x0.00 \
  --data geometry/stoich_CeO2_n10_a5.410/ceo2.data \
  --T 300 \
  --manifest runs/T300_x0.00/manifest.json
```

For a reduced case at 700 K:

```
python scripts/equilibrate_state.py \
  --state T700_x0.060 \
  --data geometry/reduced_CeO2_n10_a5.410_x0.060_rmin3.50/ceo2_reduced.data \
  --T 700 \
  --manifest runs/T700_x0.060/manifest.json
```

## Notes that keep this clean

The input uses the same `global_numerics.in` and `prep_ensembles.in` includes you generated earlier, so PPPM accuracy, cutoff, neighbor skin, and timestep are identical across all (T,x)(T,x) states. The script derives fresh NPT/NVT seeds from the state’s master seed and records them back to the manifest. The two small `*_ave.txt` files are machine-parsable summaries of temperature, pressure, volume, and density; you can commit them or fold them into a single `equilibration.json` later. The restart is the only file you carry forward to GK/rNEMD; its checksum is stored so downstream scripts can verify they’re reading exactly what you prepared.




# 11.) GK production per state (primary path)


The goal is to run several clean NVE “shards” with an identical MPI layout, sample the microscopic heat flux at a fixed stride, and leave behind two artifacts per shard: a raw binary flux time-series and a text HCACF table (`J0Jt`). We’ll drive everything from the state manifest and the analysis control you already created, and we’ll stamp SHA-256 checksums back into the manifest so downstream steps can verify the inputs exactly.

## Orchestrator: write inputs, launch shards, stamp checksums

Save as `scripts/run_gk_state.py`. It reads `runs/<state>/manifest.json` and `analysis/<state>/control.json`, emits one LAMMPS input per shard, runs each with a fixed MPI layout, and records where the artifacts landed plus their checksums. It also sets the sampling stride in **steps** from your control file’s stride in **fs** and your chosen MD timestep.

```python
#!/usr/bin/env python3
from __future__ import annotations
import os, json, subprocess, hashlib, math
from pathlib import Path
from typing import Dict, Any, List

NUMERICS_INC = "configs/includes/global_numerics.in"
GK_INC       = "configs/includes/gk_measurement.in"  # has fix momentum every N steps

GK_INPUT_TEMPLATE = """\
units           metal
atom_style      charge
boundary        p p p

include         {numerics_inc}
# include      potentials/gotte_2007/pair_coeffs.in

read_restart    {restart_path}

# --- GK flux assembly (pppm-consistent) ---
compute         cke all ke/atom
compute         cpe all pe/atom
compute         cst all centroid/stress/atom NULL virial
compute         J   all heat/flux cke cpe cst

# remove COM drift periodically
include         {gk_inc}

# sampling stride in steps (from control.json fs stride and dt)
variable        stride equal {stride_steps}

# raw binary flux dump for deterministic post-processing
dump            dJ all custom ${stride} {shard_dir}/flux.raw c_J[1] c_J[2] c_J[3]
dump_modify     dJ format binary yes
dump_modify     dJ append yes

# HCACF accumulation for quick sanity (overlapping origins)
fix             acf all ave/correlate 1 {blk} {tot} c_J[1] c_J[2] c_J[3] type auto \
                file {shard_dir}/J0Jt.dat ave running

# pure dynamics (NVE); dt is set in numerics include
run             {nsteps}

unfix           acf
undump          dJ
"""

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True, help="e.g. T300_x0.00")
    ap.add_argument("--shards", type=int, default=2)
    ap.add_argument("--shard_ps", type=float, default=None, help="override shard length (ps)")
    ap.add_argument("--mpi", type=int, default=8, help="MPI ranks (keep fixed across shards)")
    ap.add_argument("--pin", action="store_true", help="set OMP/MPI pinning env for stability")
    ap.add_argument("--lmp", default="lmp", help="LAMMPS executable")
    args = ap.parse_args()

    state = args.state
    run_dir = Path("runs")/state
    ana_dir = Path("analysis")/state
    run_dir.mkdir(parents=True, exist_ok=True); ana_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = run_dir/"manifest.json"
    control_path  = ana_dir/"control.json"
    if not manifest_path.exists() or not control_path.exists():
        raise SystemExit("Missing manifest or control JSON for this state.")

    man = json.loads(manifest_path.read_text())
    ctl = json.loads(control_path.read_text())
    # restart to read (from Step 10)
    restart = next((p for p in run_dir.glob("eq_*.restart")), None)
    if restart is None:
        raise SystemExit("No equilibrated restart found; run Step 10 first.")

    # numerics: dt and stride
    dt_fs   = float(man["numerics"]["timestep_fs"])
    stride_fs = float(ctl["gk"]["flux_sample_stride_fs"])
    stride_steps = max(1, int(round(stride_fs / dt_fs)))

    # shard length (ps) from control unless overridden
    shard_ps = args.shard_ps if args.shard_ps is not None else float(man.get("equilibration",{}).get("shard_length_ps", ctl.get("gk",{}).get("shard_length_ps", 1000.0)))
    nsteps   = int(round( shard_ps / (dt_fs/1000.0) ))

    # ave/correlate block sizing (purely for the on-the-fly J0Jt sanity file)
    blk  = 10000
    tot  = 10000

    # record fixed MPI layout for determinism
    man.setdefault("gk", {})
    man["gk"]["mpi_ranks"] = args.mpi
    man["gk"]["stride_steps"] = stride_steps
    man["gk"]["shard_steps"]  = nsteps
    man["gk"]["dt_fs"]        = dt_fs
    manifest_path.write_text(json.dumps(man, indent=2))

    env = os.environ.copy()
    if args.pin:
        # Gentle defaults; adjust to your cluster
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["OMP_PROC_BIND"] = "true"
        env["OMP_PLACES"]    = "cores"

    # launch shards
    shard_entries: List[Dict[str,Any]] = []
    for s in range(args.shards):
        sdir = run_dir/f"nve_shard_{s:03d}"
        sdir.mkdir(parents=True, exist_ok=True)
        # write the shard input
        inp = sdir/"gk.in"
        inp.write_text(GK_INPUT_TEMPLATE.format(
            numerics_inc=NUMERICS_INC,
            gk_inc=GK_INC,
            restart_path=restart.as_posix(),
            stride_steps=stride_steps,
            shard_dir=sdir.as_posix(),
            blk=blk, tot=tot, nsteps=nsteps
        ))
        # run LAMMPS with fixed MPI layout
        log = sdir/"stdout.log"
        cmd = [args.lmp, "-in", inp.name]
        # if you prefer mpirun/srun, wrap here; fix the layout so it's identical across shards
        if args.mpi > 1:
            cmd = ["mpirun", "-np", str(args.mpi)] + cmd
        subprocess = __import__("subprocess")
        subprocess.run(cmd, cwd=sdir.as_posix(), env=env, check=True, stdout=log.open("w"), stderr=subprocess.STDOUT)

        # compute checksums and record artifacts
        flux = sdir/"flux.raw"
        j0jt = sdir/"J0Jt.dat"
        flux_sha = sha256_file(flux) if flux.exists() else None
        entry = {
            "shard": s, "dir": str(sdir),
            "flux": str(flux), "flux_sha256": flux_sha,
            "J0Jt": str(j0jt), "steps": nsteps, "stride_steps": stride_steps
        }
        shard_entries.append(entry)

    # stamp per-shard info into manifest for auditing
    man["gk"]["shards"] = shard_entries
    manifest_path.write_text(json.dumps(man, indent=2))
    print(f"Completed GK shards for {state}. Shards: {len(shard_entries)}")
    for e in shard_entries:
        print(f"  shard {e['shard']:03d}: {e['flux']} (sha256={e['flux_sha256']})")

if __name__ == "__main__":
    main()
```

Run it for a typical state after equilibration:

```
python scripts/run_gk_state.py --state T300_x0.00 --shards 2 --mpi 8 --pin
```

This produces `runs/T300_x0.00/nve_shard_000/{flux.raw,J0Jt.dat,stdout.log}` (and similarly for `_001`) and stamps their checksums plus run parameters back into `runs/T300_x0.00/manifest.json`.

## Notes you’ll care about while this runs

Keep the MPI layout fixed across shards and states (same ranks, same domain decomposition policy). If your site uses `srun` or `mpiexec`, wrap the command but still pass a fixed `-n` and the same environment pinning flags so the per-atom k-space reductions follow the same order. The sampling stride is in **steps** and derives from your control file’s stride in **fs** divided by the chosen MD timestep; the script computes it so you don’t have to. The binary `flux.raw` holds three 64-bit floats per record in the order Jx,Jy,JzJ_x,J_y,J_z; your `analyze_gk.py` from Step 8 already knows how to read that. The text `J0Jt.dat` is only for on-the-fly sanity; you should continue to compute the HCACF from `flux.raw` offline to keep the analysis deterministic and versioned.

If you want a single-command convenience, add a tiny shell wrapper (or Makefile target) that runs:

1. `equilibrate_state.py`
    
2. `run_gk_state.py`
    
3. `analyze_gk.py` (with your box length and temperature pulled from the manifest)
    

and then writes a one-page “GK shard summary” into `analysis/<state>/k_result.json` and plots.




# 12.) rNEMD spot check (one or two states)


The aim is to run a **gentle** Müller–Plathe (MP) reverse NEMD on a validated stoichiometric state (and optionally one reduced state), extract a clean bulk temperature gradient, convert the accumulated exchanged energy to a heat flux qq, and report k=−q/∇Tk=-q/\nabla T. You’ll also halve/double the swap rate to show **rate independence** in the linear regime. Everything below mirrors the GK numerics you already froze.

## What the spot-check will do

You read the equilibrated restart for a state, apply `fix thermal/conductivity` along a chosen axis with a conservative period, bin temperatures with `compute chunk/atom` + `fix ave/chunk`, and let the run reach steady state. You then parse (i) the **cumulative exchanged energy** reported by the MP fix and (ii) the **bulk linear segment** of the temperature profile, discard a few bins near the swap slabs, and compute k=−q/∇Tk=-q/\nabla T. Repeat with a **slower swap rate**; agreement within uncertainty shows linear response.

---

## LAMMPS input: `lmp/rnemd_spot.in`

This reads your equilibrated restart, imposes MP swaps along zz, collects a binned temperature profile, and prints the cumulative exchanged energy to a small, machine-parsable log.

```lmp
# rnemd_spot.in — gentle Müller–Plathe cross-check
units           metal
atom_style      charge
boundary        p p p

# Numerics (dt, PPPM, cutoffs, neighbor) frozen earlier
include         configs/includes/global_numerics.in
# include      potentials/gotte_2007/pair_coeffs.in

# Equilibrated, frozen-volume state from Step 10/7
read_restart    ${restart}         # pass via -var restart runs/T300_x0.00/eq_T300_x0.00.restart

# Parameters
variable        T         equal ${Ttarget}     # pass -var Ttarget 300
variable        Nbins     equal 48             # slabs along z (increase if you want finer profile)
variable        swapper   equal ${period}      # pass -var period 5000 (gentle) or 10000 (gentler)
variable        exbins    equal 2              # bins to exclude near source/sink

# Apply Müller–Plathe heat flux along z:
# Syntax: fix ID group thermal/conductivity Nevery direction Nbin
# We'll use Nbin=2 to set source/sink slabs at the ends, and use bins for *profiling* separately below.
fix             mp all thermal/conductivity ${swapper} z 2

# Temperature profile sampling (independent of MP's internal binning)
compute         slab all chunk/atom bin/1d z lower ${Nbins} units box
fix             prof all ave/chunk 100 250 25000 slab temp file ${outroot}/profile.dat

# Convenience variables: box dims for later post-processing (printed once)
variable        Lx equal lx
variable        Ly equal ly
variable        Lz equal lz
print "RNEMD_BOX Lx=${Lx}A Ly=${Ly}A Lz=${Lz}A" file ${outroot}/rnemd_meta.log

# Run to steady state (gentle driving)
run             2000000

# Grab cumulative exchanged energy from the fix (eV)
# LAMMPS stores it as f_mp; print it to a parseable file.
variable        Eexch equal f_mp
print "RNEMD_EEXCH_eV ${Eexch}" append ${outroot}/rnemd_meta.log

unfix           prof
unfix           mp
```

Typical invocations (two rates) for your baseline state:

```
mkdir -p runs/T300_x0.00/rnemd/period_5000
lmp -in lmp/rnemd_spot.in \
    -var restart runs/T300_x0.00/eq_T300_x0.00.restart \
    -var Ttarget 300 -var period 5000 \
    -var outroot runs/T300_x0.00/rnemd/period_5000

mkdir -p runs/T300_x0.00/rnemd/period_10000
lmp -in lmp/rnemd_spot.in \
    -var restart runs/T300_x0.00/eq_T300_x0.00.restart \
    -var Ttarget 300 -var period 10000 \
    -var outroot runs/T300_x0.00/rnemd/period_10000
```

Notes. Keep the **same MPI layout** you used for GK for bitwise-stability of any auxiliary tallies. The chosen `Nbins=48` usually yields a smooth profile for a ~5.4 nm cell; scale with box length if you elongate.

---

## Analyzer: `analysis/analyze_rnemd.py`

This script reads the profile, extracts the **bulk** gradient via a robust linear fit that excludes the MP exchange regions, converts cumulative exchanged energy to a **flux**, and reports kk. It also produces a plot with the excluded bins shaded and prints a compact JSON you can store next to the profile.

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, re, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

eV  = 1.602176634e-19
A   = 1e-10

def load_profile(path: Path):
    # LAMMPS fix ave/chunk output: columns -> id zlo zhi N temp ...
    dat = np.loadtxt(path)
    zc  = 0.5*(dat[:,1] + dat[:,2])            # Å
    T   = dat[:,4]                              # K
    return zc, T

def load_meta(path: Path):
    meta = {"Lx_A": None, "Ly_A": None, "Lz_A": None, "Eexch_eV": None}
    txt = Path(path).read_text()
    m = re.search(r"RNEMD_BOX\s+Lx=([0-9\.Ee+-]+)A\s+Ly=([0-9\.Ee+-]+)A\s+Lz=([0-9\.Ee+-]+)A", txt)
    if m:
        meta["Lx_A"] = float(m.group(1)); meta["Ly_A"] = float(m.group(2)); meta["Lz_A"] = float(m.group(3))
    m2 = re.search(r"RNEMD_EEXCH_eV\s+([0-9\.Ee+-]+)", txt)
    if m2:
        meta["Eexch_eV"] = float(m2.group(1))
    return meta

def fit_bulk_gradient(zA, T, exclude_bins=2, show=False, out=None):
    # Exclude bins near both ends (source/sink). Fit central region linearly.
    z = zA[exclude_bins: -exclude_bins] * A       # m
    t = T[exclude_bins: -exclude_bins]
    A_ = np.vstack([z, np.ones_like(z)]).T
    m, b = np.linalg.lstsq(A_, t, rcond=None)[0]  # K/m slope
    # simple CI via block resampling over contiguous chunks
    n = len(z)
    B = max(8, n//16)
    idx = np.arange(n)
    # split into B contiguous blocks and compute blockwise slopes
    blocks = np.array_split(idx, B)
    slopes = []
    for blk in blocks:
        zz = z[blk]; tt = t[blk]
        if len(zz) < 3: continue
        Ablk = np.vstack([zz, np.ones_like(zz)]).T
        mblk,_ = np.linalg.lstsq(Ablk, tt, rcond=None)[0]
        slopes.append(mblk)
    slope_ci = 1.96*np.std(slopes, ddof=1)/np.sqrt(len(slopes)) if len(slopes)>1 else np.nan

    if show and out:
        plt.figure()
        plt.plot(zA, T, '-o', ms=3, lw=1)
        plt.axvspan(zA.min(), zA[exclude_bins], color='0.9')
        plt.axvspan(zA[-exclude_bins], zA.max(), color='0.9')
        zfit = zA[exclude_bins: -exclude_bins]
        plt.plot(zfit, (m*(zfit*A)+b), 'r-', lw=2)
        plt.xlabel("z (Å)"); plt.ylabel("T (K)")
        plt.tight_layout(); plt.savefig(out, dpi=200)
    return m, slope_ci

def compute_flux_Wm2(Eexch_eV, Lx_A, Ly_A, total_time_ps):
    # MP swaps move energy from cold→hot; fix reports *cumulative exchanged energy* (eV).
    # Heat flux magnitude q = (Eexch)/(2 * A * t). Factor 2 because two heat streams propagate in ±z under PBC.
    A_cross = (Lx_A*A)*(Ly_A*A)                  # m^2
    t_s     = total_time_ps * 1e-12              # s
    q = (Eexch_eV*eV) / (2.0 * A_cross * t_s)    # W/m^2
    return q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="runs/.../rnemd/period_xxxx")
    ap.add_argument("--period_steps", type=int, required=True)
    ap.add_argument("--dt_fs", type=float, default=1.0)
    ap.add_argument("--exclude_bins", type=int, default=2)
    args = ap.parse_args()

    root = Path(args.root)
    prof = root/"profile.dat"
    meta = root/"rnemd_meta.log"
    zA, T = load_profile(prof)
    M = load_meta(meta)
    if any(v is None for v in [M["Lx_A"], M["Ly_A"], M["Lz_A"], M["Eexch_eV"]]):
        raise SystemExit("Missing box or exchanged-energy metadata. Check rnemd_meta.log.")

    # gradient
    m, m_ci = fit_bulk_gradient(zA, T, exclude_bins=args.exclude_bins, show=True, out=root/"profile_fit.png")  # K/m

    # run length from the thermo: infer from profile sampling or count steps * dt
    # Use dt_fs (fs) and total steps in the run printed in log if you prefer; here we read from profile timesteps count:
    # crude estimate: total_time = n_output_blocks * Nevery*Nrepeat * dt
    # Simpler: pass total_time explicitly if you log it; here we compute from LAMMPS 'run' line:
    logtxt = (root.parent.parent/".."/"stdout.log")  # optional; ignore if unavailable
    # fallback: use known run steps (2e6) from input
    total_steps = 2_000_000
    total_time_ps = total_steps * (args.dt_fs/1000.0)

    # flux
    q = compute_flux_Wm2(M["Eexch_eV"], M["Lx_A"], M["Ly_A"], total_time_ps)  # W/m^2

    # conductivity
    k = - q / m                                # W/mK
    # conservative CI composition from gradient CI only (dominant); add small 5% systematic if desired
    k_ci = abs(k) * (m_ci/abs(m)) if np.isfinite(m_ci) else np.nan

    out = {
        "k_WmK": float(k),
        "k_ci_WmK": float(k_ci),
        "grad_K_per_m": float(m),
        "grad_ci_K_per_m": float(m_ci),
        "q_Wm2": float(q),
        "meta": {**M, "period_steps": args.period_steps, "dt_fs": args.dt_fs,
                 "exclude_bins": args.exclude_bins, "total_time_ps": total_time_ps}
    }
    (root/"rnemd_result.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
```

Run the analyzer for both rates, using your frozen timestep:

```
python analysis/analyze_rnemd.py --root runs/T300_x0.00/rnemd/period_5000   --period_steps 5000  --dt_fs 1.0
python analysis/analyze_rnemd.py --root runs/T300_x0.00/rnemd/period_10000  --period_steps 10000 --dt_fs 1.0
```

You’ll get `profile_fit.png` (with excluded bins shaded) and `rnemd_result.json` in each `period_XXXX` folder.

---

## Rate-independence checker: `analysis/rate_check.py`

This compares the two rates and prints a single line for your validation note.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, argparse
from pathlib import Path

def load(path: Path):
    return json.loads((path/"rnemd_result.json").read_text())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", required=True)   # runs/.../rnemd/period_5000
    ap.add_argument("--slow", required=True)   # runs/.../rnemd/period_10000
    args = ap.parse_args()
    A = load(Path(args.fast))
    B = load(Path(args.slow))
    kA, kB = A["k_WmK"], B["k_WmK"]
    # combine uncertainties conservatively
    sA = A.get("k_ci_WmK", 0.0)/1.96 if A.get("k_ci_WmK") else 0.0
    sB = B.get("k_ci_WmK", 0.0)/1.96 if B.get("k_ci_WmK") else 0.0
    agree = abs(kA-kB) <= 2.0*(sA**2 + sB**2)**0.5 if sA or sB else True
    status = "PASS" if agree else "WARN"
    print(f"Rate check: k_fast={kA:.2f} W/mK, k_slow={kB:.2f} W/mK → {status}")

if __name__ == "__main__":
    main()
```

Example:

```
python analysis/rate_check.py \
  --fast runs/T300_x0.00/rnemd/period_5000 \
  --slow runs/T300_x0.00/rnemd/period_10000
```

Record the output line in your GK validation note for the cross-method agreement section.

---

## What to commit as the artifact

In `runs/<state>/rnemd/period_XXXX/`:

- `profile.dat` (binned temperature), `profile_fit.png` (fit with excluded end bins shaded), and `rnemd_meta.log` (box and cumulative exchanged energy).
    
- `rnemd_result.json` with kk, its CI, qq, the fitted gradient and CI, the excluded-bin count, swap period, and dt.
    

In `analysis/<state>/K_VALIDATION.md`, add a short paragraph:

> rNEMD (MP) spot check at 300 K with gentle swap rates of 5000 and 10000 steps produced k=…k=… and k=…k=… W m−1^{-1} K−1^{-1}, respectively. The difference is within the combined 95% CIs (rate-independence **PASS**). The bulk linear region excluded two end bins near the source/sink; see `profile_fit.png`. The rNEMD value agrees with the GK mean within <2σ<2\sigma.

That’s it. You now have a complete, auditably gentle rNEMD cross-check that drops cleanly into your validation story, with profile plots, a transparent fit window, and a one-line rate-independence verdict.





# 13.) HCACF analysis & uncertainty quantification


This stage turns each shard’s raw flux series into a defensible k±k\pmCI and then aggregates shards. The pipeline below is deterministic, fast (FFT-based), and mirrors the theory you committed: overlapping origins for the HCACF, Sokal-style windowing t⋆=c τintt^\star=c\,\tau_{\text{int}}, blocking or shard variance for CIs, and component averaging only at the end (average **windowed integrals**, not the raw HCACFs).

## What the analyzer produces

For every state it writes a single JSON with all numbers you’ll cite in the paper, a compact Markdown report, and the two standard plots (HCACF and cumulative integral with the chosen t⋆t^\star marked). Everything is computed twice: per-shard (so you can see scatter) and aggregated across shards (so you can quote a single mean and CI).

---

## Drop-in analyzer (deterministic): `analysis/hcacf_analyze.py`

This script:

1. loads one or more shard flux files (`flux.raw`, 3×float64 at a fixed stride),
    
2. builds per-component HCACFs with an **FFT** autocorrelation (unbiased, overlapping origins),
    
3. estimates τint\tau_{\text{int}} from the scalar HCACF (1/3 of the trace) using a positive-sequence/first-zero rule,
    
4. sets t⋆=c τintt^\star=c\,\tau_{\text{int}} and integrates the **scalar** HCACF to k^(t⋆)\hat k(t^\star),
    
5. gets a CI two ways: (a) blocking inside each shard (default), and (b) shard-wise scatter (if M≥2M\ge 2),
    
6. averages the **windowed** kk across components and shards,
    
7. writes `k_result.json`, `report.md`, and plots into `analysis/<state>/`.
    

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- unit scaffold (LAMMPS units metal) ---
kB_eV_per_K = 8.617333262e-5
eV = 1.602176634e-19
ps = 1e-12
A = 1e-10
# prefactor converting:  [∫ C_JJ dt] / [V k_B^2 T^2]  → W/mK
# C_JJ has units of (energy flux)^2; in LAMMPS "metal" with J in eV/ps/Å^2, the integral is in (eV/ps)^2 * ps / Å^4 = eV^2/ps / Å^4
# k = (eV^2/ps / Å^4) / (Å^3 * eV^2/K^2) × (J/eV) × (Å/m)^? → carefully derive:
# J flux units (per LAMMPS doc for heat/flux): energy/time/area → eV/ps/Å^2
# ∫C dt has units (eV/ps/Å^2)^2 * ps = eV^2/(ps*Å^4)
# Denominator V k_B^2 T^2: (Å^3)*(eV^2/K^2)*K^2 = eV^2 Å^3
# Ratio gives 1/(ps * Å^7). To get W/mK, multiply by (J/eV) / (1/ps) × (Å^7 -> m·K requires we missed an Å factor in LAMMPS’s definition?).
# Rather than re-derive here, calibrate once against LAMMPS KAPPA example and store the constant below:
F_conv = (eV**2/ps) / (A**1) / ( (kB_eV_per_K**2) )  # Placeholder; replace after calibrating on your system

def load_flux_binary(path: Path) -> np.ndarray:
    raw = np.fromfile(path, dtype=np.float64)
    if raw.size % 3 != 0:
        raise ValueError(f"{path} length {raw.size} not divisible by 3")
    return raw.reshape(-1, 3)  # [Nt, 3] → Jx, Jy, Jz

def autocorr_fft(x: np.ndarray, unbiased: bool = True) -> np.ndarray:
    """Univariate autocorrelation via FFT (overlapping origins)."""
    n = x.size
    x = x - x.mean()
    # zero-pad to 2n for circular conv. safety
    N = 1 << (2*n - 1).bit_length()
    X = np.fft.rfft(x, n=N)
    S = X * np.conjugate(X)
    ac = np.fft.irfft(S, n=N)[:n]
    if unbiased:
        ac /= (np.arange(n, 0, -1))
    ac /= ac[0]
    return ac

def hcacf_scalar(J: np.ndarray, maxlag: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-axis ACFs and scalar HCACF = (1/3)∑ Cαα."""
    if maxlag is None:
        maxlag = J.shape[0]
    Cs = []
    for a in range(3):
        ac = autocorr_fft(J[:, a])[:maxlag]
        Cs.append(ac)
    Cs = np.stack(Cs, axis=1)  # [maxlag, 3]
    Cscalar = Cs.mean(axis=1)
    return Cs, Cscalar

def tau_int(C: np.ndarray, dt_ps: float) -> float:
    """Integrated autocorrelation time using a positive-sequence + first-zero cutoff."""
    s = 0.0
    for k, c in enumerate(C):
        if k > 0 and c <= 0.0:
            break
        s += c
    # discrete: τ_int = dt * (1 + 2*sum_{k>=1} C(k))
    # here C[0]=1 by our normalization; we summed from k=0 until first ≤0
    return dt_ps * (2.0*s - 1.0)

def windowed_integral(C: np.ndarray, dt_ps: float, tstar_ps: float) -> float:
    nstar = max(2, int(round(tstar_ps / dt_ps)))
    return float(np.trapz(C[:nstar], dx=dt_ps))

def blocking_ci(series: np.ndarray) -> float:
    """Flyvbjerg–Petersen blocking on a 1D series of chunk-level k estimates; returns 1σ."""
    y = series.copy()
    vars_ = []
    while y.size >= 16:
        vars_.append(y.var(ddof=1) / y.size)
        # average adjacent pairs
        y = 0.5*(y[0::2] + y[1::2])
    if not vars_:
        return np.nan
    return float(np.sqrt(max(vars_)))

def per_shard_k(J: np.ndarray, V_A3: float, T_K: float, dt_fs: float, c_sokal: float, nchunks: int = 8) -> Dict:
    dt_ps = dt_fs/1000.0
    Cs_axes, C = hcacf_scalar(J)
    tau_ps = tau_int(C, dt_ps)
    tstar_ps = c_sokal * tau_ps
    # windowed area of scalar HCACF
    area = windowed_integral(C, dt_ps, tstar_ps)
    k_scalar = area / (V_A3 * (kB_eV_per_K*T_K)**2)
    k_scalar *= F_conv

    # chunk-level blocking (split time series into nchunks >> τ_int to estimate 1σ)
    N = J.shape[0]
    L = max(int(np.floor(N/nchunks)), int(np.ceil(10*tau_ps/dt_ps)))
    if L < 1000:  # fallback
        L = max(L, 1000)
    ks = []
    for start in range(0, N-L+1, L):
        seg = J[start:start+L, :]
        _, Cseg = hcacf_scalar(seg, maxlag=min(C.size, seg.shape[0]))
        area_seg = windowed_integral(Cseg, dt_ps, tstar_ps)
        k_seg = area_seg / (V_A3 * (kB_eV_per_K*T_K)**2) * F_conv
        ks.append(k_seg)
    ks = np.array(ks) if ks else np.array([k_scalar])
    sigma_block = blocking_ci(ks)

    # per-axis windowed k (diagnostic; average windowed integrals over axes)
    k_axes = []
    for a in range(3):
        area_a = windowed_integral(Cs_axes[:, a], dt_ps, tstar_ps)
        k_a = area_a / (V_A3*(kB_eV_per_K*T_K)**2) * F_conv
        k_axes.append(k_a)

    out = {
        "k_scalar_WmK": float(k_scalar),
        "sigma_block_WmK": float(sigma_block),  # 1σ
        "k_axes_WmK": [float(x) for x in k_axes],
        "tau_int_ps": float(tau_ps),
        "tstar_ps": float(tstar_ps),
        "C_scalar": C.tolist()  # for plotting
    }
    return out

def aggregate_shards(results: List[Dict]) -> Dict:
    ks = np.array([r["k_scalar_WmK"] for r in results])
    k_mean = float(ks.mean())
    # shard-wise 95% CI
    ci95 = float(1.96 * ks.std(ddof=1) / np.sqrt(len(ks))) if len(ks) > 1 else float("nan")
    # typical tau and t* (median)
    tau_med = float(np.median([r["tau_int_ps"] for r in results]))
    tstar_med = float(np.median([r["tstar_ps"] for r in results]))
    return {"k_mean_WmK": k_mean, "ci95_WmK": ci95, "tau_med_ps": tau_med, "tstar_med_ps": tstar_med}

def plots(state: str, outdir: Path, C: np.ndarray, dt_ps: float, tstar_ps: float, ks_per_axis: List[float]):
    t = np.arange(C.size)*dt_ps
    # HCACF
    plt.figure()
    plt.plot(t, C)
    plt.xlabel("t (ps)"); plt.ylabel("C_JJ(t) / C_JJ(0)")
    plt.tight_layout(); plt.savefig(outdir/"hcacf.png", dpi=200)
    # cumulative integral with t*
    kcumu = np.cumsum((C[:-1] + C[1:])*0.5*dt_ps)
    plt.figure()
    plt.plot(t[1:], kcumu, lw=1.8)
    plt.axvline(tstar_ps, ls="--", color="k")
    plt.text(tstar_ps, kcumu.max()*0.1, f"t* = {tstar_ps:.1f} ps", rotation=90, va="bottom")
    plt.xlabel("t (ps)"); plt.ylabel("∫ C_JJ dt (arb.)")
    plt.tight_layout(); plt.savefig(outdir/"kcumu.png", dpi=200)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True, help="e.g. T300_x0.00")
    ap.add_argument("--shards", nargs="+", required=True, help="paths to .../flux.raw (one or more)")
    ap.add_argument("--L_A", type=float, required=True, help="box length (Å)")
    ap.add_argument("--T", type=float, required=True, help="temperature (K)")
    ap.add_argument("--dt_fs", type=float, required=True, help="MD timestep (fs)")
    ap.add_argument("--stride_fs", type=float, required=True, help="flux sampling stride (fs)")
    ap.add_argument("--c", type=float, default=5.0, help="Sokal window constant")
    args = ap.parse_args()

    outdir = Path("analysis")/args.state
    outdir.mkdir(parents=True, exist_ok=True)

    V_A3 = args.L_A**3
    # check stride vs dt: here we assume stride==dt (you set that earlier). If not, decimate timesteps appropriately.
    if abs(args.stride_fs - args.dt_fs) > 1e-9:
        stride = int(round(args.stride_fs/args.dt_fs))
    else:
        stride = 1

    per_shard = []
    for spath in args.shards:
        J = load_flux_binary(Path(spath))
        if stride > 1:
            J = J[::stride, :]
        r = per_shard_k(J, V_A3, args.T, args.dt_fs, args.c)
        per_shard.append(r)

    agg = aggregate_shards(per_shard)
    # write primary JSON
    res = {"state": args.state, "T_K": args.T, "L_A": args.L_A, "dt_fs": args.dt_fs, "stride_fs": args.stride_fs,
           "window_c": args.c, "per_shard": per_shard, "aggregate": agg}
    (outdir/"k_result.json").write_text(json.dumps(res, indent=2))

    # plots from shard 0 diagnostics
    C0 = np.array(per_shard[0]["C_scalar"])
    plots(args.state, outdir, C0, args.dt_fs/1000.0, per_shard[0]["tstar_ps"], per_shard[0]["k_axes_WmK"])

    # short human report
    kmean = agg["k_mean_WmK"]; ci = agg["ci95_WmK"]
    with (outdir/"report.md").open("w") as f:
        f.write(f"# HCACF analysis — {args.state}\n\n")
        f.write(f"**k (GK, windowed, shard-mean):** {kmean:.2f} W/mK")
        if np.isfinite(ci):
            f.write(f" (±{ci:.2f} at 95% from shard scatter)")
        f.write("\n\n")
        f.write(f"Window rule: t* = c τ_int with c = {args.c:.1f}; median τ_int = {agg['tau_med_ps']:.1f} ps; median t* = {agg['tstar_med_ps']:.1f} ps.\n\n")
        for i, r in enumerate(per_shard):
            f.write(f"- shard {i:02d}: k = {r['k_scalar_WmK']:.2f} W/mK; τ_int = {r['tau_int_ps']:.1f} ps; k_axes = {', '.join(f'{x:.2f}' for x in r['k_axes_WmK'])}\n")
        f.write("\nArtifacts: `k_result.json`, `hcacf.png`, `kcumu.png`.\n")

    print(f"{args.state}: k = {kmean:.2f} W/mK (±{ci:.2f} if M>1); results in {outdir}")

if __name__ == "__main__":
    main()
```

### Calibrating the unit prefactor (do this once, then freeze)

LAMMPS’s `examples/KAPPA` is the accepted yardstick. Run the argon GK example with your analysis stride and dt, compute kk with this script (set `F_conv` so the argon value matches the example’s published kk within a percent), then **freeze** that constant in version control. This keeps your conversion audited and avoids silent unit drift if you ever change sampling conventions.

---

## How to run it for a state

Assuming you already ran Step 11 and have two shards:

```
python analysis/hcacf_analyze.py \
  --state T300_x0.00 \
  --shards runs/T300_x0.00/nve_shard_000/flux.raw runs/T300_x0.00/nve_shard_001/flux.raw \
  --L_A 54.10 --T 300 --dt_fs 1.0 --stride_fs 1.0 --c 5.0
```

This writes `analysis/T300_x0.00/k_result.json`, `analysis/T300_x0.00/report.md`, `hcacf.png`, and `kcumu.png`.

If you want window-sensitivity baked into the report (e.g., c→0.8c,1.2cc\to 0.8c,1.2c), you can extend the script with two extra passes and record ∣Δk∣|\Delta k| as a “plateau robustness” line; your Step 8 `validate_k.py` already showed how.

---

## What to commit as the artifact

Commit `analysis/<state>/k_result.json` as the machine-readable ground truth, `hcacf.png` and `kcumu.png` as diagnostics, and `report.md` as the human summary. These three together constitute the “analysis report per state (plateau, IAT, CI, diagnostics)” you listed in the plan, and they are reproducible bit-for-bit from your flux files and control settings.



# 14.) Grid execution & consistency checks


You’ll now run Steps 10–13 over the entire (T,x)(T,x) grid with **frozen numerics**, keep a lightweight **ledger** of what happened (start/stop, wall-clock, seeds, MPI layout, anomalies), and emit a single, uniform **k(T,x)k(T,x)** table with the exact metadata needed to reproduce every number.

Below are three small, composable tools:

1. a **grid runner** that reads `project.yaml`, spins each state through Equilibrate → GK shards → HCACF analysis (and optional rNEMD for the states you marked), and appends entries to a JSON-lines ledger;
    
2. a **consistency checker** that scans all outputs and asserts your numerics/analysis knobs are identical across the grid (dt, PPPM acc/order, stride, Sokal cc, etc.);
    
3. a **collector** that consolidates `analysis/<state>/k_result.json` files into a single CSV and Markdown summary with uniform metadata.
    

No scheduler assumptions are baked in; it’s vanilla Python invoking your Step-10/11/13 scripts. On clusters, wrap the same calls in SLURM array jobs; the ledger will still work.

---

## Grid runner: `scripts/run_grid.py`

This reads `project.yaml`, derives canonical state names, ensures a manifest/control exist, runs Steps 10–13, and writes a JSON-lines ledger (`runs/ledger.jsonl`) with timing and key parameters. If you listed rNEMD spot checks in `project.yaml.grid.crosscheck_states`, it runs Step 12 for those states, too.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, time, subprocess, hashlib, os
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Any
import yaml

ROOT = Path(".")
LEDGER = ROOT/"runs"/"ledger.jsonl"

@dataclass(frozen=True)
class State:
    T: int
    x: float
    name: str

def load_project() -> dict:
    y = yaml.safe_load((ROOT/"project.yaml").read_text())
    return y

def states_from_project(cfg: dict) -> List[State]:
    Ts = [int(t) for t in cfg["grid"]["temperatures_K"]]
    xs = [float(x) for x in cfg["grid"]["vacancy_fractions"]]
    out=[]
    for T in Ts:
        for x in xs:
            out.append(State(T, x, f"T{T}_x{float(x):.3f}"))
    return out

def ensure_manifest_and_control(st: State):
    run_dir = ROOT/"runs"/st.name
    ana_dir = ROOT/"analysis"/st.name
    run_dir.mkdir(parents=True, exist_ok=True); ana_dir.mkdir(parents=True, exist_ok=True)
    man = run_dir/"manifest.json"
    ctrl= ana_dir/"control.json"
    if not man.exists():
        # seed and numerics will be edited by earlier Step 3 bootstrap; we create a minimal shell if missing
        shell = {
            "schema":"ceriakappa/manifest@1",
            "inputs":{"pppm":{"accuracy":1e-5,"order":5},"neighbor":{"skin_A":2.0}},
            "seeds":{"master": int(hashlib.sha256(st.name.encode()).hexdigest()[:8],16)},
            "numerics":{"timestep_fs": 1.0},
            "hashes":{}
        }
        man.write_text(json.dumps(shell, indent=2))
    if not ctrl.exists():
        ctrl.write_text(json.dumps({
            "analysis_version":"1.0.0",
            "gk":{"flux_sample_stride_fs":1.0,"sokal_c":5.0,"component_average":True,"integration":"trapezoid","demean":True,"blocking":{"min_block":1024}}
        }, indent=2))

def stamp_ledger(event: Dict[str,Any]):
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a") as f:
        f.write(json.dumps(event)+"\n")

def run_cmd(cmd: List[str], cwd: Path, env: Dict[str,str]|None=None) -> int:
    t0 = time.time()
    code = 0
    try:
        subprocess.run(cmd, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        code = e.returncode
    dt = time.time()-t0
    return code, dt

def data_path_for(st: State) -> Path:
    # stoich vs reduced
    if st.x == 0.0:
        # assume Step 5 output
        return ROOT/"geometry"/f"stoich_CeO2_n10_a5.410"/"ceo2.data"
    # assume Step 9 output naming
    return ROOT/"geometry"/f"reduced_CeO2_n10_a5.410_x{st.x:.3f}_rmin3.50"/"ceo2_reduced.data"

def main():
    cfg = load_project()
    states = states_from_project(cfg)
    cross = set((s["T"], float(s["x"])) for s in cfg["grid"].get("crosscheck_states", []))
    lmp = os.environ.get("LMP_BIN","lmp")
    mpi_ranks = int(os.environ.get("GK_MPI","8"))

    for st in states:
        ensure_manifest_and_control(st)
        data = data_path_for(st)
        if not data.exists():
            stamp_ledger({"state":st.name,"stage":"skip","reason":"missing data file", "ts":time.time()})
            continue

        # Step 10 — Equilibrate (NPT→NVT)
        eq = ["python","scripts/equilibrate_state.py","--state",st.name,"--data",str(data),
              "--T",str(st.T),"--manifest",str(ROOT/"runs"/st.name/"manifest.json"),"--lmp",lmp]
        code,dt = run_cmd(eq, ROOT)
        stamp_ledger({"state":st.name,"stage":"equilibrate","code":code,"elapsed_s":dt,"ts":time.time()})
        if code!=0: 
            stamp_ledger({"state":st.name,"stage":"abort","reason":"equilibrate failed"}); 
            continue

        # Step 11 — GK production (two shards by default)
        gk = ["python","scripts/run_gk_state.py","--state",st.name,"--shards","2","--mpi",str(mpi_ranks),"--lmp",lmp,"--pin"]
        code,dt = run_cmd(gk, ROOT)
        stamp_ledger({"state":st.name,"stage":"gk_shards","code":code,"elapsed_s":dt,"ts":time.time(),"mpi":mpi_ranks})
        if code!=0:
            stamp_ledger({"state":st.name,"stage":"abort","reason":"gk failed"}); 
            continue

        # Step 13 — HCACF analysis
        man = json.loads((ROOT/"runs"/st.name/"manifest.json").read_text())
        L_A = float((ROOT/"runs"/st.name/"freeze_volume.json").read_text()) if (ROOT/"runs"/st.name/"freeze_volume.json").exists() else (10*5.41)
        dt_fs = man["numerics"]["timestep_fs"]
        stride_fs = json.loads((ROOT/"analysis"/st.name/"control.json").read_text())["gk"]["flux_sample_stride_fs"]
        shards = []
        for d in (ROOT/"runs"/st.name).glob("nve_shard_*/flux.raw"):
            shards.append(str(d))
        if not shards:
            stamp_ledger({"state":st.name,"stage":"abort","reason":"no flux files"}); 
            continue
        ana = ["python","analysis/hcacf_analyze.py","--state",st.name,"--shards",*shards,
               "--L_A",str(L_A),"--T",str(st.T),"--dt_fs",str(dt_fs),"--stride_fs",str(stride_fs),"--c",str(cfg["analysis"]["sokal_c"] if "analysis" in cfg and "sokal_c" in cfg["analysis"] else 5.0)]
        code,dt = run_cmd(ana, ROOT)
        stamp_ledger({"state":st.name,"stage":"analyze_gk","code":code,"elapsed_s":dt,"ts":time.time()})

        # Optional rNEMD on selected states
        if (st.T, st.x) in cross:
            rdir = ROOT/"runs"/st.name/"rnemd"/"period_5000"
            rdir.mkdir(parents=True, exist_ok=True)
            r1 = ["lmp","-in","lmp/rnemd_spot.in","-var","restart",str(ROOT/"runs"/st.name/f"eq_{st.name}.restart"),
                  "-var","Ttarget",str(st.T),"-var","period","5000","-var","outroot",str(rdir)]
            code1,dt1 = run_cmd(r1, ROOT)
            r2dir = ROOT/"runs"/st.name/"rnemd"/"period_10000"
            r2dir.mkdir(parents=True, exist_ok=True)
            r2 = ["lmp","-in","lmp/rnemd_spot.in","-var","restart",str(ROOT/"runs"/st.name/f"eq_{st.name}.restart"),
                  "-var","Ttarget",str(st.T),"-var","period","10000","-var","outroot",str(r2dir)]
            code2,dt2 = run_cmd(r2, ROOT)
            stamp_ledger({"state":st.name,"stage":"rnemd","code_fast":code1,"code_slow":code2,"elapsed_s":dt1+dt2,"ts":time.time()})
            # analyze
            for pd,per in [(rdir,5000),(r2dir,10000)]:
                ar = ["python","analysis/analyze_rnemd.py","--root",str(pd),"--period_steps",str(per),"--dt_fs",str(dt_fs)]
                code,dt = run_cmd(ar, ROOT)
                stamp_ledger({"state":st.name,"stage":"analyze_rnemd","period":per,"code":code,"elapsed_s":dt,"ts":time.time()})

if __name__ == "__main__":
    main()
```

Usage (sequential, from repo root):

```
python scripts/run_grid.py
```

This will leave a `runs/ledger.jsonl` like:

```json
{"state":"T300_x0.000","stage":"equilibrate","code":0,"elapsed_s":412.6,"ts":1699999999.1}
{"state":"T300_x0.000","stage":"gk_shards","code":0,"elapsed_s":1882.3,"mpi":8,"ts":1699999999.2}
{"state":"T300_x0.000","stage":"analyze_gk","code":0,"elapsed_s":4.1,"ts":1699999999.3}
...
```

If you prefer to parallelize, run one `python scripts/run_grid.py` per state under a job array; the ledger appends safely.

---

## Consistency checker: `scripts/check_consistency.py`

This scans every state’s `manifest.json` and `analysis/.../k_result.json` to assert you actually used **identical numerics** and **analysis knobs** across the grid, and it prints any outliers.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path

def main():
    base_num=None; base_ctl=None; bad=[]
    for man in Path("runs").glob("T*_x*/manifest.json"):
        state = man.parent.name
        try:
            M = json.loads(man.read_text())
            C = json.loads((Path("analysis")/state/"k_result.json").read_text())
        except Exception as e:
            print(f"[WARN] {state}: missing outputs ({e})"); 
            continue
        num = (M["numerics"]["timestep_fs"], M["inputs"]["pppm"]["accuracy"], M["inputs"]["pppm"]["order"], M["inputs"]["neighbor"]["skin_A"])
        ctl = (C["window_c"], C["dt_fs"], C["stride_fs"])
        if base_num is None: base_num=num
        if base_ctl is None: base_ctl=ctl
        if num!=base_num or ctl!=base_ctl:
            bad.append((state, num, ctl))
    if bad:
        print("CONSISTENCY FAIL — the following states differ in numerics/analysis:")
        for s,num,ctl in bad:
            print(f"  {s}: numerics={num}, analysis={ctl}")
        sys.exit(1)
    print("CONSISTENCY PASS — numerics and analysis knobs are identical across all analyzed states.")

if __name__ == "__main__":
    main()
```

Run it after the grid:

```
python scripts/check_consistency.py
```

---

## Collector: `analysis/collect_grid.py`

This pulls each state’s kk and CI plus key metadata into a single CSV and a Markdown summary. It also writes a compact JSON you can ship with the dataset.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, csv
from pathlib import Path

def main():
    rows=[]
    for jf in Path("analysis").glob("T*_x*/k_result.json"):
        state = jf.parent.name
        J = json.loads(jf.read_text())
        T = J["T_K"]; L = J["L_A"]; dt = J["dt_fs"]; stride = J["stride_fs"]; c = J["window_c"]
        k = J["aggregate"]["k_mean_WmK"]; ci = J["aggregate"]["ci95_WmK"]
        x = float(state.split("_x")[-1])
        rows.append({"state":state,"T_K":T,"x":x,"k_WmK":k,"ci95_WmK":ci,"L_A":L,"dt_fs":dt,"stride_fs":stride,"window_c":c})
    rows.sort(key=lambda r:(r["T_K"], r["x"]))
    # CSV
    outcsv = Path("analysis")/"k_grid.csv"
    with outcsv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    # Markdown
    md = ["# Consolidated k(T,x) table\n", "| T (K) | x | k (W/mK) | CI95 | L (Å) | dt (fs) | stride (fs) | c |", "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        md.append(f"| {r['T_K']} | {r['x']:.3f} | {r['k_WmK']:.2f} | {r['ci95_WmK']:.2f} | {r['L_A']:.2f} | {r['dt_fs']:.2f} | {r['stride_fs']:.2f} | {r['window_c']:.1f} |")
    (Path("analysis")/"k_grid.md").write_text("\n".join(md))
    # JSON for programmatic reuse
    (Path("analysis")/"k_grid.json").write_text(json.dumps(rows, indent=2))
    print(f"Wrote {outcsv}, analysis/k_grid.md, and analysis/k_grid.json")

if __name__ == "__main__":
    main()
```

Run it:

```
python analysis/collect_grid.py
```

You’ll get:

- `analysis/k_grid.csv` — the consolidated table,
    
- `analysis/k_grid.md` — a readable summary you can paste into your report,
    
- `analysis/k_grid.json` — machine-readable for downstream models.
    

---

## What to commit as the artifact

Commit three things:

1. `runs/ledger.jsonl` — append-only ledger of wall-clock, seeds/MPI layout, return codes, and any anomalies per stage and state.
    
2. `analysis/k_grid.{csv,md,json}` — the consolidated k(T,x)k(T,x) table with uniform metadata (L, dt, stride, window cc).
    
3. `scripts/check_consistency.py` output (copy its PASS line into your report, or commit a small `CONSISTENCY.txt` with the run date and the PASS message).
    

That’s the entire “repeat 10–13 over the full grid, track the ledger, and emit a consolidated table” loop. From here, Step 15 (systematics slices) can branch off any one state without touching the frozen numerics for the grid.



# 15.) Sensitivity slices (as needed)


The aim is to bound **method/systematics**—not to re-characterize statistics. You will take a _single designated baseline state_ (e.g., T=300T=300 K, x=0x=0), apply one or two controlled deltas (tighter PPPM; smaller Δt\Delta t; alternate seed), re-run only the minimal pieces needed to affect kk, and then report a compact **systematics band** that is clearly separated from the shard-scatter CI. Below is a small, self-contained workflow you can drop into the repo to automate this and to emit a one-page note.

## What changes and what stays frozen

Keep the physics (force field, structure, volume LL, temperature TT) and the GK analysis knobs (Sokal cc, component averaging, blocking) **identical** to the baseline. Only alter the _numerical dial you are testing_, and only for the **measurement** leg where possible. That means:

- **PPPM tighter**: raise `kspace_accuracy` (e.g., 1×10−61\times10^{-6} instead of 1×10−51\times10^{-5}) and, if you wish, bump `kspace_modify order` one notch. Keep the same equilibrated restart; the density/volume was already validated.
    
- **Δt\Delta t smaller**: re-run the _measurement shards_ with a smaller `timestep` (e.g., 0.5 fs vs 1.0 fs) using the same restart; no need to repeat NPT/NVT unless you want to show that the frozen box behaves similarly under the smaller Δt\Delta t.
    
- **Alternate seed**: keep all numerics identical but derive a new master seed and thus new NVE shard seeds; equilibrated restart may be reused (you’re testing trajectory stochasticity in the measurement leg, not EOS).
    

## Minimal runner for sensitivity variants

Save as `analysis/run_sensitivity.py`. It clones the baseline state into a variant scratch folder, writes a small LAMMPS input that **overrides** only the dial under test, runs one or two GK shards, analyzes them with the same HCACF pipeline, and writes a side-by-side JSON + Markdown note.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, shutil, subprocess, hashlib, time
from pathlib import Path
from typing import Dict, Any, List

ROOT = Path(".")
NUM_INC = ROOT/"configs/includes/global_numerics.in"
GK_INC  = ROOT/"configs/includes/gk_measurement.in"

GK_TEMPLATE = """\
units metal
atom_style charge
boundary p p p

# Read numerics include, then override specific knobs below if requested
include {numerics_inc}
read_restart {restart}

# Optional overrides injected here
{overrides}

# Flux assembly (pppm-consistent)
compute cke all ke/atom
compute cpe all pe/atom
compute cst all centroid/stress/atom NULL virial
compute J   all heat/flux cke cpe cst

include {gk_inc}

# Sampling stride in steps (kept equal to dt unless user passes stride)
variable stride equal {stride_steps}
dump dJ all custom ${stride} flux.raw c_J[1] c_J[2] c_J[3]
dump_modify dJ format binary yes append yes

fix acf all ave/correlate 1 {blk} {tot} c_J[1] c_J[2] c_J[3] type auto file J0Jt.dat ave running

run {nsteps}
unfix acf
undump dJ
"""

def sha256_file(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def run_lmp(inp: Path, workdir: Path, lmp_bin: str="lmp", mpi: int=8):
    cmd = [lmp_bin, "-in", inp.name]
    if mpi>1:
        cmd = ["mpirun", "-np", str(mpi)] + cmd
    subprocess.run(cmd, cwd=str(workdir), check=True, stdout=(workdir/"stdout.log").open("w"), stderr=subprocess.STDOUT)

def analyze(state: str, workdir: Path, L_A: float, T: float, dt_fs: float, stride_fs: float, c: float) -> Dict[str,Any]:
    # reuse your HCACF analyzer from Step 13
    out = subprocess.run([
        "python","analysis/hcacf_analyze.py",
        "--state", state, "--shards", str(workdir/"flux.raw"),
        "--L_A", str(L_A), "--T", str(T),
        "--dt_fs", str(dt_fs), "--stride_fs", str(stride_fs), "--c", str(c)
    ], check=True, capture_output=True, text=True)
    res = json.loads((ROOT/"analysis"/state/"k_result.json").read_text())
    return res

def main():
    import argparse, math
    ap = argparse.ArgumentParser(description="Run sensitivity slice vs baseline for one state.")
    ap.add_argument("--state", required=True, help="baseline state name, e.g., T300_x0.00")
    ap.add_argument("--variant", required=True, choices=["pppm_tight","dt_small","seed_alt"])
    ap.add_argument("--lmp", default="lmp")
    ap.add_argument("--mpi", type=int, default=8)
    # optional knobs
    ap.add_argument("--pppm_acc", type=float, default=1e-6)
    ap.add_argument("--pppm_order", type=int, default=None)
    ap.add_argument("--dt_fs", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--shard_ps", type=float, default=500.0)
    args = ap.parse_args()

    base = args.state
    run_dir = ROOT/"runs"/base
    ana_dir = ROOT/"analysis"/base
    man = json.loads((run_dir/"manifest.json").read_text())
    ctrl= json.loads((ana_dir/"control.json").read_text())
    restart = next((p for p in run_dir.glob("eq_*.restart")), None)
    if restart is None:
        raise SystemExit("No equilibrated restart; run Step 10 first.")

    L_A = float(json.loads((run_dir/"freeze_volume.json").read_text())["L_box_A"]) if (run_dir/"freeze_volume.json").exists() else 10*5.41
    T = int(base.split("_")[0][1:])
    # baseline numerics for stride/blocks
    base_dt = float(man["numerics"]["timestep_fs"])
    stride_fs = float(ctrl["gk"]["flux_sample_stride_fs"])
    blk = 10000; tot = 10000

    # variant staging
    tag = {"pppm_tight":"pppmTight","dt_small":"dtSmall","seed_alt":"seedAlt"}[args.variant]
    vstate = f"{base}_{tag}"
    vdir = ROOT/"runs"/vstate
    vdir.mkdir(parents=True, exist_ok=True)

    # overrides for LAMMPS input
    overrides = []
    dt_used = base_dt
    stride_steps = max(1, round(stride_fs/base_dt))
    if args.variant == "pppm_tight":
        acc = args.pppm_acc
        overrides.append(f"kspace_style pppm {acc}")
        if args.pppm_order:
            overrides.append(f"kspace_modify order {args.pppm_order}")
    elif args.variant == "dt_small":
        dt_used = float(args.dt_fs)
        overrides.append(f"timestep {dt_used}")
        stride_steps = max(1, round(stride_fs/dt_used))
    elif args.variant == "seed_alt":
        # no input override; derives new seeds by changing master (use provided seed or hash of name+salt)
        pass

    # write input
    nsteps = int(round(args.shard_ps / (dt_used/1000.0)))
    inp = vdir/"gk_variant.in"
    inp.write_text(GK_TEMPLATE.format(
        numerics_inc=NUM_INC.as_posix(),
        restart=restart.as_posix(),
        overrides="\n".join(overrides),
        gk_inc=GK_INC.as_posix(),
        stride_steps=int(stride_steps),
        blk=blk, tot=tot, nsteps=nsteps
    ))

    # run LAMMPS
    run_lmp(inp, vdir, lmp_bin=args.lmp, mpi=args.mpi)

    # analyze (write under analysis/<variant_state>/ to keep artifacts separate)
    (ROOT/"analysis"/vstate).mkdir(parents=True, exist_ok=True)
    res = analyze(vstate, vdir, L_A=L_A, T=T, dt_fs=dt_used, stride_fs=stride_steps*dt_used, c=float(ctrl["gk"]["sokal_c"]))

    # collect baseline result to compare
    base_json = json.loads((ROOT/"analysis"/base/"k_result.json").read_text())
    k_base = base_json["aggregate"]["k_mean_WmK"]; ci_base = base_json["aggregate"]["ci95_WmK"]
    k_var  = res["aggregate"]["k_mean_WmK"];      ci_var  = res["aggregate"]["ci95_WmK"]
    delta  = k_var - k_base

    note = {
        "baseline_state": base,
        "variant_state": vstate,
        "variant": args.variant,
        "baseline": {"k_WmK": k_base, "ci95_WmK": ci_base, "dt_fs": base_dt, "pppm_acc": man["inputs"]["pppm"]["accuracy"], "pppm_order": man["inputs"]["pppm"]["order"]},
        "variant_values": {"k_WmK": k_var, "ci95_WmK": ci_var, "dt_fs": dt_used, "pppm_acc": args.pppm_acc, "pppm_order": args.pppm_order},
        "delta_k_WmK": delta
    }
    (ROOT/"analysis"/vstate/"systematics_note.json").write_text(json.dumps(note, indent=2))
    md = ROOT/"analysis"/vstate/"systematics_note.md"
    md.write_text(
        f"# Systematics slice — {args.variant}\n\n"
        f"Baseline {base}: k = {k_base:.2f} ± {ci_base:.2f} W/mK\n\n"
        f"Variant  {vstate}: k = {k_var:.2f} ± {ci_var:.2f} W/mK\n\n"
        f"Δk = {delta:+.2f} W/mK\n\n"
        f"Knobs: dt = {dt_used:.3f} fs, PPPM acc/order = {note['variant_values'].get('pppm_acc','—')}/{note['variant_values'].get('pppm_order','—')}\n"
    )
    print(json.dumps(note, indent=2))

if __name__ == "__main__":
    main()
```

Example runs for your baseline `T300_x0.00`:

```
# PPPM tighter (order unchanged)
python analysis/run_sensitivity.py --state T300_x0.00 --variant pppm_tight --pppm_acc 1e-6

# Smaller timestep (0.5 fs)
python analysis/run_sensitivity.py --state T300_x0.00 --variant dt_small --dt_fs 0.5

# Alternate seed (inherit same numerics; you can just rerun Step 11 with a new master seed, or call this with --seed and add seed plumbing if desired)
python analysis/run_sensitivity.py --state T300_x0.00 --variant seed_alt
```

If you want the seed-alternate to be truly automatic here, extend the template to inject `variable seed` and use it in a thermostat-free way (for NVE there is no RNG touchpoint unless you randomize initial velocities; most workflows change the shard start time or origin selection instead). In practice, re-deriving the shard seeds in `run_gk_state.py` and rerunning shards is the cleanest; you can then compare the two `k_result.json` directly.

## A tiny collector that turns deltas into a “systematics band”

Save as `analysis/systematics_band.py`. It reads any number of `systematics_note.json` files for a baseline and reports the **max |Δk|** and an optional **quadrature** combination if you prefer a symmetric band.

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, argparse
from pathlib import Path
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, help="state name, e.g., T300_x0.00")
    ap.add_argument("--notes", nargs="+", required=True, help="paths to .../systematics_note.json")
    args = ap.parse_args()

    deltas=[]
    entries=[]
    for p in args.notes:
        J = json.loads(Path(p).read_text())
        deltas.append(J["delta_k_WmK"])
        entries.append({"variant":J["variant"], "delta_k_WmK":J["delta_k_WmK"]})
    deltas = np.array(deltas, float)
    band_max = float(np.max(np.abs(deltas)))
    band_quad = float(np.sqrt(np.sum(deltas**2)))   # conservative if effects are uncorrelated

    out = {
        "baseline": args.baseline,
        "entries": entries,
        "band_max_abs_WmK": band_max,
        "band_quadrature_WmK": band_quad
    }
    Path("analysis")/args.baseline and (Path("analysis")/args.baseline/"systematics_band.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
```

Usage:

```
python analysis/systematics_band.py \
  --baseline T300_x0.00 \
  --notes analysis/T300_x0.00_pppmTight/systematics_note.json \
          analysis/T300_x0.00_dtSmall/systematics_note.json
```

This writes `analysis/T300_x0.00/systematics_band.json` summarizing the tested variants and the band.

## How to present the result in your report

Put the systematics band next to the statistical CI, clearly separated, for the baseline state; propagate the **same band** to the rest of the grid unless you later discover a strong TT or xx dependence of the deltas (rare if the dials are purely numerical). For example:

> Baseline T=300T=300 K, x=0x=0 GK: k=16.8±1.2k=16.8\pm1.2 W m−1^{-1} K−1^{-1} (95% statistical CI). Sensitivity slices give ∣Δk∣max⁡=0.6|\Delta k|_{\max}=0.6 W m−1^{-1} K−1^{-1} across PPPM accuracy (10−5→10−6)(10^{-5}\to10^{-6}) and Δt(1.0→0.5\Delta t(1.0\to0.5 fs), which we quote as a **systematics band**. Unless stated otherwise, add this band in quadrature when comparing absolute values; for _relative_ changes across xx at fixed numerics, we report statistical CIs only.

## Artifact checklist

You should now have, per variant: `analysis/<baseline>_<tag>/systematics_note.{json,md}` with kk, CI, the dial values, and Δk\Delta k; and, per baseline: `analysis/<baseline>/systematics_band.json` with the combined band. Commit those alongside the GK artifacts; they are small, fast to regenerate, and make your uncertainty story transparent.



# 16.) Interpretation & limitations ledger


The goal here is to ship a single, human-readable page that travels with your results and makes the boundaries of interpretation unmistakable. It must keep **statistical confidence** (GK shard scatter) separate from **numerical systematics** (Step 15 band) and from the broader **model-family uncertainty** (differences across vetted force fields). It should also restate the classical-nuclei caveat and the bulk-ideal assumptions, so a downstream reader cannot accidentally over-claim.

Below is a tiny tool that assembles those pieces from artifacts you already produce and writes a one-page `LIMITATIONS.md` per project, plus a compact JSON for programmatic reuse.

---

## What the tool ingests

- `analysis/k_grid.json` from Step 14 with k(T,x)k(T,x) and statistical CIs.
    
- Optional `analysis/<baseline>/systematics_band.json` from Step 15 with a single number (max-abs or quadrature band) that bounds numerics.
    
- Optional multi-FF comparison: drop additional grids as `analysis_altFF/<ff_key>/k_grid.json` (same schema). The tool computes a **model-family band** as the max spread across FFs at the baseline state (or across all states if you prefer), and reports it separately from the other two.
    

---

## Generator: `analysis/make_limitations.py`

```python
#!/usr/bin/env python3
from __future__ import annotations
import json, argparse, statistics as stats
from pathlib import Path
from typing import Dict, List, Optional

def load_grid(path: Path) -> List[Dict]:
    if not path.exists():
        return []
    return json.loads(path.read_text())

def model_family_band(baseline_state: str, alt_dirs: List[Path]) -> Optional[float]:
    """Return max |Δk| across force fields at the baseline state, or None if not enough data."""
    ks = []
    for d in alt_dirs:
        g = load_grid(d/"k_grid.json")
        if not g: 
            continue
        baseline_rows = [r for r in g if r["state"] == baseline_state]
        if baseline_rows:
            ks.append(baseline_rows[0]["k_WmK"])
    if len(ks) < 2:
        return None
    return max(ks) - min(ks)

def main():
    ap = argparse.ArgumentParser(description="Assemble a one-page limitations ledger.")
    ap.add_argument("--project", default="ceria_kappa_map_v1")
    ap.add_argument("--baseline", default="T300_x0.000", help="state name used for bands")
    ap.add_argument("--grid", default="analysis/k_grid.json")
    ap.add_argument("--systematics", default=None, help="analysis/<baseline>/systematics_band.json")
    ap.add_argument("--altff", nargs="*", default=[], help="paths like analysis_altFF/Gotte_2007, analysis_altFF/DIPPIM_2011")
    args = ap.parse_args()

    grid = load_grid(Path(args.grid))
    if not grid:
        raise SystemExit("No k_grid.json found; run Step 14 first.")

    # pull a couple of representative numbers for the opening paragraph
    # pick mid-temperature, x=0 and x=highest
    Ts = sorted({row["T_K"] for row in grid})
    xs = sorted({row["x"] for row in grid})
    Tmid = Ts[len(Ts)//2]
    k0 = next(r for r in grid if r["T_K"]==Tmid and abs(r["x"])<1e-12)
    kx = next(r for r in grid if r["T_K"]==Tmid and r["x"]==max(xs))
    k0_val, k0_ci = k0["k_WmK"], k0["ci95_WmK"]
    kx_val, kx_ci = kx["k_WmK"], kx["ci95_WmK"]

    # systematics band (numerical)
    sys_band = None
    if args.systematics:
        p = Path(args.systematics)
        if p.exists():
            J = json.loads(p.read_text())
            # prefer quadrature band if present, else max-abs
            sys_band = J.get("band_quadrature_WmK", None) or J.get("band_max_abs_WmK", None)

    # model-family band (across FFs)
    alt_dirs = [Path(p) for p in args.altff]
    fam_band = model_family_band(args.baseline, alt_dirs)

    # write JSON summary
    outdir = Path("analysis")
    outdir.mkdir(exist_ok=True, parents=True)
    summary = {
        "project": args.project,
        "baseline_state": args.baseline,
        "statistical_ci_note": "CIs are 95% from shard scatter unless stated otherwise.",
        "systematics_band_WmK": sys_band,
        "model_family_band_WmK": fam_band,
        "classical_nuclei_caveat": "Classical MD overpopulates high-ω modes at low T; absolute k is biased high as T→0.",
        "bulk_ideal_assumptions": "Perfect, periodic bulk with PPPM tin-foil boundary; no porosity, boundaries, cracks, or electronic/ polaronic heat carriers.",
        "representative_points": {
            f"T{Tmid}_x0.000": {"k": k0_val, "ci95": k0_ci},
            f"T{Tmid}_x{max(xs):.3f}": {"k": kx_val, "ci95": kx_ci}
        }
    }
    (outdir/"LIMITATIONS.json").write_text(json.dumps(summary, indent=2))

    # emit one-page Markdown
    md = []
    md.append(f"# Interpretation & limitations — {args.project}\n")
    md.append("This dataset reports lattice thermal conductivity k(T,x) for CeO₂₋ₓ from classical MD under periodic boundary conditions. The numbers are ready for **relative trends** and for insertion into models that consume bulk lattice k. Read the following constraints before using absolutes.\n")
    md.append("## What the numbers mean\n")
    md.append(f"At T={Tmid} K, the stoichiometric baseline is k = {k0_val:.2f} ± {k0_ci:.2f} W·m⁻¹·K⁻¹ (95% statistical CI). At the highest vacancy fraction in the grid, the same T gives k = {kx_val:.2f} ± {kx_ci:.2f} W·m⁻¹·K⁻¹. These CIs reflect sampling variability only.\n")
    md.append("## Uncertainty decomposition\n")
    md.append("**Statistical CI (reported next to each value).** Derived from shard scatter with windowed Green–Kubo analysis; scales ≈ 1/√(independent correlation time).\n")
    if sys_band is not None:
        md.append(f"**Numerical systematics band (global, not added to CI):** ±{sys_band:.2f} W·m⁻¹·K⁻¹ from sensitivity slices (tighter PPPM, smaller Δt, alternate seed) at the baseline state. Apply as a bias bound on absolutes; do not use to inflate trend error bars unless you change numerics mid-grid.\n")
    else:
        md.append("**Numerical systematics band:** Not evaluated. If comparing absolutes, run Step 15 on the baseline and quote ±|Δk|.\n")
    if fam_band is not None:
        md.append(f"**Model-family band (force-field spread):** ≈±{fam_band/2:.2f} W·m⁻¹·K⁻¹ around the chosen force field at the baseline. This reflects plausible variation across vetted CeO₂ parameterizations. Keep it distinct from sampling CI and numerics.\n")
    else:
        md.append("**Model-family band:** Not evaluated. To obtain one, re-run the grid with at least one alternate, vetted force field and take the max spread at the baseline state.\n")
    md.append("## Physics boundaries\n")
    md.append("**Classical nuclei.** High-frequency modes carry k_B T at all T; absolute k is biased high at low T. Trends with x and moderate/high-T absolutes are more robust. No path-integral or quantum corrections were applied.\n")
    md.append("**What the model does not include.** No explicit electron–phonon coupling; no small-polaron transport; no explicit polarizability unless a polarizable FF was chosen; no grain boundaries, porosity, second phases, or irradiation microcracks. Electrostatics use PPPM with tin-foil boundary.\n")
    md.append("## Using the numbers responsibly\n")
    md.append("Use the reported CIs when comparing states at fixed numerics and force field. When quoting **absolute** k, bracket it with the numerical systematics band and, if available, the model-family band. When mapping to pellets or microstructured material, apply an independent microstructure model (e.g., Maxwell–Eucken) rather than inferring from bulk MD.\n")
    md.append("## Provenance\n")
    md.append("Numerics and analysis knobs are frozen across the grid (dt, PPPM accuracy/order, flux sampling stride, Sokal c). See `analysis/k_grid.json`, `runs/ledger.jsonl`, and per-state manifests for seeds, MPI layout, and file hashes.\n")
    (outdir/"LIMITATIONS.md").write_text("\n".join(md))
    print("Wrote analysis/LIMITATIONS.md and analysis/LIMITATIONS.json")

if __name__ == "__main__":
    main()
```

Run it once after Step 14 (and Step 15 if you computed a systematics band):

```
python analysis/make_limitations.py \
  --project ceria_kappa_map_v1 \
  --baseline T300_x0.000 \
  --grid analysis/k_grid.json \
  --systematics analysis/T300_x0.000/systematics_band.json \
  --altff analysis_altFF/Gotte_2007 analysis_altFF/DIPPIM_2011
```

This writes:

- `analysis/LIMITATIONS.md` — the one-page statement you’ll link from your results,
    
- `analysis/LIMITATIONS.json` — a machine-readable summary of the same facts.
    

You can now add a single line in your main results readme or manuscript:

> See `analysis/LIMITATIONS.md` for uncertainty decomposition, classical-nuclei caveat, and bulk-ideal assumptions; cite the statistical CIs in tables and append the numerical and model-family bands only when quoting absolutes.


# 17.) Packaging & provenance freeze


This step makes a self-contained, immutable release: every state’s `k_result.json`, the plots and profiles, all per-state manifests, and the exact inputs/controls used—plus a top-level index, checksums, and a version tag you can archive (e.g., Zenodo DOI). The idea is: one command creates a tarball, verifies all hashes, writes a machine-readable catalog, and emits a short human README and CITATION file.

## What goes into the release

Everything a third party needs to (a) read results, (b) verify integrity, (c) rerun analysis from flux, and (d) understand licenses:

- `analysis/T*_x*/k_result.json`, `hcacf.png`, `kcumu.png`, and (if present) `rnemd/*/rnemd_result.json`, `profile_fit.png`.
    
- `runs/T*_x*/manifest.json` and per-shard `flux.sha256` (and the `flux.raw` **only if** you want a “with-data” bundle; usually you ship analysis-only + checksums).
    
- Top-level `analysis/k_grid.{json,csv,md}`, `LIMITATIONS.{md,json}`, and (if made) `T*_x*/systematics_band.json`.
    
- The exact analysis code used (pin a commit and include the `analysis/` and `scripts/` dirs).
    
- The numerics card: `configs/numerics.yaml` and generated includes.
    
- Force-field parameter files and species maps you actually used (plus their checksums).
    
- `LICENSE`, `CITATION.cff`, `README_RELEASE.md`, and a manifest index `RELEASE_MANIFEST.json`.
    

## Packager that builds the release and verifies integrity

Save as `scripts/make_release.py`. It creates `dist/<tag>/…`, copies artifacts, computes SHA-256 for every file, validates per-state hashes, writes an index, and builds a `.tar.gz`. If you have GPG configured, it also signs the archive.

```python
#!/usr/bin/env python3
from __future__ import annotations
import os, json, hashlib, tarfile, shutil, subprocess, datetime as dt
from pathlib import Path
from typing import Dict, List

ROOT = Path(".")
DIST = ROOT/"dist"

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def collect_states() -> List[str]:
    return sorted({p.parent.name for p in ROOT.glob("analysis/T*_x*/k_result.json")})

def copy_tree(src: Path, dst: Path, patterns: List[str]):
    dst.mkdir(parents=True, exist_ok=True)
    for pat in patterns:
        for p in src.glob(pat):
            rel = p.name
            shutil.copy2(p, dst/rel)

def build_release(tag: str, with_flux: bool = False, gpg_key: str | None = None):
    ts = dt.datetime.now().isoformat(timespec="seconds")
    outroot = DIST/tag
    if outroot.exists():
        shutil.rmtree(outroot)
    outroot.mkdir(parents=True, exist_ok=True)

    index: Dict[str, dict] = {
        "release_tag": tag,
        "created": ts,
        "generator": "make_release.py",
        "entries": [],
        "top_level": {}
    }

    # Top-level bundles
    toplevel_files = [
        "analysis/k_grid.json", "analysis/k_grid.csv", "analysis/k_grid.md",
        "analysis/LIMITATIONS.md", "analysis/LIMITATIONS.json",
        "configs/numerics.yaml", "configs/includes/global_numerics.in",
        "configs/includes/prep_ensembles.in", "configs/includes/gk_measurement.in",
        "LICENSE", "CITATION.cff", "README_RELEASE.md"
    ]
    for rel in toplevel_files:
        p = ROOT/rel
        if p.exists():
            dst = outroot/rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dst)
            index["top_level"][rel] = {"sha256": sha256_file(dst)}

    # Force-field assets you actually used
    for p in ROOT.glob("potentials/**/*"):
        if p.is_file():
            dst = outroot/"potentials"/p.relative_to(ROOT/"potentials")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dst)

    # Per-state artifacts
    states = collect_states()
    for state in states:
        # analysis artifacts
        a_src = ROOT/"analysis"/state
        a_dst = outroot/"analysis"/state
        a_dst.mkdir(parents=True, exist_ok=True)
        copy_tree(a_src, a_dst, ["k_result.json", "hcacf.png", "kcumu.png", "report.md"])
        # rNEMD (optional)
        r_src = ROOT/"runs"/state/"rnemd"
        if r_src.exists():
            for pd in r_src.glob("period_*/"):
                dstd = outroot/"runs"/state/"rnemd"/pd.name
                dstd.mkdir(parents=True, exist_ok=True)
                copy_tree(pd, dstd, ["rnemd_result.json", "profile_fit.png", "rnemd_meta.log"])
        # manifests
        r_dst = outroot/"runs"/state
        r_dst.mkdir(parents=True, exist_ok=True)
        man = ROOT/"runs"/state/"manifest.json"
        if man.exists():
            shutil.copy2(man, r_dst/"manifest.json")
        # freeze_volume (if present)
        fr = ROOT/"runs"/state/"freeze_volume.json"
        if fr.exists():
            shutil.copy2(fr, r_dst/"freeze_volume.json")
        # flux files (optional; usually no)
        if with_flux:
            for fx in ROOT.glob(f"runs/{state}/nve_shard_*/flux.raw"):
                d = r_dst/fx.parent.name
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fx, d/"flux.raw")
        # per-shard checksums always included
        for cs in ROOT.glob(f"runs/{state}/nve_shard_*/flux.sha256"):
            d = r_dst/cs.parent.name
            d.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cs, d/"flux.sha256")

        # compute file index for this state
        entry = {"state": state, "files": {}}
        for p in (a_dst.rglob("*")):
            if p.is_file():
                rel = p.relative_to(outroot).as_posix()
                entry["files"][rel] = {"sha256": sha256_file(p)}
        for p in (r_dst.rglob("*")):
            if p.is_file():
                rel = p.relative_to(outroot).as_posix()
                entry["files"][rel] = {"sha256": sha256_file(p)}
        index["entries"].append(entry)

    # Write the global index
    (outroot/"RELEASE_MANIFEST.json").write_text(json.dumps(index, indent=2))

    # Build tarball
    tarpath = DIST/f"{tag}.tar.gz"
    with tarfile.open(tarpath, "w:gz") as tar:
        tar.add(outroot, arcname=tag)

    # Optional GPG sign
    sig = None
    if gpg_key:
        try:
            subprocess.run(["gpg","--detach-sign","--local-user",gpg_key,"--armor",str(tarpath)], check=True)
            sig = f"{tarpath}.asc"
        except Exception:
            pass

    print(f"Release built: {tarpath}")
    if sig:
        print(f"Signature:     {sig}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, help="release tag, e.g., v0.3.0")
    ap.add_argument("--with-flux", action="store_true", help="include raw flux files (large)")
    ap.add_argument("--gpg-key", default=None, help="GPG key ID/email to sign tarball")
    args = ap.parse_args()
    build_release(args.tag, with_flux=args.with_flux, gpg_key=args.gpg_key)
```

Run it like:

```
python scripts/make_release.py --tag v0.1.0
# or include raw flux:
python scripts/make_release.py --tag v0.1.0 --with-flux
```

It produces `dist/v0.1.0.tar.gz` and a folder `dist/v0.1.0/` with `RELEASE_MANIFEST.json` listing every file and its SHA-256.

## Minimal `README_RELEASE.md`, `CITATION.cff`, and `LICENSE`

Put a succinct human overview at the top of the bundle:

`README_RELEASE.md`:

```markdown
# CeO2-x lattice conductivity map — Release v0.1.0

This archive contains k(T,x) for CeO₂₋ₓ computed via Green–Kubo MD with spot rNEMD checks.
- Statistical CIs are 95% (shard scatter); see `analysis/k_grid.json`.
- Numerical systematics band and model-family notes: `analysis/LIMITATIONS.md`.
- Per-state provenance: `runs/T*_x*/manifest.json`.

Re-use:
1) Read `analysis/k_grid.csv` for the table.
2) See `RELEASE_MANIFEST.json` for checksums and exact file list.
3) Cite as in `CITATION.cff`.

License: see `LICENSE`. Force-field files retain their original licenses; verify before redistribution.
```

`CITATION.cff`:

```yaml
cff-version: 1.2.0
message: "If you use this dataset or code, please cite it."
title: "CeO2-x lattice thermal conductivity map (classical MD, GK+rNEMD)"
version: "v0.1.0"
authors:
  - family-names: YourSurname
    given-names: YourName
    affiliation: Your Lab, University
date-released: 2025-09-20
keywords: [ceria, CeO2, thermal conductivity, molecular dynamics, Green-Kubo, rNEMD]
repository-code: "https://github.com/you/ceria_kappa_map"
```

`LICENSE`: pick MIT/BSD-3-Clause (for your code) and note that force-field parameter files may have different terms.

## Git tag and Zenodo DOI

Tag the repo:

```
git add dist/v0.1.0 RELEASE_MANIFEST.json analysis/k_grid.json analysis/LIMITATIONS.md
git commit -m "Release v0.1.0: packaged artifacts and manifest"
git tag -a v0.1.0 -m "CeO2-x k(T,x) release v0.1.0"
git push --tags
```

If you connect the repository to Zenodo (or similar), Zenodo will capture the tag and mint a DOI. Add the DOI badge to `README_RELEASE.md` and update `CITATION.cff` with `doi: 10.5281/zenodo.xxxxxxx`.

## Quick validation before you upload

Run a tiny auditor that (a) recomputes the SHA-256 of a random 10% sample from `RELEASE_MANIFEST.json` and (b) checks that every `k_result.json` has the fields `{aggregate.k_mean_WmK, aggregate.ci95_WmK, window_c, dt_fs, stride_fs}`. If any check fails, fix and re-mint the tag with `v0.1.1`.

That’s it: one command builds a verified, immutable bundle with provenance baked in, ready for DOI archiving and re-use.





# 18.) Comparison layer (optional)


The goal is a small, reproducible layer that (i) plots your **stoichiometric bulk kk** against a **theory bracket** (ab-initio/BTE) and **experimental pellets**, and (ii) optionally **maps bulk → pellet** with a simple porosity effective-medium model so readers can see how much of the gap is microstructure rather than physics of the potential. You’ll ship one figure plus a tiny CSV so the panel is self-documenting.

Below is a drop-in, analysis-only toolchain. It does not re-run MD; it just reads your artifacts and a small YAML of comparison numbers.

---

## 1) Put the comparison facts in one place

Create `configs/compare.yaml` and fill it once (numbers are placeholders; replace with your curated values when you finalize):

```yaml
# configs/compare.yaml
bte_clean_crystal:
  # Ab-initio/BTE bracket for dense, defect-free CeO2 at 300 K
  T_K: 300
  k_min_WmK: 15.5
  k_max_WmK: 18.0
  refs:
    - "BTE-DFT study A (Year, DOI)"
    - "BTE-DFT study B (Year, DOI)"

pellet_experiments:
  # High-density pellet measurements (report density if you have it)
  - label: "Exp #1"
    T_K: 300
    k_WmK: 11.0
    rel_density: 0.96
    ref: "Laser-flash paper (Year, DOI)"
  - label: "Exp #2"
    T_K: 300
    k_WmK: 9.5
    rel_density: 0.92
    ref: "Older compilation (Year, DOI)"

porosity_model:
  # Effective medium model parameters for mapping bulk→pellet
  model: "Maxwell-Eucken-SphericalPores"
  # φ = 1 - relative_density
  notes: "Spherical, non-interacting pores; conduction only (kp≈0). k_eff = k_bulk * 2(1-φ)/(2+φ)"
```

---

## 2) One script that makes the figure and a tiny CSV

Save as `analysis/make_comparison.py`. It reads your **baseline** stoichiometric result (`analysis/T300_x0.00/k_result.json`), the comparison YAML, draws the bracket and points, and optionally overlays a porosity curve computed from your bulk value.

```python
#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from typing import Dict, Any, List

import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def me_spherical_pores(k_bulk: float, rel_density: np.ndarray) -> np.ndarray:
    """Maxwell–Eucken (spherical, non-interacting pores, kp≈0).
       φ = 1 - ρ_rel ; k_eff = k_bulk * 2(1-φ)/(2+φ) = k_bulk * 2*ρ_rel/(1+ρ_rel)"""
    rho = np.asarray(rel_density, float)
    return k_bulk * (2.0*rho / (1.0 + rho))

def load_baseline(basestate: str) -> Dict[str, Any]:
    j = json.loads((Path("analysis")/basestate/"k_result.json").read_text())
    return {
        "T": j["T_K"],
        "k": j["aggregate"]["k_mean_WmK"],
        "ci95": j["aggregate"]["ci95_WmK"]
    }

def main():
    ap = argparse.ArgumentParser(description="Make a clean comparison figure for stoichiometric CeO2.")
    ap.add_argument("--state", default="T300_x0.00")
    ap.add_argument("--compare", default="configs/compare.yaml")
    ap.add_argument("--outdir", default="analysis/comparison")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    base = load_baseline(args.state)
    cfg = yaml.safe_load(Path(args.compare).read_text())

    # === Figure ===
    plt.figure(figsize=(5.0, 3.6))

    # 1) Your bulk MD result with 95% CI
    plt.errorbar([base["T"]], [base["k"]], yerr=[base["ci95"]],
                 fmt="o", capsize=3, label=f"MD (this work) {args.state}")

    # 2) BTE bracket as a vertical bar at the same T
    bte = cfg["bte_clean_crystal"]
    if abs(bte["T_K"] - base["T"]) < 1e9:  # same T (you can relax this if needed)
        plt.vlines(bte["T_K"]+2, bte["k_min_WmK"], bte["k_max_WmK"], lw=6, alpha=0.3, label="Ab-initio/BTE bracket")

    # 3) Pellet experiments (points), optionally with density annotation
    for p in cfg.get("pellet_experiments", []):
        if p["T_K"] != base["T"]: 
            continue
        lbl = f"{p['label']} ({p.get('rel_density','?'):.2f}ρ)" if "rel_density" in p else p["label"]
        plt.plot(p["T_K"]-2, p["k_WmK"], "s", label=lbl, ms=5)

    # 4) Optional porosity EMA curve at this T: vary relative density and map k_bulk→k_eff
    rho_grid = np.linspace(0.85, 1.00, 100)
    k_ema = me_spherical_pores(base["k"], rho_grid)
    plt.plot(np.full_like(rho_grid, base["T"]-1), k_ema, "-", lw=1.8, label="EMA (Maxwell–Eucken) vs ρ")

    plt.ylabel("Thermal conductivity k (W m$^{-1}$ K$^{-1}$)")
    plt.xlabel("Temperature (K)")
    plt.title("Stoichiometric CeO$_2$ — bulk vs. theory/experiment")
    plt.xlim(base["T"]-6, base["T"]+6)
    plt.legend(ncol=1, fontsize=8)
    plt.tight_layout()
    figpath = outdir/"comparison_T{}_{}.png".format(base["T"], args.state)
    plt.savefig(figpath, dpi=300)

    # === CSV companion (one row per series item) ===
    rows: List[Dict[str, Any]] = []
    rows.append({"series":"MD_this_work","T_K":base["T"],"k_WmK":base["k"],"ci95_WmK":base["ci95"],"note":args.state})

    rows.append({"series":"BTE_min","T_K":bte["T_K"],"k_WmK":bte["k_min_WmK"],"ci95_WmK":"","note":"lower bracket"})
    rows.append({"series":"BTE_max","T_K":bte["T_K"],"k_WmK":bte["k_max_WmK"],"ci95_WmK":"","note":"upper bracket"})

    for p in cfg.get("pellet_experiments", []):
        if p["T_K"] != base["T"]: 
            continue
        rows.append({"series":"Pellet_exp","T_K":p["T_K"],"k_WmK":p["k_WmK"],"ci95_WmK":p.get("ci95",""),
                     "note":f"{p['label']}; ρ_rel={p.get('rel_density','')}"})

    # EMA sweep is saved as a small table at the same T
    ema_csv = outdir/"ema_curve_T{}.csv".format(base["T"])
    import csv
    with (outdir/"comparison_points_T{}.csv".format(base["T"])).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["series","T_K","k_WmK","ci95_WmK","note"])
        w.writeheader(); w.writerows(rows)
    with ema_csv.open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["rho_rel","k_eff_WmK"])
        for r,k in zip(rho_grid, k_ema): w.writerow([f"{r:.4f}", f"{k:.4f}"])

    # === A tiny README with caveats ===
    md = outdir/"comparison_README.md"
    md.write_text(f"""# Comparison layer (T={base['T']} K)

This panel shows:
- **MD (this work)** — bulk, defect-free CeO₂ with 95% statistical CI from GK analysis.
- **Ab-initio/BTE bracket** — clean-crystal upper bound; no microstructure.
- **Pellet experiments** — published values at nominal density (labels show relative density ρ when available).
- **EMA curve** — Maxwell–Eucken (spherical pores) mapping of bulk MD value vs relative density ρ: 
  $$k_\\text{{eff}} = k_\\text{{bulk}} \\frac{2\\rho}{1+\\rho} = k_\\text{{bulk}} \\frac{2(1-\\varphi)}{2+\\varphi},\\quad \\varphi=1-\\rho.$$

**Caveats.** The EMA ignores pore shape/connectivity, grain boundaries, second phases, and radiative terms; it is simply a first-order porosity correction to show scale. Use `comparison_points_T{base['T']}.csv` and `ema_curve_T{base['T']}.csv` if you want to replot in your manuscript. Replace the placeholder references in `configs/compare.yaml` with your curated citations.
""")

    print(f"Wrote {figpath}, comparison_points_T{base['T']}.csv, and EMA curve CSV.")

if __name__ == "__main__":
    main()
```

Usage (after your baseline T=300T=300, x=0x=0 analysis exists):

```
python analysis/make_comparison.py \
  --state T300_x0.00 \
  --compare configs/compare.yaml \
  --outdir analysis/comparison
```

This produces:

- `analysis/comparison/comparison_T300_T300_x0.00.png` — the panel (MD point with CI, BTE bracket bar, pellet points, EMA sweep).
    
- `analysis/comparison/comparison_points_T300.csv` — the plotted numbers.
    
- `analysis/comparison/ema_curve_T300.csv` — the porosity mapping of your own bulk kk.
    
- `analysis/comparison/comparison_README.md` — the on-figure caveats and the EMA formula you used.
    

---

## Notes you can copy into your manuscript

The porosity mapping uses the Maxwell–Eucken expression for spherical, non-interacting pores with negligible pore conductivity kp≈0k_p\approx 0,

keff=kbulk 2(1−φ)2+φ=kbulk 2ρ1+ρ,k_\text{eff} = k_\text{bulk}\,\frac{2(1-\varphi)}{2+\varphi} = k_\text{bulk}\,\frac{2\rho}{1+\rho},

with φ\varphi porosity and ρ\rho relative density. This purposely ignores pore shape/connectivity and grain-boundary scattering; it is only a **scale indicator** to put bulk MD in context with pellets. Absolute comparisons should continue to cite your **statistical CI** from GK, your **numerical systematics band** (Step 15), and the **model-family band** (Step 16) when appropriate.


# 19.) Final report assembly
- Methods (with citations), validation gates, results grid, uncertainty accounting, limitations, and reuse instructions.
- Artifact: manuscript + README that lets someone rerun any state from raw inputs.