

Here’s a pragmatic, three-person plan that gets you an A in 12 weeks without drowning in “nice-to-have”s. I’ll spell out what’s **mandatory vs. optional**, how to **parallelize**, and a **week-by-week sprint map** with clear checkpoints. The through-line is: lock a lean scope; validate once, well; then replicate cleanly.

# What’s mandatory vs. optional

**Must-do (critical path to a defensible result)**  
1, 2*, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 19  
(*Force-field selection is “must-do,” but keep it lean: pick one fixed-charge FF with explicit Ce³⁺/Ce⁴⁺ and move on.)

**Recommended but trim-able (do the minimum viable version)**  
12 (rNEMD) — run on exactly **one** state (stoichiometric baseline) as a cross-check.  
15 (Sensitivity slices) — run **one** slice (tighter PPPM **or** smaller Δt) on the baseline only.  
16 (Limitations ledger) — single page; fast; **do it** (boosts credibility at low cost).  
17 (Packaging) — one tagged archive is quick and makes grading easy; **do it** once at the end.

**Truly optional (only if you have spare time)**  
18 (Comparison layer) — a one-figure appendix. Nice polish, not required for the grade.

# Dependency spine (what can run in parallel)

- **Foundations that can proceed in parallel immediately:**  
    1 (scope) ↔ 2 (FF dossiers) ↔ 3 (scaffolding).
    
- **Build & numerics:**  
    4 → 5 → 6. (While 6 is being audited, the vacancy generator in 9 can be prototyped against a placeholder structure.)
    
- **Validation fork:**  
    7 (EOS) → 8 (baseline GK). During 7, another person can finalize 9 (generator) and write 10/11 runners.
    
- **Grid phase (repeatable):**  
    10, 11, 13, 14 are a loop you can shard across people/machines once 8 passes.
    
- **Cross-checks / polish:**  
    12 and 15 (baseline only) can run while the grid loop is executing.  
    16, 17, 19 draft in parallel once the first few states are analyzed.
    

# Three-person split (roles that minimize blocking)

**Person A — Physics & validation lead**  
Owns 2, 4–8; signs off the baseline. After that, handles 12 (single rNEMD) and 15 (one sensitivity slice).

**Person B — Infrastructure & orchestration**  
Owns 3, 6 (numerics card + audit), 10–11 runners, 14 grid launcher, 17 packaging. Keeps manifests/ledgers clean.

**Person C — Defect generation & analysis**  
Owns 5 (structure scripts), 9 (vacancy + Ce³⁺ tool), 13 (HCACF analysis), 16 (limitations page), 19 (report assembly). Helps Person B execute grid.

# Scope guardrails (so you don’t bog down)

- **Grid size:** start with **2×2** (e.g., T={300,700}T=\{300,700\} K, x={0.00,0.06}x=\{0.00,0.06\}). If the baseline passes early, expand to 3×3 only if time permits.
    
- **Size checks:** elongated-axis test **only on the baseline** (already in your plan).
    
- **rNEMD:** **one** state (baseline).
    
- **Sensitivity:** **one** dial (PPPM accuracy _or_ Δt) on the baseline.
    
- **Seeds/shards:** begin with **2 shards/state**; add a third only if the plateau CI misses your target.
    
- **Force field:** pick one vetted fixed-charge model with explicit Ce³⁺/Ce⁴⁺; defer polarizable models.
    

# 12-week sprint map (concurrency baked in)

**Week 1 — Kickoff & foundations**  
A: (1) scope brief; (2) quick FF scan, pick the primary.  
B: (3) scaffolding/templates; repo skeleton; manifests.  
C: (5) stoichiometric builder; (4) sizing table for agreed xx.

**Week 2 — Numerics & EOS prep**  
B: (6) numerics card + NVE drift + PPPM audit (lock includes).  
A: (7) EOS scan input; pressure points list.  
C: (9) vacancy generator prototype (API, no Ce³⁺ yet).

**Week 3 — EOS fit & freeze; finalize generator**  
A: run EOS points; fit; freeze volume; quick NVT/RDF; NVE sanity.  
C: finish 9 (blue-noise VOV_O + nearest-neighbor Ce³⁺); write artifacts.  
B: wire 10 (equilibrate) and 11 (GK shard) runners against manifests.

**Week 4 — Stoichiometric kk baseline**  
A: (8) baseline GK shards; plateau checks; elongated-axis test.  
C: (13) HCACF analyzer; window t⋆=cτt^\star=c\tau; CI from blocking.  
B: glue: pass/FAIL gates captured in a short `K_VALIDATION.md`.

**Week 5 — rNEMD & sensitivity (baseline-only) while grid tooling finishes**  
A: (12) rNEMD on baseline; gentle rate check.  
B: (14) grid runner + ledger ready; dry-run on one reduced state.  
C: (15) one sensitivity slice (tighter PPPM **or** Δt) on baseline; record Δk “systematics band.”

**Weeks 6–9 — Grid execution loop (parallel workhorses)**  
B: launch states (10→11); keep ledger green; babysit runs.  
C: analyze as shards land (13); push `k_result.json` and plots.  
A: spot-check plateaus and occasional logs; decide if a 3×3 grid is safe—if not, freeze at 2×2.

**Week 10 — Consistency & consolidation**  
B: (14) consistency checker (numerics/analysis knobs identical across states) → PASS.  
C: assemble `k_grid.{csv,json,md}`; (16) write the one-page limitations ledger.

**Week 11 — Packaging & polish**  
B: (17) package release bundle (without raw flux); checksums; tag.  
C: optional (18) comparison figure (one panel) if time allows.  
A: sanity read of methods/results; add any missing gate screenshots.

**Week 12 — Report assembly & submission**  
C: (19) write the report: short methods, the validation gates (with figures), the grid table, uncertainty breakdown (stat CI + small systematics band), the limitations page link, and rerun instructions.  
A/B: final pass; submit archive + report.

# Minimal success checklist (what must be in the final package)

- **Validation:** EOS pass note; baseline GK plateau with elongated-axis check; rNEMD agreeing within uncertainty (single state).
    
- **Results:** k(T,x)k(T,x) table for at least **2×2** grid with 95% CIs; consistent numerics across states.
    
- **Uncertainty:** one small **systematics band** on the baseline; **limitations** page.
    
- **Reproducibility:** manifests, numerics card, analysis control; package tagged.
    
- **Report:** methods (short), validation gates, results with CIs, limitations, rerun recipe.
    

If you stick to that skeleton—lean grid, single rNEMD, one sensitivity slice—you’ll have a crisp, defensible story that hits all the grading boxes without scope creep.




---




# Roles (who owns what)

**You (PhD, comp mech/chem): Physics lead + validation gatekeeper**

- Own: 1, 2 (pick _one_ fixed-charge FF), 4–5, 7–8 (EOS + baseline GK), plus light 12 (one rNEMD) and 15 (one sensitivity slice) on the **baseline only**.
    
- Decisions: pass/fail on EOS, plateau acceptance, “freeze numerics,” and when to expand from 2×2 to 3×3 grid.
    
- Artifacts to sign off: `EOS_VALIDATION.md`, baseline `K_VALIDATION.md`.
    

**Jac (Materials Science): Potentials/defects + structure**

- Own: 2 (fill the dossier table), 5 (stoich build), 9 (vacancy + Ce³⁺ generator), 10 (equilibration protocol per state).
    
- Stretch: sanity of Ce³⁺/V_O statistics; quick RDF/peak checks after NVT.
    
- Artifacts: force-field dossier page, `vacancy_manifest.json`, separation histograms, equilibrated restarts.
    

**Akshay (CS/Stats/Physics): Tooling, numerics, analysis, automation**

- Own: 3 (repo scaffolding, manifests, seed policy), 6 (numerics card + quick PPPM/Δt audit), 11 (GK shard runner), 13 (HCACF analyzer with Sokal window + blocking), 14 (grid runner & consistency checker), 17 (packaging).
    
- Stretch: tiny comparison plot (18) if time.
    
- Artifacts: includes for numerics, `k_result.json` per state, `k_grid.{csv,json}`, `LIMITATIONS.md`, release tarball.
    

# What’s optional vs. trimmed

- **Do** rNEMD (12) only for the **baseline**.
    
- **Do** one sensitivity slice (15) on the **baseline** (either tighter PPPM **or** smaller Δt).
    
- **Skip** the full comparison layer (18) unless you have slack in Week 11.
    
- **Grid size:** start 2×2 (T={300,700}, x={0.00,0.06}); go 3×3 only if Week 8 plateaus look easy.
    

# Concurrency (who works while who’s waiting)

**Weeks 1–2 (parallel):**

- You: 1+2 decision; 4 sizing; green-light the chosen FF.
    
- Jac: 5 stoich build; 9 generator prototype.
    
- Akshay: 3 scaffolding; 6 numerics card + NVE drift + PPPM quick audit.
    

**Weeks 3–4 (parallel):**

- You: 7 EOS scan → fit → freeze; then 8 baseline GK + elongated-axis.
    
- Jac: finish 9 (blue-noise V_O + nearest-Ce³⁺) and 10 equilibration runner.
    
- Akshay: 11 GK shard runner + 13 analyzer; wire 14 grid skeleton.
    

**Week 5 (baseline polish while grid tools bake):**

- You: 12 rNEMD (one state) and 15 (one sensitivity slice) → small Δk “systematics band.”
    
- Jac: final checks on defect stats; small fixes if needed.
    
- Akshay: finalize analysis outputs and consistency checks.
    

**Weeks 6–9 (assembly line):**

- Jac: 10 equilibrate each (T,x) and hand off restarts.
    
- Akshay: 11 run GK shards; 13 analyze and push `k_result.json`.
    
- You: spot-check plateaus; decide whether to expand to 3×3 or lock 2×2.
    

**Weeks 10–12 (wrap):**

- Akshay: 14 consistency PASS, 17 package/tag; optional 18 comparison.
    
- Jac: help QA plots; make sure defect diagnostics are in the appendix.
    
- You: 16 limitations page (fast) and 19 report assembly with the gate figures.
    

# Concrete weekly checkpoints (lightweight)

- **W2 end:** FF chosen; numerics card passes drift/PPPM checks; stoich POSCAR/LAMMPS data in repo.
    
- **W4 end:** EOS PASS and baseline GK `K_VALIDATION.md` committed.
    
- **W5 end:** rNEMD single-state result within GK CI; systematics Δk note (<~5–10%).
    
- **W9 end:** all 2×2 states analyzed with `k_result.json` and plots; `k_grid.csv` exists.
    
- **W11 end:** packaging tarball + `LIMITATIONS.md`.
    
- **W12:** report submitted.
    

# Guardrails so you don’t bog down

- **One** force field (explicit Ce³⁺/Ce⁴⁺). No polarizable model this semester.
    
- **Two** shards/state to start; add a third only if a plateau misses the CI target.
    
- **One** elongated-axis test (baseline only).
    
- **One** rNEMD (baseline only).
    
- **One** sensitivity slice (baseline only).
    

# Quick team note you can drop in chat

> I’ll lead physics/validation (EOS + baseline GK + quick rNEMD/sensitivity).  
> Jac, you own potentials/defects and equilibrations; Akshay, you own tooling, numerics, GK runners, and analysis.  
> We keep a 2×2 grid unless the baseline is easy by Week 5. One rNEMD + one sensitivity slice on the baseline only. Everything else is frozen numerics and repetition. Our Week-4 gate is “EOS PASS + baseline kk plateau”; Week-9 gate is “full grid analyzed with CIs.”

If we stick to that, Akshay can flex the data/automation superpowers, Jac keeps us physically sane on defects and structures, and you make the critical calls. That’s the fastest path to an A without scope creep.