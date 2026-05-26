# Stage 3 FastRelax Remediation & Forward Workflow

**Scope:** isolate and permanently fix the apo FastRelax failure in `cc_triage.py`,
flag the other critical bugs it interacts with, and lay out the next experimental
steps (linker shortening, proline introduction). Written against the current
`cc_triage.py`, `KNOWN_ISSUES.md` (2026-05-21), and `merged_triage_results.tsv`
(40-design L10 run).

**Bottom line up front.** There is one root-cause bug and it is small. The apo and
holo poses are relaxed under *different rules* — holo is restrained, apo is not —
so the two ΔG terms in `ddG_CC = holo − apo` are not measured on the same footing.
Fix the relaxation protocol so apo and holo are scored identically and the ΔΔG
inversion you are seeing largely resolves on its own. Everything else in this doc is
either a consequence of that asymmetry or an independent bug that corrupts the same
TSV. The linker biology (shorter linkers, prolines) is real and worth doing — but
**only after** the relaxation is trustworthy, because right now you cannot rank
linker variants when 11/40 apo scores are physically impossible.

---

## 1. The data, stated plainly

From `merged_triage_results.tsv` (40 L10 designs):

- **28/40** apo ΔG fall in a plausible band (−100…+50 REU).
- **11/40** have apo ΔG **> 1000 REU** (max: `6E4H_NEGCTRL` at **+178,594**).
  These are not weak interfaces — they are nonphysical. A real α/β coiled-coil
  interface is on the order of −10 to −55 REU.
- **5/40** have apo dSASA **> 3000 Å²** — nearly half the pose treated as the A:B
  interface, exactly matching the `Interface set residues total: 232` symptom in
  the bad-run logs (job 551476).
- Holo ΔG stays sane across the board (−15 to −45 REU on the good cases) because
  holo relaxation is restrained (see §3).

The consequence the project lead flagged — **holo more stable than apo** (negative
`ddG_CC`, wrong direction for an AND-gate OFF→ON switch) — shows up even among the
28 plausible-apo designs: **10/28** have `ddG_CC < 0`. The WT reference
`7Q1R_WT_variant001` is one of them (apo −9.8, holo −20.9, ΔΔG −11.2).

> An AND-gate OFF state must be *more* stable than the antigen-bound ON state at
> the CC interface — antigen binding is what pays the thermodynamic cost of opening
> the zipper. `ddG_CC` should be **positive** (holo destabilized relative to apo).
> A systematically negative ΔΔG means either the biology is wrong or the
> measurement is. We show below it is mostly the measurement.

---

## 2. Root cause: apo and holo are relaxed under different rules

This is the single most important finding. Walk the two code paths:

**Apo path** (`process_design`, ~line 1090):
```python
make_relax(sfxn, relax_rounds).apply(apo_pose)   # unconstrained
```

**Holo path** (~line 1204–1207):
```python
n_cst = add_antigen_constraints(holo_pose, ["C", "D"], sd=0.5)  # pin antigens
make_relax(sfxn_cst, relax_rounds).apply(holo_pose)            # constrained
```

And `make_relax` (line 678) builds a MoveMap with **`mm.set_jump(True)`** —
rigid-body jumps are mobile.

Now combine that with `break_chain_in_pose` (line 303), which installs **one
FoldTree jump per chain boundary** so InterfaceAnalyzer can separate chains. In the
apo pose there is exactly one jump: the A↔B (α↔β) rigid-body transform. So apo
FastRelax is *free to translate and rotate the entire β half relative to α*, across
the cut you deliberately introduced in the GGGGS linker, with nothing holding the
two helices in register.

That is the collapse. `KNOWN_ISSUES.md::BUG-APO-IAM-RELAX` already documented the
fingerprint without naming the mechanism:

- Order matters: IAM→relax→IAM gave −51 (good); relax→IAM gave +75,497 (bad) on
  the *same PDB and same seed*. That is a deterministic basin-hopping artifact, not
  SLURM noise — confirmed by exact local repro (`dG=75496.538`).
- Bad runs show `Interface set residues total: 232`. When the β half rotates into
  the α half across the free jump, the two domains interpenetrate; IAM then counts
  half the protein as "interface," dSASA balloons to 6789 Å², and `fa_rep` from the
  steric overlap drives ΔG to +10⁴–10⁵.

In holo, the harmonic CA constraints on C/D (and the bulk of the antigen-bound
scFv arms they anchor) damp this rigid-body wandering enough that the pose stays in
a sane basin. **So the asymmetry is not cosmetic — it is the entire reason holo
looks stable and apo explodes.** The MoveMap comment in `make_relax` even says jumps
were enabled to let chains "translate apart to relieve clashes at the new termini."
That fix for the *termini-clash* problem is exactly what *creates* the
apo-collapse problem, because a single free jump between two helices with no
restraint has no reason to translate *apart* rather than *into* each other.

---

## 3. The permanent fix

The goal: **apo and holo must be relaxed identically, and neither relaxation may
let the CC helices drift out of their starting register.** This is a solved problem
in Rosetta — the standard tool is start-coordinate constraints
(`constrain_relax_to_start_coords` / `AtomCoordinateCstMover`), which the
RosettaCommons structure-prep docs recommend specifically for "if you notice that
your protein is moving too much." Right now the code half-uses this idea (only on
antigens, only in holo). Generalize it.

### 3a. Single relaxation protocol for both states (primary fix)

Replace the asymmetric relax with one constrained protocol applied to **both** apo
and holo. Constrain backbone CA atoms of the *structured* parts (both CC helices and
both scFv cores) to their start coordinates with a loose harmonic well; leave
sidechains and local backbone free to relieve clashes. Keep the per-state extra
constraints (antigen pin in holo) on top.

```python
def make_relax(sfxn, rounds: int = 1, constrain_coords: bool = True):
    """FastRelax. By default, constrain heavy-atom backbone to start coords so
    multi-domain poses cannot collapse across the single inter-helix jump.
    Sidechains and small backbone moves remain free (that is where most of the
    ref2015 energy gain comes from)."""
    from pyrosetta.rosetta.core.kinematics import MoveMap
    fr = FastRelax(sfxn, rounds)
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(False)        # <-- CRITICAL: do NOT let the α↔β jump float free
    fr.set_movemap(mm)
    if constrain_coords:
        fr.constrain_relax_to_start_coords(True)   # ramped CA coord constraints
        fr.constrain_coords(True)
    return fr
```

Two independent changes, each of which alone reduces the failure, and which together
eliminate it:

1. **`mm.set_jump(False)`** removes the rigid-body degree of freedom that lets β
   collapse into α. The original reason jumps were enabled (relieve OXT/H clashes at
   the cut) is better solved by terminus variants — which `break_chain_in_pose`
   *already* applies (`UPPER_TERMINUS_VARIANT`/`LOWER_TERMINUS_VARIANT`, line 366) —
   plus a couple of relax rounds. You do not need a free jump for that.
2. **`constrain_relax_to_start_coords(True)`** keeps backbone CA near the input CC
   geometry during relax. The score function must carry the `coordinate_constraint`
   term for this to bite — `make_scorefxn(with_coord_cst=True)` already exists; use
   it for **both** states now, not just holo.

Then in `process_design`, make the two calls symmetric:

```python
# APO  (was: make_relax(sfxn, relax_rounds).apply(apo_pose))
add_cc_constraints(apo_pose, geom, sd=0.5)             # see 3b
make_relax(sfxn_cst, relax_rounds).apply(apo_pose)
apo_pose.remove_constraints()                          # before IAM, as holo does

# HOLO  (unchanged except it now also rides the same make_relax default)
add_antigen_constraints(holo_pose, ["C","D"], sd=0.5)
add_cc_constraints(holo_pose, geom, sd=0.5)
make_relax(sfxn_cst, relax_rounds).apply(holo_pose)
holo_pose.remove_constraints()
```

Note `make_scorefxn(with_coord_cst=True)` (→ `sfxn_cst`) must be created and used in
the apo path now; today it is only built after the `skip_holo` early-return.

### 3b. Constrain the CC helices explicitly, not just antigens

Generalize `add_antigen_constraints` into a helper that pins the residues you want
held. For apo, pin both CC helices (the α residues `coil1_start..coil1_end` on chain
A and the β residues on chain B) so the interface register is preserved while
sidechains repack:

```python
def add_cc_constraints(pose, geom, sd=0.5):
    """Harmonic CA constraints on both CC helices, mapped to rechained A/B
    numbering, so the antiparallel register is preserved during relax."""
    from pyrosetta.rosetta.protocols.constraint_generator import (
        CoordinateConstraintGenerator)
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
    # chain A: coil1 occupies coil1_start..coil1_end (already 1-indexed in A)
    # chain B: coil2 occupies (coil2_start-beta_start+1)..(coil2_end-beta_start+1)
    a_lo, a_hi = geom["coil1_start"], geom["coil1_end"]
    b_lo = geom["coil2_start"] - geom["beta_start"] + 1
    b_hi = geom["coil2_end"]   - geom["beta_start"] + 1
    # ... build a ResidueIndexSelector over A:a_lo-a_hi and B:b_lo-b_hi,
    #     feed to CoordinateConstraintGenerator(ca_only=True, sd=sd), apply.
```

`sd=0.5 Å` matches the antigen pin and is loose enough to let the interface settle
without freezing it.

### 3c. Replicate + median (defends against residual basin variance)

Even constrained, FastRelax is stochastic. Per `BUG-APO-IAM-RELAX` suggestion 5,
relax the apo (and holo) pose 3× with different sub-seeds and take the **median ΔG**,
preferring replicas whose dSASA lands in the physical 600–2000 Å² window. This is a
few extra minutes per design and turns a bimodal score distribution into a stable
point estimate. Cheap insurance given what is downstream.

### 3d. Fail loud on nonphysical interfaces (stop poison entering the TSV)

Implement `BUG-APO-IAM-RELAX` suggestion 4 as a hard gate. The current apo guard
only catches dSASA≈0 (line 1096); it does **not** catch the opposite failure (dSASA
too large / ΔG exploded), which is exactly the 11/40 mode:

```python
APO_DG_MAX   = 500.0     # |REU|
APO_DSASA_MAX = 3000.0   # Å²
if abs(apo_iface["dG"]) > APO_DG_MAX or apo_iface["dSASA"] > APO_DSASA_MAX:
    row["holo_status"] = "apo_iam_failed"
    row["ddG_CC"] = float("nan")          # do NOT emit a garbage ΔΔG
    # still record the raw apo numbers for debugging, but flag the row unusable
```

A flagged row is recoverable signal; a +178,594 ΔG silently averaged into a linker
ranking is not.

### Why this is the *permanent* fix and not a patch

- It is the documented RosettaCommons protocol for multi-domain relax, not a
  bespoke workaround.
- It makes `ddG_CC` a true apples-to-apples comparison (the whole stated purpose of
  the relax step — see the Relax tutorial: "enable an apples-to-apples comparison …
  by first minimizing them … according to the same score function").
- It is deterministic-friendly: with `set_jump(False)` + start-coord constraints,
  the IAM→relax vs relax→IAM order dependence in the bug report disappears, because
  there is no longer a free rigid-body mode for the two orderings to diverge on.

### Validation gate before trusting any new run

1. Re-run the 40-design L10 set. Expect: **0** apo ΔG > 500; all dSASA in
   600–2000 Å²; WT/strong designs with **positive** `ddG_CC`.
2. Re-confirm against the one known-good reference: job 551399 gave apo −52.5 for
   `7Q1T_WT_variant001`. The fixed pipeline should reproduce a value in that
   neighborhood **every** time, not 1-in-3.
3. Spot-check that holo numbers barely move (they were already restrained; the fix
   should leave them roughly where they are, which is the proof the apo side was the
   sick one).

---

## 4. Other critical bugs that corrupt the same TSV

These are independent of the relax fix but will keep poisoning your rankings if left
alone. In priority order:

### 4a. CRITICAL — `BUG-INCOMPLETE-BACKBONE` (library is CA-only)

`KNOWN_ISSUES.md` lists this as **OPEN — library not regenerated on cluster**:
120/124 on-disk L10/L20/L30 PDBs are CA-only in the linker+CC region, and PyRosetta
silently drops those residues (`-ignore_unrecognized_res`), which mis-assigns chains
and was the original `dSASA=0` bug.

**This must be closed first.** Note the live contradiction: `cc_triage.py` line
~1022 *claims* "the production library … has been rebuilt with full N/CA/C/O on
every residue (BUG-INCOMPLETE-BACKBONE — Fixed)", but the README and KNOWN_ISSUES
scan (2026-05-21) say 120/124 still fail. `cc_triage.py` only **warns** on CA-only
input (line 1029), it does not abort — so broken PDBs flow straight into scoring.
Action:
- Run the documented rebuild (`build_cc_scaffolds.py` → `pymol_assemble.py --batch`
  for L10/L20/L30), then `validate_library_backbones.py`.
- Promote the soft warning at line 1029 to a hard `assert_full_backbone` abort
  (the function exists in `pdb_ops.py` and is already imported at line 975 — just
  call it on `apo_rechained`, not only on `design_pdb`).
- Until the rebuilt tree is validated, **do not trust any L10 numbers**, including
  the ones in `merged_triage_results.tsv`. Some of the 11 apo blow-ups may be
  compounded by dropped residues, not just the jump collapse.

### 4b. HIGH — `BUG-PYMOL-ASSEMBLY-SCAFFOLD` (4PXJ, 6Q5S helices mis-placed)

SOCKET2 batch shows 6Q5S min inter-helix CA distance **17.7 Å** and 4PXJ **11.2 Å** —
far beyond KIH packing range (7Q1T is 2.5 Å). `pymol_assemble.py` translates CC
scaffolds by N→C junction placement, which preserves each helix internally but does
**not** preserve the antiparallel dimer interface. For these two parents the CC is
not even in contact, so any apo ΔG_CC for them is meaningless regardless of the relax
fix. This explains why 6Q5S rows have near-zero dSASA_apo (60–177 Å²) and tiny ΔG —
there is no interface to score. Fix per KNOWN_ISSUES: superpose the CC dimer as a
**rigid body** from `cc_scaffolds_v2/<PARENT>.pdb`, then assert min inter-helix CA
< 8 Å before writing. **Until then, exclude 4PXJ and 6Q5S from rankings.**

### 4c. MEDIUM — `BUG-SOCKET2-PRE-RELAX` (SOCKET2 runs at the wrong time + parser/cutoff)

`cc_triage.py` runs SOCKET2 on `apo.pdb` **before** FastRelax (line 1015), but
`SOCKET2.md` is explicit that KIH packing only exists *after* relaxation — the whole
point of the negative-control test in SOCKET2.md is that pre-folding scaffolds
correctly show `socket2_detected=False`. So today's `socket2_detected=False` on every
L10 design is expected-but-uninformative, not a real signal. Three sub-fixes:
1. Move the SOCKET2 call to run on `apo_relaxed.pdb` (after the fixed relax).
2. Adaptive cutoff (7.0 Å then 8.0 Å) in `lib/socket2_filter.py` — 7Q1T only passes
   at 8.0 Å (loose AP regime). The library scanner already does this; sync the
   production wrapper to match.
3. Fix `_parse_socket_short` to parse v3.02 `-q` stdout (`coiled coil  0:`); right
   now `parsed` is always empty so even a detected CC reports `type=none`.
This is MEDIUM not CRITICAL because SOCKET2 is non-fatal and not a ranking term —
but a fixed SOCKET2 is your *independent structural check* that the relaxed CC is a
real coiled coil, which is exactly the validation you want post-fix.

### 4d. LOW — `BUG-TRIAGE-MERGE-SLURM` (stale rows in merged TSV)

`set -o pipefail` fails under `/bin/sh` on compute nodes, the merge job dies, and
part files accumulate **stale rows** (`4PXJ_MUTATED` + new `7Q1T` in the same part).
This means `merged_triage_results.tsv` may itself contain duplicate/stale
construct rows. Fix: force `#!/bin/bash`, overwrite (not append) one row per
`SLURM_ARRAY_TASK_ID`, dedupe by `construct_id` in the merge. Worth doing before you
generate the post-fix ranking, or you will be ranking a polluted table.

---

## 5. The linker work — correct hypothesis, right sequencing

The project lead's two design ideas are sound and literature-supported. The key
constraint: **do them after §3 and §4a/4b are closed**, because you cannot rank
linker variants on a scoring pipeline that produces nonphysical apo energies and
half-broken input geometry. With that ordering caveat, here is how to ground each.

### 5a. Why shorter CC↔scFv linkers help an AND-gate

The OFF-state CC must hold both scFv arms autoinhibited; antigen binding pays to
open it. A long, floppy GS linker between CC and scFv lets the arms sample
binding-competent poses *without* paying the full zipper-opening cost — it leaks
the gate. Shortening the CC↔scFv linker raises the mechanical coupling between
antigen engagement and zipper strain: the arm cannot reach its epitope unless the
zipper gives. In ΔΔG terms, a tighter linker should make holo *less* able to relax
its CC interface (more positive `ddG_CC` = better gate), which is exactly the
direction you want and the opposite of what the broken pipeline currently reports.

This is consistent with the CAR-T literature already cited in the project report
(Majzner et al.: hinge composition tunes the activation threshold across ~an order
of magnitude). Linker length here is the analogous tuning knob.

### 5b. Why prolines, and the one caveat

Proline introduces a backbone kink and lacks an amide H, so it suppresses secondary
structure and stiffens a linker — the standard rigid-linker motifs are
poly-proline `(AP)n` and `(EAAAK)n` (α-helical), versus the flexible `(GGGGS)n`. A
stiffer CC↔scFv junction further mechanically couples antigen binding to zipper
strain (the rigid-linker "guided domain separation" effect). **Caveat the lead
already anticipated by saying "step by step":** proline is also a helix-breaker, and
your inner linker sits *inside* the GGGGS that bridges the two CC helices (the α/β
cut is between residues 2 and 3 of that GGGGS, `compute_split`). Putting proline in
the *inner* linker risks disrupting the CC backbone register itself. So:

- **Outer linkers (CC↔scFv): yes**, candidates for shortening and proline
  stiffening — that is the mechanical-coupling knob.
- **Inner linker (α↔β GGGGS): leave alone** for now, or treat any change there as a
  separate, carefully-validated experiment — it defines the antiparallel turn and
  the scoring cut point.

### 5c. How to run it through the (fixed) pipeline

This maps cleanly onto **Stage 5 — Bindsweeper** in `PIPELINE.md`, which already
sweeps `linker_length` (7/14/21 aa). Extend it:

1. **Linker-length sweep first** (one variable). Generate L7/L10/L14 (and shorter:
   L4, L6) outer-linker constructs via the existing partial-diffusion / Bindsweeper
   path, holding everything else fixed. This is the secondary RFD mode in
   `PIPELINE.md` Stage 1 — "perturb and re-score known-good templates," exactly the
   `provide_seq`-free partial diffusion the report describes.
2. **Re-score with the fixed Stage 3.** Rank by `ddG_CC` (now trustworthy and
   positive-is-good), with `binderA_dG`/`binderB_dG` confirming the arms still bind
   (< −10 REU per the PIPELINE.md pass criteria) and the fixed SOCKET2 confirming
   the apo CC is still a real AP coiled coil.

(Full document continues with additional practical steps, code patches, and
validation guidance from the original FASTRELAX_WORKFLOW.md in the draft repo.)
