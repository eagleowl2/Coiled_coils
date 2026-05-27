# Coiled-Coil Logic-Gate Validation Pipeline

In-silico validation of a CC-gated bispecific antibody. The logic gate uses an
antiparallel coiled-coil (CC) to hold two scFv arms in an OFF state; antigen
binding disrupts the CC and releases both arms (ON state).

![Pipeline](docs/pipeline.svg)

Two workflows are available depending on your needs:

| Workflow | When to use |
|----------|-------------|
| [Full pipeline](#full-pipeline) | Automated end-to-end runs on a SLURM cluster |
| [Direct RFD](#direct-rfd-workflow) | Manual step-by-step, debugging, CC topology exploration |

## Chain layout convention

All downstream stages expect this canonical chain layout after `merge_scaffold.py`:

| Chain | Contents |
|-------|----------|
| Z | Antiparallel zipper / coiled-coil domain |
| A | scFv-A binder |
| B | scFv-B binder |
| C | Antigen-A target (CD103 epitope) — ON state only |
| D | Antigen-B target (CD69 epitope) — ON state only |

OFF-state constructs (Z/A/B trimers) omit chains C and D.

## PDB input prerequisites for Stage 3 (Rosetta)

Stage 3 (`scripts/cc_triage.py`) requires **full N/CA/C/O backbone atoms on
every residue**. PyRosetta silently drops residues with incomplete backbones
(the default `-ignore_unrecognized_res true` policy), which leads to wrong
chain assignment and `dSASA=0` interface scores.

RFDiffusion outputs are CA-only on the inpainted segment. The upstream
pipeline **must** run a backbone-completion step (ProteinMPNN + FastRelax,
or pdbfixer, or OpenMM Modeller) before saving to `coils_library/test_coils/`.

### Library regeneration via `pymol_assemble.py`

The original `coils_library/test_coils/L10/` PDBs were saved with CA-only
linker+CC residues — these silently fail Stage 3. To rebuild them:

```bash
# One-time: build full-backbone CC scaffolds from RCSB crystals
C:/Users/user837/miniconda3/envs/pymol_env/python.exe scripts/build_cc_scaffolds.py

# Batch-rebuild every construct in cc_sequences.tsv
C:/Users/user837/miniconda3/envs/pymol_env/python.exe scripts/pymol_assemble.py \
    --cc-seqs coils_library/test_coils/L10/cc_sequences.tsv \
    --batch \
    --out-root coils_library/test_coils_rebuilt/L10/
```

Outputs:
- `coils_library/cc_scaffolds_v2/<PARENT>.pdb` — full-backbone antiparallel
  CC scaffolds, one per parent, built from the natural dimer chains of the
  RCSB crystal joined by a fab-built GGGGS turn.
- `coils_library/test_coils_rebuilt/L10/{WT,NEGCTRL,MUTATED}/*.pdb` — full
  479-AA (or 531-AA for 4PXJ) constructs with N/CA/C/O on every residue
  and 1.33 A peptide bonds at all piece junctions.
- `*.muts.tsv` sidecar per construct — residues whose scaffold identity
  differs from the design target; Stage 3 PyRosetta applies these via
  `mutate_residue` with the Dunbrack rotamer library.

Validate a PDB locally without PyRosetta:

```bash
python3 -c "
import sys; sys.path.insert(0, 'scripts')
from pdb_ops import validate_full_backbone
print(validate_full_backbone('coils_library/test_coils/L10/WT/6E4H_WT_variant001.pdb'))
"
```

`ca_only > 0` means the PDB is not safe for Rosetta. See
[KNOWN_ISSUES.md](KNOWN_ISSUES.md) (`BUG-INCOMPLETE-BACKBONE`) for root cause, scan results, and rebuild steps.

```bash
python scripts/validate_library_backbones.py --root coils_library/test_coils/L10
```

### Coil-only libraries (Exp03)

CC dimer geometry benchmark without scFv: [`docs/COIL_CC_VARIANT_STUDY.md`](docs/COIL_CC_VARIANT_STUDY.md).

| Path | Role |
|------|------|
| `coils_library/coil_variants_unlinked/` | Crystal helix A+B (reference for NC screening) |
| `coils_library/coil_variants_sg_linked/` | GGGGS hairpin from `cc_scaffolds_v2/` |
| `experiments/exp03_coil_cc_geometry/` | Archived Rosetta scores + `REPRODUCE.sh` |

```bash
bash experiments/exp03_coil_cc_geometry/REPRODUCE.sh
```

## Project structure (scripts/)

| File | Role | PyRosetta? |
|------|------|------------|
| `pdb_ops.py` | Biopython-only PDB rewriting: `rechain_pdb`, `validate_full_backbone`, `superpose_antigen_onto_design`, `count_pdb_chain_residues`, `read_chain_layout`, `fetch_crystal` | No |
| `test_pdb_ops.py` | Synthetic + real-PDB tests for `pdb_ops` — runs locally without cluster | No |
| `cc_triage.py` | Stage 3 orchestrator: apo/holo scoring, chain-break, FoldTree, IAM | Yes |
| `diagnose_pose.py` | Standalone: load a PDB → dump pose chain state → trial IAM before/after fix | Yes |
| `test_cc_triage.py` | Smoke tests (PureTests + RosettaTests) | Partial |
| `02_rosetta.slurm` | SLURM script — single design PDB per array task | Yes |
| `submit_rosetta.sh` | Build design list + sbatch the array + dependent merge | No |

To run the pdb_ops tests locally:
```bash
python3 scripts/test_pdb_ops.py -v
```
Requires `biopython` and `numpy` only.

---

## Full Pipeline

Nextflow-orchestrated, stage-based pipeline. Runs RFDiffusion → ProteinMPNN →
Rosetta triage → AF3/Boltz geometry check → optional MD.

### Quick start

```bash
export PYROSETTA_AVAILABLE=1
# Optional — skip if not running AF3/Boltz
export AF3_BIN=/opt/alphafold3/bin/alphafold3
export AF3_MODEL_DIR=/data/af3/models
export AF3_DB_DIR=/data/af3/databases
export BOLTZ_BIN=/opt/boltz/bin/boltz

./run_pipeline.sh \
  --proteindj    /path/to/proteindj \
  --structures   /path/to/input_structures \
  --neg-controls /path/to/neg_controls \
  --outdir       ./results \
  --gpus         4 \
  --profile      slurm
```

### Required input layout

```
input_structures/
  zipper/     ← antiparallel CC scaffold PDB (chain A)
  scFv_A/     ← CD103-bound scFv complex (chain A=scFv, chain B=antigen)
  scFv_B/     ← CD69-bound  scFv complex (chain A=scFv, chain B=antigen)
```

### Stages

| Stage | Script | Description |
|-------|--------|-------------|
| 0 | `calibrate_models.py` | Negative control calibration (AF3 + Boltz + Rosetta) |
| 1 | `merge_scaffold.py` | RFDiffusion CC variants → merged Z/A/B trimers |
| 2 | `run_mpnn_batch.py` | ProteinMPNN sequence design |
| 3 | `rosetta_triage.py` | **PRIMARY ranking** — FastRelax + InterfaceAnalyzer |
| 4 | `geometry_check.py` | AF3/Boltz topology filter |
| 5 | `write_sweep_yaml.py` | Bindsweeper linker/heptad sweep (optional) |
| 6 | `run_md.py` | MD simulation — 50 ns gate → 300 ns production (optional) |
| final | `assemble_shortlist.py` | Handoff report |

Stage 3 (Rosetta ΔΔG) is the primary ranking signal.
AF3/Boltz ipTM is used only as a geometry filter, not for ranking.

### Stage shortcuts

```bash
--light   # stages 0,3,4,final — skip RFDiffusion, MPNN, sweep, MD
--merged  # stages 0,1,2,3,4,final — skip sweep and MD
--stages 3,4,final  # run specific stages only
```

### Dependencies

Required:
- `bash` (Linux)
- `python3 >= 3.10`
- `nextflow >= 24.04`
- PyRosetta (`export PYROSETTA_AVAILABLE=1`)
- `pandas`, `numpy`, `biopython`

Optional (pipeline degrades gracefully if absent):
- AlphaFold3, Boltz-2, Bindsweeper, OpenMM

### Outputs

```
OUTDIR/
  s0_nc_validation/calibration.json
  s3_rosetta_triage/rosetta_passed.csv   ← primary shortlist
  s4_geometry_check/geometry_passed.csv
  s6_md/md_passed_manifest.csv
  final_shortlist/handoff_report.tsv     ← wet-lab handoff
```

---

## Direct RFD Workflow

Manual, step-by-step workflow for CC topology exploration using RFDiffusion
directly on a SLURM GPU cluster. See [`rfd/README.md`](rfd/README.md) for
full instructions.

```
rfd/
  install_rfd_env.sh   ← one-shot environment installer
  run_rfd.slurm        ← SLURM job script
  README.md            ← step-by-step instructions
```

---

## Tests

```bash
pip install pytest biopython pandas numpy
python -m pytest tests/ -v
```

109 tests cover Kabsch correctness, chain-layout regressions, NaN handling,
calibration severity bands, and shell-syntax validation.
PyRosetta / AF3 / Boltz binaries are not required to run the tests.

## Deploying to a Linux server

```bash
git clone https://github.com/eagleowl2/Coiled_coils_validaton.git
cd Coiled_coils_validaton
# .gitattributes enforces LF line endings on clone
chmod +x run_pipeline.sh run_s0_scaffold_prep.sh
```
