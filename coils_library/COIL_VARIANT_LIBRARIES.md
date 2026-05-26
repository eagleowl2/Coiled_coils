# Coil variant libraries (index)

Controlled **CC-only** construct sets for geometry and apo-scoring experiments (Exp03).

## Libraries

| Directory | Layout | WT template | Use for |
|-----------|--------|-------------|---------|
| [`coil_variants_unlinked/`](coil_variants_unlinked/) | Chain A + B (crystal) | RCSB helices per `PARENT_SPECS` | **NC / mutation screening** (reference) |
| [`coil_variants_sg_linked/`](coil_variants_sg_linked/) | Single-chain hairpin | `cc_scaffolds_v2/<PARENT>.pdb` | Hairpin geometry stress-test |

Each library contains **40 constructs** (5 parents × [1 WT + NC + 5 MUTATED]).

## Shared metadata (per library)

| File | Description |
|------|-------------|
| `cc_sequences.tsv` | Target sequences + mutation strings |
| `construct_manifest.tsv` | Geometry indices (`alpha_end`, chain lengths) |
| `mutations_log.tsv` | Applied mutations (`A:L24A`, …) |
| `backbone_validation.tsv` | Backbone + min-Cα geometry gate |
| `WT/`, `NEGCTRL/`, `MUTATED/` | PDB files |

## Build (reproducible)

```bash
python tools/build_coil_variant_libraries.py
python tools/validate_coil_geometry_batch.py
```

## Score (apo CC)

```bash
bash pipeline/stage3_rosetta/submit_cc_only_coils.sh both
```

## Documentation

- **Study log:** [`docs/COIL_CC_VARIANT_STUDY.md`](../docs/COIL_CC_VARIANT_STUDY.md)
- **Archive:** [`experiments/exp03_coil_cc_geometry/`](../experiments/exp03_coil_cc_geometry/)
