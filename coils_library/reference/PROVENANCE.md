# Reference Linker-Frame PDB — Provenance

**File:** `linker_frame_6E4H_3H3B_7URV.pdb`

**Original filename:** `6E4H_short_two_linkers_3H3B_7URV_no_receptors_fused_final (2).pdb`

**Origin:** Manually assembled in PyMOL by fusing:
- CC scaffold: 6E4H (PALB2 antiparallel CC dimer, chains A+B)
  with short outer linkers (GGGGSGGGGS, 10 aa each)
- scFv-A arm: 3H3B crystal (anti-CD20 scFv + CD20 antigen)
- scFv-B arm: 7URV crystal (anti-CD103 scFv + CD103 antigen)
- Receptors and antigen chains removed ("no_receptors" suffix)
- Final manual clash-check in PyMOL ("final (2)" suffix)

**Role in pipeline:**  
`assemble_construct.py` Kabsch-aligns the three primary components
(scFv-A from 3H3B, CC dimer from parent crystal, scFv-B from 7URV)
onto this reference frame to place them in a mutually clash-free
spatial arrangement.  The linker residue geometries (GGGGSGGGGS at
positions 248-259 and GGGGS at positions 289-295) are also extracted
from this file for grafting.

**When committed:** 2026-05-21  
**Committed by:** user837

**Reproducibility note:**  
If this file must be rebuilt, fuse a fresh 6E4H CC dimer (chains A+B,
helix residues only) with 3H3B scFv and 7URV scFv in PyMOL using
"fuse" + "align", remove all antigen/receptor chains, verify no
steric clashes between the CC and both scFv arms, and recheck the
linker residue numbers match LINKER_TEMPLATE_SLICES in
`assemble_construct.py`.
