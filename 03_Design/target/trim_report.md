# Step 04 - Trimmed target and BindCraft configuration

Generated 2026-09-29T17:57:22+00:00 by `scripts/04_prepare_target.py`.

## What was cut

- Source: `6ARU.pdb`, chain A, 609 resolved residues
- Kept: domain III, UniProt 335–538 = PDB 311–514
- Result: 204 residues, 1571 heavy atoms
- Removed: chains B, C, D, E (the cetuximab Fab), all heteroatoms including glycans and water, hydrogens, and alternate conformations beyond the first
- No chain breaks inside the trim

Trimming is BindCraft's own recommendation: AlphaFold2 backpropagation scales badly with target size, and GPU memory is the binding constraint on the free tiers available here. Going from 609 to 204 residues is the difference between a job that fits on a T4 and one that does not.

## Numbering

The trimmed file preserves 6ARU's numbering, which is the mature protein and runs 24 below UniProt. Every hotspot number derived in steps 02 and 03 therefore remains valid without translation. Hotspots below are in PDB numbering, which is what BindCraft reads.

## Configurations

### Patch 1

- File: `bindcraft_patch1.json`
- Hotspots, PDB numbering: `403,405,406,407,408,409,410,411,434`
- Same residues, UniProt: 427, 429, 430, 431, 432, 433, 434, 435, 458
- Binder lengths: 55–100 aa
- Target accepted designs: 100

### Patch 2

- File: `bindcraft_patch2.json`
- Hotspots, PDB numbering: `407,408,409,410,411,434,436,459,463,464,465`
- Same residues, UniProt: 431, 432, 433, 434, 435, 458, 460, 483, 487, 488, 489
- Binder lengths: 55–100 aa
- Target accepted designs: 100

### Patch 5

- File: `bindcraft_patch5.json`
- Hotspots, PDB numbering: `400,403,405,406,407,431,433,434,458,459`
- Same residues, UniProt: 424, 427, 429, 430, 431, 455, 457, 458, 482, 483
- Binder lengths: 55–100 aa
- Target accepted designs: 100

## Running

Upload `EGFR_domainIII.pdb` and the JSON files to the GPU environment, then for each patch:

```bash
python -u ./bindcraft.py \
    --settings './bindcraft_patch1.json' \
    --filters './settings_filters/default_filters.json' \
    --advanced './settings_advanced/default_4stage_multimer.json'
```

Run a handful of trajectories first and time them before committing the batch, since the per-trajectory cost sets the whole budget and it is cheaper to learn it early.
