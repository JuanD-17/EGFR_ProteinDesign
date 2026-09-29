# Step 02b - Patch geometry

Generated 2026-09-29T17:07:26+00:00 by `scripts/02b_characterize_patches.py`.

Step 02 groups residues by proximity and scores them by composition. Neither operation can tell whether a patch is one usable binding face, and neither can tell whether two patches are the same face found from different seeds. This step measures both.

## Geometry

| Patch | Centre | Residues | Diameter Å | Rg Å | Planarity Å | Components | Fab | Relation |
|-------|--------|----------|------------|------|-------------|------------|-----|----------|
| 1 | 430 | 10 | 20.2 | 7.3 | 1.68 | 2 | 2.6 | direct contact |
| 2 | 460 | 11 | 17.4 | 7.5 | 1.89 | 1 | 2.6 | direct contact |
| 3 | 432 | 8 | 17.3 | 6.6 | 2.39 | 2 | 2.6 | direct contact |
| 4 | 490 | 9 | 18.8 | 7.6 | 1.8 | 1 | 2.7 | direct contact |
| 5 | 427 | 11 | 20.1 | 7.9 | 1.98 | 2 | 9.1 | same region |
| 6 | 487 | 8 | 18.9 | 7.4 | 1.44 | 3 | 2.8 | direct contact |

Planarity is the RMS deviation of the CB atoms from their best-fit plane: small values mean a flat face, large values a patch wrapped around curvature. Components counts connected groups under a 5 Å heavy-atom contact rule; more than one means the residues do not form a single continuous surface.

## Are these distinct surfaces?

| A | B | Centroid separation Å | Normal angle ° | Shared residues | Jaccard | Same face |
|---|---|----------------------|----------------|-----------------|---------|----------|
| 1 | 2 | 7.2 | 47.2 | 6 | 0.4 | no (borderline) |
| 1 | 3 | 5.3 | 55.6 | 6 | 0.5 | no |
| 1 | 4 | 19.4 | 87.8 | 0 | 0.0 | no |
| 1 | 5 | 7.2 | 135.0 | 6 | 0.4 | no |
| 1 | 6 | 12.6 | 68.5 | 1 | 0.059 | no |
| 2 | 3 | 5.5 | 48.3 | 6 | 0.462 | no (borderline) |
| 2 | 4 | 13.7 | 40.6 | 3 | 0.176 | no |
| 2 | 5 | 10.2 | 101.6 | 3 | 0.158 | no |
| 2 | 6 | 6.0 | 22.0 | 6 | 0.462 | yes |
| 3 | 4 | 16.7 | 66.3 | 0 | 0.0 | no |
| 3 | 5 | 11.4 | 79.5 | 2 | 0.118 | no |
| 3 | 6 | 11.0 | 61.2 | 2 | 0.143 | no |
| 4 | 5 | 19.5 | 69.4 | 0 | 0.0 | no |
| 4 | 6 | 9.1 | 20.4 | 4 | 0.308 | yes |
| 5 | 6 | 12.7 | 88.5 | 1 | 0.056 | no |

2 pair(s) fall within 10° of the 45.0° cutoff and their classification would flip under a modest change of threshold. They are reported as unresolved rather than decided:

- Patches 1 and 2: centroids 7.2 Å apart, normals 47.2° apart. Most likely adjacent regions of one curved surface rather than two independent faces.
- Patches 2 and 3: centroids 5.5 Å apart, normals 48.3° apart. Most likely adjacent regions of one curved surface rather than two independent faces.

Two patches are called the same face when their centroids are within 12.0 Å and their outward normals within 45.0°. Residue overlap alone cannot decide this: two patches can share few residues and still sit on one surface, which is the failure mode the step 02 deduplication has.

**4 distinct face(s) among the top 6 patches.**

- Face 1: patches 1 — best ranked is patch 1
- Face 2: patches 3 — best ranked is patch 3
- Face 3: patches 2, 4, 6 — best ranked is patch 2
- Face 4: patches 5 — best ranked is patch 5

A design campaign needs one target per distinct face. Patches on the same face are alternative framings of one surface and do not diversify the portfolio.
