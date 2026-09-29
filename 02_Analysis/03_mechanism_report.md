# Step 03 - pH-switch mechanism

Generated 2026-09-29T17:20:13+00:00 by `scripts/03_ph_mechanism.py`.

Nothing here designs a binder. It establishes what a binder would have to do, and whether each candidate surface can support it.

## How much selectivity is physically available

Wyman linkage gives the pH dependence of binding as

```
d(dG_bind)/d(pH) = 2.303 * R * T * dn(H+)
```

Over the 0.9 pH units between 7.4 and 6.5:

| Coupled protonations | Max ΔΔG kcal/mol | Max affinity change |
|---------------------:|-----------------:|--------------------:|
| 1 | 1.23 | 8× |
| 2 | 2.46 | 63× |
| 3 | 3.68 | 502× |
| 4 | 4.91 | 3987× |

A 100-fold switch requires 2.2 coupled protonation events; 1000-fold requires 3.3.

**A single histidine cannot deliver the requested behaviour.** Ideally placed and ideally titrating, one gives about one order of magnitude. The challenge asks for binding at 6.5 and none detectable at 7.4. Designs must therefore couple two to three protonation events, and this rules out the intuitive one-histidine design before any compute is spent on it.

These are ceilings. Real coupling is lower: pKa values shift at interfaces, protonation is not all-or-nothing across this narrow window, and a neutral histidine can still hydrogen bond to a carboxylate, which erodes the difference between the two states.

## Two mechanistic directions

The bound state is the acidic one, so both directions must *gain* affinity on protonation.

**Receptor-side.** A target histidine protonates at 6.5 and pairs with a carboxylate placed on the binder. Limited to the histidines the target presents, and their pKa is not ours to tune.

**Binder-side.** A histidine on the binder protonates and pairs with a carboxylate on the target. Under our control in number and placement, so this is the direction that can raise Δn, but it requires acidic residues on the target surface and is opposed by nearby basic residues, which repel the protonated histidine.

## Ionizable environment of each surface

A group counts as reachable when a partner atom could occupy a clash-free position at hydrogen-bond or salt-bridge distance from it. Relative accessibility alone is not sufficient: a group can be exposed and still have no room for a partner.

| Patch | Target His | Target acids | Target basics | Net charge | Ceiling on Δn |
|-------|-----------|--------------|---------------|-----------|---------------|
| 1 | H370, H433 | 368, 458, 460 | 427, 429, 431 | +0 | 5 |
| 2 | H433 | 368, 458, 460 | 431, 487, 489 | +0 | 4 |
| 3 | H370, H433 | 368, 458, 460 | 431 | -2 | 5 |
| 4 | — | 496 | 489 | +0 | 1 |
| 5 | — | 368, 455, 458, 460 | 427, 429, 431 | -1 | 4 |

Net charge is the count of reachable fixed-positive groups minus reachable acids. A strongly positive environment argues against the binder-side mechanism there, since a protonated histidine would be electrostatically opposed.

## Evidence tiers

Kept separate so the methods report can state what each claim rests on.

| Tier | Content |
|------|---------|
| Experimental evidence | H433, Q435 and K489 appear in published alanine scanning as positions that reduce binding. H370 and H433 are described in the literature in connection with pH dependence. |
| Structural observation | Accessibility, conservation, glycan occlusion, ionizable environment and partner room, all measured from 6ARU and UniProt in steps 01 to 03. H370 and L406 are not accessible in this structure. |
| Mechanistic hypothesis | That a given arrangement of protonatable groups at this interface will couple to binding strongly enough to produce the requested switch. Untested. |
| Computational criterion | Metrics used to rank designs at the generation stage. Not yet defined. |

## Design families

Families are organised by *how many* protonation events they attempt to couple, not by which histidine they centre on. The thermodynamics above make the count the governing variable.

- **Patch 1** (centre 430): ceiling Δn = 5, net charge +0. Supports both directions; a combined design is the only route to Δn ≥ 2 here.
- **Patch 2** (centre 460): ceiling Δn = 4, net charge +0. Supports both directions; a combined design is the only route to Δn ≥ 2 here.
- **Patch 3** (centre 432): ceiling Δn = 5, net charge -2. Supports both directions; a combined design is the only route to Δn ≥ 2 here.
- **Patch 4** (centre 490): ceiling Δn = 1, net charge +0. Binder-side only.
- **Patch 5** (centre 427): ceiling Δn = 4, net charge -1. Binder-side only.

H370 is deliberately excluded as a design target. It is buried and glycan-shadowed in this structure, so it cannot be a direct contact. That does not exclude a mechanistic role in binding or pH response, but such a role would be allosteric or conformational and is not something a first-generation design can exploit rationally.
