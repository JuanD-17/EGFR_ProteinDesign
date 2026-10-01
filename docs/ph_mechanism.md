# pH-switch mechanism

What selectivity between pH 6.5 and 7.4 is physically worth, and which groups on the target can contribute to it.

Depends on [epitope_selection.md](epitope_selection.md).

---

### `scripts/03_ph_mechanism.py` — pH-switch feasibility

Objective 1 is the highest-weighted criterion and no design pipeline
optimises for it, so the mechanism is specified before any sequence is
generated. This script does not design anything; it establishes what a
binder would have to do and whether each surface can support it.

```bash
python scripts/03_ph_mechanism.py --patches 1 2 3 4 5
```

## Findings from step 03

### One histidine cannot do this

Wyman linkage relates the pH dependence of binding to the number of protons
taken up on binding:

```
d(ΔG_bind)/d(pH) = 2.303 · R · T · Δn(H⁺)
```

Over the 0.9 pH units between 7.4 and 6.5:

| Coupled protonations | Max ΔΔG kcal/mol | Max affinity change |
|---------------------:|-----------------:|--------------------:|
| 1 | 1.23 | 8× |
| 2 | 2.46 | 63× |
| 3 | 3.68 | 502× |
| 4 | 4.91 | 3987× |

A 100-fold switch needs **2.2 coupled protonation events**; 1000-fold needs
3.3. A single histidine, ideally placed and ideally titrating, buys about
one order of magnitude — which is not "binding at 6.5, none detectable at
7.4", and may not even be resolvable in the assay.

This rules out the intuitive one-histidine design before any compute is
spent on it, and it reframes the design families: the governing variable is
*how many* protonation events couple to binding, not which histidine is
chosen.

These are ceilings, and real coupling falls short of them. pKa values shift
at interfaces, protonation is not all-or-nothing across so narrow a window,
and a neutral histidine can still hydrogen bond to a carboxylate, which
erodes the difference between the two states.

### H370 is marginal, not excluded

Step 02 reported H370 as unavailable to a binder. That was too strong, and
the correction comes from this step's own measurement.

| | Rel. SASA | Glycan | Conservation | Clash-free partner positions |
|---|---|---|---|---|
| H370 | 0.171 | **yes** | identical | 29 |
| H433 | 0.801 | no | identical | 44 |
| D368 | 0.154 | no | identical | 29 |

Relative SASA measures the whole residue against a 0.20 threshold that is a
convention rather than physics. The question that matters for a pH switch is
narrower: can a carboxylate reach the titratable nitrogen? Probing for
clash-free partner positions at hydrogen-bond distance says yes, at 29
positions.

What stands is the glycan risk: H370 lies within 12 Å of an N-glycosylation
site. That is a sequence-derived prediction, not an observation, since the
crystal did not resolve that glycan.

H370 is therefore **marginal with glycan risk** — neither the mechanistic
hotspot the literature suggests nor the excluded position claimed here
earlier. D368 sits two positions away, so H370 and D368 form an adjacent
titratable pair.

### The ionizable environment is not basic

An earlier reading of patch composition suggested a basic surface that would
oppose placing a protonated histidine on the binder. Counting only groups
with room for a partner, net charge is zero or negative:

| Patch | Target His | Target acids | Target basics | Net | Ceiling on Δn |
|-------|-----------|--------------|---------------|-----|---------------|
| 1 | H370, H433 | D368, D458, D460 | 3 | 0 | **5** |
| 3 | H370, H433 | D368, D458, D460 | 1 | −2 | **5** |
| 2 | H433 | D368, D458, D460 | 3 | 0 | 4 |
| 5 | — | D368, E455, D458, D460 | 3 | 0 | 4 |
| 4 | — | E496 | 1 | 0 | 1 |

Reachability is tested by probing for clash-free positions at hydrogen-bond
and salt-bridge distance, because a group can be solvent-exposed and still
have no room for a partner. E400 is an example: 16% relative accessibility
and only 3 clash-free positions, so it is not usable.

Patch 4 is eliminated as a pH-switch candidate on a ceiling of Δn = 1.
Patches 1 and 3 carry two target histidines each and support both
mechanistic directions. Patch 5 has no histidine but four reachable acids,
so it is a pure binder-side design.

### `scripts/03b_protonation_networks.py` — pKa and network geometry

Step 03's ceiling of 1.23 kcal/mol per coupled protonation carries a
condition that this step makes explicit and then tests.

**Selectivity comes from the pKa shift on binding, not from the pKa.** The
linkage coefficient is `Δn(pH) = f_bound(pH) − f_free(pH)`. A group whose
pKa does not change on binding contributes exactly nothing, however
favourable its interaction. The design requirement is therefore to place
each histidine so that binding stabilises its protonated form, and burying
it against a carboxylate is the mechanism that does so.

Free-state pKa values are estimated with PROPKA, and the script then
enumerates which groups can be engaged simultaneously: between 5 and 22 Å
apart, close enough for one binder face yet far enough to titrate against
the binder rather than each other, with outward normals within 80°.

```bash
python scripts/03b_protonation_networks.py --assumed-shift 2.0
```

PROPKA 3.5 does not run on Python 3.14, because it dispatches its parameter
parser on `self.__annotations__` and PEP 649 no longer resolves that through
an instance. The script carries a documented compatibility shim that reads
the class annotations instead; it changes nothing on older interpreters and
can be deleted once PROPKA supports 3.14.

## Findings from step 03b

### The two mechanisms titrate different molecules

An earlier version of this analysis scored target carboxylates as the
titrating species and concluded that D458 and D460, at pKa 1.95 and 3.70,
were too acidic to be useful. That inverted the reading.

For a binder-side switch it is the **binder's** histidine that titrates. The
target carboxylate only has to be a dependable counter-charge that upshifts
that histidine's pKa on burial. A very low pKa is therefore not a defect but
the virtue: it guarantees the anchor is fully ionised at both pH values.

| Group | Role | pKa free | kcal/mol | Verdict |
|-------|------|---------:|---------:|---------|
| H433 | receptor-side switch | 6.22 | 0.943 | usable, best switch |
| H370 | receptor-side switch | 4.93 | 0.588 | usable, weaker |
| D458 | anchor | 1.95 | 0.848 | excellent, fully ionised |
| D460 | anchor | 3.70 | 0.848 | excellent, fully ionised |
| E455 | anchor | 4.31 | 0.848 | excellent, fully ionised |
| D368 | anchor | 4.67 | 0.848 | good |
| E496 | anchor | 4.88 | 0.848 | good |

H433 at 6.22 titrates almost ideally for this window. H370 is downshifted
from its 6.50 model value, meaning its environment already destabilises the
protonated form, which is a second reason beyond glycan risk to prefer H433.

### Viable networks reach 70–90×, not 100×

Thirty-one combinations satisfy the distance and orientation constraints.
The best:

| Network | Span Å | kcal/mol | Fold | Patches |
|---------|-------:|---------:|-----:|---------|
| D368 · H433 · D458 | 15.9 | 2.64 | 86× | 1, 2, 3 |
| D368 · H433 · D460 | 15.1 | 2.64 | 86× | 1, 2, 3 |
| H433 · D458 · D460 | 15.9 | 2.64 | 86× | 1, 2, 3 |
| D368 · E455 · D458 | 17.1 | 2.54 | 73× | 5 |
| D368 · D458 · D460 | 15.1 | 2.54 | 73× | 1, 2, 3, 5 |
| D368 · H370 · H433 | 11.3 | 2.38 | 55× | 1, 3 |

Under a 2.0-unit shift assumption the best three-group networks reach about
86-fold, short of the ~100-fold that "no detectable binding" implies. Two
routes close that gap: a larger pKa shift, which deeper burial against a
carboxylate can plausibly deliver, or accepting a clearly demonstrated but
sub-maximal switch — which the competition rules explicitly value over a
high-affinity binder with no pH dependence at all.

Patch 5 supports 73× purely binder-side, with no reliance on a receptor
histidine. That makes it the family whose mechanism is entirely under our
control, and the one that does not overlap the cetuximab epitope.

### Limitations

pKa values in the complex cannot be computed before the complex exists, so
the assumed shift stands in for the quantity that actually decides the
outcome. PROPKA is empirical with errors near one pKa unit, comparable to
the effects being reasoned about. Every anchor is scored with the same
engineered-histidine model, so all anchors return identical coupling and the
ranking among anchor-only networks is driven by group count rather than
quality; in this model anchor pKa governs reliability, not the magnitude of
the shift it induces. Rigorous coupling free energies would require
constant-pH molecular dynamics, which is out of scope.

These are a screen and a geometric filter, not predictions of experimental
pH selectivity.

