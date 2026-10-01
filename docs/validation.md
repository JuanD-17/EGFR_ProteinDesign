# Histidine engineering, controls and validation

What was added to the generated binders, and the three independent tests
applied before anything was claimed about it.

Depends on [ph_mechanism.md](ph_mechanism.md) and
[design_campaign.md](design_campaign.md).

---

## Histidine placement

BindCraft optimises AlphaFold2 confidence and interface quality. It has no
notion of pH, so a design that binds well at 7.4 is exactly what it is
built to produce, which is the opposite of objective 1. The switch is
added here, on interfaces that already exist.

A binder position qualifies when its CB faces a target ionizable group
within the span a histidine can cover, the position is exposed in the
unbound binder above 0.20 relative SASA, and the residue is chemically
substitutable. Cysteine, proline and glycine are protected. Positions
already carrying a histidine are kept as existing switches rather than
discarded, since they supply a protonation event at no cost in
substitutions.

Substitutions are proposed in sets rather than singly. The thermodynamics
in [ph_mechanism.md](ph_mechanism.md) make the count the governing
variable, so variants engaging two or three handles rank above single
substitutions regardless of individual geometry.

### How far a histidine reaches

An earlier version allowed the CB to sit 4–10 Å from the handle. A rotamer
feasibility test then vetoed 33 of 33 proposals, and every veto was for
failing to reach, none for steric clash. The 10 Å ceiling was physically
impossible.

The span was re-derived by measurement rather than estimate. Across the 23
histidines resolved in 6ARU:

| Distance | Mean | Range |
|----------|-----:|-------|
| CB–ND1 | 2.53 Å | 2.52–2.54 |
| CB–NE2 | 3.67 Å | 3.67–3.68 |

With a 4.0 Å charge-assisted contact between a protonated imidazole and a
carboxylate, the ceiling is about 7.7 Å from CB. The reported selectivity
of the best variant fell from 86-fold to 55-fold once the bound was
corrected, and the earlier figure is understood as an artefact of
proposing contacts no rotamer can form.

### Rotamer feasibility

Each proposal is tested by placing an imidazole nitrogen at 200
orientations on spheres of both measured radii around CB, discarding
orientations that clash with the binder or the target, and accepting the
substitution only if a surviving orientation lands 2.6–4.0 Å from the
handle. The environment excludes the side chain being replaced: an earlier
version left it in, which asked whether two side chains could occupy one
position and vetoed everything.

Of 23 substitutions across 16 proposed variants, 18 are feasible and 11
variants survive intact.

## Negative control

The question that decides whether the placement rule is worth anything:
does it beat placing the same number of histidines at random among the
positions that were available?

The null is deliberately the hard version. It draws from exactly the pool
the method could have chosen from — interface, exposed, substitutable —
removing only the targeting. A null scattering histidines anywhere on the
binder would be easy to beat, and beating it would show only that
interface positions differ from the rest.

Across 3000 draws per design and substitution count:

| | Designed | Random |
|---|---------:|-------:|
| Mean coupling | 1.82 kcal/mol | 0.40 kcal/mol |

Between a third and four fifths of random placements engage no handle at
all, depending on the design, and that is where the separation comes from.
Five of the eight submitted engineered variants reach p < 0.05.

Two caveats are part of the result. No designed variant exceeds the
maximum of its own null distribution, so chance does occasionally find
placements as good and the claim is about frequency, not about reaching
something unreachable. And the p-values are uncorrected: variants within a
design share positions and are not independent, which makes a clean
correction awkward; under a conservative Bonferroni threshold the
strongest survive and the marginal ones do not.

Both arms are scored by the same model with the same assumed pKa shift.
This tests the placement rule under that model; it does not validate the
model.

## Structural revalidation

Introducing charged side chains into a packed interface can make it worse.
A variant with an elegant protonation network and a broken interface is
worth nothing, so the network score is never merged with the structural
scores. Three layers are reported separately: folding plausibility,
interface quality, and whether the designed network survives.

Every variant and its unmodified parent were re-predicted from sequence
against domain III under identical settings. The parent is what makes the
comparison meaningful: without it, a change in interface confidence could
not be separated from prediction noise.

| Variant | kcal/mol | i_pTM | Parent | Δ | Verdict |
|---------|---------:|------:|-------:|--:|---------|
| H7 M10H H11 | 2.38 | 0.780 | 0.778 | **+0.002** | retained |
| H7 H11 *(native)* | 1.79 | 0.778 | 0.778 | +0.000 | retained |
| W24H M25H | 1.79 | 0.742 | 0.842 | −0.100 | degraded |
| Q23H | 0.94 | 0.752 | 0.758 | −0.006 | retained |
| F59H | 0.94 | 0.872 | 0.877 | −0.005 | retained |
| W63H | 0.94 | 0.774 | 0.810 | −0.036 | retained |
| W24H | 0.94 | 0.774 | 0.842 | −0.068 | degraded |
| M25H | 0.85 | 0.806 | 0.842 | −0.036 | retained |
| M55H | 0.85 | 0.653 | 0.763 | −0.110 | degraded |
| R49H | 0.59 | 0.243 | 0.819 | **−0.576** | destroyed |

**Damage does not track the number of substitutions.** The variant
engaging three handles loses nothing, while a single substitution
elsewhere collapses the complex and takes pLDDT down 0.25 with it. What
matters is which residue is replaced. R49H is excluded.

**The cheapest variants are those built on histidines the design already
carried**, where only one position actually changes. That is a design
principle rather than a coincidence: prefer interfaces already presenting
a titratable residue, and add the minimum.

### One network was not engineered at all

Step 05 writes a pre-existing histidine as `H7` and a substitution as
`M10H`. One variant, H7/H11, consists only of pre-existing histidines: its
sequence is the parent's and it changes nothing.

It is reported as a native network rather than counted as engineered.
BindCraft optimises interface confidence and has no notion of pH, so a
design arriving with a two-histidine network engaging D368 and H433 is
incidental, and counting it as engineered would claim credit the pipeline
has not earned. It is submitted regardless, because an incidental network
is as real as a designed one.

### Caveat

These binders were generated by AlphaFold2 and are re-predicted by
AlphaFold2, so absolute confidence values are not independent evidence of
binding. The variant-versus-parent comparison is sound, since no variant
was optimised and both are predicted identically, which isolates the
effect of the substitution. Only the assay can establish binding.

## Novelty

The competition filters out designs lacking adequate sequence and
structural diversity from known proteins. These binders come from
hallucination with no starting binder, so novelty was expected, but an
expectation is not evidence.

No design returns a hit from the RCSB sequence service, which runs MMseqs2
over every polymer entity in the PDB. The highest pairwise identity within
the portfolio is 27%.

Sequence identity is a weak proxy for structural novelty — two proteins
can share a fold at 15% identity — and no structural comparison such as
Foldseek was run.

### Composition

The composition check found what the PDB search could not. Every design is
strongly acidic:

| Design | E% | D% | Net charge |
|--------|---:|---:|-----------:|
| p2g0_l66 | 15.2 | 9.1 | −9 |
| p2g0_l73 | 15.1 | 8.2 | −5 |
| p2g0_l74 | 13.5 | 5.4 | −1 |
| p2g1_l63 | 17.5 | 4.8 | −2 |
| **p2g1_l68** | **27.9** | 0.0 | −5 |

Natural proteins average about 7% glutamate. The most extreme case is the
design carrying the lead variant, so the most promising candidate is also
the most compositionally abnormal. This is characteristic of hallucinated
sequences rather than of any one trajectory, and it is an expression and
solubility liability for the whole portfolio, since these have to be made
in a laboratory. It is weighed in the submission ranking but not solved.
