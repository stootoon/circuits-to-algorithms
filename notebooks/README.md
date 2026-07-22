---
title: Notebooks
nav_order: 7
permalink: /notebooks/
---

# Companion notebooks

Three runnable models, one per flagship circuit. Each extends the code in its case study
with plots, a parameter sweep, and — the part that matters — a **null model** or an
**honest negative result**, because a model you cannot break teaches you nothing.

| Notebook | Circuit | What you build |
|---|---|---|
| [C1 — Expansion and LSH](https://github.com/stootoon/circuits-to-algorithms/blob/main/notebooks/C1-expansion-and-lsh.ipynb) | PN→KC, mushroom body | Sparse binary projection + WTA; verify the tag is locality-sensitive against $\Phi_2(t,t;\rho)/\alpha$; watch the PN layer hit Cover's bound while the KC layer sails past it; sweep the in-degree |
| [C2 — Transient synchrony](https://github.com/stootoon/circuits-to-algorithms/blob/main/notebooks/C2-transient-synchrony.ipynb) | Locust antennal lobe | PN–LN model with slow GABA and adaptation; reproduce progressive decorrelation and the Mazor–Laurent accuracy peak; compare against a spectral-radius-matched random reservoir |
| [C4 — Ring attractor](https://github.com/stootoon/circuits-to-algorithms/blob/main/notebooks/C4-ring-attractor.ipynb) | *Drosophila* ellipsoid body | Bump maintenance and angular path integration at the predicted gain $\kappa\delta/\tau$; then break the continuous symmetry and measure the drift |

Requirements: `numpy` and `matplotlib`. No scipy, no GPU, no data downloads. Each runs in
under a minute on a laptop.

```bash
pip install numpy matplotlib jupyter
jupyter lab notebooks/
```

---

## Two results worth knowing about before you run them

**C1 does not reproduce "six."** The Litwin-Kumar et al. (2017) result — that the optimal
KC in-degree matches the anatomical ~6 — is one of the best things theoretical neuroscience
has produced ([S-11](../structures/S-11-expanders-and-optimal-degree.md)). The toy sweep in
this notebook recovers the *qualitative* claim (an interior optimum exists, dimension falls
off in both directions) but peaks at $K\approx3$, because it lacks a noise model, the fly's
real PN correlation structure, and their learning-linked objective.

That gap is left in deliberately, and §4 explains it. A notebook tuned until it produced 6
would teach you the opposite of what
[S2](../part3-synthesis/S2-degeneracy-and-limits.md) is trying to teach you.

**C2's null is informative in the direction you might not expect.** A dense random recurrent
network at $g=0.95$ — strongly non-normal, maximal transient amplification, the most
favourable possible setting for the reservoir hypothesis — turns out **not** to decorrelate
the odours at all: it sits at the input correlation (~0.53) for the whole trial while the
structured PN–LN circuit drives it to ~0.16.

So decorrelation is not a generic consequence of recurrent dynamics; it requires the
structured, delayed, adapting inhibition. But note the sting in the tail: the null still
*classifies* as well or better, because it preserves input geometry rather than reshaping
it. **Decorrelation and discriminability are different objectives**, and much of the
olfaction literature slides between them without noticing.

---

## Where to take them

Each notebook ends with a "Things to try" section. The ones I would actually do:

- **C2, the slope test.** Regress asymptotic on initial correlation across a parametric
  family of odour pairs. Slope $\approx 0$ is genuine whitening; slope $\approx 1$ is mere
  rescaling. This converts "the antennal lobe decorrelates" from an observation into a
  falsifiable claim, and it is the analysis the capstone
  ([S3](../part3-synthesis/S3-capstone.md)) is built around.
- **C2, the heteroclinic competitor.** Implement the winnerless-competition model from
  C2 §2.3 and check its unique signature: dwell times growing like
  $\lambda_u^{-1}\ln(1/\varepsilon)$, linear in log noise. Neither the PN–LN circuit nor the
  reservoir null does that, which makes it the sharpest discriminator anyone has proposed.
- **C4, certify the topology.** Run persistent homology on the bump states and recover
  $\beta_0=\beta_1=1$ from unlabelled activity — the method of
  [S-05](../structures/S-05-toroidal-topology-of-grid-cells.md), on a system where you know
  the ground truth because you built it. This is the right way to learn a technique before
  trusting it on data.
