---
title: S-05 Toroidal topology of grid cells
parent: Structures
nav_order: 5
---

# S-05 — Topology as an instrument: certifying an attractor manifold
{: .no_toc }

> **The object.** Persistent homology deployed as *inference about state-space topology* — using
> Betti numbers computed from unlabelled population activity to certify that the manifold a
> theory predicts is the manifold the circuit actually uses.
> **Born from.** Twenty years of continuous-attractor theory that made a topological prediction
> nobody could test.
> **Mathematical home.** Algebraic topology (persistent homology), applied as statistics.
> **Situation.** #6, the empirical topology problem.
> **Novelty.** The topology is classical; the *use* — as a decisive test between mechanistic
> hypotheses — is new, and it is the reason this note exists.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Where it comes from

Continuous attractor network theory ([C4](../part2-case-studies/C4-continuous-attractors.md),
[C6](../part2-case-studies/C6-grid-cells.md)) makes a prediction of a very unusual kind. It does
not predict a tuning curve or a firing rate. It predicts a **shape**: the set of states the
population can occupy forms a manifold with a specific topology, determined by the symmetry of
the connectivity.

- Head-direction cells: a circle, $S^1$. Betti numbers $(\beta_0,\beta_1) = (1,1)$.
- Grid cells within one module: a twisted torus, $\mathbb{T}^2$. Betti numbers
  $(\beta_0,\beta_1,\beta_2) = (1,2,1)$.

This is a strong, falsifiable, and *purely topological* claim. And for two decades it was
effectively untestable, because testing it requires (i) many neurons recorded simultaneously
from a single module, and (ii) a way to determine the topology of a point cloud in
100-dimensional space with no access to the underlying coordinates.

The second requirement is the interesting one. You cannot just plot it. You cannot use the
animal's actual head direction or position, because that would be circular — the whole question
is whether the *intrinsic* dynamics have this structure, including when the animal is asleep and
there is no position to speak of.

You need a coordinate-free, label-free certificate of shape. That is what homology is.

---

## 2. The object, precisely

Given population activity vectors $\{x_1,\dots,x_T\} \subset \mathbb{R}^n$ (time points as
points in neuron space), the pipeline is:

1. **Denoise and reduce.** Smooth spike trains, then reduce dimension (PCA, or better, a method
   that respects the density). This step is where sins are committed; see §4.
2. **Build a filtration.** Vietoris–Rips on the point cloud: at scale $\epsilon$, include a
   simplex when all pairwise distances are below $\epsilon$.
3. **Compute persistent homology.** Track when homology classes are born and die as $\epsilon$
   grows; produce a barcode.
4. **Read the persistent Betti numbers.** Bars that survive over a long range of $\epsilon$ are
   the real topology; short bars are noise.
5. **Compare against nulls.** Shuffle in ways that destroy the hypothesized structure but
   preserve firing statistics, and confirm the long bars vanish.

For a torus you expect: one long bar in $H_0$ (connected), **two** long bars in $H_1$ (the two
independent cycles), one long bar in $H_2$ (the void). For a circle: one in $H_0$, one in $H_1$,
none in $H_2$. These signatures are different enough that the test has real power.

### The result

Gardner and colleagues (2022) recorded large populations from a single grid module with
Neuropixels probes and recovered exactly the toroidal signature $(1,2,1)$. Decisively, the
torus persisted during **REM and slow-wave sleep** and in darkness, when there is no spatial
input. The state space is intrinsic to the circuit.

Chaudhuri, Gerçek, Pandey, Peyrache & Fiete (2019) had done the corresponding job for the
head-direction system: the ring topology recovered from unlabelled activity, persisting into
sleep, with the network traversing the ring during REM.

Take a moment on what these papers accomplish. A theoretical prediction from the 1990s, about
the *shape of a state space*, was confirmed by measuring the shape of the state space. Not by a
proxy, not by a correlation — by measuring the actual predicted object. That is as close to
physics as systems neuroscience gets, and it is the model for what a mature circuit→algorithm
claim should look like.

---

## 3. Why the neuroscience forced it

Persistent homology existed. What did not exist was the practice of treating a *topological
invariant as the primary experimental observable*.

The forcing came from a specific structural feature of the hypothesis: continuous attractor
theory's core claim is invariant under everything you can easily measure. Any smooth
reparameterization of the manifold — any relabelling of which neuron prefers which direction,
any nonlinear distortion of the tuning curves — leaves the prediction unchanged. A hypothesis
whose content is invariant under a huge group demands an observable invariant under the same
group. That is Situation #1 and Situation #6 meeting: **the invariance structure of the
hypothesis dictates the mathematics of the measurement.**

Learn this. If your hypothesis is invariant under group $G$, then any analysis not invariant
under $G$ is measuring something other than your hypothesis.

---

## 4. How to recognize the pattern elsewhere — and the traps

Applicable whenever a mechanistic hypothesis predicts state-space *shape*: ring attractors,
tori, spheres (3-D head direction in bats), and — the open case —
[C2](../part2-case-studies/C2-transient-synchrony.md), where competing hypotheses about
olfactory transients predict genuinely different topologies. A stable heteroclinic sequence
through $k$ saddles predicts a cyclic structure (a loop, $\beta_1 = 1$) in the limit set; a pure
reservoir predicts no persistent topology at all, just a contractible blob whose *geometry*
carries the information; a fixed-point-per-odour code predicts $\beta_0 = $ number of odours
with all higher Betti numbers zero. **These are distinguishable by persistent homology, and as
far as I know the analysis has not been done on locust PN data.** That is a concrete, tractable
project and it is why this note is placed where it is in the course.

Now the traps, because this method is easier to misuse than to use:

- **Dimensionality reduction can create or destroy topology.** PCA to 3 dimensions can force a
  torus to self-intersect, killing $H_2$. Aggressive nonlinear embedding can manufacture loops.
  Do the homology in as high a dimension as you can afford, and verify stability across choices.
- **Density matters more than you expect.** Persistent homology of a non-uniformly sampled
  manifold produces spurious short bars and can miss real long ones. Density-normalizing
  preprocessing is often necessary and is itself a modelling choice.
- **Nulls are essential and rarely sufficient.** A shuffle that destroys temporal structure but
  preserves the covariance may still produce ring-like topology, because the covariance already
  encodes it. Design the null to break *specifically* the hypothesis, and state what it
  preserves.
- **Betti numbers are coarse.** $(1,2,1)$ is consistent with a torus, but also with other spaces;
  the inference is "consistent with," strengthened by the theory's prior prediction. This is
  fine — it is exactly how physics uses signatures — provided you do not overstate it.

---

## 5. Exercises

**Ex S5.1 ★** — Compute persistent homology of noisy samples from a circle and a torus, and
confirm the Betti signatures.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
from ripser import ripser   # pip install ripser

# circle in R^3 with noise
t = 2*np.pi*np.random.rand(400)
circ = np.c_[np.cos(t), np.sin(t), np.zeros_like(t)] + 0.05*np.random.randn(400,3)

# torus in R^3
u, v = 2*np.pi*np.random.rand(2, 1500)
R, r = 2.0, 0.8
tor = np.c_[(R+r*np.cos(v))*np.cos(u), (R+r*np.cos(v))*np.sin(u), r*np.sin(v)]
tor += 0.03*np.random.randn(*tor.shape)

for name, X in [("circle", circ), ("torus", tor)]:
    d = ripser(X, maxdim=2)['dgms']
    print(name, [sum((b[:,1]-b[:,0]) > 0.5*np.ptp(X)) for b in d])
```

Expect the circle to show one long $H_1$ bar and no $H_2$; the torus to show two long $H_1$ bars
and one long $H_2$ bar. Then reduce the torus to 2 dimensions with PCA and rerun: $H_2$ dies,
because a 2-D projection of a torus self-intersects. That single experiment teaches the main
practical lesson of the note.
</details>

**Ex S5.2 ★★** — The grid-cell manifold is a *twisted* torus, not a flat one. Explain why
persistent homology cannot distinguish these, what would, and why the inability does not damage
the inference.

<details markdown="1"><summary>Solution</summary>

Homology is a homotopy invariant, and the twisted torus (the quotient of the plane by a
hexagonal lattice) is homeomorphic to the flat torus $S^1 \times S^1$ — the "twist" refers to
the *metric/lattice* structure, i.e. the identification is by a non-rectangular lattice.
Homeomorphic spaces have identical Betti numbers, so no homological method can see the
difference.

What would see it: the *geometry*. Specifically, the ratio and relative angle of the two
generating cycles — measurable by estimating the induced metric on the recovered manifold (e.g.
by fitting a flat metric and extracting the lattice's fundamental domain), or by relating the
two toroidal coordinates to the animal's physical position and reading off the $60°$ structure.
Gardner et al. do this second step: after certifying the topology, they map the toroidal
coordinates to space and recover the hexagonal relationship.

Why it does not damage the inference: the *topological* prediction is the one that discriminates
between the competing mechanistic hypotheses. A continuous attractor with hexagonal symmetry
predicts a torus; a feedforward or oscillatory-interference account does not predict that
population activity is confined to a 2-D toroidal manifold at all, especially in sleep. The
topology carries the discriminating load; the geometry then confirms which torus. Correct
division of labour — use the coarse invariant for the decisive test, then refine.
</details>

**Ex S5.3 ★★★** — *(Situation, not object.)* Design a persistent-homology analysis to
discriminate the competing hypotheses in
[C2](../part2-case-studies/C2-transient-synchrony.md). Specify the point cloud, the null models,
the predicted signatures, and the failure modes.

<details markdown="1"><summary>Solution</summary>

*Point cloud.* PN population vectors binned at, say, 50 ms across many trials and many odours.
Two variants worth running separately: (i) pooled across odours — asks about the global
structure of the odour-response state space; (ii) per odour, pooled across trials — asks about
the structure of the transient itself.

*Predicted signatures, per hypothesis (variant ii):*

| Hypothesis | Predicted persistent topology of the per-odour trajectory bundle |
|---|---|
| Fixed-point attractor per odour | $\beta_0 = 1$, no persistent $\beta_1$; trajectories collapse |
| Stable heteroclinic sequence (cyclic) | persistent $\beta_1 = 1$: a loop through the saddle sequence |
| Stable heteroclinic sequence (open chain) | contractible, but with detectable *sequential* geometry — needs a different test |
| Generic reservoir / non-normal transient | contractible; information in geometry, not topology |

Note that two of the four predict "contractible," which is honest and important: **persistent
homology alone cannot separate the open-chain heteroclinic hypothesis from the reservoir
hypothesis.** That is a real limitation and identifying it is half the value of the exercise. To
separate those two you need a different instrument — the fixed/slow-point analysis of
[Unit 01](../part1-foundations/01-dynamical-systems.md), which asks whether the trajectory passes near
genuine saddles (heteroclinic) or through a region with no slow points (reservoir). That is
arguably the single most decisive analysis available for C2, and it is a natural capstone
deliverable.

*Nulls.* (a) Trial-shuffled PN identities preserving single-neuron statistics; (b) a
spectral-radius-matched random recurrent network driven by the same inputs — this is the
critical one, since it instantiates the reservoir hypothesis and gives the "generic dynamics"
baseline; (c) phase-randomized surrogates preserving the power spectrum, which controls for the
20–30 Hz oscillation producing spurious loops.

*Failure modes.* The oscillation is the biggest threat: a strong periodic component makes
*everything* look like a circle in $H_1$. Either analyze the oscillation-cycle-averaged
trajectory, or work in a coordinate system where the oscillation phase has been factored out —
and note that "factoring out the phase" is itself a hypothesis about the code (that phase is
gating, not content), so you must run it both ways and report both.
</details>

---

## 6. Reading

- **Gardner, Hermansen, Pachitariu, Burak, Baas, Dunn, Moser & Moser (2022)**, *Toroidal
  topology of population activity in grid cells* (Nature) — read it for: the result, and the
  methods section, which is a masterclass in doing this carefully.
- **Chaudhuri, Gerçek, Pandey, Peyrache & Fiete (2019)**, on the intrinsic attractor manifold of
  the head-direction circuit across wake and sleep (Nature Neuroscience) — read it for: the
  simpler ring case done first, and the sleep analysis that makes the intrinsic-dynamics
  argument.
- **Rybakken, Baas & Dunn (2019)**, on decoding of neural data using topological methods — read
  it for: the technical machinery of getting from spikes to a certified manifold.
- **Ghrist (2008)**, *Barcodes: the persistent topology of data* (Bull. AMS) — read it for: the
  best short introduction to persistence.
- **Burak & Fiete (2009)**, on accurate path integration in continuous attractor grid-cell
  models (PLOS Comput. Biol.) — read it for: the theory whose prediction was being tested. Read
  it *before* Gardner et al., so the confirmation lands properly.
