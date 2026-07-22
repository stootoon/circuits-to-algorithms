---
title: S-04 Manifold capacity
parent: Structures
nav_order: 4
---

# S-04 — Manifold capacity: statistical mechanics when the objects are not points
{: .no_toc }

> **The object.** The capacity $\alpha_M$ for linearly separating $P$ *manifolds* in $N$
> dimensions, together with two emergent order parameters — the **effective radius** $R_M$ and
> **effective dimension** $D_M$ — and the **anchor points** that generate them.
> **Born from.** The observation that a neural representation of an "object" is not a vector but
> a manifold swept out by nuisance variation (pose, concentration, speed, plume dynamics).
> **Mathematical home.** Replica theory of disordered systems; high-dimensional convex geometry.
> **Situation.** #3, the capacity problem.
> **Novelty.** **Genuinely new.** A real extension of Gardner's calculation, with new order
> parameters that had no analogue in the point case.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Where it comes from

Gardner's calculation is one of the jewels of statistical mechanics: how many random points can
a perceptron in $N$ dimensions separate with random labels? The answer is $P/N \to \alpha = 2$
(no margin, points in general position), obtained by a replica computation of the volume of
weight space.

But nothing in the brain is a point. A given odour at a given concentration in a given plume
produces not a PN activity vector but a *set* of vectors — a manifold parameterized by nuisance
variables you do not care about and cannot control. The computational problem the animal faces
is not "separate these points" but "separate these manifolds," and the geometry of the manifolds
— how big, how curved, how high-dimensional — is precisely what a sensory circuit is in the
business of changing.

So the question "how many objects can this representation support?" needs a theory for
manifolds. Gardner does not give you one, and neither does anything else on the shelf.

---

## 2. The object, precisely

### Setup

Let $\mathcal{M}^1,\dots,\mathcal{M}^P \subset \mathbb{R}^N$ be manifolds with random labels
$y^\mu = \pm 1$. We ask for a weight vector $w$ (norm 1) achieving margin $\kappa$ on *every
point of every manifold*:

$$y^\mu \, (w \cdot x) \geq \kappa \qquad \forall x \in \mathcal{M}^\mu, \ \forall \mu.$$

The capacity is the largest $\alpha_M = P/N$ for which a solution exists with probability 1 as
$N \to \infty$.

### The point limit, for calibration

For points ($\mathcal{M}^\mu = \{x^\mu\}$, i.i.d. Gaussian), the classical Gardner result is

$$\alpha^{-1}(\kappa) = \int_{-\kappa}^{\infty} Dt\ (t+\kappa)^2, \qquad
Dt = \frac{e^{-t^2/2}}{\sqrt{2\pi}}\,dt,$$

giving $\alpha(0) = 2$. Memorize this; it is the yardstick against which every manifold result
is read.

### The new ingredient: anchor points

The technical obstruction to extending this is that the constraint is now a *continuum* of
inequalities per manifold. The resolution is the key idea:

> For a given $w$, only one point of each manifold matters — the point that minimizes the
> margin. Call it the **anchor point** $\tilde{x}^\mu(w)$.

The anchor point is a support-vector-like object, but it is $w$-dependent, and in the replica
calculation $w$ is itself a fluctuating quantity being integrated over. So the anchor point
becomes a *self-consistently determined* random variable, and the theory must average over its
distribution. That is where the new structure lives.

Carrying this through gives self-consistent equations whose solution is summarized by two
scalars:

$$R_M = \text{effective radius}, \qquad D_M = \text{effective dimension},$$

defined as moments of the anchor-point distribution rather than as naive geometric properties of
the manifold. This distinction is the substance of the result. $D_M$ is *not* the manifold's
topological or PCA dimension, and $R_M$ is *not* its diameter — they are the dimension and
radius *as seen by a random linear classifier*, which weights the manifold's extent by how
often each direction is the binding constraint. A manifold with a few long thin spikes has small
$D_M$ despite large ambient extent; a fat isotropic ball has $D_M$ close to its true dimension.

The capacity is then, to good approximation, a function of these two numbers alone:

$$\alpha_M \approx \alpha\big(R_M, D_M; \kappa\big),$$

decreasing in both, and recovering $\alpha \to 2$ as $R_M \to 0$ (the point limit). The
compression of "all the geometry of $P$ manifolds" into two effective scalars is the practical
payoff, because it says what a circuit must do to increase capacity: **shrink $R_M$, shrink
$D_M$, or both.** That is a level-1 objective stated in geometric terms, and it is directly
measurable from data.

### Why this is usable

Given recorded population responses to many objects under many nuisance conditions, you can
estimate $R_M$, $D_M$ and $\alpha_M$ empirically and ask how they change across a processing
hierarchy. Applied to deep networks and the ventral visual stream, the finding is that
successive stages progressively reduce manifold radius and dimension — an untangling, made
quantitative. This turned a metaphor ("untangling," DiCarlo & Cox) into a measurement.

---

## 3. Why the neuroscience forced it

Machine learning had margins and support vectors, and statistical physics had Gardner. Neither
had the case where each *class* is a continuous set whose geometry is the object of study.
Machine learning did not need it because in ML the data are points; the manifold structure is
implicit in the sample. Neuroscience needed it because the manifold structure is the *thing the
circuit is manipulating* — sensory processing is largely the business of reshaping object
manifolds — and so it must appear explicitly in the theory.

This is a common and underrated way that new math appears from biology: the *dependent variable*
in biology is the *nuisance parameter* elsewhere.

---

## 4. How to recognize the pattern elsewhere

You are in Situation #3 when: **an existing theory answers your question for objects of type A,
and your objects are of type B.** The instinct is to approximate B by A (sample points from the
manifold and use Gardner). That instinct is wrong when the difference between A and B is the
scientifically interesting part — and here it is, because the whole point is how the manifold's
shape changes.

The diagnostic question: *if I approximate my objects as points, does the phenomenon I care
about disappear?* If yes, you need new theory, not a better approximation.

For olfaction, the application is sharp and I think underexploited. A concentration series of a
single odorant traces a curve in PN space; adding plume dynamics and background odours makes it
a manifold. Concentration invariance — one of the central unsolved problems in olfactory coding
— is precisely the statement that the circuit shrinks $R_M$ and $D_M$ along the concentration
direction. Manifold capacity gives you the metric to ask *how much* invariance the antennal lobe
achieves, and to compare PN and KC layers on a common scale. Given that
[C1](../part2-case-studies/C1-expansion-and-sparsening.md) says the KC layer expands
dimensionality, and capacity theory says expansion helps only if $D_M$ does not grow with it,
this is a genuinely non-trivial prediction with an unclear answer.

---

## 5. Exercises

**Ex S4.1 ★** — Verify $\alpha(0) = 2$ from the Gardner formula.

<details markdown="1"><summary>Solution</summary>

$$\alpha^{-1}(0) = \int_0^\infty Dt\ t^2 = \frac{1}{2}\int_{-\infty}^{\infty} Dt\ t^2
= \frac{1}{2}\,\mathrm{Var}(t) = \frac{1}{2},$$

using symmetry of the Gaussian and $\mathbb{E}[t^2]=1$. So $\alpha(0) = 2$. $\square$

Sanity check against the classical Cover count: the number of dichotomies of $P$ points in
general position in $\mathbb{R}^N$ realizable by a homogeneous hyperplane is
$2\sum_{k=0}^{N-1}\binom{P-1}{k}$, which equals $2^P$ (all dichotomies) iff $P \leq N$, and
drops through $1/2$ of all dichotomies exactly at $P = 2N$. Same answer, different route —
Cover (1965) is the combinatorial version of Gardner's statistical-mechanical one.
</details>

**Ex S4.2 ★★** — Show that $D_M$ as "dimension seen by a random classifier" differs from PCA
dimension, using a concrete manifold: a $D$-dimensional ellipsoid with one long axis of length
$L$ and $D-1$ short axes of length $\ell \ll L$.

<details markdown="1"><summary>Solution</summary>

Take $w$ uniform on the sphere in $\mathbb{R}^N$ and consider the manifold's extent along $w$.
For an ellipsoid with semi-axes $a_1 = L, a_2=\dots=a_D=\ell$, the support function in direction
$w$ is $h(w) = \sqrt{\sum_i a_i^2 w_i^2}$, so the anchor point is at the ellipsoid boundary in
the direction that maximizes projection, and

$$\mathbb{E}_w[h(w)^2] = \frac{1}{N}\sum_i a_i^2 = \frac{L^2 + (D-1)\ell^2}{N}.$$

If $L \gg \ell\sqrt{D}$, this is dominated by the single long axis: the manifold "looks
one-dimensional" to a random classifier, and the anchor point lies at one of the two ends of the
long axis with probability approaching 1. Effective dimension $\to 1$, effective radius set by
$L$.

PCA dimension, by contrast: the participation ratio $\big(\sum a_i^2\big)^2 / \sum a_i^4$ also
$\to 1$ in this limit — so here they agree. To separate them, take $L = \ell\sqrt{D}$ so all
axes contribute equally to variance: PCA participation ratio is $\approx D$, but the anchor-point
distribution concentrates on the *extremes*, and the capacity-relevant dimension is governed by
the tail behaviour of $\max_i$ rather than the bulk variance. Curvature makes this sharper still:
a spiky manifold (high curvature at extremes) has anchor points pinned to the spikes, so $D_M$
counts *spikes*, not variance directions. A hypercube's corners and its faces give the same
covariance but very different $D_M$.

The moral: **capacity cares about extremes, PCA cares about variance.** Any time you summarize a
neural manifold by explained variance, you are answering a different question from the one a
downstream linear readout is asking.
</details>

**Ex S4.3 ★★★** — *(Situation, not object.)* Identify a theory in your area that solves your
problem for the wrong type of object. State the type mismatch precisely, and say what the new
order parameters would have to describe.

<details markdown="1"><summary>Solution</summary>

Open. A worked example in the spirit of the note, for temporal codes — which is where I would
push it for [C2](../part2-case-studies/C2-transient-synchrony.md):

*Existing theory for wrong object.* Manifold capacity itself assumes the readout is an
instantaneous linear classifier acting on a static population vector. But a synchrony-based
olfactory code is read out by *coincidence detectors integrating over a window* — the object
being classified is a trajectory segment, not a point or even a static manifold.

*Type mismatch.* Class = a set of *trajectories* $\{x(t)\}$ over a window $[0,T]$, with nuisance
variation in both the trajectory's shape and its *timing*. The readout is a spatiotemporal
filter, so the relevant separability question is about a set in function space under a group
action that includes time-warping.

*What new order parameters would describe.* Plausibly three: an effective radius (as before), an
effective dimension of the *trajectory bundle*, and a new one with no static analogue — an
effective *temporal extent* or "warp tolerance," measuring how much of the margin is consumed by
timing jitter rather than by shape variation. The interesting physics would be the tradeoff
between the last two: a code can buy timing robustness by spending trajectory dimension, and the
optimal allocation would depend on the jitter statistics of the biology.

I do not know whether this has been done. If it has not, it is a thesis.
</details>

---

## 6. Reading

- **Gardner (1988)**, *The space of interactions in neural network models* (J. Phys. A) — read
  it for: the original calculation. You do not need to reproduce the replica algebra, but you
  should understand what quantity is being computed (the volume of solution space) and why it
  vanishes at capacity.
- **Cover (1965)**, on geometrical and statistical properties of linear inequalities — read it
  for: the combinatorial route to $\alpha = 2$, which is more intuitive than the replica route.
- **Chung, Lee & Sompolinsky (2018)**, *Classification and geometry of general perceptual
  manifolds* (Phys. Rev. X) — read it for: anchor points, the self-consistent equations, and
  $R_M$/$D_M$. This is the paper; it is long and worth it.
- **Cohen, Chung, Lee & Sompolinsky (2020)**, on separability and geometry of object manifolds
  in deep networks (Nature Communications) — read it for: the empirical method — how to
  actually estimate $\alpha_M$, $R_M$, $D_M$ from recorded or simulated activity — and the
  untangling-across-a-hierarchy result.
- **DiCarlo & Cox (2007)**, *Untangling invariant object recognition* (Trends Cogn. Sci.) —
  read it for: the metaphor that this theory made quantitative. Read it *after*, to appreciate
  the difference between a picture and a measurement.
