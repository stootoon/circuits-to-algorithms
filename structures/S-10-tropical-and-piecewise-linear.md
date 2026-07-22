---
title: S-10 Tropical & piecewise-linear geometry
parent: Structures
nav_order: 10
---

# S-10 — Hyperplane arrangements and tropical geometry: the combinatorics hiding inside ReLU

> **The object.** The chamber structure of the hyperplane arrangement induced by threshold
> nonlinearities — and its algebraic counterpart, the identification of ReLU networks with
> **tropical rational functions**.
> **Born from.** Threshold-linear dynamics, where "which units are active" is a combinatorial
> variable governing a locally linear system.
> **Mathematical home.** Polyhedral and tropical geometry; matroid theory.
> **Situation.** #4 and #8 — it is the geometric substrate underneath
> [S-03](S-03-combinatorial-threshold-linear-networks.md).
> **Novelty.** Tropical geometry is classical; the synthesis with neural dynamics and the
> resulting counting results are new, and they explain *why* expansion works.

---

## 1. Where it comes from

Rectification is the most consequential nonlinearity in both neuroscience and machine learning,
and it has a special property: it is **piecewise linear**. That is usually treated as an
analytical inconvenience. It is actually a gift, because piecewise linearity converts a
nonlinear problem into a *combinatorial* one plus a family of *linear* ones.

Once you see this, several apparently unrelated things become the same thing: why CTLN graph
rules exist ([S-03](S-03-combinatorial-threshold-linear-networks.md)), why depth beats width in
deep networks, and why expansion-plus-threshold increases separability
([C1](../part2-case-studies/C1-expansion-and-sparsening.md)).

---

## 2. The object, precisely

### Chambers

Consider $N$ threshold-linear units with input $u_i(x) = \sum_j W_{ij}x_j + \theta_i$ and output
$[u_i]_+$. Each unit defines a hyperplane $H_i = \{x : u_i(x) = 0\}$ in state space. The
arrangement $\mathcal{A} = \{H_1,\dots,H_N\}$ partitions $\mathbb{R}^d$ into **chambers**, and
within each chamber the *activity pattern* $\sigma = \{i : u_i(x) > 0\}$ is constant — so the
dynamics are **affine**.

Nonlinear dynamics thus decompose into: a combinatorial object (which chamber) plus linear
algebra (the affine flow within it). Fixed points are found by solving the linear system in each
chamber and checking the solution actually lies in that chamber — precisely the *positivity +
consistency* conditions of [S-03](S-03-combinatorial-threshold-linear-networks.md). The graph
rules are theorems about which chambers can contain their own linear solution.

**Counting chambers.** Not all $2^N$ sign patterns are realizable. For $N$ hyperplanes in general
position in $\mathbb{R}^d$, the number of chambers is

$$R(N,d) = \sum_{k=0}^{d} \binom{N}{k},$$

Zaslavsky's theorem in its simplest form. For $N \gg d$ this is $\Theta(N^d)$ — polynomial in the
number of units, exponential in the dimension.

This one formula does a lot of work. It is the same sum that appears in **Cover's function
counting theorem** for linear dichotomies, and it is the reason expansion helps: mapping a
$d$-dimensional input through $N \gg d$ thresholded units gives access to $\sum_{k\le d}
\binom{N}{k}$ distinct activity patterns, and a downstream linear readout can assign each
chamber its own output. Expansion buys you chambers; chambers buy you separability.
[C1](../part2-case-studies/C1-expansion-and-sparsening.md) is this fact wearing a biological hat.

### The tropical identification

Now the algebra. Work in the **tropical semiring** $(\mathbb{R}\cup\{-\infty\}, \oplus, \odot)$
with

$$a \oplus b = \max(a,b), \qquad a \odot b = a + b.$$

A **tropical polynomial** in $x_1,\dots,x_d$ is a max of finitely many affine functions,

$$p(x) = \max_{\alpha} \big(\langle \alpha, x\rangle + c_\alpha\big),$$

which is exactly a convex piecewise-linear function. A **tropical rational function** is a
difference $p - q$ of two tropical polynomials, which is exactly a general (not necessarily
convex) piecewise-linear function.

**Theorem (Zhang, Naitzat & Lim).** A feedforward ReLU network with integer-or-rational weights
computes a tropical rational function; conversely every tropical rational function is
representable by such a network. The correspondence is exact.

What this buys you is that the machinery of tropical geometry applies. The **Newton polytope**
of $p$ encodes the linear regions; the number of linear regions relates to the number of
vertices of the Newton polytope; composition of layers corresponds to *Minkowski sums* of
polytopes, and Minkowski sums multiply vertex counts. That is the clean explanation of the
depth result:

**Depth is exponentially more efficient than width at carving space.** A network of depth $L$ and
width $n$ can realize $\Omega\big((n/d)^{d(L-1)} n^d\big)$ linear regions, exponential in depth,
whereas a shallow network needs exponentially many units to match. Montúfar, Pascanu, Cho &
Bengio proved the counting; the tropical picture explains it — each layer takes a Minkowski sum,
and vertex counts multiply.

---

## 3. Why the neuroscience forced it

Partly it did not — the deep-learning community drove the tropical connection. But the
*dynamical* half is neuroscience's, and it is the half that matters here: recurrent
threshold-linear networks, where the chamber structure interacts with the flow. The question
"which chambers can host a fixed point of their own affine dynamics?" is not a machine-learning
question, and it is the question [S-03](S-03-combinatorial-threshold-linear-networks.md) answers.

The general lesson is worth extracting: **a nonlinearity that is piecewise linear hands you a
combinatorial skeleton for free.** Most theorists reach for smooth approximations (sigmoid,
softplus) because smoothness is familiar. That reflex throws away the skeleton. Sometimes the
right move is to keep the kink.

---

## 4. How to recognize the pattern elsewhere

Look for **piecewise structure with a finite index set**: rectification, hard thresholds,
winner-take-all, max-pooling, spiking (which is a threshold), and any "regime-switching" model.
In every such case there is a combinatorial variable — which piece — and a family of tractable
problems indexed by it.

The recipe: (1) name the combinatorial variable; (2) count its realizable values; (3) characterize
which values are consistent with the dynamics; (4) look for local rules that decide consistency.

Note that $k$-winner-take-all — the mushroom body operation of
[C1](../part2-case-studies/C1-expansion-and-sparsening.md) and the $k$-cap of
[S-09](S-09-assembly-calculus.md) — has its own arrangement: the regions where a *fixed* set of
$k$ units is the top-$k$ are the chambers of the arrangement of the $\binom{N}{2}$ hyperplanes
$\{u_i = u_j\}$, which is the **braid arrangement**, whose combinatorics is classical and rich
(its chambers are the $N!$ orderings; its intersection lattice is the partition lattice). So the
combinatorics of a sparse olfactory code is the combinatorics of the braid arrangement
restricted to a $d$-dimensional subspace. I have not seen this exploited and it looks like it
should be.

---

## 5. Exercises

**Ex S10.1 ★** — Verify Zaslavsky's count numerically for random hyperplanes in $\mathbb{R}^2$
and $\mathbb{R}^3$.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
from math import comb

def n_regions(N, d, n_samples=400000, rng=np.random.default_rng(0)):
    W = rng.standard_normal((N, d)); b = rng.standard_normal(N)
    X = rng.standard_normal((n_samples, d)) * 3
    S = (X @ W.T + b > 0)
    return len(np.unique(S, axis=0))

for d in (2,3):
    for N in (3,5,8):
        print(d, N, n_regions(N,d), sum(comb(N,k) for k in range(d+1)))
```

Sampled counts approach $\sum_{k\le d}\binom{N}{k}$ (undercounting slightly, since Monte Carlo
misses tiny regions — increase samples or restrict the sampling range to see convergence).
For $d=2, N=5$: $1+5+10 = 16$ regions. For $d=3, N=8$: $1+8+28+56 = 93$.

The point to feel: with $N=8$ units in $d=3$ you get 93 regions, not 256. The gap between $2^N$
and $\sum_{k\le d}\binom{N}{k}$ is exactly the constraint that the input is low-dimensional, and
it is the reason a high-dimensional code driven by low-dimensional input has far less capacity
than its dimensionality suggests. This is a fact worth carrying into
[C1](../part2-case-studies/C1-expansion-and-sparsening.md).
</details>

**Ex S10.2 ★★** — Show that $\max(0, x)$ is a tropical polynomial and that a one-hidden-layer
ReLU network computes a tropical rational function.

<details markdown="1"><summary>Solution</summary>

$\max(0,x) = 0 \oplus x$ — a tropical polynomial with two terms (the constant $0$, and the
monomial $x$ with coefficient $0$). Convex piecewise linear, as required.

For a one-hidden-layer network $f(x) = \sum_i c_i [\,w_i \cdot x + b_i\,]_+$ with $c_i \in
\mathbb{R}$: split the sum into positive and negative coefficients,

$$f(x) = \underbrace{\sum_{i : c_i > 0} c_i [w_i\cdot x + b_i]_+}_{=:P(x)}
\;-\; \underbrace{\sum_{i : c_i < 0} |c_i| [w_i \cdot x + b_i]_+}_{=:Q(x)} .$$

Each of $P$ and $Q$ is a nonnegative combination of convex piecewise-linear functions, hence
convex piecewise linear, hence a tropical polynomial (tropical addition is $\max$, and a sum of
maxes of affine functions is again a max of affine functions — this is where the tropical
*product* $\odot = +$ does the work, since $\sum_i \max_j(\cdot)$ expands into a max over
tuples). So $f = P - Q$ is a tropical rational function. $\square$

The corollary that makes it useful: the linear regions of $f$ are determined by the common
refinement of the regions of $P$ and $Q$, and the Newton polytope of the product is the
Minkowski sum of the Newton polytopes. Iterating over layers multiplies complexity, which is
the depth theorem.
</details>

**Ex S10.3 ★★★** — *(Situation, not object.)* The set of inputs for which a *fixed* set of $k$
Kenyon cells is the top-$k$ is a chamber of a restricted braid arrangement. Work out the
implications for the coding capacity of the mushroom body, and compare with the LSH picture in
[C1](../part2-case-studies/C1-expansion-and-sparsening.md).

<details markdown="1"><summary>Solution</summary>

Open-ended; here is the intended line of attack.

Let PN activity be $y \in \mathbb{R}^{d}$ ($d \approx 50$), KC input $u_i = a_i \cdot y$ with
$a_i$ the sparse random sampling vector of KC $i$ ($N \approx 2000$). The top-$k$ set is constant
on chambers of the arrangement $\{u_i = u_j\}_{i<j}$ restricted to the image of PN space — i.e.
the braid arrangement $\mathcal{A}_{N-1}$ pulled back along a linear map $\mathbb{R}^d \to
\mathbb{R}^N$.

*Capacity.* The number of distinct top-$k$ tags is at most the number of chambers of an
arrangement of $\binom{N}{2}$ hyperplanes in $\mathbb{R}^d$, i.e. at most $\sum_{j \le d}
\binom{\binom{N}{2}}{j} = O(N^{2d})$ — enormous, and *not* the binding constraint. The binding
constraint is instead noise: chambers must be large enough that jitter does not move you across
a boundary. So the right question is not how many chambers exist but the **distribution of
chamber sizes** under the relevant input distribution, which is a geometric probability question
about the induced arrangement.

*Comparison with LSH.* The LSH picture says the tag is locality-preserving: nearby $y$ give
overlapping tags. In the arrangement picture, "nearby $y$ give the same tag" is exactly "nearby
points lie in the same chamber," and the *degree* of locality preservation is governed by the
density of hyperplanes near a typical point, i.e. by how close the top-$k$ competition is. These
are the same statement, but the arrangement version is quantitative and gives you the tail
behaviour — the probability that a small perturbation flips a KC is the probability of being
within $\epsilon$ of a boundary, computable from the density of the $u_{(k)} - u_{(k+1)}$ gap.

*The payoff.* That order-statistic gap is directly measurable and directly manipulable
(it depends on $k$, on the sampling degree, and on PN correlation structure), so this route turns
"the mushroom body is an LSH" from a metaphor into a formula for the code's noise robustness —
and it predicts an optimal $k$ trading tag specificity against gap size. Compare with the optimal
*degree* result in [S-11](S-11-expanders-and-optimal-degree.md), which is the same style of
argument applied to the other free parameter.
</details>

---

## 6. Reading

- **Zaslavsky (1975)**, on facing up to arrangements — read it for: the chamber-counting theorem.
  You need the statement, not the proof.
- **Cover (1965)**, on geometrical and statistical properties of systems of linear inequalities
  — read it for: the same count in the guise of linear separability, and the cleanest
  explanation of why expansion helps.
- **Montúfar, Pascanu, Cho & Bengio (2014)**, *On the number of linear regions of deep neural
  networks* (NeurIPS) — read it for: the depth-vs-width counting result.
- **Zhang, Naitzat & Lim (2018)**, *Tropical geometry of deep neural networks* (ICML) — read it
  for: the exact correspondence and the Newton-polytope machinery.
- **Curto, Geneson & Morrison (2019)** — read it for: the dynamical half, where chambers meet
  fixed points. See [S-03](S-03-combinatorial-threshold-linear-networks.md).
- **Stanley**, *An Introduction to Hyperplane Arrangements* (lecture notes) — read it for: the
  braid arrangement, the intersection lattice, and the characteristic polynomial, if Ex S10.3
  appeals.
