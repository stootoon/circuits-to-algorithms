---
title: S-02 Clique topology
parent: Structures
nav_order: 2
---

# S-02 — Clique topology: an invariant of matrix *order*
{: .no_toc }

> **The object.** The Betti curves of the order complex of a symmetric matrix — a topological
> invariant that depends only on the *rank ordering* of the matrix entries, hence is invariant
> under every monotone transformation of them.
> **Born from.** Hippocampal pairwise correlations, where the map from correlation to underlying
> distance is unknown and certainly nonlinear.
> **Mathematical home.** Persistent homology; random matrix theory.
> **Situation.** #1, the invariant problem.
> **Novelty.** **Genuinely new.** Persistent homology is classical; using it as an invariant of
> the *order type* of a matrix, to detect geometric organization through an unknown monotone
> nuisance, is not.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Where it comes from

You record $n$ neurons and compute the pairwise correlation matrix $C$. You suspect the neurons
have receptive fields sitting in some low-dimensional space, so that $C_{ij}$ is large when
fields $i$ and $j$ are close. But you have no idea what the function relating correlation to
distance is. It is monotone decreasing, presumably. Beyond that: unknown, nonlinear, and
dependent on firing rates, recording quality, behavioural state, and the definition of
correlation you happened to use.

So the *values* $C_{ij}$ mean nothing. Any analysis that uses them — spectral methods, MDS,
thresholding at a fixed value — is measuring your nuisance function as much as your biology.
This is Situation #1 in its purest form: you are working modulo the group of monotone increasing
transformations, and you need an invariant of that group.

What survives a monotone transformation? Exactly one thing: **the ordering of the entries.** So
the invariant you need is a function of the order type of $C$ and nothing else. Clique topology
is such a function, and it is a substantive one.

---

## 2. The object, precisely

Let $C$ be symmetric with distinct off-diagonal entries. Order them
$C_{i_1j_1} < C_{i_2j_2} < \cdots$ and build a nested sequence of graphs:

$$G_0 \subset G_1 \subset \cdots \subset G_{\binom{n}{2}}, \qquad
G_t = \big([n],\ \{\text{the } t \text{ smallest-valued edges}\}\big).$$

(If $C$ is a similarity rather than a dissimilarity, reverse the order.) Now replace each graph
by its **clique complex** (flag complex) $X(G_t)$: the simplicial complex whose $k$-simplices are
the $(k{+}1)$-cliques of $G_t$. This is the *order complex* filtration.

Compute persistent homology of the filtration and record the **Betti curves**

$$\beta_k(t) = \dim H_k\big(X(G_t)\big), \qquad k = 1,2,3,\dots$$

as functions of edge density $\rho = t/\binom{n}{2}$.

**The key property, which is the whole design:** $G_t$ depends only on which entries are among
the $t$ smallest. Applying any strictly increasing $f$ to every entry of $C$ leaves that
unchanged. Therefore

$$\beta_k^{f(C)}(t) = \beta_k^{C}(t) \quad \text{for all monotone increasing } f.$$

The Betti curves are an invariant of the order type. That is the trick, and it is worth pausing
on how cheap and how powerful it is.

### What the curves discriminate

The invariant is only useful if it separates hypotheses. It does:

- **Random (i.i.d.) symmetric matrix.** The order type is a uniformly random ordering. Betti
  curves are large, with a characteristic shape and known asymptotics from the theory of random
  clique complexes (Kahle).
- **Geometric matrix.** $C_{ij} = f(\|z_i - z_j\|)$ for points $z_i$ drawn i.i.d. in
  $\mathbb{R}^d$, $f$ monotone. Betti curves are *far smaller* and have a distinct shape,
  because the underlying space is contractible — the cycles that appear are finite-sampling
  artifacts, not real topology.

The gap between these two families is large and robust, so "geometric vs. not" is decidable from
Betti curves. Dimension $d$ is discriminable too, but much more weakly — the curves depend on
$d$ only mildly, which is a real limitation.

### The result

Applied to hippocampal place-cell correlations, the Betti curves fall firmly in the geometric
family rather than the random family. Crucially this holds not only during spatial exploration
— where it would be unsurprising, since the animal is moving through a real 2-D environment —
but also during *sleep and rest*, when there is no spatial input at all. The geometric
organization is intrinsic to the network, not imposed by the stimulus.

That is a genuinely strong inference, and note the form it takes: a conclusion about the
*existence of an underlying geometry* drawn from data in which the geometry was never measured
and the measurement scale was never calibrated.

---

## 3. Why the neuroscience forced it

Every standard method for detecting geometric structure — MDS, Isomap, spectral embedding,
UMAP — requires you to commit to a metric. In neural data that commitment is unjustifiable, and
worse, it is *consequential*: the answer changes with the commitment. The field's usual response
is to try several and report the one that looks best, which is not a method.

The honest alternative is to use only what you actually know, which is the ordering. Once you
frame it that way the object almost designs itself — but nobody framed it that way until a
concrete dataset made the nuisance function impossible to ignore.

---

## 4. How to recognize the pattern elsewhere

The signature phrase is **"up to an unknown monotone transformation."** Related, equally
diagnostic phrases: "the units are arbitrary," "only relative values are meaningful," "the gain
is unknown."

Whenever you catch yourself saying one of these, you are working with a group action, and the
right move is to construct an invariant rather than to normalize and hope. Clique topology is
the invariant for the monotone group acting on a symmetric matrix. Other groups want other
invariants — the general recipe is: identify the group, identify what is preserved, build a
functor out of the preserved structure.

For olfaction, the direct application is to PN or KC pairwise correlation matrices across an
odour panel. The claim "odour representations have geometric structure" is usually supported by
an MDS plot, which presupposes a metric. Clique topology tests it without one. The negative
result would also be informative: if the Betti curves look *random*, that is evidence against
any low-dimensional geometric odour space, which would be a substantive and surprising finding.
See also [S-06](S-06-hyperbolic-odor-space.md) — and note that clique topology as usually
deployed compares against *Euclidean* geometric nulls, so extending the null family to
hyperbolic point processes is a straightforward and unpublished piece of work.

---

## 5. Exercises

**Ex S2.1 ★** — Verify the invariance claim concretely. Take a random $10 \times 10$ symmetric
matrix, apply $f(x) = x^3$ and $f(x) = \tanh(5x)$, and confirm the order complex filtration is
identical. Then apply $f(x) = -x$ and explain what happens.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
from itertools import combinations

n = 10
A = np.random.randn(n, n); C = (A + A.T) / 2; np.fill_diagonal(C, 0)

def order_type(M):
    pairs = list(combinations(range(M.shape[0]), 2))
    return tuple(sorted(pairs, key=lambda p: M[p]))

print(order_type(C) == order_type(C**3))          # True: x^3 strictly increasing
print(order_type(C) == order_type(np.tanh(5*C)))  # True: tanh strictly increasing
print(order_type(C) == order_type(-C))            # False: order reversed
```

$x^3$ and $\tanh$ are strictly increasing on $\mathbb{R}$, so they preserve the ordering and
hence every $G_t$, hence every Betti curve. $f(x) = -x$ is strictly *decreasing*: it reverses
the filtration, turning the "build up from smallest" filtration into "build up from largest."
The Betti curves are not preserved — they become the curves of the complementary filtration.

This is why the similarity/dissimilarity convention must be fixed before you start: the
invariance is under the *increasing* monotone group, not the full monotone group. In practice
for a correlation matrix you filter by decreasing correlation, i.e. increasing dissimilarity.
</details>

**Ex S2.2 ★★** — Explain why the *clique* complex is used rather than, say, the Vietoris–Rips
complex of some metric or a Čech complex. What would go wrong with the alternatives?

<details markdown="1"><summary>Solution</summary>

Vietoris–Rips of a metric space at scale $\epsilon$ is defined by: include a simplex when all
pairwise distances are $< \epsilon$. That is *already* the clique complex of the
$\epsilon$-threshold graph. So Rips and clique complex agree — but Rips is normally described in
terms of a *metric and a scale*, which reintroduces exactly the nuisance we are trying to
eliminate. Framing it as the clique complex of the order filtration makes manifest that no
metric or scale is used, only the ordering. The mathematics is the same; the epistemics are not,
and the epistemics are the contribution.

Čech is worse for our purposes: it requires actual balls in an actual ambient space, i.e. it
requires the embedding we do not have. (Čech has better theoretical guarantees — the nerve lemma
gives it the correct homotopy type — but it is uncomputable here.)

There is also a computational reason. The clique complex is determined by the graph alone, so
the filtration is indexed by the $\binom{n}{2}$ edges and is completely combinatorial. No
geometric predicates, no numerical tolerances.
</details>

**Ex S2.3 ★★★** — *(Situation, not object.)* Your data are related to the quantity of interest
by an unknown transformation from some group $G$. For each of the following, name $G$ and
propose an invariant: (a) calcium fluorescence vs. spike rate; (b) population vectors recorded
on different days with different subsets of neurons; (c) inter-areal spike-count correlations
where each area has its own unknown gain.

<details markdown="1"><summary>Solution</summary>

(a) **$G$ = monotone increasing maps** (fluorescence is a saturating, kernel-smoothed, monotone
function of rate, with unknown parameters). Invariants: anything built from rank order —
Spearman rather than Pearson correlation; clique topology directly; rank-based mutual
information. Note that the temporal kernel also introduces a *convolution*, which is not in the
monotone group, so this is honestly only approximate — a good illustration that identifying $G$
correctly is the hard part, not building the invariant.

(b) **$G$ = the symmetric group acting on neuron identity, composed with subsampling.** Absolute
neuron labels are meaningless across days. Invariants: the spectrum of the covariance matrix
(permutation invariant), participation ratio, the *distribution* of pairwise correlations, and
topological summaries of the population manifold. This is exactly why
[S-05](S-05-toroidal-topology-of-grid-cells.md) uses persistent homology — it is a label-free
invariant, so it can compare wake to sleep with different sortable units. Note that
subsampling is not a group action (no inverse), which is why stability results here are
statistical rather than exact.

(c) **$G = (\mathbb{R}^+)^2$, independent positive scalings per area.** Correlation coefficients
are already invariant to per-neuron scaling, which is why people use them; but if the gain is
*nonlinear*, you are back in case (a). Genuine invariants: canonical correlation coefficients
between areas (invariant under any invertible linear map on either side, a much bigger group);
or mutual information (invariant under any invertible reparameterization). The lesson is that
each larger group buys robustness at the cost of discarding more structure, and the art is
choosing the smallest group you actually need.
</details>

---

## 6. Reading

- **Giusti, Pastalkova, Curto & Itskov (2015)**, *Clique topology reveals intrinsic geometric
  structure in neural correlations* (PNAS) — read it for: the method, the random-vs-geometric
  discrimination, and the hippocampal result including the sleep condition. This is the paper.
- **Curto (2017)**, *What can topology tell us about the neural code?* (Bull. AMS) — read it
  for: context, and the relationship to [S-01](S-01-convex-codes-and-the-neural-ideal.md).
- **Kahle**, work on the topology of random clique complexes — read it for: the asymptotics of
  Betti numbers of random flag complexes, which is what makes the "random" null quantitative
  rather than empirical.
- **Ghrist (2008)**, *Barcodes: the persistent topology of data* (Bull. AMS) — read it for: the
  cleanest short introduction to persistent homology if it is unfamiliar.
- **Carlsson (2009)**, *Topology and data* (Bull. AMS) — read it for: the broader programme, and
  a sober account of what topological methods can and cannot deliver on real data.
