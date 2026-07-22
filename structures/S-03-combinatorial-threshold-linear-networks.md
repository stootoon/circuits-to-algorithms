---
title: S-03 Combinatorial threshold-linear networks
parent: Structures
nav_order: 3
---

# S-03 — Combinatorial threshold-linear networks: a dictionary from graphs to attractors
{: .no_toc }

> **The object.** A correspondence $G \mapsto \mathrm{FP}(G)$ sending a directed graph to the
> set of fixed-point supports of an associated threshold-linear dynamical system, together with
> a growing list of **graph rules** — combinatorial theorems that determine the dynamics without
> integrating anything.
> **Born from.** The desire to know which recurrent circuits produce which attractors and
> sequences, without simulating parameter by parameter.
> **Mathematical home.** Graph theory meets nonlinear dynamics; polyhedral geometry.
> **Situation.** #4, the structure→dynamics dictionary.
> **Novelty.** **Genuinely new.** A correspondence between directed graph combinatorics and
> attractor structure, with theorems, that did not previously exist.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Where it comes from

Here is a predicament every modeller has been in. You build a recurrent network, you run it, you
get sequential activity, and you write a paper saying the connectivity "supports sequences."
Then a referee asks *why*, and the honest answer is: because I ran it and that is what came out.
Change $\epsilon$ and it might not.

The gap between "I simulated it" and "this graph *must* produce this behaviour" is a missing
dictionary. Curto, Morrison and collaborators built one, for a restricted but rich class.

The restriction is the trick. Rather than allowing arbitrary weights — where nothing general is
provable — let the weights be determined *combinatorially* by a graph, leaving only two global
parameters. Then ask what is true for all parameters in a legal range. The answer turns out to be:
a great deal.

---

## 2. The object, precisely

A **combinatorial threshold-linear network** (CTLN) on a directed graph $G$ with $n$ nodes:

$$\frac{dx_i}{dt} = -x_i + \Big[\ \sum_{j} W_{ij}x_j + \theta\ \Big]_+, \qquad [z]_+ = \max(z,0),$$

with weights determined entirely by $G$:

$$W_{ij} = \begin{cases}
0 & i = j,\\
-1 + \varepsilon & j \to i \text{ in } G,\\
-1 - \delta & j \not\to i \text{ in } G.
\end{cases}$$

Everything is inhibitory; the graph only decides *how strongly*. Edges are "less inhibited"
connections. The parameters must lie in the **legal range**

$$\delta > 0, \qquad 0 < \varepsilon < \frac{\delta}{\delta + 1}, \qquad \theta > 0,$$

and the remarkable empirical and theoretical fact is that within this range the qualitative
dynamics depend on $G$ *only*.

### Fixed points and their supports

For a fixed point $x^*$, let $\sigma = \mathrm{supp}(x^*) = \{i : x^*_i > 0\}$. On $\sigma$ the
dynamics are linear, so $x^*|_\sigma$ solves a linear system; for $x^*$ to be a genuine fixed
point it must additionally satisfy (i) *positivity* — the solution is strictly positive on
$\sigma$ — and (ii) *consistency* — every $j \notin \sigma$ receives net input $\leq 0$, so it
stays off. Define

$$\mathrm{FP}(G) = \{\sigma \subseteq [n] : \sigma \text{ is the support of a fixed point}\}.$$

$\mathrm{FP}(G)$ is a combinatorial object — a set system — and the central discovery is that it
is determined by $G$ alone across the legal range. The state space is carved by a hyperplane
arrangement into chambers, one per activity pattern, and within each chamber the flow is
affine; see [S-10](S-10-tropical-and-piecewise-linear.md) for that geometric picture, which is
the reason the combinatorics is available at all.

### Graph rules

The theorems are the payoff. A representative sample of the flavour (statements simplified;
see the reading for precise hypotheses):

- **Cliques.** A clique of $G$ that is not "dominated" by any outside node supports a stable
  fixed point. Cliques are the CTLN analogue of stored memories.
- **Domination.** If node $j$ *graphically dominates* node $k$ with respect to $\sigma$ —
  roughly, $j$ receives everything $k$ receives from $\sigma$ and more, and sends at least as
  much — then certain supports containing $k$ are ruled out of $\mathrm{FP}(G)$. This is a
  purely local, checkable condition that kills fixed points.
- **Sinks and sources.** Nodes with no outgoing edges, and other degenerate structures, admit
  clean rules for what survives.
- **Cyclic unions.** If $G$ is built by arranging component subgraphs in a directed cycle, with
  all forward edges present between consecutive components and no backward edges, then the
  network has **no stable fixed points** and instead produces a *sequential attractor* cycling
  through the components in order. The sequence is a theorem about the graph, not an observation
  about a simulation.
- **Core motifs.** The small graphs $\sigma$ for which $\mathrm{FP}(G|_\sigma) = \{\sigma\}$ —
  a unique full-support fixed point — are the irreducible units; they correspond to attractors,
  and larger networks' attractors are built from their core motifs. All core motifs on $\leq 5$
  nodes have been catalogued.

This is exactly the dictionary that was missing. Given a graph you can, in many cases, *read
off* whether the network is a memory network or a sequence generator.

---

## 3. Why the neuroscience forced it

Nothing in dynamical systems theory answers "which graphs give sequences." Dynamical systems
theory is largely local and analytic; this question is global and combinatorial. And nothing in
graph theory concerns attractors. The correspondence had to be built, and it was built because a
neuroscientific question — how does connectivity determine the repertoire of a circuit? — has no
useful answer at the level of individual simulations.

Note also the methodological move that made it tractable: **deliberately shrinking the parameter
space to expose combinatorial structure.** The CTLN weight scheme is biologically crude, and
that crudeness is the point. Fully general threshold-linear networks admit no such theory. This
is a trade every theorist should be able to make consciously: give up parametric generality to
buy structural theorems.

---

## 4. How to recognize the pattern elsewhere

You are in Situation #4 when you find yourself running parameter sweeps to establish a
qualitative claim. Parameter sweeps are evidence that a dictionary is missing.

The recipe that worked here, and that generalizes:

1. Identify the qualitative behaviours you care about (fixed points, sequences, oscillations).
2. Find a coarsening of the parameter space that preserves those behaviours — here, "which sign
   pattern of connections" rather than "which weights."
3. Prove that the coarsening is *sufficient*: the behaviour depends only on the coarse data.
4. Then hunt for local, checkable rules on the coarse data.

Step 3 is the hard one and the one that turns a modelling convenience into a theorem.

For [C2](../part2-case-studies/C2-transient-synchrony.md) this is directly relevant: "stable
heteroclinic sequences" and "sequential attractors" are competing formalizations of the same
intuition about transient olfactory dynamics, and CTLN cyclic unions give a third, fully
combinatorial one. Knowing all three, and knowing that they make different predictions about
what happens when you perturb connectivity, is exactly the kind of leverage the course is
trying to build.

---

## 5. Exercises

**Ex S3.1 ★** — For $n=2$ with $G$ the graph $1 \to 2$ only, write down the linear system for a
full-support fixed point and determine when it exists in the legal range.

<details markdown="1"><summary>Solution</summary>

With $\sigma = \{1,2\}$ both active, the fixed-point equations are $x_i = \sum_j W_{ij}x_j +
\theta$, i.e. $(I - W)x = \theta \mathbf{1}$ restricted to $\sigma$. Here

$$W = \begin{pmatrix} 0 & -1-\delta \\ -1+\varepsilon & 0 \end{pmatrix}$$

(row $i$, column $j$; $W_{21} = -1+\varepsilon$ because $1 \to 2$, and $W_{12} = -1-\delta$
because $2 \not\to 1$).

$$(I-W) = \begin{pmatrix} 1 & 1+\delta \\ 1-\varepsilon & 1\end{pmatrix}, \qquad
\det = 1 - (1+\delta)(1-\varepsilon).$$

Solving, $x_1 = \theta\,\dfrac{1 - (1+\delta)}{\det} = \dfrac{-\theta\delta}{\det}$ and
$x_2 = \theta\,\dfrac{1-(1-\varepsilon)}{\det} = \dfrac{\theta\varepsilon}{\det}$.

For both to be positive we need $x_1 > 0$, requiring $\det < 0$, and then $x_2 > 0$ requires
$\det > 0$. Contradiction: **there is no full-support fixed point**, for any legal parameters.

Sanity check against the graph rules: node 1 is a source that inhibits node 2 strongly
($-1-\delta$ from 2 back to 1 means 2 barely inhibits 1... careful with the direction). The
surviving fixed points are the singletons; you can verify $\{1\} \in \mathrm{FP}(G)$ (node 1
alone at $x_1 = \theta$, with node 2 receiving $(-1+\varepsilon)\theta + \theta = \varepsilon
\theta > 0$ — so actually $\{1\}$ fails consistency and node 2 turns on). Working through it,
$\{2\}$ survives: $x_2 = \theta$, and node 1 receives $(-1-\delta)\theta + \theta = -\delta\theta
< 0$, so it stays off. Good — the network settles with node 2 on, node 1 off, which is the
sensible reading of "1 excites 2 (relatively) and 2 strongly inhibits 1."
</details>

**Ex S3.2 ★★** — Implement a CTLN and verify the cyclic-union prediction. Take $n = 6$ arranged
as three components $\{1,2\},\{3,4\},\{5,6\}$ in a directed cycle (all edges forward between
consecutive components, all edges present within components). Integrate and confirm sequential
activity with no stable fixed point.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
from scipy.integrate import solve_ivp

n, eps, delta, theta = 6, 0.25, 0.5, 1.0
comp = [[0,1],[2,3],[4,5]]

A = np.zeros((n,n), dtype=bool)            # A[i,j] = True means j -> i
for c in range(3):
    nxt = comp[(c+1) % 3]
    for j in comp[c]:
        for i in nxt: A[i,j] = True        # forward edges between components
    for j in comp[c]:
        for i in comp[c]:
            if i != j: A[i,j] = True       # within-component edges (clique)

W = np.where(A, -1+eps, -1-delta); np.fill_diagonal(W, 0)
f = lambda t,x: -x + np.maximum(W @ x + theta, 0)
sol = solve_ivp(f, [0,150], 0.1*np.random.rand(n), dense_output=True, max_step=0.05)
```

Plot `sol.y`. You should see a robust periodic orbit in which components $\{1,2\}$,
$\{3,4\}$, $\{5,6\}$ activate in cyclic order, with each component's two neurons co-active
(they form a clique). Verify no stable fixed point by starting from many random initial
conditions — all converge to the same limit cycle — and by checking numerically that no
candidate support satisfies positivity + consistency + stability.

The point to internalize: you did not tune anything. Any legal $(\varepsilon, \delta, \theta)$
gives the same qualitative answer, because the theorem is about the graph.
</details>

**Ex S3.3 ★★★** — *(Situation, not object.)* Take a claim from your own modelling where you
established a qualitative behaviour by parameter sweep. Identify the coarsening of parameter
space that you believe preserves the behaviour, and state precisely the theorem you would need
in order to stop sweeping. Is it plausibly true?

<details markdown="1"><summary>Solution</summary>

Open; this is the transferable skill. A worked template for olfaction:

*Claim established by sweep.* "A PN–LN network with this connectivity produces odour-specific
transient trajectories that decorrelate over time."

*Proposed coarsening.* The behaviour depends on the connectivity only through (i) the sign
pattern of PN→LN and LN→PN connections, (ii) the spectral radius of the effective recurrent
matrix, and (iii) the degree of non-normality (measured, say, by
$\|W^\top W - WW^\top\|_F / \|W\|_F^2$).

*Theorem needed.* "For all connectivity matrices with the given sign pattern, spectral radius in
$(a,b)$, and non-normality index above $\nu$, the pairwise correlation between responses to
distinct inputs is monotonically decreasing over $t \in (0, T)$."

*Plausibility.* Partially. The non-normality→transient-amplification link is real and provable
(see [S-07](S-07-random-matrices-and-chaos.md); Murphy & Miller's balanced amplification,
Hennequin et al.'s stability-optimized circuits). The *decorrelation* part is harder because it
is a statement about a pair of trajectories, not one — you would likely need to restrict to
inputs in general position and would get a probabilistic rather than universal statement. Which
is itself informative: it tells you the honest claim is "generic decorrelation," and that a
paper claiming universal decorrelation is over-claiming.

That reframing — from "I swept and it worked" to "here is the class over which it holds, and
here is why it can only be generic, not universal" — is the deliverable.
</details>

---

## 6. Reading

- **Curto, Degeratu & Itskov (2013)**, on flexible memory networks and permitted sets of
  threshold-linear networks — read it for: the origin of the "which supports are permitted"
  framing, and the connection to stable symmetric matrices.
- **Morrison, Degeratu, Itskov & Curto**, on diversity of emergent dynamics in competitive
  threshold-linear networks — read it for: the CTLN definition and the first systematic survey
  of what the graphs produce.
- **Curto, Geneson & Morrison (2019)**, *Fixed points of competitive threshold-linear networks*
  (Neural Computation) — read it for: the core theory of $\mathrm{FP}(G)$ and the domination
  rules. This is the technical centre.
- **Parmelee, Moore, Morrison & Curto**, on sequential attractors and cyclic unions — read it
  for: the sequence-generation theorems, which are the most directly relevant to
  [C2](../part2-case-studies/C2-transient-synchrony.md).
- **Curto & Morrison (2019)**, review of relating network connectivity to dynamics (Curr. Opin.
  Neurobiol.) — read it for: the short version, and the best entry point if you want the shape
  before the proofs.
