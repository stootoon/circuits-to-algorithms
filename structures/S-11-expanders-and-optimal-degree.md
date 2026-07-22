---
title: S-11 Expanders & optimal degree
parent: Structures
nav_order: 11
---

# S-11 — Sparse random expansion, and a theory that predicted an anatomical number

> **The object.** Sparse random bipartite connectivity analysed as an expander / sparse
> random-projection matrix — and the derivation of the **optimal in-degree** for such a
> projection, which turned out to match the measured anatomy of two unrelated circuits.
> **Born from.** The discovery that Kenyon cells sample projection neurons at random, with a
> conspicuously small and stereotyped number of inputs.
> **Mathematical home.** Expander graphs, random projections, high-dimensional probability.
> **Situation.** #2 and #9 — the realizability and optimal-parameter problems.
> **Novelty.** Expanders and sparse Johnson–Lindenstrauss are classical. The *result* — a
> derived optimum that matches biology to within one synapse — is new, and it is one of the
> best things theoretical neuroscience has produced.

---

## 1. Where it comes from

In the *Drosophila* mushroom body, each Kenyon cell receives input from roughly **6** projection
neurons out of about 50. Caron, Ruta, Abbott & Axel traced this connectivity and found it to be,
as far as could be detected, **random** — no discernible structure in which PNs a given KC
samples, beyond a bias reflecting glomerular size.

Two questions follow, and they have very different characters. The first is qualitative: what is
random sparse connectivity *for*? That is [C1](../part2-case-studies/C1-expansion-and-sparsening.md)
— it is a random projection, a hash, a random-feature expansion. The second is quantitative and
much sharper: **why six?**

The second question is the one worth sitting with. It is rare for biology to hand you a specific
small integer and rarer still for theory to predict it. When both happen, you are looking at a
genuine explanation rather than a story.

---

## 2. The object, precisely

### The tradeoff

Let there be $M$ input units (PNs) and $N \gg M$ output units (KCs), each KC sampling $K$ inputs
at random. We want the KC representation to be **high-dimensional**, because dimension is what
buys linear separability
([S-10](S-10-tropical-and-piecewise-linear.md), [C1](../part2-case-studies/C1-expansion-and-sparsening.md)).
Measure it with the participation ratio of the KC covariance $C$:

$$\mathrm{PR} = \frac{\big(\sum_i \lambda_i\big)^2}{\sum_i \lambda_i^2}
= \frac{\big(\mathrm{tr}\,C\big)^2}{\mathrm{tr}\,C^2}.$$

Now the two failure modes, which pull in opposite directions:

**$K$ too small.** Each KC is essentially a copy of one or two PNs. It inherits their tuning and
their correlations; the KC layer is a redundant re-representation of the PN layer and adds no
dimensions. In the limit $K=1$, KC dimension cannot exceed PN dimension.

**$K$ too large.** Every KC samples a large fraction of the PNs, so any two KCs share many
inputs. Shared input means correlated responses. In the limit $K = M$, all KCs receive the *same*
input (up to weights) and the representation collapses to low dimension again.

Between these, an optimum. The calculation — carried out by Litwin-Kumar, Harris, Axel,
Sompolinsky & Abbott — expresses $\mathrm{PR}$ in terms of $K$, the PN correlation structure, and
the sparsening threshold, and maximizes.

### The result

The optimum is **$K \approx 6$–$7$** for the fly's parameters. The measured anatomical value is
about 6.

The theory was then applied to the cerebellar granule layer, an entirely different circuit in an
entirely different phylum: granule cells receive input from a small number of mossy fibres. The
predicted optimum is about **4**. The measured number of dendritic claws on a cerebellar granule
cell is about **4**.

Two independent circuits, separated by 600 million years of evolution, both sitting at the
theoretical optimum of the same quantity. That is as good as it gets. And note the epistemic
form: this is not "the model can be fit to the data," it is "the model has no free parameters
left and predicts a number, and the number is right." Hold every circuit→algorithm claim you
meet against this standard — see [S2](../part3-synthesis/S2-degeneracy-and-limits.md).

### The expander connection

Sparse random bipartite graphs with small left-degree are, with high probability, **expanders**:
every not-too-large subset of the left vertices has a neighbourhood substantially larger than
itself. Expansion is exactly the property that makes a sparse binary matrix behave like a good
measurement matrix — it is the basis of expander codes and of sparse compressed sensing
(Berinde & Indyk), and it is why sparse binary Johnson–Lindenstrauss constructions work (Kane &
Nelson).

The connection to make: **the mushroom body's connectivity is a sparse JL matrix, and the reason
$K$ can be as small as 6 is the same reason sparse JL works** — you do not need dense random
projections to preserve geometry, only enough expansion. Biology gets the metric-preservation
guarantee of a random projection at a sixth of the wiring cost.

---

## 3. Why the neuroscience forced it

Nothing here required new mathematics. What it required was asking a *quantitative* question of
a structure everyone had already described qualitatively.

The literature had "the mushroom body performs a random expansion" for years. That statement is
compatible with $K = 2$ and with $K = 30$; it forbids nothing; by the criteria of
[00 §2.1](../00-orientation/README.md) it is barely a level-2 claim at all. Turning it into a
prediction required specifying the objective (dimension of the representation), the constraint
(fixed $K$), and doing the optimization.

**This is the most transferable lesson in the structures thread.** The move from "this circuit
does X" to "this circuit does X, and therefore parameter $p$ should equal $6$" is the move that
converts a story into science, and it is almost always available if you are willing to write down
the objective explicitly.

---

## 4. How to recognize the pattern elsewhere

You are in Situation #9 when: **your qualitative story has a free parameter that biology has
already measured.**

Look for the numbers biology hands you and that everyone treats as incidental: degrees,
convergence and divergence ratios, layer sizes, time constants, sparsity levels, module counts,
scale ratios. Each is an opportunity. The recipe:

1. Write your functional story as an explicit objective $\mathcal{J}(p)$.
2. Identify the tradeoff — there must be one, or there is no interior optimum, and if your
   objective is monotone in $p$ you have mis-specified it (a very common error: "more dimensions
   is better" has no optimum; "more dimensions but shared input hurts" does).
3. Optimize; predict the number.
4. Compare to anatomy *without fitting*.

Worked precedents beyond this one: the optimal grid-cell scale ratio ($\approx\sqrt{2}$ to $e$,
depending on the objective — Wei, Prentice & Balasubramanian; see
[C6](../part2-case-studies/C6-grid-cells.md)); optimal KC coding sparsity; optimal numbers of
receptor types.

Open in your own area, and worth thinking about: the locust mushroom body has vastly more Kenyon
cells than the fly (∼50,000 vs ∼2,000) with a larger PN population — does the same optimization
predict the locust's in-degree? If the theory transfers across two phyla it should transfer
across two insects, and if it does not, that is more interesting still.

---

## 5. Exercises

**Ex S11.1 ★★** — Derive the qualitative form of the tradeoff: compute how the correlation
between two KCs depends on $K$, and show it produces an interior optimum for dimension.

<details markdown="1"><summary>Solution</summary>

Two KCs $a,b$ sample sets $S_a, S_b \subset [M]$ of size $K$ uniformly at random. The expected
overlap is

$$\mathbb{E}\,|S_a \cap S_b| = \frac{K^2}{M}.$$

Take PN activity $y$ with covariance $\Sigma$; for intuition take $\Sigma = I$ (uncorrelated PNs,
unit variance) and unit weights. Then the KC *input* covariance is

$$\mathrm{Cov}(u_a, u_b) = |S_a \cap S_b|, \qquad \mathrm{Var}(u_a) = K,$$

so the input correlation is

$$\rho_{ab} = \frac{|S_a \cap S_b|}{K} \approx \frac{K}{M}.$$

Now the two regimes:

- **Correlation grows linearly in $K$.** As $K \to M$, $\rho \to 1$: all KCs identical, and
  $\mathrm{PR} \to 1$. For a population with uniform pairwise correlation $\rho$, the covariance
  has one large eigenvalue $\approx N\rho$ and $N-1$ eigenvalues $\approx 1-\rho$, giving
  $$\mathrm{PR} \approx \frac{N}{1 + (N-1)\rho^2} \xrightarrow[\rho \to 1]{} 1 .$$
  So large $K$ kills dimension, at rate set by $\rho^2 \approx K^2/M^2$.

- **Small $K$ limits the achievable dimension differently.** With $K$ small the KC responses are
  nearly copies of individual PNs, so the KC representation cannot exceed the rank of the PN
  representation: $\mathrm{PR} \lesssim \mathrm{PR}_{\text{PN}}$. More carefully, with realistic
  *correlated* PNs, small $K$ means each KC inherits a single PN's correlation structure and no
  mixing occurs — mixing is what allows the expansion to create new dimensions, and mixing
  requires $K \geq 2$ and improves with $K$.

The two effects — mixing benefit rising with $K$, shared-input cost rising as $K^2/M^2$ — give an
interior maximum. Putting in the sparsening threshold (which suppresses the shared-input
correlations somewhat, shifting the optimum upward) and realistic PN correlations, the optimum
lands at $K \approx 6$–$7$ for $M \approx 50$.

Note the scaling intuition worth carrying away: the optimum grows *sublinearly* in $M$ — roughly
like $\sqrt{M}$ from balancing $K$-mixing against $K^2/M$ correlation — which is why a circuit
with 50 inputs wants ~6 and not ~25.
</details>

**Ex S11.2 ★★** — Verify numerically that dimension is maximized at intermediate $K$.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
rng = np.random.default_rng(0)
M, N, P = 50, 2000, 500           # PNs, KCs, odours

L = rng.standard_normal((M, 8)); Sigma = L@L.T + 0.5*np.eye(M)   # correlated PNs
Y = rng.multivariate_normal(np.zeros(M), Sigma, size=P).T        # M x P

def pr(K, sparsity=0.05):
    A = np.zeros((N, M))
    for i in range(N):
        A[i, rng.choice(M, K, replace=False)] = 1
    U = A @ Y
    thr = np.quantile(U, 1-sparsity, axis=1, keepdims=True)      # per-KC threshold
    R = np.maximum(U - thr, 0)
    C = np.cov(R)
    ev = np.linalg.eigvalsh(C)
    return ev.sum()**2 / (ev**2).sum()

for K in [1,2,3,4,6,8,12,20,35,50]:
    print(K, round(pr(K), 1))
```

You should see PR rise, peak at small-to-moderate $K$, and fall. The exact peak depends on the
PN correlation structure and the sparsity — which is the honest caveat: the "6" is not
universal, it is the optimum *for the fly's measured PN statistics and KC sparsity*. Vary
`sparsity` and the rank of `L` and watch the optimum move. Reproducing that sensitivity is the
real content of the exercise, because it tells you what the prediction is actually conditioned
on.
</details>

**Ex S11.3 ★★★** — *(Situation, not object.)* Find a number in your own system that theory
should predict. State the objective, the tradeoff, and the prediction — then look up the
measured value.

<details markdown="1"><summary>Solution</summary>

Open. Candidates from olfaction, roughly in order of tractability:

- **Locust KC in-degree.** Direct transfer of the calculation above with locust parameters.
  Requires locust PN correlation statistics and KC sparsity; both are approximately available.
  Tractable and, as far as I know, undone.
- **Number of LNs relative to PNs** in the antennal lobe. Objective: decorrelation or
  normalization quality ([C8](../part2-case-studies/C8-divisive-normalization.md)); tradeoff:
  more inhibitory neurons cost energy and add noise.
- **The 20–30 Hz oscillation frequency.** Objective: maximize information transmitted per unit
  time subject to downstream coincidence-detection windows and the membrane time constant of
  KCs. Tradeoff: faster cycles pack more coding epochs into the sniff but shrink the number of
  spikes per cycle below what a coincidence detector can resolve. This one is genuinely
  interesting because the constraint — the KC integration time constant — is measurable, and the
  prediction would be a *frequency*, which is about as sharp as it gets. I do not know of anyone
  having done it.
- **The duration of the transient** (~1–1.5 s in locust). Objective: enough time for
  decorrelation to converge; cost: behavioural latency. A speed–accuracy tradeoff, formally
  analogous to [C10](../part2-case-studies/C10-evidence-accumulation.md), which gives you the
  machinery for free.

The discipline: state the objective *before* you look up the number. Otherwise you are fitting,
and fitting a single number is not evidence of anything.
</details>

---

## 6. Reading

- **Caron, Ruta, Abbott & Axel (2013)**, on random convergence of olfactory inputs in the
  *Drosophila* mushroom body (Nature) — read it for: the anatomical result that random sparse
  sampling is what the circuit actually does.
- **Litwin-Kumar, Harris, Axel, Sompolinsky & Abbott (2017)**, *Optimal degrees of synaptic
  connectivity* (Neuron) — read it for: the derivation and the two matching predictions. This is
  the paper; it is short, and it is the best single example in the course of theory earning its
  keep.
- **Babadi & Sompolinsky (2014)**, on sparseness and expansion in sensory representations
  (Neuron) — read it for: the complementary analysis of what expansion does to signal and noise,
  and why more dimensions is not automatically better.
- **Berinde & Indyk**, on sparse recovery using sparse random matrices — read it for: the
  expander view of sparse measurement matrices, and the guarantees.
- **Kane & Nelson (2014)**, *Sparser Johnson–Lindenstrauss transforms* (J. ACM) — read it for:
  how sparse a random projection can be while preserving geometry. The theoretical ceiling that
  biology appears to be operating near.
- **Dasgupta, Stevens & Navlakha (2017)**, *A neural algorithm for a fundamental computing
  problem* (Science) — read it for: the LSH reading of the same circuit. See
  [C1](../part2-case-studies/C1-expansion-and-sparsening.md).
