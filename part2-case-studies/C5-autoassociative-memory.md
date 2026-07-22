---
title: C5 · Autoassociative memory
nav_order: 45
---

# C5 — Recurrent autoassociation → content-addressable memory

> **Circuit.** Hippocampal CA3: ~$10^4$ recurrent collateral synapses per pyramidal cell onto a population of the same cells, fed by a sparse, expanded, detonator-synapse input from dentate gyrus and a diffuse one from entorhinal cortex.
> **Primitive extracted.** Store patterns as attractors of a symmetric recurrent network built by Hebbian outer products; retrieve by relaxation from a partial cue. The memory is addressed by *content*, not location, and capacity is a computable function of connectivity, sparseness, and the form of the energy.
> **Status.** The best-developed *theory* in computational neuroscience and one of the weakest *conversions*. The mathematics is exact and beautiful; the mapping to CA3 is suggestive, partially tested, and has never survived a decisive experiment. Meanwhile the object escaped and conquered machine learning.
> **Structures thread.** [`../structures/S-04-manifold-capacity.md`](../structures/S-04-manifold-capacity.md) (Gardner-style capacity calculations, of which this unit's $\alpha_c$ is the ancestor); [`../structures/S-03-combinatorial-threshold-linear-networks.md`](../structures/S-03-combinatorial-threshold-linear-networks.md) (reading attractors off connectivity without a Lyapunov function); [`../structures/README.md`](../structures/README.md), situations 3 and 4.

---

## 1. The phenomenon

David Marr, in 1971, looked at the archicortex and asked a question nobody had asked in that form: *what is a piece of cortex with dense recurrent collaterals for?* His answer — "simple memory," a device that stores events as patterns of activity and completes a stored pattern from a fragment of it — was posed with an explicit capacity calculation, an explicit division of labour between input pathways, and an explicit prediction about sparseness. It was eleven years before Hopfield's paper, forty years before anyone could test it, and it is still the single most impressive act of circuit→algorithm conversion in the field's history.

The anatomy that motivated it: CA3 pyramidal cells receive recurrent collateral input from other CA3 cells with high divergence and low probability — roughly $C\approx1.2\times10^4$ recurrent contacts per cell in rat, in a population of a few hundred thousand. This is the only place in the mammalian brain with a large, genuinely *associational* recurrent network onto a single homogeneous cell population. It also receives ~50 enormous mossy-fibre synapses from dentate granule cells, of which about 1% are active at any time, and each of which can single-handedly fire the postsynaptic cell.

Marr and, much more sharply, Treves & Rolls (1992) read this as a two-input architecture with two jobs: **mossy fibres write** (sparse, powerful, forcing a new pattern onto CA3 during storage) and **perforant path reads** (diffuse, weak, supplying a cue that recurrent dynamics completes). That is an algorithmic claim about pathways, and it predicts a dissociation: impair the mossy fibres and you lose encoding; impair CA3 recurrence and you lose completion from partial cues (Nakazawa et al. 2002).

Upstream, the dentate gyrus does the opposite job: ~$10^6$ granule cells fed by ~$2\times10^5$ entorhinal layer II cells, firing at 1–2%, strongly inhibited. This is **expansion followed by sparsification** — precisely the motif of [C1](C1-expansion-and-sparsening.md). Its job is *pattern separation*: make similar inputs dissimilar before storage, so the autoassociator does not confuse them.

Three parts, three primitives — separate, store, complete — and the middle one is where the mathematics is.

---

## 2. The conversion

### 2.1 The Hopfield model and its energy

Take $N$ binary units $\sigma_i\in\{-1,+1\}$, $p$ random patterns $\xi^\mu \in\{-1,+1\}^N$, $\mu=1,\dots,p$, and the Hebbian outer-product rule

$$J_{ij} = \frac{1}{N}\sum_{\mu=1}^{p}\xi_i^\mu\xi_j^\mu,\qquad J_{ii}=0 .$$

Dynamics: pick a unit $k$ at random and set $\sigma_k \leftarrow \mathrm{sgn}(h_k)$, $h_k=\sum_j J_{kj}\sigma_j$.

Hopfield's move — the one that turned a heuristic into a theory — was to notice that $J$ is symmetric and therefore admits a Lyapunov function:

$$E(\boldsymbol\sigma) = -\tfrac12\sum_{i\ne j}J_{ij}\sigma_i\sigma_j .$$

Flip unit $k$. Only terms containing $\sigma_k$ change, and by symmetry there are two of them, so

$$\Delta E = -\tfrac12\Big[\sum_{j\ne k}J_{kj}(-2\sigma_k)\sigma_j + \sum_{i\ne k}J_{ik}\sigma_i(-2\sigma_k)\Big] = 2\sigma_k h_k .$$

The update flips $\sigma_k$ only when $\mathrm{sgn}(h_k)\ne\sigma_k$, i.e. only when $\sigma_kh_k<0$, so $\Delta E<0$ on every flip. $E$ is bounded below and decreases strictly, hence the dynamics converges to a fixed point in finite time. **The network is a descent algorithm on a landscape whose minima are the memories.** Two hypotheses are load-bearing and both are biologically false: symmetry of $J$, and asynchronous updating. Neither is fatal — asymmetric networks still retrieve, they just may cycle — but you should notice exactly where the mathematics stops being about brains.

### 2.2 Stability of a stored pattern, and the crosstalk term

Set $\boldsymbol\sigma=\xi^\nu$ and compute the local field:

$$h_i = \frac1N\sum_{j\ne i}\sum_{\mu}\xi_i^\mu\xi_j^\mu\xi_j^\nu = \underbrace{\frac{N-1}{N}\,\xi_i^\nu}_{\text{signal}} \;+\; \underbrace{\frac1N\sum_{\mu\ne\nu}\xi_i^\mu\sum_{j\ne i}\xi_j^\mu\xi_j^\nu}_{\text{crosstalk } C_i^\nu}.$$

Unit $i$ is stable iff $\xi_i^\nu h_i>0$. Define the *aligned* crosstalk $ C \equiv -\xi_i^\nu C_i^\nu$; then instability means $C>1$ (dropping the $1/N$ correction). Now count: $C$ is a sum of $(p-1)(N-1)$ terms each equal to $\pm1/N$ with equal probability and, for random patterns, essentially independently. So

$$\mathbb E[C]=0,\qquad \mathrm{Var}[C] = \frac{(p-1)(N-1)}{N^2}\;\approx\;\frac pN \equiv \alpha .$$

By the CLT, $C\sim\mathcal N(0,\alpha)$ and the per-bit error probability at the moment you present the stored pattern is

$$\boxed{\;P_{\rm err} = \Pr[C>1] = \tfrac12\,\mathrm{erfc}\!\Big(\frac{1}{\sqrt{2\alpha}}\Big) = Q\!\Big(\frac{1}{\sqrt\alpha}\Big)\;}$$

with $Q$ the Gaussian tail. Invert: $\alpha = 1/z^2$ where $Q(z)=P_{\rm err}$.

| $P_{\rm err}$ | $z$ | $\alpha = p/N$ |
|---|---|---|
| 0.05 | 1.645 | 0.370 |
| 0.01 | 2.326 | 0.185 |
| **0.0036** | **2.69** | **0.138** |
| 0.001 | 3.09 | 0.105 |

**Two things to take from this table.** First, signal-to-noise gets you the right order of magnitude — $\alpha\sim0.1$–0.2 — in three lines and no statistical mechanics. Second, and this is the part usually skipped: the true capacity $\alpha_c=0.138$ sits at a per-bit initial error rate of $0.36\%$. That is not a coincidence. It is where the *avalanche* becomes supercritical: each wrong bit corrupts the field at every other unit, and the question is whether the corrections outrun the corruptions.

### 2.3 Where $\alpha_c=0.138$ actually comes from

The naive calculation treats the crosstalk as external noise. It is not: the network *responds* to the crosstalk, and its response feeds back into the crosstalk. Amit, Gutfreund & Sompolinsky (1985, 1987) did the calculation properly with replicas. The structure of the answer is worth knowing even if you never do the algebra.

Define the overlap with the retrieved pattern $m=\frac1N\sum_i\xi_i^1\sigma_i$, the Edwards–Anderson order parameter $q=\frac1N\sum_i\langle\sigma_i\rangle^2$, and $r$ = the variance of the non-condensed pattern noise *in units of $\alpha$*. The replica-symmetric mean-field equations at temperature $T=1/\beta$ are

$$m=\int\! Dz\;\tanh\big[\beta(m+\sqrt{\alpha r}\,z)\big],\qquad q=\int\! Dz\;\tanh^2\big[\beta(m+\sqrt{\alpha r}\,z)\big],\qquad r=\frac{q}{\big[1-\beta(1-q)\big]^2},$$

with $Dz$ the standard Gaussian measure. Look at the third equation. If there were no feedback we would have $r=1$ and the crosstalk variance would be exactly $\alpha$, reproducing §2.2. Instead $r=q/[1-\beta(1-q)]^2 > 1$: the denominator is an **Onsager reaction term**, the network's own susceptibility amplifying the noise it is trying to suppress. That single factor is the whole difference between $0.185$ and $0.138$.

At $T=0$ ($q\to1$, $\beta(1-q)\to C$ finite) these collapse to

$$m=\mathrm{erf}\!\Big(\frac{m}{\sqrt{2\alpha r}}\Big),\qquad r=\frac{1}{(1-C)^2},\qquad C=\sqrt{\frac{2}{\pi\alpha r}}\;e^{-m^2/2\alpha r}.$$

Solve numerically: a retrieval solution with $m$ close to 1 exists for $\alpha<\alpha_c$ and vanishes *discontinuously* at

$$\alpha_c = 0.138,\qquad m_c\approx0.967 .$$

The transition is first-order in $m$: the memory does not degrade gracefully, it collapses. Push past $\alpha_c$ and the network does not retrieve a slightly noisy version of $\xi^1$ — it falls into the spin-glass phase, where an exponential number of spurious minima exist and none of them is a memory. Replica symmetry breaking corrects $\alpha_c$ in the third decimal place, which tells you the retrieval states are not strongly hierarchical.

**The opinionated part.** Every textbook quotes $0.138$; almost none says what is bad about it. $0.138N$ patterns of $N$ bits is $0.138N^2$ bits stored in $N^2/2$ synapses: about $0.28$ bits per synapse, with an all-or-nothing failure mode and a requirement that patterns be uncorrelated. For a brain, all three properties are wrong. Sparseness fixes the first two; the modern-Hopfield story of §2.5 fixes the first spectacularly.

### 2.4 Sparse coding: the capacity that actually applies to CA3

Real hippocampal representations are sparse: a fraction $f\approx0.02$–$0.05$ of CA3 cells are active for a given environment/event. Take $\xi_i^\mu\in\{0,1\}$ with $\Pr[\xi=1]=f$ and use the covariance rule

$$J_{ij}=\frac{1}{Nf(1-f)}\sum_\mu(\xi_i^\mu-f)(\xi_j^\mu-f),$$

with a threshold set (equivalently, a global inhibition) so that exactly $fN$ units are active at retrieval. Redo the signal-to-noise argument. Signal $\propto (1-f)$ per active unit; crosstalk variance now scales as $\alpha f/(1-f)$ rather than $\alpha$; and, crucially, the *threshold* need only separate active from inactive rather than $+1$ from $-1$. The result (Tsodyks & Feigel'man 1988; Treves & Rolls 1991) is

$$\boxed{\;p_{\max}\;\approx\;\frac{k\,C}{f\ln(1/f)}\;},\qquad k\approx0.2\text{–}0.3,$$

with $C$ the number of *recurrent connections per neuron* (not $N$ — an important detail, since CA3 is not fully connected).

There is a beautiful information-theoretic shortcut to the scaling. A sparse binary pattern carries $\approx N f\log_2(1/f)$ bits. If the network stores $O(1)$ bits per synapse, the total is $O(NC)$ bits, so

$$p_{\max}\cdot Nf\log_2(1/f) \sim NC \quad\Longrightarrow\quad p_{\max}\sim \frac{C}{f\log(1/f)} .$$

Same answer, no replicas, and it tells you *why*: sparseness helps because each pattern is cheap in bits, not because interference magically vanishes.

**Numbers for rat CA3.** $C\approx1.2\times10^4$, $f\approx0.03$, $k=0.25$:

$$p_{\max}\approx\frac{0.25\times1.2\times10^4}{0.03\times\ln(33)}\approx\frac{3000}{0.105}\approx 3\times10^4 .$$

Tens of thousands of storable episodes, set by *connections per cell* — which is why $C$, not the number of neurons, is the parameter evolution should care about, and why CA3 recurrent divergence is so conspicuously high. That inference is, to me, the strongest piece of evidence that the conversion is right: a quantitative anatomical peculiarity explained by an algorithmic requirement.

### 2.5 Dense associative memory, and the escape into machine learning

Rewrite the Hopfield energy in terms of pattern overlaps:

$$E(\boldsymbol\sigma) = -\frac{1}{2N}\sum_\mu\big(\boldsymbol\xi^\mu\!\cdot\!\boldsymbol\sigma\big)^2 + \text{const}.$$

The energy is a *quadratic* function of the overlaps. Krotov & Hopfield (2016) asked the obvious question nobody had asked: why quadratic? Take

$$E = -\sum_\mu F\big(\boldsymbol\xi^\mu\!\cdot\!\boldsymbol\sigma\big),\qquad F(x)=x^n .$$

Redo the stability analysis. For $\boldsymbol\sigma=\boldsymbol\xi^1$, the energy difference for flipping unit $i$ is $\Delta_i = \sum_\mu\big[F(h^\mu+\xi_i^\mu)-F(h^\mu-\xi_i^\mu)\big]$ with $h^\mu=\sum_{j\ne i}\xi_j^\mu\sigma_j$. For the target, $h^1\approx N$ and the contribution is $\approx 2\xi_i^1 F'(N)=2n\xi_i^1N^{n-1}$. For each other pattern, $h^\mu\sim\sqrt N$ and the contribution is $\sim 2nN^{(n-1)/2}$ with a random sign; summing $p$ of them gives $\sim n\sqrt p\,N^{(n-1)/2}$. Stability requires signal $\gg$ noise:

$$n N^{n-1}\gg n\sqrt p\,N^{(n-1)/2} \quad\Longrightarrow\quad \boxed{\;p_{\max}\sim N^{\,n-1}\;}$$

For $n=2$ this recovers the classical linear-in-$N$ capacity. For $n=3$ it is quadratic. Demircigil et al. (2017) pushed to $F(x)=e^{x}$ and obtained capacity exponential in $N$ ($\sim 2^{N/2}$).

Then Ramsauer et al. (2020) made the states continuous. Let $X\in\mathbb R^{d\times M}$ hold $M$ stored patterns as columns and let $\boldsymbol\xi\in\mathbb R^d$ be the state. Define

$$E(\boldsymbol\xi) = -\beta^{-1}\log\sum_{\mu=1}^{M}e^{\beta\,\boldsymbol x_\mu\cdot\boldsymbol\xi} + \tfrac12\|\boldsymbol\xi\|^2 + \text{const}.$$

Minimizing by the concave–convex procedure gives the update

$$\boxed{\;\boldsymbol\xi^{\rm new} = X\,\mathrm{softmax}\big(\beta X^\top\boldsymbol\xi\big)\;}$$

Stare at that. With $\boldsymbol\xi$ playing the role of a query, $X^\top$ the keys, and $X$ the values, this is *exactly* the transformer attention operation $\mathrm{softmax}(QK^\top/\sqrt{d})V$. One step of retrieval in a modern Hopfield network **is** one attention head. And Ramsauer et al. prove the retrieval is one-step-exact when patterns are well separated, with exponential-in-$d$ storage.

This deserves to be flagged loudly, because it is the course's cleanest example of the arrow running backwards. The chain is: *archicortex anatomy* (Marr 1971) → *energy-based associative memory* (Hopfield 1982) → *statistical mechanics of storage* (AGS 1987) → *higher-order energies* (Krotov & Hopfield 2016) → *the attention mechanism of every large language model*. Attention was invented independently in NLP, but the Hopfield derivation explains *why it works*: it is the retrieval dynamics of an associative memory with exponential capacity, and the softmax temperature $\beta$ is the sharpness knob between "average all memories" and "retrieve one."

**The biological catch, and my position on it.** $F(x)=x^n$ requires $n$-th order synaptic interactions, which is not a thing. Krotov & Hopfield (2021) resolve this by introducing an explicit layer of "memory neurons," one per stored pattern, connected to the "feature neurons" by ordinary pairwise synapses; the many-body energy emerges from eliminating the hidden layer. But note what that architecture *is*: an expanded layer of highly selective units, each responding to one pattern, feeding back onto a smaller feature layer. That is not CA3's recurrent collateral network. That is **DG→CA3**, or **PN→Kenyon cell→MBON**, or **mossy fibre→granule cell→Purkinje**. The dense associative memory looks far more like the expansion architectures of [C1](C1-expansion-and-sparsening.md) than like the recurrent autoassociator it is named after. I think this is the most interesting unexploited hypothesis in the area: *the high-capacity memory in a real brain may live in the expansion, not the recurrence*, with the recurrence contributing completion and the expansion contributing capacity.

### 2.6 Pattern separation as an expansion operation

The autoassociator's capacity formula assumes uncorrelated patterns. Real inputs are not: two similar rooms produce similar entorhinal activity, and storing them in the same network guarantees interference. The dentate gyrus's job is to fix this before storage.

Model the operation as a random expansion with thresholding: $M$ inputs $\to N\gg M$ outputs, $y_a = \Theta\big(\sum_b W_{ab}x_b - \theta\big)$ with random $W$. For jointly Gaussian inputs with correlation $\rho$, the pre-threshold outputs of a given cell for two stimuli are jointly Gaussian with the same $\rho$, so the correlation of the binary outputs is

$$\rho_{\rm out} = \frac{\Phi_2(\theta,\theta;\rho)-f^2}{f(1-f)},\qquad f=Q(\theta),$$

with $\Phi_2$ the bivariate Gaussian tail. For $\rho<1$ and $f\to0$, $\rho_{\rm out}\to 0$ faster than $\rho$: **thresholding at high sparseness decorrelates.** The expansion does not by itself decorrelate (a random linear map preserves correlations up to sampling noise); the *nonlinearity at high threshold* does, and the expansion is what makes the sparse thresholded code still carry enough information to be worth storing. This is exactly the analysis of [C1](C1-expansion-and-sparsening.md), where the same three-part motif (expand, randomize, sparsify) appears in the insect mushroom body and the cerebellar granule layer, and where Babadi & Sompolinsky (2014) and Litwin-Kumar et al. (2017) work out the optimal expansion ratio and connection degree.

So the hippocampal trio reads, algorithmically:

| Stage | Operation | Level-2 primitive |
|---|---|---|
| EC → DG | expand $\times5$, sparsify to $f\approx0.01$ | decorrelation / random features |
| DG → CA3 | few, powerful, forcing synapses | *write* a specific pattern |
| CA3 → CA3 | Hebbian recurrent collaterals | store & complete: energy descent |
| EC → CA3 | diffuse weak | *read*: supply a cue |

Empirically: Leutgeb et al. (2007) showed DG representations diverge for similar environments while CA3 representations converge — separation upstream, completion downstream, in the same animal; Neunuebel & Knierim (2014) sharpened it with simultaneous DG/CA3 recording. McHugh et al. (2007) deleted NMDA receptors in DG granule cells and impaired separation specifically; Nakazawa et al. (2002) did CA3 and impaired completion specifically. About as good as double dissociation gets in this business.

---

## 3. Worked example / model to build

Two experiments: the classical capacity cliff, and one-step attention-style retrieval far beyond it.

```python
import numpy as np
rng = np.random.default_rng(1)

# ---------- A. classical Hopfield: find the capacity cliff ----------
def retrieve(N, alpha, flip=0.10, sweeps=15):
    p  = max(1, int(round(alpha*N)))
    xi = rng.choice([-1.0, 1.0], size=(p, N))
    J  = xi.T @ xi / N
    np.fill_diagonal(J, 0.0)
    s  = xi[0].copy()
    s[rng.random(N) < flip] *= -1                 # corrupt 10% of the bits
    for _ in range(sweeps):
        for i in rng.permutation(N):              # asynchronous updates
            s[i] = 1.0 if J[i] @ s >= 0 else -1.0
    return abs(s @ xi[0]) / N                     # final overlap with the target

N = 400
for alpha in [0.05, 0.08, 0.11, 0.13, 0.138, 0.15, 0.18, 0.25]:
    m = [retrieve(N, alpha) for _ in range(8)]
    print(f"alpha={alpha:6.3f}   <m>={np.mean(m):.3f}   min={np.min(m):.3f}")

# ---------- B. modern (dense) associative memory = one attention step ----------
def modern(N, p, beta, flip=0.30):
    X = rng.choice([-1.0, 1.0], size=(p, N))      # p patterns, p >> N allowed
    q = X[0].copy(); q[rng.random(N) < flip] *= -1
    z = beta * (X @ q)
    a = np.exp(z - z.max()); a /= a.sum()         # softmax(beta * K q)
    return ((a @ X) @ X[0]) / N                   # overlap after ONE step

for p in [50, 500, 5000, 50000]:
    print(f"p={p:6d}  alpha={p/200:8.1f}   m_after_one_step={modern(200, p, beta=0.4):.3f}")
```

**What to look for.**
1. In part A the mean overlap sits at $0.99$–$1.00$ for $\alpha\lesssim0.12$, wobbles around $\alpha\approx0.138$, and falls to $\approx0.5$ by $\alpha=0.18$ and $\approx0.4$ at $\alpha=0.25$ (the residue is the spin-glass state's incidental correlation with the target, not retrieval). The interesting column is `min`: at $\alpha=0.138$ you should see mean $\approx0.92$ but individual runs as low as $0.5$ — the distribution is **bimodal**, runs either retrieve or fail, and it is that bimodality, not the mean, that is the signature of AGS's *discontinuous* transition.
2. Increase $N$ from 200 to 800 and watch the cliff sharpen around $0.138$; finite-size rounding is why $N=400$ looks soft.
3. In part B the classical capacity at $N=200$ is $\approx27$ patterns. You will retrieve at $p=50\,000$ — three orders of magnitude past it — with overlap $1.000$, in a *single* step, no iteration, no energy descent. Tune $\beta$: too small gives the average of all patterns ($m\to0$); too large makes the softmax a hard max and you lose the noise-averaging that made retrieval robust to the 30% corrupted bits.
4. Modify part A to sparse patterns ($\xi\in\{0,1\}$, $f=0.05$, covariance rule, global inhibition holding $fN$ units active) and verify $p_{\max}\propto 1/(f\ln(1/f))$ by scanning $f$.

---

## 4. Exercises

**★ E5.1 — The Lyapunov function and what breaks it.** Prove $\Delta E=2\sigma_kh_k$ for an asynchronous flip, and then show by explicit construction of a 2-unit network that an antisymmetric $J$ produces a limit cycle rather than convergence.

<details markdown="1"><summary>Solution</summary>

**Proof.** $E=-\frac12\sum_{i\ne j}J_{ij}\sigma_i\sigma_j$. Only terms with $i=k$ or $j=k$ change when $\sigma_k\to-\sigma_k$:
$$\Delta E = -\tfrac12\Big[\sum_{j\ne k}J_{kj}(-\sigma_k-\sigma_k)\sigma_j + \sum_{i\ne k}J_{ik}\sigma_i(-\sigma_k-\sigma_k)\Big] = \sigma_k\sum_{j\ne k}J_{kj}\sigma_j + \sigma_k\sum_{i\ne k}J_{ik}\sigma_i.$$
With $J_{ij}=J_{ji}$ both sums equal $h_k$, giving $\Delta E=2\sigma_kh_k$. A flip happens only if $\mathrm{sgn}(h_k)\ne\sigma_k$, i.e. $\sigma_kh_k<0$, hence $\Delta E<0$. Since $E$ takes finitely many values and strictly decreases, the dynamics halts.

**Counterexample.** $J_{12}=+1$, $J_{21}=-1$, $J_{11}=J_{22}=0$. Then $h_1=\sigma_2$, $h_2=-\sigma_1$. Start at $(+,+)$ and update alternately: $\sigma_1\leftarrow\mathrm{sgn}(\sigma_2)=+$ (no change); $\sigma_2\leftarrow\mathrm{sgn}(-\sigma_1)=-$, giving $(+,-)$; then $\sigma_1\leftarrow\mathrm{sgn}(-1)=-$, giving $(-,-)$; then $\sigma_2\leftarrow\mathrm{sgn}(+1)=+$, giving $(-,+)$; then $\sigma_1\leftarrow+$, giving $(+,+)$. Period-4 cycle. No $E$ exists because the "force" field is not a gradient: the circulation around this 4-cycle is nonzero.

**Why it matters.** Cortical connectivity is emphatically not symmetric, and this counterexample is the reason the Hopfield framework does not simply transfer. The honest statement is that symmetry is a *sufficient* condition for convergence, not a necessary one; and networks with $J = J_{\rm sym}+\epsilon J_{\rm asym}$ retain fixed-point retrieval up to a critical $\epsilon$, beyond which retrieval states become slow limit cycles or chaotic wandering among memories — which some people, not unreasonably, regard as a feature.
</details>

**★★ E5.2 — Capacity from signal-to-noise.** Derive $P_{\rm err}=Q(1/\sqrt\alpha)$ and compute the maximum $\alpha$ for which fewer than 1 bit in 1000 is initially wrong. Then explain, quantitatively, why the true $\alpha_c$ is smaller than the value this criterion suggests.

<details markdown="1"><summary>Solution</summary>

**Derivation** as in §2.2: $C=-\xi_i^\nu C_i^\nu = -\frac1N\sum_{\mu\ne\nu}\sum_{j\ne i}\xi_i^\nu\xi_i^\mu\xi_j^\mu\xi_j^\nu$. Each of the $(p-1)(N-1)$ summands is a product of four independent $\pm1$'s, hence $\pm1$ with equal probability; they are pairwise uncorrelated. So $\mathbb E C=0$, $\mathrm{Var}\,C=(p-1)(N-1)/N^2\approx\alpha$. CLT gives $C\sim\mathcal N(0,\alpha)$, and instability of bit $i$ is $C>1$, so $P_{\rm err}=\Pr[\mathcal N(0,1)>1/\sqrt\alpha]=Q(1/\sqrt\alpha)$.

**$P_{\rm err}=10^{-3}$:** $Q(z)=10^{-3}\Rightarrow z=3.090$, so $\alpha=1/z^2=0.105$. So $p\le0.105N$.

**Why the truth is smaller.** Two effects, both missed by the calculation.
1. *Dynamics.* The $P_{\rm err}$ above is the error rate at $t=0$, presenting the exact stored pattern. Those wrong bits then corrupt the fields of all other units at the next step. Let $x$ be the fraction of wrong bits. Each wrong bit reduces the signal from $1$ to $1-2x$ and adds to the noise. Self-consistency, $x = Q\big((1-2x)/\sqrt{\alpha}\big)$, has a stable small-$x$ root only up to some $\alpha$; beyond it the only root is $x\approx1/2$. This "avalanche" analysis already lowers the threshold and reproduces the first-order character.
2. *Reaction.* The replica calculation shows the residual noise variance is $\alpha r$, not $\alpha$, with $r=q/[1-\beta(1-q)]^2>1$: the network's own susceptibility to the noise re-amplifies it (Onsager term). At $T=0$ near threshold $r\approx1.3$–$1.5$, which converts $\alpha\approx0.185$ into $\alpha\approx0.138$: $0.185/1.34\approx0.138$. The naive criterion is off by exactly this reaction factor.

Note that $\alpha_c=0.138$ corresponds to $P_{\rm err}=Q(1/\sqrt{0.138})=Q(2.69)=0.0036$ — the capacity is reached when about one bit in 280 is initially wrong. In an $N=1000$ network that is 3.6 wrong bits, and the avalanche from 3.6 wrong bits is enough to destroy the memory. Fragile.
</details>

**★★ E5.3 — Spurious states.** Show that the "odd mixture" $\sigma_i=\mathrm{sgn}(\xi_i^1+\xi_i^2+\xi_i^3)$ is a stable fixed point of the Hopfield dynamics in the low-loading limit, and that the corresponding 2-pattern mixture is not.

<details markdown="1"><summary>Solution</summary>

**Overlaps.** Let $\sigma=\mathrm{sgn}(\xi^1+\xi^2+\xi^3)$. Compute $m^1=\langle\xi_i^1\sigma_i\rangle$ over the 8 equiprobable sign patterns of $(\xi_i^1,\xi_i^2,\xi_i^3)$. Condition on $\xi^1_i=+1$; the other two give sums $\xi_i^2+\xi_i^3\in\{+2,0,0,-2\}$, so the total is $\{3,1,1,-1\}$ and $\sigma_i=\{+,+,+,-\}$. Hence $\xi_i^1\sigma_i$ averages to $(1+1+1-1)/4=\tfrac12$. By symmetry $m^1=m^2=m^3=\tfrac12$, and $m^\mu=0$ for $\mu>3$.

**Stability.** In the low-loading limit ($p/N\to0$) the field is $h_i=\sum_\mu m^\mu\xi_i^\mu = \tfrac12(\xi_i^1+\xi_i^2+\xi_i^3)$. Then $\mathrm{sgn}(h_i)=\mathrm{sgn}(\xi_i^1+\xi_i^2+\xi_i^3)=\sigma_i$. Self-consistent: stable fixed point. (The sum is never zero for three $\pm1$'s, so $\mathrm{sgn}$ is well defined.)

**Two patterns.** $\xi_i^1+\xi_i^2\in\{+2,0,0,-2\}$, and $\mathrm{sgn}(0)$ is undefined for half the units — the state is not even well defined. If you break the tie arbitrarily, the resulting state has $m^1=m^2=\tfrac12$ from the aligned units plus $\pm$ contributions from the tied ones; the field on a tied unit is $\tfrac12(\xi_i^1+\xi_i^2)=0$, so those units are marginal, not stable. Any noise destroys the state. More generally, only *odd* mixtures of $2k+1$ patterns are stable; even mixtures are not, which is a genuinely surprising and completely elementary result.

**Consequence.** The energy landscape contains $\binom{p}{3}$ 3-mixtures, $\binom p5$ 5-mixtures, etc. — a combinatorially large set of attractors that are not memories. At finite loading these plus the spin-glass states dominate. A brain that used this scheme would confabulate blends of episodes, which is either a bug or an unusually good model of human memory, depending on your mood.
</details>

**★★ E5.4 — Sparse capacity, two ways.** Derive $p_{\max}\sim C/(f\log(1/f))$ (a) by information counting and (b) by redoing the signal-to-noise calculation with the covariance rule. Then evaluate for rat CA3 and for a fly mushroom body ($C\approx6$ claw inputs per Kenyon cell — what goes wrong?).

<details markdown="1"><summary>Solution</summary>

**(a) Information counting.** A pattern with $fN$ of $N$ units active carries $H=N\,H_2(f)\approx Nf\log_2(1/f)$ bits for small $f$ (using $H_2(f)=-f\log_2f-(1-f)\log_2(1-f)\approx f\log_2(1/f)+f\log_2 e$). If the network's synapses store $O(1)$ bits each and each neuron has $C$ of them, total capacity is $O(NC)$ bits. Then $p_{\max}Nf\log_2(1/f)\lesssim \kappa NC$, giving $p_{\max}\lesssim \kappa C/(f\log_2(1/f))$.

**(b) Signal-to-noise.** With $J_{ij}=\frac{1}{Nf(1-f)}\sum_\mu(\xi_i^\mu-f)(\xi_j^\mu-f)$ and $\sigma=\xi^\nu$, the field on unit $i$ is
$$h_i = \frac{1}{Nf(1-f)}\sum_\mu(\xi_i^\mu-f)\sum_{j}(\xi_j^\mu-f)\xi_j^\nu .$$
Signal ($\mu=\nu$): $\sum_j(\xi_j^\nu-f)\xi_j^\nu = fN(1-f)$, giving signal $=(\xi_i^\nu-f)$ — separation $1$ between active and inactive units.
Noise ($\mu\ne\nu$): $\sum_j(\xi_j^\mu-f)\xi_j^\nu$ has mean 0 and variance $fN\cdot f(1-f)$ (only the $fN$ units with $\xi_j^\nu=1$ contribute, each with variance $f(1-f)$). Then the $\mu$-sum over $p-1$ patterns, weighted by $(\xi_i^\mu-f)$ of variance $f(1-f)$, gives
$$\mathrm{Var}[h_i^{\rm noise}] = \frac{p\,f(1-f)\cdot f^2N(1-f)}{N^2f^2(1-f)^2} = \frac{p f}{N}\cdot\frac{1}{1}\cdot\frac{1}{1}\;\sim\;\alpha f\; \text{(up to }O(1)).$$
So the SNR is $1/\sqrt{\alpha f}$ — better than the dense case by $1/\sqrt f$, giving $p\sim N/f$. The extra $\log(1/f)$ comes from requiring a *fixed number of errors per pattern* rather than a fixed per-bit error rate: with only $fN$ active units, tolerating a fixed error *fraction of the active set* demands $Q(1/\sqrt{\alpha f})\lesssim \epsilon f$, and inverting the Gaussian tail supplies the logarithm. The careful version is Tsodyks & Feigel'man (1988) / Treves & Rolls (1991).

**Rat CA3:** $C=1.2\times10^4$, $f=0.03$, $k=0.25$ gives $p\approx3\times10^4$ (as in §2.4).

**Mushroom body:** $C\approx6$ (claws per Kenyon cell) with $f\approx0.05$ gives $p\approx 0.25\times6/(0.05\times3.0)\approx10$. Ten memories. The formula is telling you, correctly, that **the mushroom body is not an autoassociative memory** — its Kenyon cell layer has essentially no recurrence, and its capacity as a *classifier* (which is what it is: KC→MBON with plastic output synapses) is governed by an entirely different calculation, the perceptron capacity of [C1](C1-expansion-and-sparsening.md). Applying the CA3 formula to a mushroom body is a category error, and noticing that is the point of the exercise.
</details>

**★★★ E5.5 — Dense associative memory exponent.** Derive $p_{\max}\sim N^{n-1}$ for $E=-\sum_\mu(\boldsymbol\xi^\mu\!\cdot\!\boldsymbol\sigma)^n$, and find the capacity for $F(x)=e^{\beta x}$.

<details markdown="1"><summary>Solution</summary>

**Polynomial case.** Let $\boldsymbol\sigma=\boldsymbol\xi^1$ and consider flipping unit $i$. Write $h^\mu = \sum_{j\ne i}\xi_j^\mu\sigma_j$. The energy gap is
$$\Delta_i=\sum_\mu\Big[F(h^\mu+\xi_i^\mu)-F(h^\mu-\xi_i^\mu)\Big]\approx 2\sum_\mu \xi_i^\mu F'(h^\mu),$$
valid when $F$ varies slowly on scale 1. Stability of unit $i$ requires $\xi_i^1\Delta_i>0$.
- $\mu=1$: $h^1=N-1\approx N$, contribution $2\xi_i^1F'(N)=2n\xi_i^1N^{n-1}$. Multiply by $\xi_i^1$: signal $=2nN^{n-1}$.
- $\mu\ne1$: $h^\mu$ is a sum of $N$ random signs, so $h^\mu\sim\mathcal N(0,N)$, $|h^\mu|\sim\sqrt N$, and $F'(h^\mu)=n(h^\mu)^{n-1}\sim nN^{(n-1)/2}$ with random sign (for even $n-1$ the sign comes from $\xi_i^\mu$; either way the terms are zero-mean and independent). Summing $p-1$ of them: standard deviation $\approx 2n\sqrt{p}\,N^{(n-1)/2}c_n$ with $c_n=\sqrt{\mathbb E[Z^{2(n-1)}]}=\sqrt{(2n-3)!!}$ for $Z\sim\mathcal N(0,1)$.

Stability with high probability needs signal $\gtrsim$ a few noise standard deviations:
$$N^{n-1}\gtrsim \sqrt{p}\,N^{(n-1)/2}\;\Longrightarrow\; p\lesssim N^{n-1}.$$
More carefully, requiring per-bit error $\le1/N$ costs a factor $\log N$: $p_{\max}\asymp N^{n-1}/\big(c_n^2\log N\big)$. For $n=2$: $p\sim N/\log N$, which is the known rigorous result for the classical model with all patterns stable.

**Exponential case.** $F(x)=e^{\beta x}$, $F'=\beta e^{\beta x}$. Signal $\sim \beta e^{\beta N}$. Noise: $\sqrt{p}\cdot\beta\sqrt{\mathbb E[e^{2\beta h}]}=\sqrt p\,\beta e^{2\beta^2N}$ using $\mathbb E[e^{2\beta h}]=e^{2\beta^2N}$ for $h\sim\mathcal N(0,N)$. Stability:
$$e^{\beta N}\gtrsim \sqrt p\,e^{2\beta^2N}\;\Longrightarrow\; \log p \lesssim 2N(\beta-2\beta^2).$$
Maximize $\beta-2\beta^2$ at $\beta=1/4$, giving $\log p\lesssim N/4$, i.e. $p\lesssim e^{N/4}$. This reproduces the exponential scaling; the sharp constant ($p\sim2^{N/2}$) requires the more careful large-deviations treatment of Demircigil et al. (2017), because the Gaussian moment bound is loose in the tail that matters.

**The moral.** Capacity is not a property of "associative memory" — it is a property of the *shape of the energy as a function of overlap*. Sharper functions separate the target from the distractors more aggressively. Softmax is the sharpest useful choice, and that is why attention works.
</details>

**★★ E5.6 — Attention is retrieval.** Starting from $E(\boldsymbol\xi)=-\beta^{-1}\mathrm{lse}(\beta,X^\top\boldsymbol\xi)+\tfrac12\|\boldsymbol\xi\|^2$, derive the update $\boldsymbol\xi^{\rm new}=X\,\mathrm{softmax}(\beta X^\top\boldsymbol\xi)$ and state the condition for one-step exact retrieval.

<details markdown="1"><summary>Solution</summary>

Write $\mathrm{lse}(\beta,\mathbf z)=\beta^{-1}\log\sum_\mu e^{\beta z_\mu}$ so $E(\boldsymbol\xi)=-\mathrm{lse}(\beta,X^\top\boldsymbol\xi)+\tfrac12\|\boldsymbol\xi\|^2$. Split $E=E_{\rm cave}+E_{\rm vex}$ with $E_{\rm cave}=-\mathrm{lse}(\cdot)$ (concave, since lse is convex) and $E_{\rm vex}=\tfrac12\|\boldsymbol\xi\|^2$ (convex). The concave–convex procedure iterates
$$\nabla E_{\rm vex}(\boldsymbol\xi^{t+1}) = -\nabla E_{\rm cave}(\boldsymbol\xi^{t}).$$
Now $\nabla E_{\rm vex}(\boldsymbol\xi)=\boldsymbol\xi$, and
$$\nabla_{\boldsymbol\xi}\,\mathrm{lse}(\beta,X^\top\boldsymbol\xi) = X\,\mathrm{softmax}(\beta X^\top\boldsymbol\xi),$$
since $\partial_{z_\mu}\mathrm{lse}=\mathrm{softmax}(\beta\mathbf z)_\mu$ and $\partial\mathbf z/\partial\boldsymbol\xi=X^\top$. Hence
$$\boxed{\boldsymbol\xi^{t+1}=X\,\mathrm{softmax}(\beta X^\top\boldsymbol\xi^{t})}$$
CCCP guarantees monotone decrease of $E$, so this is a genuine retrieval dynamics, not just a formal resemblance.

**Attention.** Given a query matrix $Q$, keys $K=XW_K$, values $V=XW_V$: $\mathrm{softmax}(QK^\top/\sqrt{d})V$ is the same operation with $\beta=1/\sqrt d$ and separate linear maps on keys and values. Ramsauer et al.'s point is that a transformer attention layer *is* one CCCP step of a modern Hopfield network with $\beta=1/\sqrt{d}$.

**One-step exact retrieval.** Define the separation of pattern $\mu$ as $\Delta_\mu = \boldsymbol x_\mu\!\cdot\!\boldsymbol x_\mu - \max_{\nu\ne\mu}\boldsymbol x_\mu\!\cdot\!\boldsymbol x_\nu$. If the query lies within a ball around $\boldsymbol x_\mu$ and $\Delta_\mu$ is large enough that
$$\beta\Delta_\mu \gtrsim \log\big(2(M-1)M\beta\,R^2\big)\quad (R=\max_\mu\|\boldsymbol x_\mu\|),$$
then the softmax is dominated by component $\mu$ to within exponentially small error, and $\boldsymbol\xi^{1}=\boldsymbol x_\mu + O(e^{-\beta\Delta_\mu})$. Retrieval is one step, with exponentially small error — nothing like the many-sweep relaxation of the classical model. For random patterns on the sphere in $\mathbb R^d$, $\Delta_\mu\sim d - \sqrt{d\log M}$, so $M$ can be exponential in $d$ while $\beta\Delta_\mu$ stays large. That is the exponential capacity, restated geometrically.
</details>

**★★ E5.7 — Computational: separation before storage.** Using the Section 3 code, store $p$ patterns with pairwise correlation $c$ (rather than independent ones) and measure how capacity degrades with $c$. Then insert a random expansion + threshold stage and show it restores capacity.

<details markdown="1"><summary>Solution</summary>

Generate correlated patterns as $\xi^\mu = \mathrm{sgn}(\sqrt{c}\,\mathbf g_0 + \sqrt{1-c}\,\mathbf g_\mu)$ with $\mathbf g$'s i.i.d. standard normal; then $\mathbb E[\xi_i^\mu\xi_i^\nu]=\frac{2}{\pi}\arcsin c$ for $\mu\ne\nu$.

**Result without separation.** The crosstalk analysis changes: the cross terms $\sum_j\xi_j^\mu\xi_j^\nu$ no longer have mean zero but mean $N\rho$ with $\rho=\frac2\pi\arcsin c$. Then $\mathrm{Var}[C]\approx\alpha(1+N\rho^2)$, so the effective load is $\alpha_{\rm eff}=\alpha(1+N\rho^2)$ and capacity collapses as
$$\alpha_c^{\rm eff}\approx \frac{0.138}{1+N\rho^2}.$$
For $N=400$ and even $\rho=0.05$ this is $0.138/2 = 0.069$ — half the capacity from a 5% correlation. Correlated patterns are catastrophic, and this is precisely why a separation stage is not optional.

**With separation.** Map $\mathbf x\mapsto \mathbf y=\Theta(W\mathbf x-\theta)$ with $W\in\mathbb R^{N'\times N}$ random, $N'=5N$, and $\theta$ set so that $f=0.05$ of the $y$'s are active. Then store $\mathbf y$ instead. Measured output correlation drops from $\rho\approx0.05$ to $\rho_{\rm out}\lesssim0.005$ (use the $\Phi_2$ formula of §2.6 to predict it), and the recovered capacity in units of the *expanded* network is close to the sparse bound $k/(f\ln(1/f))\approx 1.7$ patterns per neuron — vastly more than $0.138$.

**What you have just shown.** DG is not decoration. In a network fed by naturalistic, correlated cortical input, the autoassociator's capacity without a separation stage is a small fraction of its nominal value, and the expansion recovers it. The two anatomically weird facts about the hippocampus — a $5\times$ expansion into a hyper-sparse layer, and a hugely recurrent layer downstream of it — are two halves of a single algorithm.
</details>

---

## 5. Reading path

- **Marr (1971)**, *Simple memory: a theory for archicortex* (Phil. Trans. R. Soc. B) — read it for: the conversion done before any of the tools existed. Read the actual paper, not a summary of it.
- **Hopfield (1982)**, *Neural networks and physical systems with emergent collective computational abilities* (PNAS) — read it for: the energy function, and the four pages that created a field.
- **Amit, Gutfreund & Sompolinsky (1985 PRL; 1987 Ann. Phys.)**, on the storage capacity and statistical mechanics of the Hopfield model near saturation — read them for: $\alpha_c=0.138$, the phase diagram, and the reaction term.
- **Hertz, Krogh & Palmer (1991)**, *Introduction to the Theory of Neural Computation*, ch. 2–3 — read it for: the cleanest exposition of §2.2–2.3, including the $P_{\rm err}$ table.
- **Tsodyks & Feigel'man (1988)**, on storage capacity of networks with low activity levels — read it for: the sparse-coding capacity result in its original form.
- **Treves & Rolls (1991)**, *What determines the capacity of autoassociative memories in the brain?* (Network) — read it for: capacity as a function of *connections per neuron*, the number that matters anatomically.
- **Treves & Rolls (1992)**, on why CA3 needs two distinct input systems (Hippocampus) — read it for: the sharpest algorithmic argument here — mossy fibres write, perforant path reads.
- **Nakazawa et al. (2002)** and **McHugh et al. (2007)** (both Science) — read together for: the double dissociation of CA3 completion and DG separation by region-specific NMDA receptor deletion.
- **Leutgeb, Leutgeb, Moser & Moser (2007)** (Science) and **Neunuebel & Knierim (2014)** (Neuron) — read for: separation and completion measured in the same animals.
- **Krotov & Hopfield (2016)**, *Dense associative memory for pattern recognition* (NeurIPS) — read it for: the energy's *shape* sets the capacity; $p\sim N^{n-1}$.
- **Demircigil, Heusel, Löwe, Upgang & Vermet (2017)** (J. Stat. Phys.) — read it for: the exponential-capacity result made rigorous.
- **Ramsauer et al. (2020)**, *Hopfield networks is all you need* — read it for: the continuous-state energy, the CCCP update, and the identification with transformer attention. Section 2 and Appendix A are the substance.
- **Krotov & Hopfield (2021)** (ICLR) — read it for: higher-order energies realized with pairwise synapses and a hidden memory layer — then ask yourself what circuit that looks like.
- **Babadi & Sompolinsky (2014)** (Neuron) and **Litwin-Kumar, Harris, Axel, Sompolinsky & Abbott (2017)** (Neuron) — read for: the theory of the expansion stage; the formal link back to [C1](C1-expansion-and-sparsening.md).

---

## 6. Open problems and what would settle them

**1. Is CA3 actually an attractor network, or does it just look like one?** Everything usually cited as evidence — pattern completion from partial cues, attractor-like transitions in morphed environments (Wills et al. 2005 versus Leutgeb et al. 2005, which famously disagreed) — is consistent with a feedforward network with the right tuning. **What would settle it:** the fly ring-attractor playbook of [C4](C4-continuous-attractors.md). Perturb the population off the putative attractor state optogenetically and measure the *return trajectory*. An attractor has a characteristic anisotropic relaxation: fast along transverse directions, slow or absent along the manifold. Nobody has done this in CA3 with adequate spatial precision. It is the decisive experiment and it is now technically feasible.

**2. Is the learning rule outer-product Hebbian?** The capacity theory assumes $\Delta J_{ij}\propto\xi_i\xi_j$ summed over patterns with no interference correction. Real CA3 plasticity is spike-timing dependent, asymmetric, and *saturating*, and palimpsest schemes with bounded synapses have capacity $O(\sqrt N)$ or worse rather than $O(N)$ (Fusi and colleagues on cascade models). **What would settle it:** measure $\Delta J$ as a function of pre/post activity in CA3 recurrent collaterals *in vivo* during encoding, and ask whether the covariance rule's mean subtraction — essential for sparse capacity — has any physiological counterpart.

**3. Where does high capacity live?** §2.5 argues that dense associative memories look architecturally like expansion layers, not recurrent ones. If so, hippocampal capacity may reside in DG's million granule cells rather than CA3's recurrence, with CA3 doing completion only. **What would settle it:** compare information-theoretic capacity estimates of DG and CA3 representations, and test whether DG lesions reduce the *number* of storable memories rather than their *discriminability*.

**4. Do brains implement softmax retrieval?** The Hopfield↔attention correspondence is currently a fact about mathematics, not about brains. A softmax over stored patterns requires population-wide normalization, which is exactly what divisive normalization does ([C8](C8-divisive-normalization.md)). **What would settle it:** find a circuit where normalized, temperature-controlled competition among pattern-selective units produces one-step retrieval. My candidate is the mushroom body output layer, where KC→MBON plasticity plus MBON–MBON competition has roughly this form. Speculative, and the most interesting thing in the unit.

**5. What is the level-1 problem, really?** "Store and retrieve episodes" is a description, not a specification. The modern answer — complementary learning systems, and the successor-representation view of [C6](C6-grid-cells.md) — is that the hippocampus supports rapid learning of a predictive model the neocortex slowly distils. If that is right, the algorithm is closer to experience replay for a model-based RL agent, and the capacity calculations answer a question nobody should have asked. **What would settle it:** whether hippocampal representations read out during replay are organized by *predictive utility* rather than recency and similarity.
