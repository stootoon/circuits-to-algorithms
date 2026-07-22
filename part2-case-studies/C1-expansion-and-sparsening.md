---
title: "C1 · Expansion & sparsening"
parent: Case studies
nav_order: 1
---

# C1 — Expansion + sparsening → random projection, LSH, and kernels

> **Circuit.** Insect antennal lobe projection neurons fan out onto a vastly larger population of Kenyon cells (~50 PN channels → ~2,000 KCs in *Drosophila*; ~830 PNs → ~50,000 KCs in locust), each KC sampling a handful of PNs, with a single giant inhibitory interneuron (GGN in locust, APL in *Drosophila*) closing a global negative-feedback loop.
> **Primitive extracted.** A sparse binary random projection into a high-dimensional space followed by a winner-take-all threshold — i.e. a locality-sensitive hash whose output is simultaneously a random-feature map that linearises downstream classification.
> **Status.** Well established as a *sufficient* reading; the LSH framing is a genuine contribution of theory to biology and is largely uncontested. What remains contested is whether the projection is truly random (structured sampling has been reported) and where on the sparseness–noise tradeoff real circuits actually sit.
> **Structures thread.** `structures/S-04-manifold-capacity.md` (Cover counting and its manifold generalisation — the capacity story below is the scalar special case), `structures/S-06-hyperbolic-odor-space.md` (the hash is only useful if the metric it preserves is the metric that matters — if odour space is hyperbolic, Euclidean-cosine LSH is preserving the wrong thing), `structures/S-03-combinatorial-threshold-linear-networks.md` (the GGN/APL loop is a threshold-linear network with uniform global inhibition; its fixed-point combinatorics *is* the winner-take-all). Start from `structures/README.md`.

## 1. The phenomenon

Three facts, all robust, all measured.

**Expansion.** The insect olfactory pathway performs one of the most extreme dimensional expansions in any nervous system. In *Drosophila* roughly 50 glomerular channels project to ~2,000 Kenyon cells; each KC extends ~6–7 dendritic "claws," each claw receiving input from a single PN bouton. In locust the numbers are ~830 PNs to ~50,000 KCs, with each KC sampling ~10–20 PNs. The fan-out ratio is 40× in fly, 60× in locust.

**Randomness (approximately).** Caron, Ruta, Abbott & Axel (2013) traced the PN inputs to individual *Drosophila* KCs and found no glomerular structure beyond what you would expect from independent sampling with unequal glomerular "sampling probabilities." Later work using the hemibrain connectome (Zheng et al., 2022) found statistically detectable structure — some glomerular pairs co-occur more or less often than chance, and there is a food-odour-biased community — so "random" is a good null model that is not exactly true. Hold onto that: it is the sharpest empirical wedge into the algorithmic claim.

**Sparsening.** Perez-Orive et al. (2002) showed that locust KCs respond to odours with an average of ~1 spike, in a narrow window locked to the 20–30 Hz local field potential, and that only a few percent of KCs respond to any given odour. Turner, Bazhenov & Laurent (2008) found the same in fly. The sparseness is actively enforced: Papadopoulou et al. (2011) showed that the locust GGN, a single non-spiking wide-field neuron, receives excitation from the whole KC population and returns inhibition to all of it, and that hyperpolarising GGN *destroys* KC sparseness and odour specificity. Lin et al. (2014) showed the *Drosophila* homologue APL does the same, and that blocking it degrades learned fine discrimination.

So: expand, then subtract off a term proportional to total activity and keep only what survives. That is the circuit. Now extract the algorithm.

## 2. The conversion

### 2.1 The map

Write the PN population vector as $x \in \mathbb{R}^N_{\ge 0}$, $N \approx 50$. Let $A \in \{0,1\}^{M\times N}$ be the PN→KC adjacency matrix with $M \gg N$ and each row containing exactly $s$ ones ($s \approx 6$), placed uniformly at random. The KC drive is $m = Ax$. The GGN/APL loop subtracts a term proportional to $\mathbf{1}^\top m$ and the KC threshold keeps a fixed fraction $\alpha$ of the population:

$$
h(x) \;=\; \mathrm{WTA}_{k}\!\big(Ax\big), \qquad k = \alpha M, \quad h \in \{0,1\}^M .
$$

The **tag** of $x$ is the index set $T(x) = \{ j : h_j(x) = 1\}$. Dasgupta, Stevens & Navlakha (2017) proposed exactly this as a locality-sensitive hash. Three things need proving: (i) that the nonlinearity is doing the work, not the expansion; (ii) that the tag is locality-sensitive, with a quantitative similarity-preservation law; (iii) that the map buys linear separability downstream.

### 2.2 Expansion alone is worthless — a one-line proof

Suppose the KC layer were linear, $y = Ax$. Take any labelling of patterns $\{x^\mu\}$ that is linearly separable in KC space: there exist $w \in \mathbb{R}^M, b$ with $\mathrm{sign}(w^\top A x^\mu + b) = \ell^\mu$. But $w^\top A x = (A^\top w)^\top x$, so $u = A^\top w \in \mathbb{R}^N$ separates the same labelling in *PN* space. Contrapositive: anything not separable at the PN layer is not separable at the KC layer either.

$$
\boxed{\text{A linear map cannot increase the number of linearly separable dichotomies.}}
$$

Every drop of computational value in the expansion comes from the threshold. This is worth internalising because it kills a whole family of hand-waving about "high-dimensional representations" that never says where the nonlinearity is.

### 2.3 The KC layer as a random feature map / kernel

Because $A$ has i.i.d. rows, the KC layer is a Monte-Carlo estimator of a kernel. For the simpler sign nonlinearity $\phi_w(x) = \mathrm{sign}(w^\top x)$ with $w \sim \mathcal{N}(0,I)$,

$$
\mathbb{E}_w\big[\phi_w(x)\,\phi_w(y)\big] \;=\; 1 - \frac{2}{\pi}\arccos\!\left(\frac{x^\top y}{\|x\|\|y\|}\right),
$$

which is the Goemans–Williamson / Charikar SimHash identity: the probability that a random hyperplane separates $x$ and $y$ is the angle between them over $\pi$. With a nonzero threshold you get Cho & Saul's arc-cosine kernels; with a rectifier instead of a step you get the $n=1$ arc-cosine kernel, which is the infinite-width ReLU kernel. So the KC layer is a Rahimi–Recht random-features approximation to a fixed nonlinear kernel, and $M$ controls the Monte-Carlo error, not the kernel. That is a level-2 statement with real content: **the mushroom body implements a fixed kernel; learning at KC→MBON is kernel regression.**

### 2.4 Why the fly's tag is locality-sensitive: do the calculation

Now the sparse binary case. Let $A_{ji} \sim \mathrm{Bernoulli}(p)$ i.i.d. with $p = s/N$ (fixed row sums change nothing to leading order). For two odours $x, y$, consider the pair of drives $(m_j(x), m_j(y))$ across the KC index $j$. Then

$$
\mathbb{E}_j[m_j(x)] = p\,\textstyle\sum_i x_i, \qquad
\mathrm{Cov}_j\big(m(x), m(y)\big) = \sum_i \mathrm{Var}(A_{ji})\,x_i y_i = p(1-p)\,x^\top y ,
$$

and $\mathrm{Var}_j(m(x)) = p(1-p)\|x\|^2$. Hence the correlation coefficient of the two drives across the KC population is

$$
\rho \;=\; \frac{x^\top y}{\|x\|\,\|y\|},
$$

**independent of $p$**. Sparse binary sampling preserves cosine similarity exactly in expectation, with no dependence on how sparse the sampling is. That is the whole "random projection" content, and it is why an anatomy that samples 6 out of 50 inputs is not throwing information away.

One caveat with teeth. The GGN/APL loop removes the population mean before thresholding. After mean removal the relevant $\rho$ is not the raw cosine but the **Pearson correlation across PN channels**, $\rho = \mathrm{corr}(x,y)$. For non-negative firing rates these differ enormously: two statistically independent rate vectors have raw cosine $\approx (\mathbb{E}|z|)^2/\mathbb{E}[z^2] \approx 0.64$ but Pearson correlation $\approx 0$. Normalisation is what makes the hash discriminative at all — see C8, where this is the same computation seen from the other side.

Now the winner-take-all. By the CLT over $N$ terms, $(m_j(x), m_j(y))$ is approximately bivariate normal with correlation $\rho$. Standardise; the top-$\alpha$ threshold is $t = \Phi^{-1}(1-\alpha)$. Then

$$
\Pr\big[j \in T(x) \ \wedge\ j \in T(y)\big] \;=\; \Phi_2(t,t;\rho),
$$

the bivariate Gaussian orthant probability, so the expected fractional tag overlap is

$$
\boxed{\ \mathcal{O}(\rho) \;=\; \frac{|T(x)\cap T(y)|}{\alpha M} \;=\; \frac{\Phi_2(t,t;\rho)}{\alpha}\ }
$$

$\Phi_2(t,t;\rho)$ is strictly increasing in $\rho$ for fixed $t$ (Slepian's inequality), with $\mathcal{O}(0)=\alpha$ and $\mathcal{O}(1)=1$. **That is the locality-sensitivity proof**, and it is quantitative rather than asymptotic hand-waving.

The small-$\alpha$ asymptotics are more memorable. Using $\Pr[X>t, Y>t] \approx \exp\!\big(-t^2/(1+\rho)\big)$ and $\alpha \approx \exp(-t^2/2)$,

$$
\mathcal{O}(\rho) \;\approx\; \alpha^{\frac{2}{1+\rho} - 1} \;=\; \alpha^{\frac{1-\rho}{1+\rho}} .
$$

The exponent $\frac{1-\rho}{1+\rho}$ is exactly the LSH quality exponent that appears in the Indyk–Motwani complexity of approximate nearest-neighbour search. It is not an analogy: the mushroom body's sparseness parameter $\alpha$ is the LSH amplification knob, and sparser codes give sharper similarity discrimination at the cost of more false negatives. (Numerically, at $\alpha=0.05$ this asymptotic overshoots the exact orthant probability by 20–50% at intermediate $\rho$; the code in §3 shows both.)

Convert to correlations: the tags are binary with density $\alpha$, so the *centred* correlation between tags is $(\mathcal{O}-\alpha)/(1-\alpha)$. At $\rho = 0.43$, $\alpha=0.05$, this is $(0.24-0.05)/0.95 \approx 0.20$. **The KC layer more than halves pairwise correlations.** Expansion-plus-WTA is a decorrelator. Keep that in view for C2 and C8: three different circuits, one recurring objective.

### 2.5 Contrast with classical LSH — the fly is doing something genuinely different

Classical SimHash: draw $M \ll N$ dense Gaussian vectors, keep $\mathrm{sign}(w_j^\top x)$, get a compact $M$-bit key; the whole point is Johnson–Lindenstrauss **dimensionality reduction**, and the design target is a short key for a hash table. The fly does the opposite on two axes and the same on a third:

| | classical LSH | fly |
|---|---|---|
| dimension | $M \ll N$ (reduce) | $M \gg N$ (expand) |
| projection | dense Gaussian, $\pm$ | sparse binary, non-negative |
| output | $M$ dense bits | $\alpha M$ set bits out of $M$ |
| bits in the tag | $M$ | $\log_2\binom{M}{\alpha M} \approx \alpha M \log_2(e/\alpha)$ |

Both are LSH. The fly's tag has comparable *information content* to a short dense key ($\alpha M \log_2(e/\alpha) \approx 100 \times 5.6 \approx 560$ bits for $M=2000, \alpha=0.05$) but pays for it with a huge, mostly-silent population. Why?

My position: because the fly is not building a hash table. It is building a substrate for **one-shot, local, Hebbian association**. A dense $M$-bit key requires a downstream reader that can compute Hamming distance; a sparse index set can be associated to a valence by depressing $\sim 100$ synapses onto a single MBON, in one trial, with no credit-assignment problem. Sparse-and-expanded is what you choose when the readout has to be a *linear* function learnable by local plasticity, not when memory is the constraint. Dimensionality reduction and dimensionality expansion solve different downstream problems; only the middle step — the random projection — is shared.

### 2.6 Capacity: Marr, Albus, and Cover

Marr (1969) and Albus (1971) proposed the same architecture for cerebellum: mossy fibres → an enormous granule cell layer with ~4 inputs each and Golgi-cell-enforced sparseness → Purkinje cells learning a linear readout. The mushroom body is the same three-layer motif. The formal justification is Cover's function-counting theorem: for $P$ points in general position in $\mathbb{R}^M$, the number of linearly separable dichotomies is

$$
C(P,M) \;=\; 2\sum_{k=0}^{M-1}\binom{P-1}{k},
$$

so the fraction of random dichotomies that are separable is $C(P,M)/2^P$, which is $\approx 1$ for $P < 2M$ and falls off sharply above it. The capacity is $P^\star = 2M$. Expanding $N=50 \to M=2000$ raises capacity from ~100 patterns to ~4000. The §3 code shows this transition crisply.

But "general position" is doing work. Correlated patterns are not in general position in the relevant sense — the margin collapses long before $P = 2M$. This is why §2.4 matters: the WTA decorrelation is what pushes the KC representation toward the regime where Cover's bound is achievable.

Litwin-Kumar, Harris, Axel, Sompolinsky & Abbott (2017) made this quantitative. Define the dimension of the KC representation by the participation ratio $D = (\sum_i \lambda_i)^2 / \sum_i \lambda_i^2$ of its covariance eigenvalues. For sparse random connectivity with in-degree $K$, $D$ is non-monotonic in $K$: small $K$ makes KCs too similar to individual PNs (low $D$ from input correlations passing through), large $K$ makes all KCs sample nearly the same average (low $D$ from convergence to the mean). The optimum sits near $K \approx 7$ for realistic parameters. *Drosophila* KCs have 6–7 claws. That is one of the cleanest quantitative level-1↔level-2 contacts in systems neuroscience, and you should treat it as the standard to which other conversions in this course are held.

### 2.7 The catch: Babadi & Sompolinsky's tradeoff

Expansion improves separability of *noiseless* patterns. Real PN responses are noisy, and the threshold amplifies noise in a specific, derivable way.

Let KC $j$ have standardised drive $u_j$ and let trial-to-trial input noise contribute a perturbation with standard deviation $\epsilon$ in the same units. KC $j$ flips state iff $|u_j - t| \lesssim |\delta u_j|$, so to first order the probability of a flip is

$$
q \;\approx\; \phi(t)\,\epsilon\,\sqrt{2/\pi}\ \propto\ \phi(t)\,\epsilon .
$$

The number of *active* KCs is $\alpha M = M\,\Phi(-t)$. The corruption that matters is the fraction of the tag that is wrong, i.e. the ratio

$$
\frac{q}{\alpha} \;\approx\; \epsilon\,\frac{\phi(t)}{\Phi(-t)} \;\xrightarrow[t \to \infty]{}\; \epsilon\, t \;=\; \epsilon\sqrt{2\ln(1/\alpha)} .
$$

using the Mills ratio $\phi(t)/\Phi(-t) \to t$. So:

$$
\boxed{\ \text{relative noise in a sparse tag grows as } \sqrt{2\ln(1/\alpha)} \text{ as the code sparsens.}\ }
$$

Meanwhile the benefit — decorrelation — improves only as $\log \mathcal{O} \propto \frac{1-\rho}{1+\rho}\log\alpha$, i.e. also logarithmically. Two logarithms fighting, with different prefactors set by $\epsilon$ and $\rho$: hence an interior optimum in $\alpha$, and hence Babadi & Sompolinsky's (2014) central result that expansive sparse recoding improves discriminability up to a point and then degrades it, with the optimum depending on the input SNR. Their fuller treatment shows that the achievable signal-to-noise ratio in the expanded layer *saturates* as $M \to \infty$: you cannot manufacture information by expansion, you can only redistribute it into a form a linear readout can use.

The corollary you should carry: **the observed sparseness of a code is a measurement of the input SNR**, given the theory. A circuit operating at $\alpha = 5\%$ is telling you something about $\epsilon$.

### 2.8 Mixed selectivity and shattering dimensionality

Rigotti et al. (2013) made the same argument for prefrontal cortex without a random-projection anatomy in hand. Their observable is **shattering dimensionality**: the fraction of the $2^P$ possible dichotomies of $P$ experimental conditions that a linear decoder can implement. Pure ("diagonal") selectivity — neurons tuned to single task variables — supports only the dichotomies that are linear in those variables. Nonlinear mixed selectivity — neurons responding to conjunctions — makes almost all dichotomies decodable, which is exactly Cover's count with $M$ effectively large.

KCs are the extreme case of nonlinear mixed selectivity: each is a random conjunction of ~6 glomerular channels, thresholded. The mushroom body is what a mixed-selectivity layer looks like when evolution has optimised it and you can count the wires. Bernardi et al. (2020) later showed the tension: maximal shattering dimensionality destroys *abstraction* (the ability to generalise across conditions), and real cortex sits at a compromise. The insect answer to that tension is architectural — the mushroom body shatters, while the lateral horn, receiving the same PNs through stereotyped wiring, preserves innate categories. Two readouts, two points on the abstraction–shattering curve, one input.

## 3. Worked example / model to build

Two experiments in one script: the LSH similarity-preservation law against theory, and the Cover capacity transition.

```python
import numpy as np
from math import erf, sqrt, pi
rng = np.random.default_rng(0)
N, M, s, alpha = 50, 2000, 6, 0.05        # PNs, KCs, PNs per KC, WTA fraction

A = np.zeros((M, N))                       # sparse binary projection
for j in range(M):
    A[j, rng.choice(N, s, replace=False)] = 1.0

def tag(X):                                # expansion -> normalise (APL) -> top-k WTA
    Y = X @ A.T
    Y = (Y - Y.mean(1, keepdims=True)) / (Y.std(1, keepdims=True) + 1e-12)
    k = int(alpha * M)
    return (Y >= np.partition(Y, -k, axis=1)[:, -k][:, None]).astype(float)

Phi  = np.vectorize(lambda z: 0.5*(1 + erf(z/sqrt(2))))
def Phi2(t, r, n=20000):                   # P(X>t, Y>t) for standard bivariate normal
    u = np.linspace(t, t+12, n)
    return float(np.trapezoid(np.exp(-u**2/2)/sqrt(2*pi) *
                              Phi((r*u - t)/sqrt(max(1-r*r, 1e-15))), u))
t = 1.6449                                 # Phi^{-1}(1-alpha), alpha = 0.05

print(" corr(x,y)  tag overlap   Phi2(t,t;r)/alpha   alpha^((1-r)/(1+r))")
for rho in [0.0, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99]:
    x = rng.standard_normal((400, N)); z = rng.standard_normal((400, N))
    y = np.abs(rho*x + sqrt(1-rho**2)*z) + 0.5
    x = np.abs(x) + 0.5                    # non-negative PN firing rates
    xc, yc = x - x.mean(1, keepdims=True), y - y.mean(1, keepdims=True)
    r = float(np.mean((xc*yc).sum(1) /
              (np.linalg.norm(xc,axis=1)*np.linalg.norm(yc,axis=1))))
    ov = float(((tag(x)*tag(y)).sum(1)/(alpha*M)).mean())
    print("   %6.3f     %6.3f          %6.3f              %6.3f"
          % (r, ov, Phi2(t, r)/alpha, alpha**((1-r)/(1+r))))

def sep_frac(F, P, trials=20):             # fraction of random dichotomies separable
    F = np.hstack([F[:P], np.ones((P,1))]); ok = 0
    for _ in range(trials):
        lab = rng.choice([-1., 1.], P)
        w = np.linalg.lstsq(F, lab, rcond=None)[0]
        for _ in range(400):
            bad = lab*(F@w) < 1e-9
            if not bad.any(): break
            w += 0.05*(F[bad]*lab[bad, None]).mean(0)
        ok += bool((lab*(F@w) > 0).all())
    return ok/trials

X = np.abs(rng.standard_normal((600, N))); KC = tag(X)
print("\n  P    PN layer (N=50)   KC layer (M=2000)")
for P in [40, 80, 150, 300, 600]:
    print("%4d       %.2f              %.2f" % (P, sep_frac(X, P), sep_frac(KC, P)))
```

**What to look for.** The measured tag overlap should track $\Phi_2(t,t;\rho)/\alpha$ within a few percent across the whole range (you should see $\approx 0.049, 0.073, 0.125, 0.239, 0.422, 0.643, 0.834$ against theory $0.048, 0.067, 0.105, 0.206, 0.370, 0.599, 0.812$) — the small excess is real and comes from KCs sharing PNs, which breaks the independence assumption. The $\alpha^{(1-\rho)/(1+\rho)}$ column visibly overshoots: that is what an asymptotic in $\alpha \to 0$ looks like at $\alpha=0.05$. In the capacity table, the PN layer collapses between $P=40$ and $P=80$ — Cover predicts the transition at $P = 2N = 100$ — while the KC layer is still perfectly separable at $P = 600$. Then push $P$ toward 4000 and watch the KC layer fail at $2M$.

Three modifications worth doing yourself: (i) set the WTA aside and use a fixed absolute threshold — watch the overlap law fall apart when total drive varies across stimuli (this is C8's job); (ii) replace the sparse binary $A$ with a dense Gaussian one — the overlap law is *unchanged*, which is the point of §2.4; (iii) add trial noise to $x$ and sweep $\alpha$ to find the Babadi–Sompolinsky optimum.

## 4. Exercises

**E1 (★).** Show that the expected tag overlap $\mathcal{O}(\rho)$ is strictly increasing in $\rho$ without invoking Slepian's inequality.

<details markdown="1"><summary>Solution</summary>

Write $\Phi_2(t,t;\rho) = \Pr[X>t, Y>t]$ for standard bivariate normal with correlation $\rho$. Plackett's identity states
$$\frac{\partial}{\partial \rho}\Phi_2(a,b;\rho) = \varphi_2(a,b;\rho),$$
the bivariate normal density at $(a,b)$. (Proof: the bivariate normal density satisfies the diffusion-like identity $\partial \varphi_2/\partial\rho = \partial^2\varphi_2/\partial a\, \partial b$; integrate over $[a,\infty)\times[b,\infty)$ and the boundary terms give $\varphi_2(a,b;\rho)$.)

Hence
$$\frac{d\mathcal{O}}{d\rho} = \frac{1}{\alpha}\,\varphi_2(t,t;\rho) = \frac{1}{\alpha}\cdot\frac{1}{2\pi\sqrt{1-\rho^2}}\exp\!\left(-\frac{t^2}{1+\rho}\right) > 0$$
for all $|\rho|<1$. Strictly increasing, and the derivative is largest near $\rho = 1$, which is why the code has good *resolution* for near-identical odours — a desirable property for a discrimination system, and the opposite of what a linear similarity measure gives you.
</details>

**E2 (★★).** A KC has $s$ claws sampling PNs uniformly at random. Suppose one glomerulus is "hot" (high firing) for odour $x$ and all others are silent. Compute the probability that a KC is driven at all, and use it to explain why $s=1$ would be a bad design even though it maximally preserves input sparseness.

<details markdown="1"><summary>Solution</summary>

With $s$ claws drawn without replacement from $N$ PNs, $\Pr[\text{KC touches the hot glomerulus}] = 1 - \binom{N-1}{s}/\binom{N}{s} = s/N$. For $N=50$, $s=6$: $0.12$, so $\approx 240$ of 2000 KCs receive the hot input; the WTA at $\alpha=5\%$ then selects 100 of those 240 by the tie-breaking heterogeneity of synaptic weights.

Now $s=1$. Each KC is a relabelled copy of a single PN, so the KC representation is a *permutation-with-multiplicity* of the PN representation: its covariance has the same eigenvalues as the input covariance (up to multiplicities), the participation-ratio dimension is unchanged, and — by §2.2, since $A$ is now a scaled selection matrix — no dichotomy becomes separable that was not separable at the PN layer, thresholding notwithstanding, because thresholding a single input is a monotone function of a single coordinate. You have paid 2000 neurons for exactly nothing.

The nonlinearity only creates new features when the threshold acts on a *conjunction*, which requires $s \ge 2$. And $s$ cannot be too large either: as $s \to N$ every KC computes (nearly) $\mathbf{1}^\top x$ and all KCs become identical, so dimension collapses again. Hence a maximum at intermediate $s$ — Litwin-Kumar et al.'s $K \approx 7$.
</details>

**E3 (★★).** Derive the number of bits in a fly tag and compare it to a dense SimHash of the same information content. Then argue quantitatively which is cheaper *for a downstream reader implementing a linear readout with one synapse per KC*.

<details markdown="1"><summary>Solution</summary>

Tag entropy: choosing $k=\alpha M$ active units out of $M$ gives at most
$$\log_2\binom{M}{\alpha M} = M H_2(\alpha)\ \text{bits}, \qquad H_2(\alpha) = -\alpha\log_2\alpha - (1-\alpha)\log_2(1-\alpha).$$
For $M=2000, \alpha=0.05$: $H_2(0.05)=0.286$, so $572$ bits. A dense SimHash achieving 572 bits needs 572 projections.

Reader cost. A linear readout needs one weight per input unit: 2000 synapses for the fly, 572 for SimHash. So SimHash wins on synapse count by 3.5×. But now count the *plasticity events per learning trial*: the fly changes only the $\alpha M = 100$ synapses that were active (Hebbian coincidence requires presynaptic activity), whereas the dense code requires updating all 572. And count the *interference*: two random fly tags share $\alpha^2 M = 5$ units on average, so associating a new odour perturbs an old association in only 5 of its 100 relevant synapses — a 5% crosstalk. Two random dense codes share half their bits, so every new association perturbs half of every old one, and you need explicit error-driven learning to undo it.

Quantitatively, for $P$ stored associations the sparse readout's signal-to-interference ratio scales as $\sqrt{\alpha M / (P\alpha^2 M)} = \sqrt{1/(P\alpha)}$, whereas the dense one scales as $\sqrt{1/P}$ — a factor $1/\sqrt\alpha \approx 4.5$ improvement. That is the trade: 3.5× more synapses buys 4.5× less crosstalk and 5.7× fewer plasticity events. For a one-shot learner, sparse-and-expanded wins.
</details>

**E4 (★★).** Cover's theorem gives $P^\star = 2M$ for points in general position. Show that if all $P$ patterns have pairwise correlation exactly $c > 0$ and are otherwise generic, the effective capacity is not reduced at all — and then explain why correlated inputs nevertheless hurt in practice.

<details markdown="1"><summary>Solution</summary>

Let $G = (1-c)I_P + c\,\mathbf{1}\mathbf{1}^\top$ be the Gram matrix of unit-norm patterns. As long as $c<1$, $G$ is full rank, so the patterns are linearly independent whenever $P \le M$ and are in general position for generic perturbations. Cover's count depends only on general position, not on the Gram matrix, so the *count* of separable dichotomies is unchanged.

What changes is the **margin**. For a labelling $\ell$, the max-margin solution has $\kappa = 1/\|w\|$ with $w = \sum_\mu a_\mu \ell^\mu x^\mu$ and $\|w\|^2 = a^\top \tilde G a$ where $\tilde G_{\mu\nu} = \ell^\mu\ell^\nu G_{\mu\nu}$. For a balanced random labelling the rank-one term $c\,\mathbf{1}\mathbf{1}^\top$ contributes $c(\sum_\mu a_\mu\ell^\mu)^2$, which is $O(c\,P)$ for the uniform solution $a_\mu = 1/P$, versus $(1-c)/P$ for the identity part. So $\|w\|^2 \approx (1-c)/P + c\,\varepsilon^2$ and the margin degrades as $c \to 1$ like $\kappa \sim \sqrt{(1-c)/P}$.

In practice you never have exact separability with zero noise: what matters is $\kappa/\sigma_{\text{noise}}$. Correlation shrinks $\kappa$ toward zero while leaving noise untouched, so the *robustly* separable capacity collapses roughly as $(1-c)$ even though the *nominal* capacity does not. This is why the WTA decorrelation in §2.4 is not a luxury — it is what converts nominal capacity into usable capacity.
</details>

**E5 (★★, computational).** Using the §3 code, add i.i.d. Gaussian trial noise of standard deviation $\epsilon\|x\|/\sqrt{N}$ to each PN vector and measure, as a function of $\alpha \in \{0.5, 0.2, 0.1, 0.05, 0.02, 0.01\}$, the *discriminability index* $d = (\mathcal{O}_{\text{same odour, two trials}} - \mathcal{O}_{\text{different odours}})/\sqrt{\mathrm{Var}}$. Find the optimum and check it against §2.7.

<details markdown="1"><summary>Solution</summary>

Implementation: for each $\alpha$, generate 200 odours; for each, two noisy trials; compute within-odour tag overlap $\mathcal{O}_{\text{w}}$ and across-odour overlap $\mathcal{O}_{\text{a}}$; report $d = (\bar{\mathcal{O}}_{\text w} - \bar{\mathcal{O}}_{\text a})/\sqrt{\tfrac12(\mathrm{Var}_{\text w}+\mathrm{Var}_{\text a})}$.

Expected shape. As $\alpha$ decreases from $0.5$: $\mathcal{O}_{\text a} \to \alpha \to 0$, so the *separation* between the two distributions grows — that is the decorrelation benefit, and it is the dominant effect down to moderate sparseness. But $\mathcal{O}_{\text w}$ also falls, because §2.7's relative corruption $\epsilon\sqrt{2\ln(1/\alpha)}$ grows, and simultaneously the *variance* of the overlap estimate grows because it is a binomial fraction over only $\alpha M$ trials: $\mathrm{Var} \sim \mathcal{O}(1-\mathcal{O})/(\alpha M)$, i.e. $\propto 1/\alpha$.

Combining, $d \propto (1 - \epsilon\sqrt{2\ln(1/\alpha)} - \alpha)\sqrt{\alpha M}$ to leading order. Differentiating with respect to $\ln\alpha$ and setting to zero gives the condition
$$\tfrac12\left(1 - \epsilon\sqrt{2\ln(1/\alpha^\star)}\right) = \frac{\epsilon}{\sqrt{2\ln(1/\alpha^\star)}},$$
i.e. writing $L=\sqrt{2\ln(1/\alpha^\star)}$, $\ \epsilon L^2 + 2\epsilon - L = 0$, so $L = (1-\sqrt{1-8\epsilon^2})/(2\epsilon)$. For $\epsilon = 0.2$: $L \approx 2.19$, $\alpha^\star \approx e^{-2.4} \approx 0.09$. For $\epsilon = 0.05$: $L \approx 2.03$, $\alpha^\star \approx 0.13$ — but the optimum is very flat, and the numerics will show a broad plateau from $\alpha \approx 0.02$ to $0.2$. Two conclusions: (i) an interior optimum genuinely exists and sits within a factor of a few of the biological $\alpha \approx 0.05$; (ii) the flatness means you should *not* claim the fly is finely tuned to it. Report the plateau, not the argmax.
</details>

**E6 (★★★).** Zheng et al. (2022) reported that PN→KC connectivity in the hemibrain is not uniformly random: some glomerular pairs converge onto shared KCs more often than chance. Design an analysis that decides whether this structure *improves* or *degrades* the LSH property, and state what result would falsify the random-projection reading.

<details markdown="1"><summary>Solution</summary>

Formalise. Let $\Pi \in [0,1]^{N\times N}$ be the measured co-sampling matrix, $\Pi_{ik} = \Pr[\text{a KC touches both } i \text{ and } k]$, versus the null $\Pi^0_{ik} = p_ip_k$ with $p_i$ the measured marginal glomerular sampling rates (this is Caron et al.'s null, and it already accounts for unequal $p_i$). Define $\Delta = \Pi - \Pi^0$.

Redo §2.4 with correlated rows. The KC-drive covariance becomes
$$\mathrm{Cov}_j(m(x),m(y)) = \sum_i p_i(1-p_i)x_iy_i + \sum_{i\ne k}\Delta_{ik}\,x_iy_k,$$
so the preserved similarity is no longer the plain cosine but $\rho_\Delta = \frac{x^\top (D + \Delta) y}{\sqrt{x^\top(D+\Delta)x}\sqrt{y^\top(D+\Delta)y}}$ with $D = \mathrm{diag}(p_i(1-p_i))$. The hash is still an LSH — but for the **Mahalanobis metric induced by $D+\Delta$**, not the Euclidean one.

The analysis. Take a large panel of odour-evoked PN responses (Drosophila hallem-style ORN data pushed through the C8 normalisation gives you $x$). Compute (a) the behavioural or ecological similarity structure you believe the animal cares about, or failing that the natural-odour covariance $C_{\text{nat}}$; (b) test whether $D+\Delta$ is closer to $C_{\text{nat}}^{-1}$ than $D$ alone is. If $\Delta$ *whitens* — i.e. $\rho_\Delta$ under natural statistics is closer to the ideal-observer discriminability than $\rho_0$ — then the structure is an improvement and the correct level-2 statement upgrades from "random projection" to "projection matched to natural odour statistics."

Falsification of the random reading: measure tag overlap directly. Optogenetically or thermogenetically activate defined glomerular combinations, image KC populations (calcium or, better, voltage), and compute $\mathcal{O}$ as a function of input Pearson correlation. Random projection predicts a single curve $\Phi_2(t,t;\rho)/\alpha$ collapsing all glomerular combinations. Structured sampling predicts systematic, *reproducible-across-animals* residuals organised by glomerular identity. Since the wiring is not stereotyped across individuals in the fly, cross-animal reproducibility of the residual is the decisive test: it would show the structure is genetically specified and therefore functional, rather than a developmental byproduct.
</details>

**E7 (★★★).** The mushroom body and the lateral horn read the same PNs. Formalise the shattering–abstraction tradeoff (Bernardi et al. 2020) and show that a two-readout architecture strictly dominates any single readout on a task set containing both innate and learned categories.

<details markdown="1"><summary>Solution</summary>

Setup. Let odours be $x \in \mathbb{R}^N$ with a latent innate category $g(x) \in \{\pm 1\}$ that is a *smooth, low-dimensional* function (e.g. $g = \mathrm{sign}(v^\top x)$ for a genetically specified $v$). Let learned tasks be arbitrary dichotomies $\ell$ of a training set of $P$ odours, drawn uniformly from the $2^P$ possibilities.

Single readout from a representation $\phi$. Generalisation error on $g$ for a linear readout trained on $P$ samples scales, in the kernel view, with the alignment between the kernel $K_\phi$ and the target: $\varepsilon_g \sim 1 - \langle g, K_\phi g\rangle/(\|g\|\|K_\phi\|)$. A shattering representation has $K_\phi \approx I$ (all patterns near-orthogonal), so $\langle g, K_\phi g\rangle$ carries no information about the smoothness of $g$ and generalisation to *new* odours is at chance — memorisation without abstraction. Conversely a smooth, low-rank representation ($\phi(x)=v^\top x$) generalises perfectly on $g$ but shatters only $O(P)$ of $2^P$ dichotomies (Cover with $M=1$).

Formally: shattering dimensionality $S(\phi) = \frac{1}{2^P}|\{\ell : \ell \text{ linearly decodable from } \phi\}|$ and abstraction $A(\phi) = 1-\varepsilon_g$. Cover gives $S \le C(P,M_{\text{eff}})/2^P$ with $M_{\text{eff}}$ the effective rank of $K_\phi$, so $S$ is increasing in effective rank; kernel-alignment bounds give $A$ decreasing in effective rank once the rank exceeds that of $g$. Hence $S$ and $A$ trade off along the rank axis for any *single* representation: no $\phi$ simultaneously attains $S \approx 1$ and $A \approx 1$.

Two readouts. Let $\phi_{\text{MB}}$ be the KC tag (rank $\approx \alpha M$, $S\approx 1$, $A\approx 0$) and $\phi_{\text{LH}} = Vx$ a stereotyped low-rank projection with $V$ containing $v$ ($S$ small, $A \approx 1$). The concatenated readout $[\phi_{\text{LH}}, \phi_{\text{MB}}]$ has $S \ge S(\phi_{\text{MB}})$ and $A \ge A(\phi_{\text{LH}})$ because a linear readout can zero out either block. Therefore the pair strictly dominates. The cost is that the two blocks must be *architecturally* separated — you cannot get this by tuning one population's sparseness, because the tradeoff is along a single scalar (rank) for a single population.

Prediction: lateral horn representations should be low-dimensional, stereotyped across animals, and support cross-odour generalisation of valence; KC representations should be high-dimensional, individual-specific, and support arbitrary but non-generalising associations. Both halves have empirical support, and the second half — that MBON learning generalises poorly to novel odours in exactly the way a random tag predicts — is the more testable and less tested.
</details>

## 5. Reading path

- **Marr (1969)**, *A theory of cerebellar cortex* (J. Physiol.) — read it for: the original expansion-recoding argument, and for how much of modern theory was already there in 1969.
- **Albus (1971)**, *A theory of cerebellar function* (Math. Biosci.) — read it for: the explicit perceptron formulation and the case for sparse granule codes.
- **Cover (1965)**, *Geometrical and statistical properties of systems of linear inequalities with applications in pattern recognition* — read it for: the function-counting theorem, which is the load-bearing mathematics under every "expansion helps" claim.
- **Perez-Orive, Mazor, Turner, Cassenaer, Wilson & Laurent (2002)**, *Oscillations and sparsening of odor representations in the mushroom body* (Science) — read it for: the primary measurement of KC sparseness and its oscillatory gating; also the bridge to C2.
- **Caron, Ruta, Abbott & Axel (2013)**, *Random convergence of olfactory inputs in the Drosophila mushroom body* (Nature) — read it for: the anatomical claim of randomness and, importantly, exactly what null model it tests.
- **Papadopoulou, Cassenaer, Nowotny & Laurent (2011)**, *Normalization for sparse encoding of odors by a wide-field interneuron* (Science) — read it for: the GGN loop as the physical winner-take-all, and the causal test.
- **Lin, Bygrave, de Calignon, Lee & Miesenböck (2014)**, *Sparse, decorrelated odor coding in the mushroom body enhances learned odor discrimination* (Nat. Neurosci.) — read it for: APL, and for the behavioural consequence of losing sparseness.
- **Dasgupta, Stevens & Navlakha (2017)**, *A neural algorithm for a fundamental computing problem* (Science) — read it for: the LSH conversion itself; then reread §2.4 above and check their similarity-preservation claim against the exact orthant probability.
- **Babadi & Sompolinsky (2014)**, *Sparseness and expansion in sensory representations* (Neuron) — read it for: the noise/decorrelation tradeoff, and for the discipline of asking what expansion costs.
- **Litwin-Kumar, Harris, Axel, Sompolinsky & Abbott (2017)**, *Optimal degrees of synaptic connectivity* (Neuron) — read it for: the derivation that predicts the number of claws. This is the standard.
- **Rigotti, Barak, Warden, Wang, Daw, Miller & Fusi (2013)**, *The importance of mixed selectivity in complex cognitive tasks* (Nature) — read it for: shattering dimensionality as an experimentally measurable quantity.
- **Bernardi, Benna, Rigotti, Munuera, Fusi & Salzman (2020)**, *The geometry of abstraction in the hippocampus and prefrontal cortex* (Cell) — read it for: the counterweight — why you should not want maximal dimensionality.
- **Rahimi & Recht (2007)**, *Random features for large-scale kernel machines* (NIPS) — read it for: the machine-learning statement of §2.3, and for how to think about $M$ as a Monte-Carlo budget.
- **Zheng et al. (2022)**, *Structured sampling of olfactory input by the fly mushroom body* (Curr. Biol.) — read it for: the connectomic evidence that "random" is an approximation, and for E6.

## 6. Open problems and what would settle them

**Is the projection random, matched, or learned?** The Caron null is "independent sampling with unequal glomerular rates," and the deviations from it (Zheng et al.) are small but real. Three hypotheses are live: developmental noise, a genetically specified bias matched to natural odour statistics, and activity-dependent refinement. *Settling it:* rear flies in odour-controlled environments and re-measure the co-sampling matrix $\Delta$; a genetic bias is environment-invariant, refinement is not. In parallel, test whether $\Delta$ predicts anything — see E6 — because a structure that improves nothing is developmental noise regardless of its statistical significance.

**Where on the sparseness–noise curve does the animal sit, and is it adaptive?** §2.7 and E5 predict an interior optimum $\alpha^\star(\epsilon)$. GGN/APL gain is under neuromodulatory control, so $\alpha$ is a state variable, not a constant. *Settling it:* measure PN trial-to-trial noise $\epsilon$ and KC sparseness $\alpha$ in the same preparation across conditions that change $\epsilon$ (odour concentration, background, arousal, satiety). The theory predicts $\alpha$ should *increase* (code should densify) when $\epsilon$ increases. Nobody has done this cleanly, and it is a completely tractable experiment.

**Does the KC code preserve the right metric?** §2.4 proves the tag preserves Pearson correlation in PN space. Whether that is the metric of perceptual similarity is an entirely separate question, and there is no reason to assume it. *Settling it:* generate odour pairs matched for PN-space correlation but differing in behavioural generalisation, or vice versa. If behaviour tracks PN correlation, the conversion is complete; if it does not, the interesting computation is somewhere else (upstream in receptor tuning, or downstream in MBON compartments). This is where `structures/S-06-hyperbolic-odor-space.md` bites: if the natural geometry of odour space is hyperbolic, a hash that preserves Euclidean cosine is systematically distorting exactly the hierarchical relations that matter.

**Does the locust do the same thing?** Almost everything quantitative above is *Drosophila*. The locust has 20× more KCs, ~10–20 claws each, and — critically — an oscillatory gating regime (C2) that the fly's mushroom body does not obviously share. A KC that integrates over one 50 ms oscillation cycle is not computing $\mathrm{WTA}(Ax)$ on a static $x$; it is computing $\mathrm{WTA}$ on a *synchronous subset* of PNs, which is a different projection matrix every cycle. *Settling it:* measure KC tag overlap as a function of PN-population correlation cycle by cycle in locust, and test whether the effective $A$ is cycle-invariant. If it is not, C1's static LSH is a special case of C2's dynamic one, and the sequence of tags carries more than any single tag — which is the thesis of the next unit.
