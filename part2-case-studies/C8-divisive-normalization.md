---
title: "C8 · Divisive normalization"
parent: Case studies
nav_order: 8
---

# C8 — Divisive normalization → a canonical computation with several algorithmic readings

> **Circuit.** A driven neuron whose response is divided by a pooled measure of the activity of its neighbours: V1 simple cells suppressed by a broadly tuned normalisation pool (Heeger), *Drosophila* antennal lobe PNs suppressed by lateral presynaptic inhibition proportional to total ORN activity (Olsen, Bhandawat & Wilson), locust Kenyon cells suppressed by the GGN in proportion to total KC activity (Papadopoulou et al.).
> **Primitive extracted.** $r_i = g\,x_i^n/(\sigma^n + \sum_j w_{ij}x_j^n)$ — and, depending on which parameter you look at, this single expression *is* automatic gain control, contrast-invariant tuning, redundancy reduction, softmax, winner-take-all, attention, and approximate marginalisation.
> **Status.** The equation is extremely well established empirically across systems and species. Every *algorithmic reading* of it is contested, and this unit argues that the contest is largely a category error — which makes it the best available case study in level-2 non-uniqueness.
> **Structures thread.** `structures/S-03-combinatorial-threshold-linear-networks.md` — a threshold-linear network with uniform global inhibition is the subtractive cousin of this motif, and its fixed-point combinatorics tell you exactly when normalisation collapses to winner-take-all. `structures/S-04-manifold-capacity.md` — normalisation reshapes object manifolds by fixing their radius while leaving their dimension alone, which is exactly the quantity that theory says controls capacity. See `structures/README.md`.

## 1. The phenomenon

Heeger (1992) proposed that V1 simple cells are not linear filters followed by a static nonlinearity, but linear filters whose output is *divided* by the pooled activity of a large population of nearby cells. The model accounts for a list of otherwise awkward facts: response saturation with contrast; cross-orientation suppression by a mask that on its own evokes no response; surround suppression; and — the crown jewel — **contrast-invariant orientation tuning**, the fact that a V1 cell's tuning width barely changes as contrast varies over two log units even though its response amplitude changes by an order of magnitude. Carandini, Heeger & Movshon (1997) put it on a firm quantitative footing, and Carandini & Heeger (2012) surveyed the same computation in retina, V1, MT, IT, olfactory systems, auditory cortex, parietal value coding, and multisensory integration, and gave it the name "canonical."

The olfactory instance is the one to lead with, because the circuit is fully identified. Olsen, Bhandawat & Wilson (2010) measured ORN→PN transfer in *Drosophila* and found that a PN's response to its own glomerulus's ORN input is suppressed in proportion to the *total* ORN activity across all glomeruli. The transformation is well fit by

$$
\mathrm{PN}_i \;=\; R_{\max}\,\frac{\mathrm{ORN}_i^{\,n}}{\mathrm{ORN}_i^{\,n} + \big(\sigma + m\sum_j \mathrm{ORN}_j\big)^{n}},\qquad n \approx 1.5 ,
$$

with the lateral term mediated by GABA_B receptors on ORN presynaptic terminals (Olsen & Wilson, 2008; Root et al., 2008). Bhandawat et al. (2007) had already shown that AL processing increases the reliability and separability of ensemble odour representations. And one synapse later, the mushroom body does it again: the GGN/APL loop divides KC drive by total KC activity (Papadopoulou et al., 2011; Lin et al., 2014).

So the same motif appears twice in three synapses of the insect olfactory pathway, at the input to the expansion (C1) and at its output. That is not an accident, and §2.7 says why.

## 2. The conversion

### 2.1 The equation and its parameters

$$
\boxed{\ r_i \;=\; g\,\frac{x_i^{\,n}}{\sigma^n + \sum_j w_{ij}x_j^{\,n}}\ }
$$

Four knobs, and essentially every algorithmic reading below corresponds to taking a limit in one of them:

| knob | limit | what you get |
|---|---|---|
| $\sigma$ | $\sigma \to \infty$ | linear (or power-law) response, no normalisation |
| $\sigma$ | $\sigma \to 0$ | pure ratio; scale invariance |
| $n$ | $n \to \infty$ | winner-take-all |
| $n$ | finite, on log inputs | softmax at inverse temperature $n$ |
| $w_{ij}$ | $w_{ij}=1$ (uniform) | pure gain control; tuning unchanged |
| $w_{ij}$ | structured | whitening / redundancy reduction |

Keep this table in view. Almost every dispute in the literature is two people taking different limits and each believing they have identified "the" computation.

### 2.2 Reading 1 — gain control and dynamic range

Write the input as a contrast times a pattern, $x = c\,\hat{x}$ with $\|\hat x\|$ fixed. Then

$$
r_i = g\,\frac{c^n \hat x_i^{\,n}}{\sigma^n + c^n\sum_j w_{ij}\hat x_j^{\,n}}
= g\,\hat x_i^{\,n}\,\frac{c^n}{\sigma^n + c^n S_i},\qquad S_i = \sum_j w_{ij}\hat x_j^{\,n},
$$

which is the Naka–Rushton function of contrast with half-maximum at

$$
c_{50}^{(i)} = \sigma\, S_i^{-1/n} .
$$

Two consequences. First, the response saturates: a neuron with a limited output range covers an unbounded input range. Second, and more interesting, $c_{50}$ *depends on the stimulus* through $S_i$ — the operating point slides to wherever the input distribution currently is. That is adaptive gain control with no explicit adaptation mechanism and no time constant; it is instantaneous, a property of the divisive form itself.

The olfactory numbers make the point brutally. In §3, total ORN drive across four log units of concentration spans a factor of 57; total PN drive spans a factor of 6.6. An eightfold compression of dynamic range, obtained by a single presynaptic inhibitory pathway.

### 2.3 Reading 2 — contrast invariance of tuning

Take $w_{ij} = w$ uniform, so $S_i = S$ is the same for all $i$. Then in the high-contrast limit $c^n S \gg \sigma^n$,

$$
r_i \;\longrightarrow\; g\,\frac{\hat x_i^{\,n}}{S} ,
$$

**independent of $c$.** The response *pattern* — hence the tuning curve shape, hence the decoded stimulus identity — is invariant to contrast, while the amplitude saturates. This is Heeger's original headline result for orientation tuning.

Notice what it costs: contrast invariance requires $w_{ij}$ to be *uniform*, i.e. the normalisation pool must be untuned. A tuned normalisation pool ($w_{ij}$ depending on the similarity of $i$ and $j$) *breaks* contrast invariance and instead sharpens tuning. These are different computations and they need different $w$. Any claim that a circuit "does normalisation" without specifying $w_{ij}$ has not said what the circuit computes.

An important non-consequence: if $w_{ij}$ is exactly uniform, the divisor is the same scalar for every $i$, so the operation is a pure gain change and the *cosine similarity* between any two population patterns is exactly unchanged. Uniform normalisation cannot decorrelate anything. Everything in Reading 3 depends on $w_{ij}$ being structured or on the self-term $x_i^n$ appearing in the denominator.

### 2.4 Reading 3 — redundancy reduction and whitening

This is the normative reading, and it is the one that makes falsifiable predictions about $w_{ij}$.

Natural signals are well described by a **Gaussian scale mixture**: a filter-response vector $x$ is generated as

$$
x = \sqrt{z}\;u,\qquad u\sim\mathcal N(0,C),\quad z>0\ \text{a scalar "local contrast" variable},
$$

with $z$ independent of $u$. This model captures the dominant statistical dependency in natural images and sounds: filter responses are *uncorrelated* but strongly dependent, because they share $z$ — the variance-dependency structure, visible as the bowtie-shaped conditional histograms Schwartz & Simoncelli (2001) plotted.

Divisive normalisation is the exact inverse of this generative model. Take $n=2$, $\sigma=0$, and $w_{ij}$ chosen so that $\sum_j w_{ij}x_j^2$ estimates $z$ times a constant. Then

$$
r_i \;=\; \frac{x_i^2}{\sum_j w_{ij}x_j^2}
\;=\; \frac{z\,u_i^2}{z\sum_j w_{ij}u_j^2}
\;=\; \frac{u_i^2}{\sum_j w_{ij}u_j^2},
$$

**and $z$ has cancelled exactly.** The normalised response is independent of the shared variance variable, which is precisely what "removing the dependency" means. Schwartz & Simoncelli's contribution was to *fit* $w_{ij}$ from natural image statistics and show that the resulting weights predict the measured pattern of surround suppression in V1 and of two-tone suppression in auditory nerve. That is the discipline this unit is arguing for: a normative reading earns its keep when it predicts circuit parameters it was not fit to.

Note the relationship to whitening. Linear whitening removes *second-order* correlations, $\langle x_ix_j\rangle \to \delta_{ij}$. That is what centre-surround receptive fields do (Atick & Redlich). Divisive normalisation removes the *residual higher-order* dependency that linear whitening leaves behind. They are complementary stages, and in the retina you can see both — the linear one in the receptive field, the divisive one in contrast gain control.

The engineering vindication is worth knowing: Ballé, Laparra & Simoncelli's generalised divisive normalisation, used as a learned invertible transform, remains competitive as a Gaussianising layer in learned image compression, and batch/layer/group normalisation in deep networks are the same motif with the statistics estimated over different axes.

### 2.5 Reading 4 — softmax, max, and winner-take-all

Set $w_{ij}=1$ and $\sigma \to 0$:

$$
r_i \;=\; \frac{x_i^{\,n}}{\sum_j x_j^{\,n}} .
$$

Substitute $x_i = e^{u_i}$:

$$
r_i \;=\; \frac{e^{n u_i}}{\sum_j e^{n u_j}} \;=\; \mathrm{softmax}(n\mathbf u)_i .
$$

**Divisive normalisation with exponent $n$ acting on exponentiated inputs is exactly softmax at inverse temperature $n$.** This is not an analogy; it is an identity. The attention operation in a transformer is a divisive normalisation, and $n$ is the temperature.

Take $n\to\infty$ with a unique maximum $x_{i^\star}$:

$$
r_i = \frac{1}{\sum_j (x_j/x_i)^n} \;\longrightarrow\; \mathbb 1[i = i^\star],
$$

**winner-take-all.** So the single expression interpolates continuously between "pass everything through, normalised" ($n\to0$: $r_i\to1/N$, complete information destruction), "divide by the mean" ($n=1$), "soft competition" ($n\sim2$), and "select the maximum" ($n\to\infty$). Section 3 measures this by the entropy of the normalised pattern, which falls from 4.4 bits at $n=0.5$ to 0.4 bits at $n=40$ for a 24-channel population.

Biology sits at $n \approx 1.5$–$2$ almost everywhere it has been measured. That is a *soft* competition: enough to compress dynamic range and equalise, not enough to discard the pattern. If you believe the WTA reading of C1, note that the KC layer needs $n$ effectively large — and it gets there not by raising the exponent but by adding a *threshold* on top of the normalisation, which is a cheaper way to achieve the same selectivity.

### 2.6 Reading 5 — attention, and Reading 6 — marginalisation

**Attention (Reynolds & Heeger, 2009).** Model attention as an *attention field* $A_i$ that multiplies the stimulus drive before normalisation:

$$
r_i = \frac{A_i E_i}{\sigma + \sum_j w_{ij}A_jE_j}.
$$

Two limits, and they explain a decade of apparently contradictory experiments:
- **Narrow attention field** (attention smaller than the stimulus): $A$ boosts the numerator for the attended neuron but barely changes the pooled denominator, giving multiplicative **response gain** — the tuning curve scales up.
- **Broad attention field** (attention larger than the stimulus): $A$ scales numerator and denominator together for a high-contrast stimulus, so the effect appears only where the denominator is dominated by $\sigma$, i.e. at low contrast — giving **contrast gain**, a leftward shift of the contrast-response function.

The model therefore predicts *which* form of gain modulation you will observe from the relative sizes of the attention field and the stimulus, which is a genuinely risky prediction and has held up.

**Marginalisation (Beck, Latham & Pouget, 2011).** Suppose a population has Poisson-like variability, so the log-likelihood over a stimulus $s$ is linear in the spike counts: $\log p(r|s) = h(s)^\top r + \text{const}$. Then multiplying likelihoods (cue combination) is *addition* of activity — trivially neural. But **marginalising** a nuisance variable, $p(s|r) = \int p(s,\nu|r)d\nu$, is not linear, and generically requires a quadratic nonlinearity followed by a division. Beck et al. showed that quadratic pooling plus divisive normalisation performs near-optimal marginalisation for a broad class of problems.

For an olfactory reader this is the most interesting reading in the unit, because concentration is exactly a nuisance variable for identity. If the AL's normalisation is marginalising concentration, then the PN population should support a concentration-invariant identity readout that is *close to optimal*, not merely improved. §3 measures the improvement (cross-concentration decoding rises from 0.42 to 0.83 in the model); whether it approaches the ideal-observer bound is a question nobody has answered for a real olfactory system, and it is answerable.

### 2.7 The olfactory hook: why normalisation is the precondition for C1

Here is the argument I most want you to take from this unit.

C1's locality-sensitive hash requires a threshold that selects a fixed fraction of Kenyon cells. A *fixed absolute* threshold cannot do this if total PN drive varies with concentration: at low concentration nothing crosses threshold, at high concentration everything does, and the "hash" degenerates to a concentration meter. §3 measures exactly this: with unnormalised ORN drive into a fixed-threshold KC layer, the active fraction goes from $0.000$ at the lowest concentration to $1.000$ at the highest. The code is destroyed.

With AL normalisation in place, total PN drive is compressed eightfold and the active KC fraction is far better behaved across the middle three log units. And the second normalisation stage — GGN/APL — finishes the job, because subtractive/divisive global feedback makes the *active set* scale-invariant by construction (E4).

So the correct architectural statement is: **the insect olfactory pathway normalises twice, and it has to, because the random-projection-plus-threshold code of C1 is only well-defined on a normalised input.** Normalisation is not a separate computation bolted on; it is the precondition that makes the expansion code meaningful. This is a general moral. Thresholds are only meaningful relative to a scale, and any circuit that uses a threshold must contain, somewhere upstream, the machinery that sets that scale.

### 2.8 The methodological point: what to do with a computation that has six readings

Six readings of one equation. Are five of them wrong?

No — and the question is malformed. Sort them into two categories.

**Descriptive readings** — gain control (2.2), contrast invariance (2.3), softmax/WTA (2.5) — are *restatements of the input–output function in different regimes*. They are not competing hypotheses; they are true by algebra once you accept the equation and specify the parameters. Arguing about whether normalisation "is really" gain control or "is really" softmax is like arguing whether $y=1/x$ is really a hyperbola or really an inverse: it depends on which coordinate you are varying.

**Normative readings** — redundancy reduction (2.4), marginalisation (2.6), and the attention account insofar as it claims an objective — are genuinely competing, because each claims the *parameters* $(n, \sigma, w_{ij})$ take particular values *for a reason*, and those reasons make different predictions:

- Redundancy reduction predicts $w_{ij}$ from the covariance structure of natural inputs, and predicts that $w_{ij}$ should *change* if the input statistics change (over development, or adaptively).
- Marginalisation predicts $w_{ij}$ from the structure of the nuisance variable — for concentration invariance in olfaction, $w_{ij}$ should be uniform (since concentration scales all channels), whereas for background invariance it should be tuned to the background's glomerular profile.
- The attention account predicts $w$ from a resource-allocation objective and says nothing about input statistics.

These make *different, incompatible* predictions about the same measurable matrix. That is what a real dispute looks like.

**The discipline.** A level-2 reading is idle unless it predicts a level-1 parameter it was not fitted to. Schwartz & Simoncelli did this (natural image statistics → surround suppression weights). Litwin-Kumar et al. did it in C1 (dimensionality objective → number of claws). Harper & McAlpine did it in C3 (Fisher information → best-delay distribution). When you evaluate an algorithmic reading, the question is never "is this consistent with the data" — with six readings available, everything is consistent with the data. The question is: **what number does it predict, and was that number measured afterwards?**

## 3. Worked example / model to build

An ORN→PN→KC olfactory model with Olsen–Bhandawat–Wilson normalisation, measuring dynamic-range compression, cross-concentration decoding, the exponent sweep from softmax to WTA, and the effect on a downstream fixed-threshold KC layer.

```python
import numpy as np
rng = np.random.default_rng(1)
NG, NO = 24, 12                                    # glomeruli, odours
K  = np.exp(rng.normal(-1., 2.0, (NO, NG)))        # affinity of odour a for ORN i
C  = np.logspace(-2, 2, 5)                         # 4 decades of concentration

ORN = lambda a, c, R=200.: R*(c*K[a])/(1. + c*K[a])          # saturating binding
PN  = lambda o, n=1.5, s0=12., m=0.05, R=165.: \
          R*o**n/(o**n + (s0 + m*o.sum())**n)                # Olsen-Bhandawat-Wilson

O = np.array([[ORN(a, c) for c in C] for a in range(NO)])    # (NO, nC, NG)
P = np.array([[PN(O[a, k]) for k in range(len(C))] for a in range(NO)])

print(" conc      sum(ORN)   sum(PN)   max(PN)")
for j, c in enumerate(C):
    print("%7.2f   %8.0f  %8.0f  %8.0f" % (c, O[:,j].sum(1).mean(),
                                           P[:,j].sum(1).mean(), P[:,j].max(1).mean()))

r2 = np.random.default_rng(5); A = np.zeros((2000, NG))      # KC layer (see C1)
for j in range(2000): A[j, r2.choice(NG, 6, replace=False)] = 1.
thO = np.quantile(O[:,2] @ A.T, 0.95); thP = np.quantile(P[:,2] @ A.T, 0.95)
print("\n conc    KC active frac from ORN | from PN")
for j, c in enumerate(C):
    print("%7.2f      %.3f              %.3f"
          % (c, ((O[:,j]@A.T) > thO).mean(), ((P[:,j]@A.T) > thP).mean()))

def acc(V, tr, te):                                          # cross-concentration ID
    Cn = V[:,tr]/np.linalg.norm(V[:,tr], axis=1, keepdims=True)
    X  = V[:,te]/np.linalg.norm(V[:,te], axis=1, keepdims=True)
    return float((np.argmax(X @ Cn.T, 1) == np.arange(NO)).mean())
print()
for tr, te in [(1,3), (2,4), (0,4)]:
    print("train c=%.2f test c=%.2f :  ORN %.2f   PN %.2f"
          % (C[tr], C[te], acc(O,tr,te), acc(P,tr,te)))

print("\n  n     entropy of normalised pattern (bits, max %.2f)" % np.log2(NG))
for n in [0.5, 1., 1.5, 3., 8., 40.]:                        # softmax -> WTA
    Q = np.array([PN(O[a,2], n=n, s0=1e-6, m=1.0) for a in range(NO)])
    Q = Q/Q.sum(1, keepdims=True)
    print("%5.1f          %.2f" % (n, -np.mean((Q*np.log2(Q+1e-300)).sum(1))))
```

**What to look for.**

*Dynamic-range compression (§2.2).* Total ORN drive rises $77 \to 4363$ over four decades; total PN drive rises only $245 \to 1621$. Peak PN response actually *falls* at the highest concentrations, from 137 to 74 — normalisation is trading peak response for population-level stability.

*Cross-concentration identity (§2.6).* Decoding odour identity with centroids from one concentration and test data from another: for a two-decade gap, ORN $0.42 \to$ PN $0.83$; for a different two-decade gap at higher concentration, ORN $0.17 \to$ PN $0.33$. Both large improvements. But for the full four-decade gap the normalisation *hurts*: ORN $0.17 \to$ PN $0.08$. Do not skip past this. The marginalisation is approximate and it fails at the extremes, for the same reason the KC table below fails there — at the top of the range the ORN pattern itself has saturated into near-uniformity, and no amount of division recovers a pattern that the receptors have already destroyed. Normalisation marginalises concentration over the range where the *input* still carries identity, and not one log unit further. This is exactly the behaviour Beck et al. would predict for a circuit doing quadratic-plus-divisive pooling rather than exact integration, and it is a useful antidote to reading "invariance" as an absolute.

*The precondition argument (§2.7).* This is the important table. A fixed-threshold KC layer fed by raw ORN activity gives active fractions $0.000, 0.000, 0.050, 0.955, 1.000$ across four decades — catastrophic. Fed by normalised PN activity: $0.000, 0.002, 0.050, 0.074, 0.000$. Far flatter over the middle range, and still failing at the extremes, which is precisely why a *second* normalisation stage exists (E4).

*Softmax to WTA (§2.5).* Entropy of the normalised pattern falls monotonically $4.43 \to 4.09 \to 3.74 \to 2.94 \to 1.63 \to 0.40$ bits for $n = 0.5, 1, 1.5, 3, 8, 40$. At the biological $n=1.5$ the population retains 3.7 of a possible 4.6 bits — soft competition, not selection.

## 4. Exercises

**E1 (★).** Show that with uniform $w_{ij}=w$ and $\sigma=0$, divisive normalisation leaves the cosine similarity between any two population patterns exactly unchanged, and identify the two ingredients that let it decorrelate anyway.

<details markdown="1"><summary>Solution</summary>

With $w_{ij}=w$ and $\sigma=0$, the divisor $D(x) = w\sum_j x_j^n$ does not depend on $i$, so $r(x) = g\,x^{\circ n}/D(x)$, where $x^{\circ n}$ is the elementwise power. For two inputs $x,y$,

$$\cos(r(x), r(y)) = \frac{r(x)^\top r(y)}{\|r(x)\|\|r(y)\|}
= \frac{[g/D(x)][g/D(y)]\,(x^{\circ n})^\top y^{\circ n}}{[g/D(x)]\|x^{\circ n}\|\,[g/D(y)]\|y^{\circ n}\|}
= \cos(x^{\circ n}, y^{\circ n}).$$

The scalar divisors cancel identically. So uniform divisive normalisation changes the *pattern geometry* only through the pointwise power $n$ — and if $n=1$ it changes nothing at all.

Two ingredients recover decorrelation:
1. **A self-term in the denominator.** The Olsen form has $x_i^n$ in the denominator as well: $r_i = x_i^n/(x_i^n + D)$. This is a *compressive, saturating* function applied per channel, which shrinks strong channels relative to weak ones and hence flattens the pattern — a nonlinear reshaping that does change cosine similarity. In §3 this shows up as $\max(\mathrm{PN})$ falling while $\sum(\mathrm{PN})$ rises.
2. **Structured $w_{ij}$.** If $w_{ij}$ is larger for channels that are correlated in natural stimuli, then each channel is divided by a quantity that predicts it, which is the divisive analogue of predictive coding and is exactly Reading 3.

The moral: "normalisation decorrelates" is *false* for the textbook uniform-pool model and true only for the two variants above. Be specific about which one you mean.
</details>

**E2 (★★).** Derive the contrast-response function and show that divisive normalisation makes the *Fisher information about stimulus identity* approximately contrast-invariant at high contrast, while the information about contrast itself saturates.

<details markdown="1"><summary>Solution</summary>

Let $x = c\hat x$, uniform pool, $r_i = g c^n\hat x_i^n/(\sigma^n + c^nS)$ with $S = \sum_j\hat x_j^n$. Write $\theta$ for a stimulus-identity parameter that changes $\hat x$ but not its norm (e.g. orientation, or odour identity at fixed total drive), so $\partial_\theta S = 0$.

$$\frac{\partial r_i}{\partial\theta} = \frac{gc^n}{\sigma^n + c^nS}\,\frac{\partial \hat x_i^n}{\partial\theta} = \underbrace{\frac{gc^n}{\sigma^n+c^nS}}_{\textstyle \kappa(c)}\ \partial_\theta \hat x_i^n .$$

With Poisson-like variability $\sigma_i^2 = \phi\, r_i$, Fisher information about $\theta$ is

$$I_\theta = \sum_i\frac{(\partial_\theta r_i)^2}{\phi r_i}
= \frac{\kappa(c)^2}{\phi}\sum_i\frac{(\partial_\theta\hat x_i^n)^2}{\kappa(c)\hat x_i^n}
= \frac{\kappa(c)}{\phi}\sum_i\frac{(\partial_\theta\hat x_i^n)^2}{\hat x_i^n}.$$

So $I_\theta \propto \kappa(c) = gc^n/(\sigma^n + c^nS)$, which **saturates at $g/S$** as $c\to\infty$. Above $c_{50}$, identity information is essentially contrast-independent: the code becomes contrast-invariant not only in its mean but in its precision.

Now contrast itself. $\partial_c r_i = g\hat x_i^n\, n c^{n-1}\sigma^n/(\sigma^n+c^nS)^2$, so

$$I_c = \sum_i \frac{(\partial_c r_i)^2}{\phi r_i}
= \frac{n^2\sigma^{2n}c^{n-2}}{\phi\,(\sigma^n+c^nS)^{3}}\,g\sum_i \hat x_i^{\,n} \;\propto\; \frac{c^{n-2}\sigma^{2n}}{(\sigma^n+c^nS)^3},$$

which for $c \gg c_{50}$ falls as $c^{-2n-2}$ — **contrast information is destroyed at high contrast, fast.**

The system has therefore made a choice: it spends its dynamic range on identity and throws away intensity above the saturation point. In olfaction that is exactly the marginalisation reading — concentration is the nuisance variable, and the circuit is *deliberately* discarding it. It also predicts that concentration discrimination should be poor precisely where identity discrimination is best, which is a testable psychophysical prediction and, as far as I know, an open one in insects.
</details>

**E3 (★★).** Show that the Gaussian scale mixture cancellation of §2.4 is exact only for a particular choice of $w_{ij}$, and derive that choice.

<details markdown="1"><summary>Solution</summary>

Model: $x = \sqrt z\,u$, $u \sim \mathcal N(0,C)$, $z \perp u$. We want $r_i = x_i^2/\sum_j w_{ij}x_j^2$ to be independent of $z$.

Substituting, $r_i = u_i^2/\sum_j w_{ij}u_j^2$ for *any* $w$: the $z$ cancels identically because it multiplies numerator and denominator. So $z$-independence is automatic and imposes no constraint. The real objective is stronger: we want the $r_i$ to be **mutually independent**, or at least uncorrelated with the local variance.

Impose $\mathbb{E}[r_i] = $ const and, more usefully, that the divisor be an unbiased estimate of the *conditional* variance of $x_i$:

$$\sum_j w_{ij}x_j^2 \;\approx\; \mathbb{E}[x_i^2 \mid \{x_j\}_{j\ne i}] .$$

For a jointly Gaussian $u$ this conditional expectation is available in closed form. Writing $C$ for the covariance, the linear MMSE predictor of $x_i^2$ from $\{x_j^2\}$ under the GSM has coefficients proportional to the *squared* correlations:

$$w_{ij} \;\propto\; \frac{C_{ij}^2}{C_{jj}} \quad\text{(up to a normalisation making } \textstyle\sum_j w_{ij}C_{jj} = C_{ii}).$$

Derivation sketch: for zero-mean jointly Gaussian $(x_i,x_j)$ with correlation $\rho_{ij}$, $\mathrm{Cov}(x_i^2, x_j^2) = 2C_{ij}^2$ and $\mathrm{Var}(x_j^2) = 2C_{jj}^2$. The optimal linear predictor of $x_i^2$ from the vector $\{x_j^2\}$ is $w = \Sigma^{-1}\gamma$ with $\Sigma_{jk}=2C_{jk}^2$ and $\gamma_j = 2C_{ij}^2$; in the diagonal-$\Sigma$ approximation this gives $w_{ij}\propto C_{ij}^2/C_{jj}^2 \cdot C_{jj}$.

**The prediction that matters:** $w_{ij}$ should scale with the *squared* correlation between channels $i$ and $j$ under natural stimulus statistics — not the correlation itself, and not the tuning-curve similarity. This is measurable. In olfaction it says the strength of lateral inhibition between two glomeruli should track $C_{ij}^2$ computed over natural odour ensembles. Since *Drosophila* lateral inhibition is largely global rather than glomerulus-specific, the honest conclusion is that the fly AL is *not* implementing this reading; it is implementing the uniform-pool version, which is gain control and marginalisation, not redundancy reduction. That is a genuine discrimination between readings using existing anatomy, and it is why the olfactory case is worth more here than the visual one.
</details>

**E4 (★★, computational).** In §3 the normalised PN input still produces a KC active fraction that collapses at the highest concentration. Add a GGN-style second normalisation at the KC layer and show that the sparseness becomes far more concentration-stable. Then derive why.

<details markdown="1"><summary>Solution</summary>

Implementation: replace the fixed threshold with subtractive global feedback, solved to self-consistency,
```python
def ggn(d, beta=0.05, it=500):
    G = 0.
    for _ in range(it):
        y = np.maximum(d - beta*G, 0.); G = 0.8*G + 0.2*y.sum()
    return y
```
and report `(ggn(V @ A.T) > 0).mean()`. Running this on the §3 model gives active fractions of roughly $0.06, 0.08, 0.08, 0.14, 0.16$ across the four decades for ORN drive and $0.07, 0.08, 0.05, 0.08, 0.17$ for PN drive — compared with $0.000 \to 1.000$ for a fixed threshold.

**Why.** At the self-consistent solution, let $k$ units be active with total drive $S_k = \sum_{i \in \text{active}} d_i$. Then
$$G = \sum_{i\in\text{active}}(d_i - \beta G) = S_k - k\beta G \quad\Longrightarrow\quad G = \frac{S_k}{1+k\beta},$$
and the effective cut level is
$$c = \beta G = \frac{\beta S_k}{1 + k\beta}.$$
Now scale all drives by $\lambda$: $S_k \to \lambda S_k$, so $c \to \lambda c$ *with the same $k$* — the self-consistency equation is homogeneous of degree one. Therefore **the active set is exactly invariant to a uniform gain change**, for any $\beta$. Global feedback inhibition is a scale-invariant thresholding device; that is its entire computational content, and it is why C1's top-$k$ idealisation is a good model of the GGN loop rather than a convenient fiction.

The residual drift you do observe comes from the ORN *pattern* changing shape with concentration (receptor saturation flattens it), not from gain. That is exactly the part that AL normalisation cannot fix and that receptor-level adaptation must handle. Three stages, three different scale problems.
</details>

**E5 (★★).** Show that the stabilised supralinear network (SSN) produces divisive normalisation dynamically, and identify what it predicts that the static equation does not.

<details markdown="1"><summary>Solution</summary>

The SSN (Ahmadian, Rubin & Miller, 2013; Rubin, Van Hooser & Miller, 2015) is a recurrent E–I network with a *supralinear* power-law transfer function $r = k[V]_+^{\,p}$, $p>1$. Its key property is that the effective gain $\partial r/\partial V = pk[V]_+^{p-1}$ *increases* with activity.

At low input, gains are small, the effective recurrent coupling is weak, and the network behaves like a feedforward supralinear unit: $r \propto x^p$, supralinear summation. As input grows, gains grow, and the recurrent inhibitory loop becomes strong enough to dominate — the network transitions to an inhibition-stabilised regime where the effective gain is set by the *inhibitory* feedback, and the steady-state response becomes

$$r \;\approx\; \frac{x}{1 + \text{(recurrent inhibitory gain)}} \;\propto\; \frac{x}{\text{pooled activity}},$$

i.e. divisive normalisation, with the divisor emerging from the fixed point rather than being written down. Linearising the E–I system about its fixed point and eliminating $I$ gives $r_E = (I + W_{EI}\,g_I\,W_{IE})^{-1}x$, and since $g_I$ grows with activity, the effective divisor grows with the pooled drive: the Heeger form, derived rather than assumed.

**What SSN predicts that the static equation cannot:**
1. **A supralinear-to-sublinear transition** in summation as contrast increases — additive at low contrast, normalising at high contrast. Observed.
2. **Dynamics.** Normalisation should not be instantaneous: the network takes a few time constants to reach the new fixed point, and the transient should show a characteristic overshoot whose size grows with the strength of the normalisation. The static equation predicts zero latency.
3. **Variability quenching.** In the inhibition-stabilised regime, the network suppresses shared fluctuations, so stimulus onset should *reduce* trial-to-trial variability — one of the most robust and initially puzzling observations in cortex, which the static equation says nothing about.
4. **Surround suppression with a specific spatial dependence** determined by the recurrent connectivity footprint rather than by a freely fitted $w_{ij}$.

This is the template for upgrading a descriptive level-2 equation into a mechanistic one: the mechanism must reproduce the equation *and* predict its failures. Ask the same question of the olfactory case — does the *Drosophila* AL normalisation have a latency, and does it quench PN variability at odour onset? Bhandawat et al. reported increased reliability, which is suggestive. It has not, to my knowledge, been analysed as an SSN.
</details>

**E6 (★★★).** Devise a single experiment in the *Drosophila* antennal lobe that discriminates the redundancy-reduction reading from the marginalisation reading of the same normalisation circuit.

<details markdown="1"><summary>Solution</summary>

The two readings make incompatible predictions about $w_{ij}$, and — crucially — about how it should *change*.

**Predictions.**
- *Redundancy reduction* (§2.4, E3): $w_{ij} \propto C_{ij}^2/C_{jj}$ where $C$ is the ORN covariance over the animal's natural odour ensemble. Therefore $w$ is **structured**, and it must **track the environment**: an animal reared in an odour environment with different second-order statistics should develop different lateral weights.
- *Marginalisation of concentration* (§2.6): concentration multiplies all channels by a common factor, so the sufficient statistic for the nuisance variable is $\sum_j x_j$, and the optimal $w_{ij}$ is **uniform**. It should be **environment-independent**, because the nuisance variable's structure is a physical fact about dilution, not about the odour ensemble.

**The experiment.** Two arms.

*Arm 1 — measure $w_{ij}$ directly.* Use paired glomerular stimulation: optogenetically activate ORN class $j$ alone at graded intensities while recording PN $i$, for all ordered pairs $(i,j)$ across a substantial subset of glomeruli. Fit the Olsen model per pair to extract $w_{ij}$. Then test: is $\hat w$ statistically indistinguishable from uniform (marginalisation), or does it correlate with $C_{ij}^2$ estimated from a natural odour panel (redundancy reduction)? Existing evidence (largely global GABA_B presynaptic inhibition with modest glomerulus-specific structure) leans uniform, but nobody has done the full pairwise matrix with the $C_{ij}^2$ prediction pre-registered.

*Arm 2 — the causal test.* Rear flies from eclosion in two odour environments engineered to have the *same* marginal glomerular activation statistics but *different* covariance structure (achievable by delivering the same set of monomolecular odorants either independently or in fixed correlated blends). Redundancy reduction predicts $w$ diverges between the groups, in the direction predicted by $C_{ij}^2$. Marginalisation predicts $w$ is identical.

**Why this is decisive and why it is worth doing.** The two hypotheses agree on every static input–output measurement — both produce dynamic-range compression, both improve cross-concentration decoding, both are consistent with the fitted $n\approx1.5$. They disagree only on the *matrix* and on its *plasticity*. That is the general shape of a good level-2 discrimination: find the parameter each theory claims to explain, and vary the thing each theory says the parameter depends on.

A prediction I will commit to: the fly AL will come out close to uniform and environment-independent, i.e. **marginalisation, not redundancy reduction** — because concentration variation in a turbulent plume spans four log units within a single sniff, whereas the covariance structure of the odour ensemble varies slowly and can be handled downstream by plasticity in the mushroom body. Fast nuisance variables get handled by fixed circuitry; slow statistical structure gets handled by learning. If that principle is right, it also tells you where to look for the redundancy-reduction machinery: not in the AL, but in KC→MBON plasticity.
</details>

## 5. Reading path

- **Heeger (1992)**, *Normalization of cell responses in cat striate cortex* (Visual Neurosci.) — read it for: the original proposal, and for how much explanatory work one equation does.
- **Carandini, Heeger & Movshon (1997)**, *Linearity and normalization in simple cells of the macaque primary visual cortex* (J. Neurosci.) — read it for: the quantitative test, and for the careful separation of what is linear from what is not.
- **Schwartz & Simoncelli (2001)**, *Natural signal statistics and sensory gain control* (Nat. Neurosci.) — read it for: the normative derivation of §2.4 and, more importantly, for the method of predicting circuit parameters from stimulus statistics. This is the paper to imitate.
- **Olsen & Wilson (2008)**, *Lateral presynaptic inhibition mediates gain control in an olfactory circuit* (Nature) — read it for: the mechanism (GABA_B on ORN terminals) in a circuit you can enumerate.
- **Olsen, Bhandawat & Wilson (2010)**, *Divisive normalization in olfactory population codes* (Neuron) — read it for: the fitted olfactory normalisation model used in §3, and for the demonstration that the pool is the *total* ORN activity.
- **Bhandawat, Olsen, Gouwens, Schlief & Wilson (2007)**, *Sensory processing in the Drosophila antennal lobe increases reliability and separability of ensemble odor representations* (Nat. Neurosci.) — read it for: the population-level consequence, and as the closest thing to a direct test of the marginalisation reading.
- **Papadopoulou, Cassenaer, Nowotny & Laurent (2011)**, *Normalization for sparse encoding of odors by a wide-field interneuron* (Science) — read it for: the second normalisation stage, and for §2.7's argument.
- **Reynolds & Heeger (2009)**, *The normalization model of attention* (Neuron) — read it for: the response-gain versus contrast-gain resolution, which is the model's most impressive risky prediction.
- **Beck, Latham & Pouget (2011)**, *Marginalization in neural circuits with divisive normalization* (J. Neurosci.) — read it for: the inference reading, which is the one an olfaction person should care most about.
- **Carandini & Heeger (2012)**, *Normalization as a canonical neural computation* (Nat. Rev. Neurosci.) — read it for: the breadth of the phenomenon; read it *after* the primary papers, or you will absorb the slogan without the constraints.
- **Ohshiro, Angelaki & DeAngelis (2011)**, *A normalization model of multisensory integration* (Nat. Neurosci.) — read it for: normalisation explaining inverse effectiveness, which no other account handles as cleanly.
- **Louie, Grattan & Glimcher (2011)**, *Reward value-based gain control: divisive normalization in parietal cortex* (J. Neurosci.) — read it for: the same motif in value coding, and for the behavioural consequence (context-dependent choice), which is where normalisation stops being obviously beneficial.
- **Ahmadian, Rubin & Miller (2013)**, *Analysis of the stabilized supralinear network* (Neural Comput.), and **Rubin, Van Hooser & Miller (2015)**, *The stabilized supralinear network: a unifying circuit motif underlying multi-input integration in sensory cortex* (Neuron) — read them for: normalisation as an emergent property of recurrent dynamics rather than a written-down equation. This is the mechanistic upgrade.
- **Wilson (2013)**, *Early olfactory processing in Drosophila: mechanisms and principles* (Annu. Rev. Neurosci.) — read it for: the integrated olfactory account, and for a sober assessment of what the normalisation buys.
- **Ballé, Laparra & Simoncelli (2016)**, *Density modeling of images using a generalized normalization transformation* (ICLR) — read it for: the engineering vindication of Reading 3, and for what a learned $w_{ij}$ looks like.

## 6. Open problems and what would settle them

**Is the normalisation pool uniform or structured, and does it adapt?** E6 lays out the decisive experiment for the fly. The general version applies everywhere the motif appears: nobody has measured a full $w_{ij}$ matrix in any system and tested it against a pre-registered normative prediction. Schwartz & Simoncelli came closest, in V1, with a fitted rather than predicted matrix. *Settling it:* full pairwise measurement plus environmental manipulation, with the prediction stated in advance.

**What is the exponent $n$ for?** It is consistently $1.5$–$2$ and consistently unexplained. Candidate accounts: it is whatever the biophysics gives you (a power-law arising from the threshold nonlinearity of a noisy spiking neuron, which yields $p\approx2$); it is set by a redundancy-reduction objective on natural statistics; it sets the softmax temperature to trade selectivity against information. These make different predictions about whether $n$ varies with input statistics or with noise level. *Settling it:* measure $n$ in the same neurons under manipulations that change trial-to-trial variability without changing mean drive. The biophysical account predicts $n$ falls as noise rises (the effective transfer function linearises); the normative accounts predict it does not.

**Is normalisation instantaneous or dynamic, and does the transient matter?** The static equation has no time in it, but the SSN says normalisation has a latency and a transient. In olfaction, that transient would land squarely on the timescale of C2's decorrelation trajectory. *Settling it:* measure the time course of the normalisation strength itself — deliver a step change in the pooled drive (activate other glomeruli) at various delays after odour onset and measure how quickly the PN response is suppressed. If normalisation takes 100–300 ms to develop, then a substantial part of what C2 calls "progressive decorrelation" may simply be normalisation coming online. That would be a deflationary but important result, and it is directly testable with existing tools.

**Does normalisation help or hurt at the behavioural level?** Every reading in §2 is a story about a benefit. But normalisation makes responses *context-dependent* — Louie et al. showed that value normalisation produces violations of independence-from-irrelevant-alternatives in choice, and the visual analogue is that a stimulus's representation depends on what surrounds it. In olfaction, the prediction is that identifying odour A should be *impaired* by the presence of an unrelated background odour B, in proportion to B's total ORN drive — even though B carries no information about A. The counter-observation is that locusts and moths track odours against backgrounds remarkably well. *Settling it:* measure the psychophysical cost of background as a function of background *total drive* at fixed background *identity overlap*. Normalisation predicts a cost that scales with total drive; a background-subtraction mechanism predicts a cost that scales with overlap. This is the cleanest behavioural discriminator available, and it connects directly back to C2, where background-invariant odour recognition is one of the strongest arguments that the antennal lobe is doing something more than gain control.
