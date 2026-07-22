---
title: S-07 Random matrices & chaos
parent: Structures
nav_order: 7
---

# S-07 — Random matrices, chaos, and the theory you need to build a null

> **The object.** Dynamic mean-field theory for randomly connected networks, the transition to
> chaos at $g=1$, and the spectral theory of *structured* non-Hermitian random matrices —
> including the non-normality that generates transient amplification.
> **Born from.** The need to know what a generic recurrent circuit does, before claiming a
> specific one does something special.
> **Mathematical home.** Random matrix theory; statistical field theory.
> **Situation.** #8, the generic-null problem.
> **Novelty.** New results in an old field, driven by neural questions — structured spectra,
> non-normal transient theory, and the DMFT machinery in the form the field now uses.

---

## 1. Where it comes from

You have a claim: "the antennal lobe's recurrent dynamics decorrelate odour representations."
Before that means anything you must answer: what would a *random* network with the same
statistics do? If a random network decorrelates too, you have discovered a property of recurrent
dynamics in general, not of the antennal lobe.

Building that null requires knowing the generic behaviour of a randomly connected nonlinear
network. In 1988 nobody did. Answering it produced a body of theory that has long outlived the
original motivation, and which you now need as a working tool.

---

## 2. The object, precisely

### The model and the transition

$$\dot{x}_i = -x_i + \sum_{j=1}^N J_{ij}\,\phi(x_j), \qquad J_{ij} \sim \mathcal{N}(0, g^2/N),
\qquad \phi = \tanh.$$

**Linear stability of $x=0$.** The Jacobian is $-I + J\phi'(0) = -I + J$. By **Girko's circular
law**, the eigenvalues of $J$ fill the disk of radius $g$ in the complex plane, uniformly. So the
Jacobian's spectrum is that disk shifted to be centred at $-1$, and the origin loses stability
exactly when the disk's right edge crosses zero:

$$\boxed{g = 1.}$$

Above it, the network is chaotic. This is one of the cleanest results in theoretical
neuroscience: a sharp, parameter-free transition, derived in two lines from a random matrix law.

### Dynamic mean-field theory

Beyond the transition you need the full nonlinear theory. The DMFT argument: in the $N \to
\infty$ limit each neuron receives input $\eta_i(t) = \sum_j J_{ij}\phi(x_j)$ which, by a
central-limit argument over the disorder, is a **Gaussian process** with autocorrelation

$$\langle \eta(t)\eta(t+\tau)\rangle = g^2\,\big\langle \phi(x(t))\phi(x(t+\tau))\big\rangle.$$

Self-consistency closes the loop: the statistics of the input determine the statistics of the
output, which determine the input. Writing $C(\tau) = \langle x(t)x(t+\tau)\rangle$, the
self-consistent equation takes the form of a Newtonian particle in a potential,

$$\frac{d^2 C}{d\tau^2} = -\frac{\partial V(C)}{\partial C}, \qquad
V(C) = -\frac{C^2}{2} + g^2\,\Phi(C),$$

where $\Phi$ is built from Gaussian averages of $\phi$. The chaotic solution is the trajectory in
this effective potential; its decay time is the correlation time of the chaos, and it
**diverges as $g \to 1^+$**. That divergence is the reason "edge of chaos" is a real phenomenon
and not a slogan: the network's memory timescale becomes arbitrarily long at the transition.

### Structure: the part that matters for real circuits

Real connectivity is not i.i.d. Two extensions you will actually use:

**Excitatory/inhibitory structure.** With $E$ and $I$ populations having different means and
variances, the spectrum is no longer a symmetric disk. Rajan & Abbott showed that when each
column is *balanced* (the summed input from each neuron is fixed), the outlier eigenvalue
associated with the mean connectivity is eliminated and the bulk spectrum's radius is set by the
variances — with the radius depending on the E/I composition. This matters because the outlier,
not the bulk, usually controls stability in unbalanced networks.

**Non-normality.** $J$ random is strongly **non-normal**: $JJ^\top \neq J^\top J$. For non-normal
matrices, eigenvalues do not control transient behaviour. Even with every $\mathrm{Re}\,\lambda_i
< 0$, $\|e^{tA}\|$ can grow by a large factor before decaying — governed by the
**pseudospectrum**, not the spectrum. The relevant quantity is the numerical abscissa
$\omega(A) = \lambda_{\max}\big(\tfrac{1}{2}(A + A^\top)\big)$: transient growth occurs whenever
$\omega(A) > 0$ even though $\max\mathrm{Re}\,\lambda(A) < 0$.

**This is the single most important technical fact in this note for your purposes.** A stable,
non-normal network produces large, structured, input-specific transients — which is precisely
the phenomenology of the locust antennal lobe
([C2](../part2-case-studies/C2-transient-synchrony.md)). "Balanced amplification" (Murphy &
Miller) and "stability-optimized circuits" (Hennequin, Vogels & Gerstner) are the two developed
versions of this idea, and the second explicitly constructs networks that are stable yet produce
long, rich, amplified transients — a mechanistic competitor to both the chaotic-reservoir and
heteroclinic accounts.

---

## 3. Why the neuroscience forced it

Random matrix theory was mature; Girko's law predates the neural application. What was missing
was (i) the *nonlinear dynamical* theory (DMFT in this form), (ii) spectra of matrices with the
particular structures biology imposes — E/I sign constraints, column balance, low-rank additions
— and (iii) taking non-normality seriously as a computational resource rather than a nuisance.
All three were driven by neural questions and all three are now used outside neuroscience.

---

## 4. How to recognize the pattern elsewhere

Situation #8 announces itself as: **"I need to know what a generic instance does, and there is
no closed form."** The temptation is to simulate a few random instances and eyeball it. That is
fine for reassurance and useless as theory, because you cannot tell which features are generic
and which are finite-size accidents.

The payoff for doing it properly is usually larger than the original question, because a null
model theory is reusable across every claim in your subfield.

Practical rule for your own work: **every claim of the form "this circuit does X" needs a
matched random control, and matched means matched in the statistics that plausibly matter** —
spectral radius, E/I balance, degree distribution, and non-normality index. Matching only the
spectral radius while leaving non-normality free is the most common way to build a null that is
secretly much weaker than the real network.

---

## 5. Exercises

**Ex S7.1 ★** — Verify the circular law and the $g=1$ transition numerically.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
N = 2000
for g in [0.8, 1.2]:
    J = g*np.random.randn(N,N)/np.sqrt(N)
    ev = np.linalg.eigvals(J)
    print(g, "max |lambda|:", np.abs(ev).max(), " max Re(-1+lambda):", (ev.real-1).max())
```

Eigenvalues fill a disk of radius $\approx g$ (plot them). The Jacobian $-I+J$ has
$\max\mathrm{Re} \approx g-1$, so the origin is stable for $g<1$ and unstable for $g>1$.
Integrate the full nonlinear system at $g=0.8$ and $g=1.5$ from a small random initial
condition: the first decays to zero, the second sustains irregular activity indefinitely.
</details>

**Ex S7.2 ★★** — Demonstrate transient amplification in a stable non-normal network, and show
it is invisible in the eigenvalues.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
from scipy.linalg import expm

N = 200
J = 0.9*np.random.randn(N,N)/np.sqrt(N)
A = -np.eye(N) + J
print("max Re lambda:", np.linalg.eigvals(A).real.max())          # < 0: stable
print("numerical abscissa:", np.linalg.eigvalsh((A+A.T)/2).max()) # > 0: transient growth

ts = np.linspace(0, 5, 200)
norms = [np.linalg.norm(expm(A*t), 2) for t in ts]
print("peak amplification:", max(norms))                          # >> 1
```

Every eigenvalue has negative real part, so the eigen-analysis says "decays." Yet the operator
norm $\|e^{At}\|_2$ rises well above 1 before decaying — the network amplifies certain input
directions substantially on the way to zero.

To find *which* directions: the top left singular vector of $e^{At}$ at the peak time is the
optimally-amplified initial condition. In an E/I network these turn out to be
difference-like patterns amplified into sum-like patterns, which is the content of balanced
amplification.

Why this matters for [C2](../part2-case-studies/C2-transient-synchrony.md): a stable non-normal
network produces long, input-specific, high-dimensional transients *without being chaotic and
without heteroclinic saddles*. Any claim that observed olfactory transients require chaos or
heteroclinic channels must first rule out plain non-normal amplification. It usually is not
ruled out.
</details>

**Ex S7.3 ★★★** — *(Situation, not object.)* For a claim in your own work, specify the null
model completely: what statistics are matched, what are free, and what property of the real
network you are asserting is *not* generic. Then estimate whether the null is actually weaker
than the data in that specific respect.

<details markdown="1"><summary>Solution</summary>

Open. The discipline to impose: write the null as an explicit generative model with named matched
statistics, and name the one number that separates data from null. If you cannot name that
number, you do not have a claim, you have an observation.

For "the antennal lobe decorrelates," a defensible specification: match $N$, E/I ratio, mean and
variance of connection strength, spectral radius, degree distribution, and the input statistics.
Leave free: the specific realization. Assert as non-generic: that the *asymptotic* correlation
between odour representations is approximately independent of their initial correlation. As
noted in [00 Ex 0.5](../00-orientation/README.md), a generic network decorrelates by roughly a
common factor — preserving the ordering and rough spacing of initial correlations — whereas true
whitening compresses them toward a common value. So the discriminating statistic is the *slope*
of asymptotic vs. initial correlation, and the prediction is slope $\approx 0$ for whitening,
slope $\approx 1$ for generic expansion.

That is a well-posed claim. Note how much work went into getting from "it decorrelates" to a
single measurable slope, and that essentially all of that work was figuring out what the null
does.
</details>

---

## 6. Reading

- **Sompolinsky, Crisanti & Sommers (1988)**, *Chaos in random neural networks* (Phys. Rev.
  Lett.) — read it for: the transition and the mean-field derivation. Four pages; read all of
  them.
- **Rajan & Abbott (2006)**, *Eigenvalue spectra of random matrices for neural networks* (Phys.
  Rev. Lett.) — read it for: what E/I structure and column balance do to the spectrum.
- **Ahmadian, Fumarola & Miller (2015)**, on properties of networks with partially structured
  and partially random connectivity — read it for: the general machinery for spectra of
  structured-plus-random matrices, which you will need for
  [S-08](S-08-low-rank-connectivity.md).
- **Murphy & Miller (2009)**, *Balanced amplification: a new mechanism of selective
  amplification of neural activity patterns* (Neuron) — read it for: non-normality as a
  computational mechanism, in the clearest possible setting.
- **Hennequin, Vogels & Gerstner (2014)**, on optimal control of transient dynamics in balanced
  networks (Neuron) — read it for: networks engineered to be stable yet produce long rich
  transients. The most direct mechanistic competitor to reservoir and heteroclinic accounts of
  [C2](../part2-case-studies/C2-transient-synchrony.md).
- **Trefethen & Embree (2005)**, *Spectra and Pseudospectra* — read it for: the definitive
  treatment of why eigenvalues mislead. Chapters 14–17 if you are in a hurry.
- **Kadmon & Sompolinsky (2015)**, on the transition to chaos in random neuronal networks
  (Phys. Rev. X) — read it for: the modern, more general treatment of the transition, including
  spiking and more realistic transfer functions.
