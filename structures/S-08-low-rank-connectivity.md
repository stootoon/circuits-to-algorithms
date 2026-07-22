---
title: S-08 Low-rank connectivity
parent: Structures
nav_order: 8
---

# S-08 — Connectivity space: when the population's *statistics* are the algorithm
{: .no_toc }

> **The object.** The joint distribution over neurons of their connectivity loadings — each
> neuron a point in a low-dimensional "connectivity space," the population a cloud, and the
> cloud's statistics the carrier of the computation. Plus **rank** as a resource with a
> minimum required per task.
> **Born from.** Wanting to read a trained RNN's algorithm off its weights.
> **Mathematical home.** Mean-field theory of structured random matrices.
> **Situation.** #4, the structure→dynamics dictionary — in the continuous-statistics regime,
> where [S-03](S-03-combinatorial-threshold-linear-networks.md) is the combinatorial regime.
> **Novelty.** **Genuinely new.** "The distribution over connectivity space" is an object that
> did not exist as a modelling primitive before this line of work.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Where it comes from

You train an RNN on a task. It works. You now have an $N \times N$ weight matrix and no idea
what it is doing. Fixed-point analysis
([Unit 01](../part1-foundations/01-dynamical-systems.md), [S1](../part3-synthesis/S1-rnn-as-model-organism.md))
tells you the *dynamical* structure, but not how the connectivity produces it.

The insight is that trained networks are typically **low-rank plus random**, and that low-rank
structure is legible. Better: it is legible in a specific, transferable way — through the
statistics of a cloud of points, one per neuron.

---

## 2. The object, precisely

Write

$$J = \underbrace{\frac{g}{\sqrt{N}}\,\chi}_{\text{random}}
\;+\; \underbrace{\frac{1}{N}\sum_{r=1}^{R} m^{(r)}\big(n^{(r)}\big)^{\!\top}}_{\text{rank } R},
\qquad \chi_{ij}\sim\mathcal{N}(0,1),$$

with dynamics $\dot{x} = -x + J\phi(x) + I$.

**The reduction.** Project the dynamics onto the $m^{(r)}$:

$$x_i(t) \;=\; \sum_{r=1}^R \kappa_r(t)\, m_i^{(r)} \;+\; I_i \;+\; \big(\text{random component}\big),
\qquad \kappa_r(t) = \frac{1}{N}\sum_i n_i^{(r)}\,\phi\big(x_i(t)\big).$$

So the $N$-dimensional dynamics collapse onto **$R$ latent variables** $\kappa_r$, whose
evolution is closed under mean-field. The random part does not vanish — it sets the effective
gain and can drive chaos exactly as in [S-07](S-07-random-matrices-and-chaos.md) — but it does
not carry the computation.

**The new object.** Each neuron $i$ is now a point

$$\big(m_i^{(1)},\dots,m_i^{(R)},\, n_i^{(1)},\dots,n_i^{(R)},\, I_i\big) \in \mathbb{R}^{2R+n_{\text{in}}},$$

and the population is a cloud of $N$ such points. The mean-field equations depend on the
population **only through the distribution of this cloud**. That distribution is the object.
Two networks with completely different weight matrices but the same connectivity-space
distribution compute the same thing. This is a precise, non-trivial statement of what "same
algorithm" means at the level of connectivity — a genuine equivalence relation, which is exactly
what [00 §1](../00-orientation/README.md) said an algorithm is.

**Consequences that make it a theory rather than a reformulation:**

- **Gaussian is not enough.** If the cloud is a single multivariate Gaussian, the achievable
  computations are limited — the latent dynamics inherit strong restrictions. Dubreuil,
  Valente, Beiran, Mastrogiuseppe & Ostojic showed that many tasks require the cloud to be a
  **mixture** of Gaussians, i.e. genuinely distinct sub-populations. The number of mixture
  components is a meaningful complexity measure of the task, and it makes contact with
  experimentally identified cell types — a rare case where a theoretical quantity has an
  obvious biological referent.
- **Rank is a resource.** Beiran and colleagues characterized the minimal rank needed to perform
  given classes of tasks. Rank $1$ gives you a one-dimensional latent: bistability, integration
  along one direction. Rank $2$ buys oscillations and ring-like structures. This is a
  *complexity theory for recurrent computation* measured in rank rather than in neurons.
- **It is a reverse-engineering method.** Given a trained RNN, fit a low-rank approximation, plot
  the neurons in connectivity space, look for clusters, and read the algorithm off the resulting
  reduced system.

---

## 3. Why the neuroscience forced it

Physics had structured random matrices and knew about rank-one perturbations (the BBP
transition). What it did not have was the framing in which **the low-rank part is the signal and
the random part is the substrate**, with the population distribution over loadings as the
primary object of study.

That framing comes from a specifically neuroscientific need: to explain a computation in terms
of an *interpretable, biologically-meaningful* decomposition of connectivity, rather than to
compute a spectrum. The mathematics that resulted — mean-field theory parameterized by a
distribution over neuron-level loadings — is genuinely new machinery.

---

## 4. How to recognize the pattern elsewhere

You are here when: **the high-dimensional system's behaviour is controlled by a low-dimensional
structure embedded in it, and the right description is the joint distribution of how each unit
participates in that structure.**

The transferable move is subtle and worth naming: *stop describing the system by its parameters,
and start describing it by the empirical distribution of its units' parameter-tuples.* This
converts an $N$-dependent description into an $N$-independent one — a distribution — which is
what makes a thermodynamic limit and hence a theory possible. It is the same move as going from
a list of particle positions to a density, and it applies far outside RNNs.

For [C2](../part2-case-studies/C2-transient-synchrony.md) the application is concrete: if PN–LN
connectivity has low-rank structure, the odour-evoked transient is a trajectory in a
low-dimensional latent space whose geometry is set by the connectivity-space cloud. Fitting a
low-rank model to a data-constrained antennal lobe model, and examining whether the loading
cloud is unimodal or clustered, would say directly whether the transient dynamics require
distinct PN/LN sub-populations — a question with an anatomical answer available to check
against.

---

## 5. Exercises

**Ex S8.1 ★★** — Build a rank-one network that implements bistability, and verify the latent
reduction.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
N, g = 1500, 0.4
m = np.random.randn(N); n = np.random.randn(N)
n += 2.5*m                                     # correlate n with m -> positive feedback
J = g*np.random.randn(N,N)/np.sqrt(N) + np.outer(m,n)/N

def run(x0, T=60, dt=0.05):
    x = x0.copy(); traj=[]
    for _ in range(int(T/dt)):
        x += dt*(-x + J@np.tanh(x))
        traj.append(n@np.tanh(x)/N)            # the latent kappa
    return np.array(traj)

print(run(+0.5*m)[-1], run(-0.5*m)[-1])        # two opposite-sign fixed points
print(run(1e-3*np.random.randn(N))[-1])        # near zero -> stays near a branch
```

You should see two stable values of $\kappa$ of opposite sign: bistability. The whole
$1500$-dimensional system is described by one number. Verify the reduction by checking that
$x_i \approx \kappa m_i$ at the fixed point (regress $x$ on $m$; $R^2$ should be high, with the
residual set by $g$).

Then vary the correlation $\langle m_i n_i\rangle$: below a threshold the network is monostable
at $\kappa = 0$; above it, bistable. That threshold is computable in closed form from the
mean-field equations, and finding it is the natural follow-up.
</details>

**Ex S8.2 ★★** — Show why a rank-one network with a *single Gaussian* connectivity-space
distribution cannot produce a stable oscillation, and what rank-two buys you.

<details markdown="1"><summary>Solution</summary>

With rank one the latent dynamics are one-dimensional: $\dot{\kappa} = F(\kappa)$ for a scalar
$\kappa$ (after the mean-field average over the Gaussian cloud). A scalar autonomous ODE on
$\mathbb{R}$ cannot oscillate — trajectories are monotone between fixed points, by the
intermediate value theorem. So the only attractors available are fixed points: monostability,
bistability, or multistability.

Rank two gives $(\kappa_1,\kappa_2) \in \mathbb{R}^2$, and planar autonomous systems admit limit
cycles (Poincaré–Bendixson permits them; they are excluded in 1-D). Concretely, choosing the
loadings so that the effective $2\times 2$ latent Jacobian has a complex-conjugate pair crossing
the imaginary axis gives a Hopf bifurcation and a stable oscillation.

This is the rank-as-resource statement in miniature, and note its form: **a topological
constraint on low-dimensional dynamics becomes a constraint on connectivity rank.** Rank $\geq 2$
for oscillation; rank $\geq 2$ for a ring attractor (which needs a 2-D latent plane to host the
circle — the [C4](../part2-case-studies/C4-continuous-attractors.md) connection); higher rank for
richer repertoires. Combined with the mixture result, you get a genuine complexity theory:
tasks are graded by (rank, number of populations).
</details>

**Ex S8.3 ★★★** — *(Situation, not object.)* Take a high-dimensional model you use and rewrite
its description as a distribution over per-unit parameter tuples. What becomes visible?

<details markdown="1"><summary>Solution</summary>

Open. Worked template for a data-constrained antennal lobe model:

Each PN $i$ has a tuple: (its odour input loading, its loading onto the dominant PC of the
population transient, its LN input weight, its LN output weight). Plot the cloud in these
coordinates.

What becomes visible: (a) whether PNs form clusters — i.e. whether the model implicitly relies
on distinct sub-populations, which is checkable against glomerular identity; (b) whether the
input loading and the transient loading are *correlated*, which is the signature of the input
directly driving the latent rather than the recurrent dynamics generating it; (c) whether the
LN loadings are aligned with the transient direction, which distinguishes "inhibition shapes the
trajectory" from "inhibition provides uniform gain control" — and that distinction is exactly
the C8-vs-C2 question about whether lateral inhibition is normalization or dynamics.

The general payoff: you stop asking "what are the weights" and start asking "what is the shape
of the cloud," and the second question has answers that transfer between models.
</details>

---

## 6. Reading

- **Mastrogiuseppe & Ostojic (2018)**, *Linking connectivity, dynamics, and computations in
  low-rank recurrent neural networks* (Neuron) — read it for: the framework. This is the founding
  paper; the mean-field derivation in the supplement is worth working through.
- **Dubreuil, Valente, Beiran, Mastrogiuseppe & Ostojic (2022)**, on the role of population
  structure in computations through neural dynamics (Nature Neuroscience) — read it for: the
  Gaussian-mixture result and the link to cell types. This is where the theory becomes
  biologically pointed.
- **Beiran, Dubreuil, Valente, Mastrogiuseppe & Ostojic (2021)**, on shaping dynamics with
  multiple populations in low-rank recurrent networks (Neural Computation) — read it for: rank
  as a resource, and the minimal-rank results.
- **Schuessler, Dubreuil, Mastrogiuseppe, Ostojic & Barak (2020)**, on the dynamics of random
  recurrent networks with correlated low-rank structure — read it for: what training actually
  does to connectivity, i.e. why the low-rank-plus-random ansatz is the right one.
- **Ahmadian, Fumarola & Miller (2015)** — read it for: the spectral theory underlying
  structured-plus-random matrices. The technical prerequisite.
- **Sussillo & Barak (2013)** — read it for: the complementary dynamical (rather than
  structural) reverse-engineering method. Use both; they answer different questions.
