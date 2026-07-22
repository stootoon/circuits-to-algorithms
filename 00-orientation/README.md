---
title: 00 · Orientation
nav_order: 2
permalink: /orientation/
---

# 00 — Marr's levels, multiple realizability, and the identifiability problem
{: .no_toc }

> **Primitive.** The abstraction ladder: what kind of object an "algorithm" is, why it is not
> determined by the circuit, and what therefore counts as evidence for one.
> **Prerequisites.** None.
> **Structures thread.** [`structures/README.md`](../structures/README.md) — the methodology
> note there is the Thread-B counterpart to this unit.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Why this matters for circuit→algorithm conversion

Almost everything that goes wrong in this field goes wrong in the first ten minutes, when
someone decides what kind of thing they are looking for. So we start there.

You already know Marr's three levels. What is usually taught badly is the *logical relationship*
between them, so let us be precise about it, because the entire course depends on getting this
right.

- **Level 1 (computational).** What problem is being solved, and why is it the right problem?
  Formally: a specification. A map $f : \mathcal{X} \to \mathcal{Y}$, or — more usefully — an
  objective $\mathcal{L}$ together with constraints, such that the system's behaviour is
  $\arg\min \mathcal{L}$ subject to those constraints.
- **Level 2 (algorithmic).** The procedure and representation used to compute $f$. Formally:
  an equivalence class of procedures under some notion of "same algorithm."
- **Level 3 (implementational).** The physical substrate. Formally: the circuit, with its
  conductances, delays, and wiring.

The relationships are **many-to-one going up** and — this is the part that bites — **not
functionally determined going down**. A level-2 algorithm has many level-3 implementations
(multiple realizability). But it is also true that a level-3 circuit description does not
uniquely determine a level-2 algorithm, because "same algorithm" is an equivalence relation
*you choose*, and different choices carve the same circuit differently. C6 (grid cells) and C8
(normalization) are the course's two cleanest demonstrations that one circuit can carry two
incompatible, both-defensible level-2 readings.

So "extracting the algorithm" is not measurement. It is **model selection over descriptions**,
and it needs the epistemology that goes with that.

---

## 2. The core theory

### 2.1 What makes a level-2 claim non-vacuous

Here is the failure mode. Take any circuit with dynamics $\dot{x} = F(x, u)$. I can always say:
"the algorithm is: integrate $F$." That is true, and it is worthless. It is worthless for a
reason we can make precise.

**Compression.** A good level-2 description is *shorter* than the level-3 description it
replaces, and remains predictive. This is Kolmogorov-flavoured but usable in practice: if your
algorithmic story requires as many free parameters as the circuit has synapses, you have not
abstracted, you have re-parameterized. The ring attractor story for the ellipsoid body (C4)
replaces $\sim\!50$ neurons and their connectivity with *one* number — the bump's angular
position $\theta(t)$ — plus a one-line update rule $\dot{\theta} = \omega$. That is a
compression of enormous ratio, and it still predicts.

**Transfer.** A good level-2 description applies to systems that share no level-3 properties. TD
learning (C7) was extracted from midbrain dopamine and it describes a chess engine. Cross-
correlation (C3) describes a barn owl and a radar receiver. If your algorithm only ever
describes the one circuit it was extracted from, suspect that you have described the circuit.

**Falsifiability with quantitative teeth.** The dopamine/TD conversion is the field's crown
jewel not because it is elegant but because it said, in advance, that the response to a fully
predicted reward must be *zero* and the response to an omitted reward must go *below baseline* —
and those are strange, non-obvious, and were then observed. Contrast a story of the form "the
circuit performs Bayesian inference," which, without a specified likelihood, prior, and
approximation scheme, forbids nothing.

Formally, a useful criterion: your level-2 claim must induce a constraint on the level-3
behaviour that a generic circuit of the same architecture would violate. If a random network
with the same connectivity statistics also satisfies your prediction, the prediction is testing
the architecture, not your algorithm. **Always build that null.**

### 2.2 The identifiability problem, stated properly

Let $\mathcal{C}$ be a circuit and let $\mathcal{A}$ be a candidate algorithm. Define the
implementation relation $\mathcal{C} \models \mathcal{A}$ to mean: there is a map $\phi$ from
circuit states to algorithm states such that $\phi$ commutes with the dynamics —

$$\phi\big(\Phi^{\mathcal{C}}_t(x)\big) \;=\; \Phi^{\mathcal{A}}_t\big(\phi(x)\big)$$

i.e. $\phi$ is a semiconjugacy from the circuit's flow to the algorithm's. This is the honest
formalisation of "the circuit implements the algorithm," and it makes the problem obvious:
**$\phi$ is not unique, and it is not observable.** Worse, if $\phi$ is allowed to be arbitrarily
complicated, *every* sufficiently rich circuit implements *every* algorithm — this is the
computational-triviality objection (Putnam, Searle), and the standard repair is to demand that
$\phi$ be simple, causal-structure-respecting, and fixed in advance rather than fitted post hoc.

That repair is not philosophy pedantry; it is your experimental design. "Fixed in advance"
means: *state the decoder before you see the data.*

### 2.3 The microprocessor result

Jonas & Kording (2017) took a MOS 6502 — a system where the ground-truth algorithm is fully
known — ran it on *Donkey Kong*, *Space Invaders*, and *Pitfall*, and applied the standard
neuroscience toolkit: lesions (transistor knockouts), tuning curves, spike-triggered averages,
connectomics, dimensionality reduction, Granger causality. The analyses produce results that
*look* like neuroscience results. Single transistors appear selective for individual games.
Lesion studies identify "game-specific" transistors. None of it recovers the algorithm, or
even the concept of an instruction, a register, or a clock.

Read this early and take the correct lesson, which is *not* "neuroscience is hopeless." It is:
**the standard toolkit is a level-3 toolkit, and level-3 tools do not spontaneously yield
level-2 descriptions.** You need an independent source of level-2 hypotheses. In this course
that source is almost always either (a) a normative derivation — Units 03, 04, 06 — or (b) a trained
model you can reverse-engineer — S1. There is no third option that has ever worked.

### 2.4 The two master strategies

Everything in Phase 1 is one of these:

**Bottom-up / mechanistic.** Reduce the circuit's dynamics until a recognizable primitive falls
out. Find the fixed points and classify them (Unit 01). Find the low-dimensional latent and the
connectivity structure generating it (Unit 02). The output is a *dynamical* characterization: this
circuit is an integrator, a switch, a filter, a sequence generator.

**Top-down / normative.** Posit that the circuit solves some problem optimally, derive the
optimal solution analytically, and check the derived solution against the circuit (Units 03, 04, 06).
The output is a *functional* characterization: this circuit whitens, this circuit computes a
posterior, this circuit cancels a prediction.

Each alone is weak. Bottom-up gives you mechanism without meaning — you know it is an
integrator, not what it integrates or why. Top-down gives you meaning without mechanism, and is
where just-so stories breed. **The results that have lasted are the ones where both were run
independently and met in the middle.** Retinal whitening: Barlow's normative argument met
Kuffler's centre-surround receptive fields. Dopamine: Sutton & Barto's TD algorithm met
Schultz's recordings. Fly heading: Zhang's ring-attractor theory (1996) met Seelig &
Jayaraman's imaging (2015), nineteen years later, with the theory published first and
essentially correct.

That last one is the template. Write the theory so sharply that when the data arrive two decades
later, they can falsify it.

---

## 3. Worked example: two readings of one circuit

Consider a two-population mutual-inhibition circuit — a motif that appears in perceptual
decision-making (C10), in olfactory lateral inhibition (C2), and in half the models in the
literature:

$$\dot{r}_1 = -r_1 + [\,I_1 - w r_2 + \beta r_1\,]_+, \qquad
\dot{r}_2 = -r_2 + [\,I_2 - w r_1 + \beta r_2\,]_+.$$

**Reading A (level 2 = winner-take-all).** For $w$ large, the symmetric fixed point is unstable;
the system has two stable fixed points, one per population. The algorithm is $\arg\max(I_1,
I_2)$ — a decision, a selection, a max operation. This reading is correct.

**Reading B (level 2 = evidence accumulation / SPRT).** Change coordinates to
$u = r_1 - r_2$, $s = r_1 + r_2$. Near the symmetric fixed point, when $\beta - 1 \approx w$,
the $u$ dynamics have a near-zero eigenvalue: $\dot{u} \approx \lambda u + (I_1 - I_2)$ with
$\lambda \approx 0$. So $u(t) \approx \int_0^t (I_1 - I_2)\,dt'$ — the circuit *integrates the
difference in evidence* until $|u|$ hits the basin boundary. That is a drift-diffusion process,
which is the continuous limit of the sequential probability ratio test, which is Wald-optimal.
This reading is also correct.

These are not the same algorithm. One is a memoryless comparator; the other is an optimal
sequential statistical test with a speed–accuracy tradeoff. They are readings of the *same
equations* in different parameter regimes and on different timescales, and the circuit does not
come labelled with which regime it is in. Deciding between them is an empirical question about
$\beta - 1 - w$ and about the timescale of the input — i.e. **the level-2 question reduces to a
quantitative level-3 measurement, but only once you have both candidate algorithms in hand.**

You could not have discovered Reading B by staring at the connectivity. You discovered it by
knowing what SPRT looks like and going to check.

That is the whole course, in one example.

---

## 4. Exercises

**Ex 0.1 ★** — State the compression criterion quantitatively. For the ring attractor (C4), the
circuit has $N$ neurons with $N^2$ synaptic weights; the algorithm has one state variable. Under
what condition on the connectivity is the compression *legitimate* rather than a lossy summary
that discards computationally relevant structure?

<details markdown="1"><summary>Solution</summary>

The compression is legitimate when the discarded degrees of freedom are (i) not driven by the
inputs of interest and (ii) contracting, so trajectories collapse onto the retained manifold
faster than the retained variable evolves.

Concretely: linearize about a point on the attractor manifold $\mathcal{M}$. The Jacobian has
one eigenvalue $\lambda_0 \approx 0$ (the Goldstone mode tangent to $\mathcal{M}$, arising from
the continuous rotational symmetry of $W(\theta - \theta')$) and $N-1$ eigenvalues with
$\mathrm{Re}\,\lambda_i \le -\kappa < 0$. The reduction $x \mapsto \theta$ is valid on timescales
$t \gg 1/\kappa$, with error $O(\|u\|/\kappa)$ for input magnitude $\|u\|$.

So the condition is a **spectral gap**: $\kappa \gg |\lambda_0|$ and $\kappa \gg$ the input
timescale. If a second eigenvalue is also near zero — e.g. the bump *amplitude* mode is slow —
then amplitude is computationally relevant, the true latent is 2-D, and the 1-D story is lossy.
This is not hypothetical: amplitude modulation on the fly heading bump carries information
about confidence/stimulus reliability, so the honest reduction there is 2-D.

The general statement is normal hyperbolicity of $\mathcal{M}$, which is what guarantees the
reduced description persists under perturbation.
</details>

**Ex 0.2 ★** — Give an example, from any field, of two systems that implement the same algorithm
with no shared level-3 properties. Then give an example of two systems with near-identical
level-3 descriptions implementing different algorithms.

<details markdown="1"><summary>Solution</summary>

*Same algorithm, disjoint implementation:* the barn owl's nucleus laminaris and a phased-array
radar receiver both compute cross-correlation to localize a source. Shared level-3 properties:
none. Shared level-2 structure: delay-and-multiply-and-integrate. (C3.)

*Same implementation, different algorithms:* the mutual-inhibition circuit in §3, in two
parameter regimes — WTA vs. SPRT. Or, more strikingly and biologically: the crustacean
stomatogastric ganglion's pyloric circuit versus its gastric mill circuit share neurons and
much connectivity but generate different rhythms serving different functions, with the
difference set by neuromodulatory state rather than wiring. Neuromodulation is precisely a
mechanism for making one circuit run different algorithms — which is why "the connectome
determines the computation" is false as stated. (S2.)
</details>

**Ex 0.3 ★★** — Read Jonas & Kording (2017). Identify one analysis in the paper that *would*
have recovered part of the 6502's algorithm if it had been applied with a correct level-2
hypothesis in hand. Argue for what that hypothesis would have to be.

<details markdown="1"><summary>Solution</summary>

The most promising is the dimensionality-reduction / whole-system-activity analysis. It fails as
applied because it is unsupervised and looks for low-dimensional continuous structure, which is
the wrong prior for a digital machine.

With the right level-2 hypothesis — "this system executes a fetch–decode–execute cycle with a
periodic clock and a small set of discrete internal states" — the corresponding analysis is:
(1) recover the clock by finding the dominant periodicity in aggregate switching activity;
(2) *bin transistor states by clock phase*; (3) look for a small number of discrete recurring
state-vectors within each phase, i.e. cluster in the binary state space rather than seeking
continuous latents.

That analysis would recover the instruction cycle and plausibly the opcode structure. The point
generalizes and is the main lesson of the unit: the analysis method must be *chosen by* the
level-2 hypothesis. Unsupervised dimensionality reduction encodes a hypothesis — "the
computation lives on a smooth low-dimensional manifold" — that is often smuggled in unexamined.
For continuous attractors (C4) that hypothesis is right; for a 6502 it is wrong; for a
transient, synchrony-gated olfactory code (C2) it is *exactly the open question*.
</details>

**Ex 0.4 ★★** — Construct the semiconjugacy $\phi$ explicitly for the $N$-neuron ring attractor
and the algorithm $\dot{\theta} = \omega$. Where does the commutation fail, and by how much?

<details markdown="1"><summary>Solution</summary>

Take $\phi(x) = \arg\big(\sum_j x_j e^{i\theta_j}\big)$, the population vector angle — the
standard decoder, and note it is fixed in advance, not fitted.

On the attractor manifold with an exact bump profile $x_j = f(\theta_j - \theta)$ and symmetric
$W(\theta-\theta')$, $\phi$ commutes exactly: rotational equivariance of $W$ means the flow maps
bumps to bumps and $\phi$ tracks the shift.

Commutation fails in three places, in increasing order of severity:

1. **Off-manifold transients.** For $x \notin \mathcal{M}$, $\phi(x)$ is defined but the flow is
   not tangent to $\mathcal{M}$. Error decays as $e^{-\kappa t}$ (Ex 0.1).
2. **Finite $N$.** Discrete sampling breaks continuous rotational symmetry down to the cyclic
   group $\mathbb{Z}_N$, replacing the attractor circle with $N$ (or $2N$) discrete fixed points
   separated by saddles. The bump then has preferred angles and drifts to the nearest one; the
   residual "pinning" force scales roughly as $e^{-cN\sigma}$ for bump width $\sigma$ — 
   exponentially small in the number of neurons per bump width, which is why moderate $N$
   already looks continuous.
3. **Heterogeneity.** Quenched disorder in $W$ of magnitude $\epsilon$ pins the bump with force
   $O(\epsilon)$, producing systematic drift $\dot{\theta} = \omega + g(\theta)$ with
   $\|g\| = O(\epsilon)$ and a nonzero winding-number condition for whether full rotation
   survives at all. This is the dominant error in real circuits and is the subject of C4 §2.

The moral: the algorithm $\dot\theta = \omega$ is an *idealization* whose fidelity is
quantifiable. That is what a good level-2 claim looks like — not "the circuit is a ring
attractor," but "the circuit is a ring attractor to within drift of $O(\epsilon)$ rad/s, and
here is the measurement of $\epsilon$."
</details>

**Ex 0.5 ★★★** — Take a circuit→algorithm claim from your own field that you believe. Write
down (a) the claim as a semiconjugacy, specifying $\phi$; (b) the prediction it makes that a
generic circuit of the same architecture would violate; (c) the experiment that would kill it.
If you cannot do (b), the claim is architectural, not algorithmic — diagnose which.

<details markdown="1"><summary>Solution</summary>

No canonical answer; this is the first entry in your conversion log. A worked template, using
"the locust antennal lobe decorrelates odour representations":

(a) $\phi$: the map from the $\sim\!800$-dimensional PN population state to the trajectory of
pairwise correlations $\rho_{ab}(t)$ between the representations of odours $a$ and $b$. The
"algorithm" is a flow on correlation space with $\dot{\rho} < 0$.

(b) Prediction distinguishing it from generic recurrent dynamics: decorrelation must be
*specific* — it should be strongest for the most similar odour pairs, and the asymptotic
correlation should be approximately independent of the initial correlation (a genuine
whitening-like normalization, not just a generic transient expansion which would decorrelate
everything by roughly a common factor). A random recurrent network with matched spectral radius
decorrelates too; it does not produce initial-condition-independent asymptotics.

(c) Killer experiment: present a parametric family of odour mixtures spanning a wide range of
initial PN correlations. If asymptotic correlation is a monotone function of initial
correlation with slope $\approx 1$ on log axes, the circuit is a generic expander, not a
decorrelator, and the normative story is dead. Compare against a spectral-radius-matched random
network null in the same analysis pipeline.

Note that step (b) is the one that took real work, and that it required knowing what generic
recurrent dynamics do — which is [S-07](../structures/S-07-random-matrices-and-chaos.md).
</details>

---

## 5. Reading path

- **Marr (1982)**, *Vision*, ch. 1 — read it for: the three levels stated by the person who
  meant them, including the parts everyone forgets (the primacy of level 1, and that levels are
  levels of *explanation*, not of *description*).
- **Jonas & Kording (2017)**, *Could a neuroscientist understand a microprocessor?* (PLOS
  Comput Biol) — read it for: the demonstration that level-3 tools do not yield level-2
  answers, on a system with known ground truth.
- **Krakauer, Ghazanfar, Gomez-Marin, MacIver & Poeppel (2017)**, *Neuroscience needs
  behavior: correcting a reductionist bias* (Neuron) — read it for: the argument that level-1
  specification requires studying behaviour, not just circuits; the necessary complement to
  Jonas & Kording.
- **Marder & Goaillard (2006)**, review of variability and degeneracy in neuronal and network
  parameters (Nat Rev Neurosci) — read it for: hard experimental evidence that the level-3
  parameters are not identifiable even in principle. Revisited properly in
  [S2](../part3-synthesis/S2-degeneracy-and-limits.md).
- **Kording, Blohm, Schrater & Kay (2020)**, *Appreciating the variety of goals in
  computational neuroscience* — read it for: a taxonomy of what different modelling efforts are
  actually *for*, which prevents a lot of pointless argument.

Optional but clarifying if the semiconjugacy framing appealed: **Chalmers (1996)**, *Does a rock
implement every finite-state automaton?* — the cleanest treatment of the triviality objection
and the standard repair.

---

## 6. Where it breaks

Marr's hierarchy is a useful ladder and a misleading picture of science. Three honest problems:

**Levels are not cleanly separable in biology.** Marr's own framing assumes near-decomposability
— that you can specify the algorithm without knowing the hardware. Evolution had no reason to
respect this. Dendritic nonlinearities, spike-timing, and neuromodulation are level-3 facts that
*change what algorithm is cheap*, and therefore leak into level 2. The clean cases in this
course (C3, C4, C7) are clean partly because they were selected for being clean.

**Level 1 is often unavailable.** For vision, "recover scene structure" is a defensible level-1
statement. For the locust antennal lobe, what *is* the level-1 problem? "Classify odours"
under-determines everything; the animal must also generalize across concentration, segment
mixtures, detect novelty, and do it in 300 ms. Absent a defensible level-1 specification, the
normative strategy has nothing to optimize, and you are reduced to bottom-up characterization
plus taste. Much of C2 is about this predicament, and it is the reason the capstone exists.

**"The algorithm" may not be unique or even well-posed.** A circuit shaped by evolution is not
a designed artifact with an intended function. It may implement several partially-overlapping
computations, none cleanly separable, all of them contingent on a phylogenetic accident. The
assumption that a crisp algorithm exists and is waiting to be found is itself a hypothesis, and
it is not always true. Hold it as a working assumption that has paid off spectacularly in the
handful of cases in Phase 2 — while remembering that Phase 2 is a survey of successes, which is
the most misleading sample in science.
