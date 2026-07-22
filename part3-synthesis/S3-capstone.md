---
title: "S3 · Capstone"
parent: Synthesis
nav_order: 3
---

# S3 — Capstone: a level-2 specification for transient-synchrony olfactory coding
{: .no_toc }

> **The task.** Produce a complete, defensible level-2 specification of what the locust antennal lobe and mushroom body are computing during the odour-evoked transient — one that would score ≥20/24 on the checklist in `./S2-degeneracy-and-limits.md` §4.
> **Why this system.** Because the data are extraordinary and the algorithm is genuinely unsettled. Thirty years of recordings have produced a phenomenology of surpassing richness — oscillatory synchronisation, transient trajectories, sparse readout, phase-gated one-shot plasticity — and at least four incompatible algorithmic readings, none of which has been made to bind hard enough to fail.
> **How to use this note.** As a staged project with deliverables, not as reading. Each milestone ends in an artefact you can show someone. Expect M0–M3 to take a few weeks and M4–M7 to take a few months.
> **Structures thread.** `../structures/S-11-expanders-and-optimal-degree.md` fixes the PN→KC matrix $A$ in M0; `../structures/S-03-combinatorial-threshold-linear-networks.md` gives a second, discrete route to the sequence dynamics of H3; `../structures/S-06-hyperbolic-odor-space.md` bears on the input geometry you assume in M1. Index in `../structures/README.md`.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## M0 — State the computational problem (level 1)

**Deliverable: a two-page problem statement containing (a) inputs, (b) outputs and who reads them, (c) an objective functional, (d) a constraint table, and (e) an explicit list of what you are *not* explaining.**

Most work in this area skips straight to mechanism. Do not. Per `./S2-degeneracy-and-limits.md` §5, level 1 is the load-bearing part.

**Inputs.** An odour is a time-varying mixture concentration $c(t)\in\mathbb{R}^{N_o}_{\ge0}$, delivered in the wild as an intermittent plume: filaments with sharp onsets, log-normally distributed concentrations, and inter-filament intervals of tens to hundreds of milliseconds. Transduction gives ORN activity $u(t)\in\mathbb{R}^{N_r}$, $N_r\approx 90{,}000$ receptor neurons of $\sim$tens of types.

**Outputs and their consumer.** The antennal lobe produces PN activity $x(t)\in\mathbb{R}^{N_p}$, $N_p\approx830$. This is read by $\sim50{,}000$ Kenyon cells, each sampling roughly half the PN population with a high spike threshold (Jortner, Farivar & Laurent, 2007), under global feedback normalisation by the wide-field GGN (Papadopoulou, Cassenaer, Nowotny & Laurent, 2011), with an effective integration window set by oscillation-locked feedforward inhibition (Perez-Orive et al., 2002). KC output converges on a few tens of mushroom-body extrinsic neurons whose synapses undergo one-shot, phase-gated, neuromodulator-conditional STDP (Cassenaer & Laurent, 2007, 2012).

**Objective.** Write it down, commit to it, and note that this is the step at which the whole enterprise can be a just-so story (`./S2-degeneracy-and-limits.md` §3.2). My proposal, which you should feel free to reject in writing:

$$\mathcal J[E] = \underbrace{I\big(\,\iota\,;\,\Theta(A\,x_{[0,T]})\,\big)}_{\text{identity decodable by a \emph{fixed random sparse} readout}} \;-\;\lambda\,\underbrace{\mathbb E\big[\|\Theta(Ax)\|_0\big]}_{\text{metabolic/sparsity cost}} \;-\;\mu\,\underbrace{T}_{\text{decision latency}},$$

where $\iota$ is odour identity, $A$ is a *random, fixed* PN→KC matrix with the measured statistics, and $\Theta$ is a high threshold with GGN normalisation. The critical modelling commitment here is that **the readout is not optimised**: it is random and fixed, so the encoder must do all the work. That single commitment turns a vague "the AL makes odours separable" into a constrained optimisation with a specific answer, and it is what makes the specification falsifiable.

**Constraints table** (fill this in with citations and numbers): $N_p$, $N_{KC}$, connectivity fraction $\approx0.5$, KC threshold, $\tau_m^{KC}$, LFP frequency $\approx20$ Hz, PN spike jitter relative to LFP phase, GGN gain, plasticity window and its phase dependence, behavioural discrimination latency.

**What you are not explaining.** Concentration invariance across four log units. Mixture segmentation in plumes. Learning of the readout. Write these down; a specification that quietly claims everything explains nothing.

---

## M1 — Build the substrate

**Deliverable: a simulator plus a benchmark odour ensemble, validated against at least five published summary statistics.**

Two levels, both worth having:

1. **Biophysical.** The Bazhenov, Stopfer, Rabinovich, Abarbanel, Sejnowski & Laurent (2001) AL model — Hodgkin–Huxley PNs and LNs with fast GABA$_A$ (synchronisation) and slow GABA$_B$-like (slow patterning) inhibition. Its central and correct claim is that *synchronisation and slow temporal patterning are produced by different mechanisms*, which MacLeod & Laurent (1996) established experimentally by dissociating them with picrotoxin.
2. **Reduced rate model.** $N_p=90$ PNs, $N_l=30$ LNs, current-based, with an explicit 20 Hz oscillatory drive and a slow inhibitory variable. This is what you will actually run 10,000 times.

**Validation targets** (each is a published, quotable number or curve): PN response transience and the fact that trajectories are more discriminable during the transient than at the fixed point (Mazor & Laurent, 2005); oscillatory phase-locking of a subset of PN spikes; loss of synchrony but preservation of slow patterning under picrotoxin (MacLeod & Laurent, 1996); KC sparseness and specificity (Perez-Orive et al., 2002); concentration-dependent but identity-preserving trajectory families (Stopfer, Jayaraman & Laurent, 2003).

**On data.** Be honest about this: the primary locust datasets are largely not public. Three viable routes, in order of preference: (i) collaborate for raw data — this project is worth an email; (ii) fit to *published summary statistics* rather than raw traces, which is scientifically weaker but perfectly respectable if you say so; (iii) use *Drosophila* analogues where data are public — Hallem & Carlson's (2006) ORN response matrix as a realistic input ensemble, and the hemibrain connectome (Scheffer et al., 2020) for AL and MB connectivity statistics — while stating clearly which claims transfer. Route (iii) plus route (ii) is enough to complete M0–M6.

---

## M2 — Formalise the competing hypotheses

**Deliverable: four model classes with a common interface `encode(u(t)) -> x(t)`, each with an explicit parameter count and prior.**

These are the candidates enumerated in `../part2-case-studies/C2-transient-synchrony.md`. Formalising them is the intellectual core of the project, because three of the four are usually stated too loosely to be wrong.

### H1 — Reservoir / liquid state
$$\tau\dot x = -x + \phi\big(g\,Wx + B u(t)\big),\qquad y = V\rho(x),$$
$W$ random with spectral radius $g$. The claim is that the AL is a generic high-dimensional nonlinear expansion with fading memory, and that separability is a property of the *ensemble* of such systems, not of any particular connectivity (Maass, Natschläger & Markram, 2002; Buonomano & Maass, 2009). Quantify with the separation property, $\mathrm{Sep}(u,v;t)=\|x_u(t)-x_v(t)\|$ against within-odour variability, and with the rank-based kernel-quality and generalisation measures of Legenstein & Maass (2007). Free parameters: $g,\tau$, input gain, noise — four.

**The signature that makes H1 falsifiable:** genericity. Rewire $W$ at random, preserving only its statistics, and performance is unchanged. No other hypothesis survives that manipulation.

### H2 — Dynamic decorrelation
$$\tau\dot x = -x + u(t) - W x, \qquad \dot W_{ij} \propto \eta\,x_ix_j\;(i\ne j),$$
whose fixed point drives $\langle xx^\top\rangle$ toward diagonal — i.e. the anti-Hebbian decorrelation network. Note this is the *same* mathematics as `../part2-case-studies/C9-adaptive-filters.md`: each PN subtracts a learned prediction of itself formed from its neighbours, and the residual is the innovation. At steady state $x=(I+W)^{-1}u$; during the transient, correlations decay with time constants set by the eigenvalues of $(I+W)/\tau$. The zebrafish evidence for exactly this is Friedrich & Laurent (2001). Free parameters: $\tau$, $\eta$, gain — three.

**Signature:** pairwise pattern correlation $\rho_{uv}(t)$ decreases monotonically over the transient, *most* for initially most-similar odours, and the time course of decorrelation quantitatively predicts the time course of discriminability with no additional parameters.

### H3 — Stable heteroclinic sequences (winnerless competition)
Generalised Lotka–Volterra on PN-group amplitudes,
$$\dot a_i = a_i\Big(\sigma_i(u) - \sum_j\rho_{ij}a_j\Big) + \epsilon\,\xi_i(t),$$
with asymmetric $\rho$. For suitable $\rho$ (May–Leonard condition: cyclic with $\alpha<1<\beta$, $\alpha+\beta>2$) the attractor is a *heteroclinic channel* — a chain of saddles, each with a one-dimensional unstable manifold, joined by connecting orbits. Odour identity is the channel; reproducibility comes from transverse stability; the trajectory is structurally the message (Rabinovich et al., 2001; Afraimovich, Zhigulin & Rabinovich, 2004). Free parameters: $\rho$ structure (constrained to two numbers $\alpha,\beta$ if cyclic), $\sigma(u)$ mapping, $\epsilon$.

**Signature — and this is the best discriminating prediction in the whole project:** the dwell time at saddle $k$ obeys

$$T_k \;\simeq\; \frac{1}{\lambda_k^{u}}\,\ln\frac{1}{\epsilon},$$

because escape from a saddle requires the noise-seeded perturbation $\epsilon$ to grow to $O(1)$ along the unstable direction at rate $\lambda^u_k$. **Logarithmic** in noise amplitude. Change the noise by two decades and dwell times change by an *additive* $\ln(100)/\lambda^u \approx 4.6/\lambda^u$, not by a factor. Nothing else on this list predicts that.

### H4 — Synchrony as gating
The AL's job is to make PN spikes *coincident* within an oscillation cycle so that KCs, as high-threshold coincidence detectors with a short effective integration window, respond selectively. Model KC input in a cycle as $n$ PN spikes with phase jitter $\sigma_j$, each contributing a kernel of width $\Delta$:
$$V(t) = g\sum_{i=1}^{n}\mathcal N(t;t_i,\Delta^2),\quad t_i\sim\mathcal N(\mu,\sigma_j^2) \;\Longrightarrow\; \mathbb E[V]_{\max} = \frac{n\,g}{\sqrt{2\pi(\Delta^2+\sigma_j^2)}} .$$

The peak depolarisation falls as $(\Delta^2+\sigma_j^2)^{-1/2}$ **at constant spike count**. With threshold $\theta$ near the peak, KC output collapses under desynchronisation even though every PN fires exactly as before. This is precisely the picrotoxin manipulation (MacLeod & Laurent, 1996), its downstream consequence (MacLeod, Bäcker & Laurent, 1998), and its behavioural consequence in bees (Stopfer, Bhagavan, Smith & Laurent, 1997). Free parameters: $\Delta,\theta,g,n$ — but $n$ and $g$ are constrained by Jortner et al. (2007) and $\Delta$ by KC membrane properties, so effectively one or two.

**A position, since this note is supposed to have one.** H4 is not a rival to H1–H3; it is a statement about the *readout*, and it can be true simultaneously with any of them. The genuine three-way contest is H1 vs. H2 vs. H3 about what the AL dynamics are *for*. My prior, stated so it can be held against me: H2 and H4 are both substantially true and mechanistically established, H3 is an elegant description whose distinctive predictions have never been tested properly, and H1 is nearly unfalsifiable as usually stated and should be forced to make the rewiring prediction or withdraw.

---

## M3 — Fit, and account for every parameter

**Deliverable: a table with one row per hypothesis giving $k_{\text{fitted here}}$, $k_{\text{fixed from literature}}$, $n_{\text{constraints}}$, and the resulting $\tfrac{k}{2}\log n$ penalty (see `./S2-degeneracy-and-limits.md` Ex. 5).**

Fit all four model classes to the *same* M1 validation targets. Report which targets each can and cannot match. Then compute prior predictive spread: sample parameters from your priors and plot the range of predicted curves. If a model can produce every curve, its fit is worthless, and you should say so at this stage rather than after you have grown attached to it.

---

## M4 — Derive discriminating predictions

**Deliverable: a prediction matrix with effect sizes and a power analysis.**

| Manipulation / measurement | H1 reservoir | H2 decorrelation | H3 heteroclinic | H4 gating |
|---|---|---|---|---|
| Desynchronise at constant PN rates (picrotoxin) | mild loss | mild loss | loss only if channel disrupted | **catastrophic KC loss** |
| Random rewiring of AL recurrence, statistics matched | **no change** | loss of decorrelation | loss of channel structure | no change (if PN synchrony survives) |
| Noise amplitude $\times100$ | dwell/timing scales as power law or not at all | correlation time constants unchanged | **dwell times shift by additive $\ln 100/\lambda^u$** | jitter ↑, KC output ↓ |
| Pairwise correlation over the transient | no required trend | **monotone decrease, largest for similar pairs** | decrease as trajectories separate, not correlated with similarity | n/a |
| Stimulus intensity ↑ | new trajectory | faster decorrelation | **same channel, dilated/compressed timing** | more spikes, more synchrony |
| Odour B presented during A's transient | strong path dependence (fading memory) | additive in the decorrelated basis | return to A's channel or switch to a new one; hysteresis | n/a |
| Trial-to-trial dwell-time distribution | ~Gaussian | ~Gaussian | **skewed, Gumbel-like (log of an exponential-tailed escape)** | n/a |

Rows 2, 3 and 7 are the money. Row 3 in particular is a *qualitative* signature — logarithmic versus power-law — which survives almost any amount of parameter uncertainty, and is therefore the kind of prediction that `./S2-degeneracy-and-limits.md` says you should be hunting for. Note that Broome, Jayaraman & Laurent (2006) already performed something close to row 6, and Stopfer, Jayaraman & Laurent (2003) close to row 5; go and check whether the published results discriminate, before you design anything new.

For each row, compute the effect size in your M1 simulator and the number of trials/neurons needed for 80% power. A prediction without a power analysis is not an experiment.

---

## M5 — Simulate, including model recovery

**Deliverable: a model-recovery confusion matrix.**

This is the milestone people skip, and it is the one that `./S2-degeneracy-and-limits.md` §2 says is mandatory. Generate synthetic data from each of H1–H4 (with realistic recording constraints: 50–200 simultaneously recorded PNs, 20–50 trials, spike sorting noise), run your full analysis pipeline blind, and build the confusion matrix. **If you cannot recover the generating model from data you generated, you cannot adjudicate on real data, and every subsequent claim is void.**

Two starter programs.

**(a) The heteroclinic dwell-time law.**

```python
import numpy as np
rng = np.random.default_rng(0)

n, alpha, beta = 3, 0.8, 1.4                      # May-Leonard: a<1<b, a+b>2
rho = np.eye(n) + alpha*np.roll(np.eye(n), 1, 1) + beta*np.roll(np.eye(n), -1, 1)
sigma = np.ones(n)

def run(eps, T=200000, dt=0.005):
    a = np.array([0.6, 0.2, 0.05]); lead = np.empty(T, int)
    for t in range(T):
        a += dt*a*(sigma - rho @ a) + eps*np.sqrt(dt)*rng.standard_normal(n)
        a = np.clip(a, 1e-14, None); lead[t] = a.argmax()
    sw = np.flatnonzero(np.diff(lead) != 0)
    return np.diff(sw)*dt                          # dwell times, in order

for eps in [1e-8, 1e-6, 1e-4, 1e-2]:
    d = run(eps)
    print("eps=%.0e  dwell times: %s" % (eps, np.round(d[2:8], 2)))
```

**What to look for.** Successive dwell times *grow* (the cycle is attracting) and then plateau at a level set by the noise. Plot plateau dwell time against $\log_{10}\epsilon$: a **straight line of slope $-\ln(10)/\lambda^u$**, where $\lambda^u$ is the unstable eigenvalue of the saddle. Compute $\lambda^u$ independently from the Jacobian at the single-species equilibrium $a=(1,0,0)$ — for these parameters the eigenvalues are $\{-1,\,-\beta+1,\,1-\alpha\} = \{-1,\,-0.4,\,+0.2\}$, so $\lambda^u=0.2$ and the predicted slope is $2.303/0.2 \approx 11.5$ time units per decade of noise. The simulation gives plateaux of roughly $90,\,66,\,46,\,24$ at $\epsilon = 10^{-8},10^{-6},10^{-4},10^{-2}$ — about $11$ per decade. **Parameter-free, and it works.** Six decades of noise change dwell times by less than a factor of four. (At $\epsilon=10^{-2}$ you will also see spurious near-zero dwells from flicker between near-equal amplitudes; filter them.)

Then run the same measurement on H1: a reservoir's transient timescales are set by $\tau$ and $g$ and are independent of noise amplitude. The two hypotheses are trivially distinguishable *if you vary the noise* — which, remarkably, nobody has done in a locust. If you do only one experiment from this note, do this one.

**(b) Synchrony gating of a Kenyon cell.**

```python
def kc_dprime(n_pn=50, g=1.0, delta=5.0, jitter=np.linspace(0.1, 40, 12), n_trial=4000):
    out = []
    for sj in jitter:
        peak = []
        for _ in range(2):                          # odour A (n_pn) vs odour B (0.7*n_pn)
            k = int(n_pn if not peak else 0.7*n_pn)
            ts = rng.normal(0.0, sj, size=(n_trial, k))
            grid = np.linspace(-60, 60, 241)
            v = np.exp(-0.5*((grid[None, :, None] - ts[:, None, :])/delta)**2).sum(2)
            peak.append(g*v.max(1)/np.sqrt(2*np.pi*delta**2))
        a, b = peak
        out.append((a.mean()-b.mean())/np.sqrt(0.5*(a.var()+b.var())))
    return jitter, np.array(out)
```

**What to look for.** The mean peak depolarisation tracks the closed form $n g/\sqrt{2\pi(\Delta^2+\sigma_j^2)}$ to within a percent — check it directly, it is the whole of H4 in one line. Consequently $d'$ for the 50-PN versus 35-PN odours falls steeply with jitter (roughly $8.8 \to 1.7$ as $\sigma_j$ goes from 4 ms to 40 ms with $\Delta=5$ ms), while the total spike count is *identical* in every condition. Ignore the $\sigma_j\to0$ point: with no jitter the peak is deterministic and $d'$ diverges.

This is the quantitative content of "synchrony carries information that rate does not", and it converts a slogan into a curve you can lay against KC recordings under picrotoxin. Note that it also gives you the free parameter you most want: fitting the measured $d'(\sigma_j)$ recovers $\Delta$, the KC integration window, from the *psychophysics* rather than from a patch pipette — and the two estimates had better agree.

---

## M6 — Pre-register the analysis for real data

**Deliverable: a pre-registration document plus analysis code that already runs end-to-end on M5's synthetic data.**

Specify in advance: the decoding analysis (cross-validated, with the readout constrained to be random and sparse per M0 — not an optimised linear decoder, which would answer a different question); the correlation time-course analysis with its null; the dwell-time estimator and how you will identify saddles from data (e.g. via low-velocity epochs in a smoothed state-space trajectory, with the velocity threshold fixed in advance); the model-comparison statistic; and the falsification conditions, one per hypothesis, written as sentences of the form *"If X, I will abandon H_i."*

Write those sentences before you look at data. If you cannot write one for a hypothesis, that hypothesis is not a scientific claim and should be dropped from the specification.

---

## M7 — Write the specification

**Deliverable: a document that scores ≥20/24 on the `./S2-degeneracy-and-limits.md` checklist, with the scoring shown.**

Structure it exactly as the case studies in `../part2-case-studies/`: phenomenon, conversion, model, predictions, open problems. It should contain (i) the level-1 statement from M0; (ii) the winning algorithm stated as mathematics compact enough that someone could implement it from the text; (iii) the mapping onto AL/MB anatomy with the specific mechanistic commitments (which cell class does what); (iv) at least one zero-parameter quantitative prediction; (v) the causal experiment that would kill it; (vi) the compression argument from `./S2-degeneracy-and-limits.md` §3.1, with numbers.

A note on the likely outcome, which you should be prepared for: the honest answer may be a *composite* — H2 for what the AL dynamics achieve, H4 for how the MB reads it, H3 as a description of the geometry rather than a mechanism, H1 rejected as unfalsifiable or accepted as a null model that the others must beat. A composite specification is a real result provided each component is separately falsifiable. What is not a result is a specification that accommodates every finding.

---

## Rubric: what "done" looks like

| Dimension | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| Level-1 statement | absent or circular | prose only | objective written, constraints listed | objective written, constraints sourced with numbers, non-goals stated |
| Formalisation | verbal hypotheses | one model implemented | all four implemented, shared interface | all four with priors, parameter counts, and prior predictive plots |
| Discrimination | "consistent with" | qualitative differences | prediction matrix with effect sizes | matrix + power analysis + at least one qualitative (shape-level) discriminator |
| Model recovery | not attempted | attempted informally | confusion matrix on ideal data | confusion matrix under realistic recording constraints |
| Contact with data | none | qualitative comparison to figures | quantitative fit to published statistics | pre-registered analysis run on raw data |
| Compression | not addressed | asserted | free parameters counted | parameter→rule relocation demonstrated with description lengths |
| Falsifiability | none stated | one vague condition | one per hypothesis | one per hypothesis, pre-registered, at least one already tested |

Total ≥17/21 with no dimension at 0 or 1: you have a level-2 specification. Below that, you have a proposal — which is fine, but call it one.

---

## Reading path

- **Laurent & Davidowitz (1994)**, *Encoding of olfactory information with oscillating neural assemblies* — read it for: the founding observation, and the framing of assemblies as dynamic.
- **Wehr & Laurent (1996)**, *Odour encoding by temporal sequences of firing in oscillating neural assemblies* — read it for: the sequences that H3 is a theory of.
- **MacLeod & Laurent (1996)**, *Distinct mechanisms for synchronization and temporal patterning of odor-encoding neural assemblies* — read it for: the dissociation that makes H4 testable; the single most important experiment for this project.
- **Stopfer, Bhagavan, Smith & Laurent (1997)**, *Impaired odour discrimination on desynchronization of odour-encoding neural assemblies* — read it for: the behavioural consequence, i.e. the level-1 relevance of synchrony.
- **MacLeod, Bäcker & Laurent (1998)**, *Who reads temporal information contained across synchronized and oscillatory spike trains?* — read it for: the downstream readout, and the correct question in its title.
- **Rabinovich et al. (2001)**, *Dynamical encoding by networks of competing neuron groups: winnerless competition* — read it for: H3, stated properly, with the dynamical-systems machinery.
- **Bazhenov, Stopfer, Rabinovich, Abarbanel, Sejnowski & Laurent (2001)** (two companion papers) — read them for: the biophysical substrate you will build in M1.
- **Friedrich & Laurent (2001)**, *Dynamic optimization of odor representations by slow temporal patterning of mitral cell activity* — read it for: H2, in a different animal, with the decorrelation measured.
- **Perez-Orive, Mazor, Wilson, Cassenaer, Stopfer & Laurent (2002)**, *Oscillations and sparsening of odor representations in the mushroom body* — read it for: the coincidence-detection window and the sparse readout — the empirical heart of H4.
- **Stopfer, Jayaraman & Laurent (2003)**, *Intensity versus identity coding in an olfactory system* — read it for: the concentration manipulation in row 5 of the M4 matrix.
- **Mazor & Laurent (2005)**, *Transient dynamics versus fixed points in odor representations by locust antennal lobe projection neurons* — read it for: the demonstration that the transient, not the steady state, carries the information. This is the paper the whole capstone is about.
- **Broome, Jayaraman & Laurent (2006)**, *Encoding and decoding of overlapping odor sequences* — read it for: the history-dependence manipulation in row 6.
- **Cassenaer & Laurent (2007)**, *Hebbian STDP in mushroom body efferent neurons in the insect brain* — read it for: the plasticity rule at the readout, which any level-2 spec must be compatible with.
- **Jortner, Farivar & Laurent (2007)**, *A simple connectivity scheme for sparse coding in an olfactory system* — read it for: the PN→KC connectivity numbers that fix $A$ in M0.
- **Legenstein & Maass (2007)**, *Edge of chaos and prediction of computational performance for neural circuit models* — read it for: how to measure kernel quality and generalisation, i.e. how to make H1 quantitative.
- **Papadopoulou, Cassenaer, Nowotny & Laurent (2011)**, *Normalization for sparse encoding of odors by a wide-field interneuron* — read it for: the GGN, and the normalisation term in your readout model.
- **Cassenaer & Laurent (2012)**, *Conditional modulation of spike-timing-dependent plasticity for olfactory learning* — read it for: phase-gated, neuromodulator-conditional plasticity — the strongest evidence that oscillation phase is used, not epiphenomenal.
- **Tootoonian & Lengyel (2014)**, *A Dual Algorithm for Olfactory Computation in the Locust Brain* — read it for: a worked example of exactly the kind of level-2 specification M7 asks for, in this system.

Also revisit, in this order: `../part2-case-studies/C2-transient-synchrony.md` (the phenomenology), `../part2-case-studies/C9-adaptive-filters.md` (H2 is an adaptive filter), `../part2-case-studies/C10-evidence-accumulation.md` (the latency term in M0's objective is an accumulation problem), `./S1-rnn-as-model-organism.md` (train an RNN on the same level-1 problem and reverse-engineer it — this is an excellent, and cheap, fifth hypothesis), and `./S2-degeneracy-and-limits.md` (the checklist you must pass).

---

## Milestone summary

| M | Deliverable | Done when |
|---|---|---|
| M0 | Level-1 problem statement | Objective, constraints with numbers, and non-goals written down |
| M1 | Simulator + benchmark ensemble | Reproduces ≥5 published summary statistics |
| M2 | Four formalised hypotheses | Common interface; parameter counts and priors stated |
| M3 | Fits + parameter accounting | Table of $k$, $n$, penalties, prior predictive spread |
| M4 | Prediction matrix | Every cell has an effect size and required $n$ |
| M5 | Model-recovery confusion matrix | Off-diagonal mass < 20% under realistic constraints |
| M6 | Pre-registration + pipeline | Falsification sentence written for each hypothesis |
| M7 | The specification | ≥20/24 on the S2 checklist, scoring shown |

The point of the capstone is not to be right about the locust. It is to find out, by doing it, exactly how much work a defensible level-2 claim costs — and then to notice how few published claims have paid it.
