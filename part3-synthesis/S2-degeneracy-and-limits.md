---
title: "S2 · Degeneracy & limits"
parent: Synthesis
nav_order: 2
---

# S2 — When there is no clean algorithm
{: .no_toc }

> **The worry.** Nervous systems are massively degenerate: many parameter sets produce identical behaviour. Standard analyses applied to a system whose algorithm we know exactly (a microprocessor) fail to recover it. And normative stories are cheap — for almost any behaviour, *some* objective makes it optimal.
> **The claim of this note.** Degeneracy at the parameter level is not an obstacle to level-2 claims; it is an argument *for* them. But the enterprise only works if "level-2 claim" is given teeth. Here are the teeth.
> **Structures thread.** Identifiability and Fisher information; convexity and level sets; minimum description length; convex cones and inverse problems. See `../structures/README.md`.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. Degeneracy, stated properly

Let $\Theta\subseteq\mathbb{R}^p$ be the parameter space of a circuit model (maximal conductances, synaptic strengths, time constants) and let $\Phi:\Theta\to\mathcal B$ map parameters to behaviour, where $\mathcal B$ is whatever we can measure. Marder's crustacean stomatogastric ganglion is the definitive case: Prinz, Bucher & Marder (2004) ground through ~20 million versions of a three-cell pyloric model and found on the order of $10^5$–$10^6$ that produce pyloric-like output, spread across the full range of every maximal conductance, with individual conductances varying several-fold between good models. Biology agrees: Schulz, Goaillard & Marder measured mRNA levels for conductances across identified neurons and found large cell-to-cell variability with *correlations* between conductances — the parameters wander over a manifold rather than sitting at a point.

The formal content is about the **fibres** of $\Phi$. Define $G_b = \Phi^{-1}(b)$, the set of parameters producing behaviour $b$. The empirical fact is that $\dim G_b$ is large — most directions in $\Theta$ are behaviourally neutral. Three consequences follow immediately, and each has been observed:

**(i) Parameter identification is hopeless, and no amount of data fixes it.** This is not measurement noise; it is structural non-identifiability. The Fisher information matrix
$$\mathcal I_{ij} = \sum_k \frac{1}{\sigma_k^2}\frac{\partial y_k}{\partial\theta_i}\frac{\partial y_k}{\partial\theta_j}$$
is rank-deficient (structural) or has eigenvalues spanning many decades (practical). Gutenkunst et al. (2007) showed the latter is generic across systems-biology models: eigenvalues roughly log-uniform over six or more decades — "sloppiness". Parameter uncertainty along eigen-direction $i$ scales as $1/\sqrt{\lambda_i}$, so the sloppy directions are unconstrained by *any* feasible experiment. Transtrum et al. (2015) argue this is what it looks like when a microscopic model has a low-dimensional effective theory — which is, note, exactly the level-2 claim in different words.

**(ii) Averaging over good models is invalid.** Golowasch, Goldman, Abbott & Marder (2002) made this concrete: average the conductances of several models that each produce correct bursting and the average model does not burst. The reason is trivial once stated — $G_b$ is not convex, so $\theta_1,\theta_2\in G_b \not\Rightarrow \tfrac12(\theta_1+\theta_2)\in G_b$ (Exercise 1). Every "average model of the neuron" in the literature is potentially a model of nothing.

**(iii) Therefore the interesting invariant lives downstream of $\Phi$.** If $\theta$ is not recoverable but $b$ is, then any well-posed claim must be a claim about $\mathcal B$, or about a quotient $\Theta/\!\sim$ where $\theta\sim\theta'$ iff $\Phi(\theta)=\Phi(\theta')$. **The algorithm is a function on the image of $\Phi$, not on its domain.** Degeneracy is thus precisely the condition under which level-2 description is the *right* level of description: it says there is a coarse-graining that loses nothing behaviourally relevant. O'Leary et al. (2014) went further and showed how homeostatic regulation *generates* the degenerate manifold — the correlations are the trace of a controller that steers parameters onto $G_b$ from anywhere. The circuit is engineered to make its own parameters uninteresting.

That is the optimistic reading, and I think it is right. Now the pessimistic case.

## 2. The microprocessor

Jonas & Kording (2017) took a MOS 6502 — a system whose algorithm is documented to the bit — ran three games on a simulator, and applied the standard toolkit: lesions (knock out each transistor, see which games break), tuning curves, spike-triggered-style analyses, local field potentials, dimensionality reduction, Granger causality. The analyses produce beautiful, publishable, wrong results. Single-transistor lesions identify transistors "necessary for Donkey Kong but not Space Invaders" — which is true and completely uninformative about the algorithm. Dimensionality reduction finds low-dimensional structure that reflects the clock, not the computation. Granger causality recovers connections that do not exist.

Two objections are usually raised, and both are half right. *"A processor is not brain-like — no redundancy, no robustness."* True, and it means the failure of lesioning is unsurprising there and might be less bad in brains; but the failure of dimensionality reduction and of causal inference has nothing to do with redundancy. *"They did not ask the questions a neuroscientist would ask."* Also partly true — but the deeper point stands, and it is a methodological one: **we have essentially no analysis method that has been validated against a known ground-truth algorithm.** We validate methods on simulated data generated by the models the methods presuppose. That is circular, and the 6502 is the cheapest available escape.

The constructive response is not despair; it is benchmarking. Any new analysis method for extracting computation from recordings should be run on (a) the 6502, (b) trained RNNs with reverse-engineered mechanisms (see `./S1-rnn-as-model-organism.md`), and (c) a biophysical circuit model whose algorithm you designed. If it fails on all three, its success on your data is not evidence.

## 3. Fitting is not explaining

A model fits when it assigns high probability to held-out data. It explains when it *compresses*, *transfers*, and *constrains counterfactuals*. These come apart constantly.

Roberts & Pashler (2000) put it most sharply: a good fit is persuasive only in proportion to how much of the possible data space the model *excluded in advance*. A model that could have fitted anything, fitting your data, is worth nothing. Formalise it as prior predictive risk. Let $M$ be your model with parameters $\psi$ and prior $p(\psi)$; the prior predictive is $p(y\mid M)=\int p(y\mid\psi)p(\psi)\,d\psi$. Evidence for $M$ over $M'$ is $\log p(y\mid M) - \log p(y\mid M')$, and the first term is small whenever $M$ spreads its predictive mass widely. Fitting after the fact, with parameters free, systematically inflates apparent support because it replaces $p(y\mid M)$ with $\max_\psi p(y\mid\psi)$. Pitt & Myung (2002) make the same point in terms of functional-form complexity: a model can be more flexible than another with the same number of parameters.

The practical version of this that I would enforce: **count the parameters you spent, and where you spent them.** A conversion that fits $k$ parameters on dataset A and then predicts dataset B with *zero* refitted parameters is worth more than one that refits on every dataset, and the difference is quantifiable, not rhetorical.

### 3.1 Compression, made precise

The intuition "the algorithm should be shorter than the circuit description" is right and can be made rigorous with the two-part MDL code. For a model class $M$,
$$\mathrm{DL}(y) = \mathrm{DL}(M) + \mathrm{DL}(y\mid M),$$
and we prefer the model minimising the total. Applied to a circuit→algorithm conversion, the comparison is:

$$\underbrace{\mathrm{DL}(\text{algorithm}) + \mathrm{DL}(\text{data}\mid\text{algorithm})}_{\text{level-2 account}} \quad\text{vs.}\quad \underbrace{\mathrm{DL}(\text{parameters}) + \mathrm{DL}(\text{data}\mid\text{parameters})}_{\text{circuit account}} .$$

The subtle and important point is *where the compression comes from*. It is not that the algorithm has fewer symbols than the circuit has synapses — that comparison is unfair, because a bare algorithm cannot predict a specific negative image while the synaptic weights can. The compression comes from the algorithm **deriving the parameters from something else you already had to specify**: in `../part2-case-studies/C9-adaptive-filters.md`, the $\sim10^4$ parallel-fibre weights are not free — they are $w^*=-C^{-1}q$, determined by the input statistics. The algorithm converts $10^4$ free parameters into two ($\eta,\beta$) plus a description of the environment you needed anyway. *That* is the compression, and it is exactly what "explanation" buys over "fit" (Exercise 4).

### 3.2 Just-so stories, and the geometry of underdetermination

Normative explanation has a specific, provable underdetermination problem. Bowers & Davis (2012) complained that Bayesian models in psychology can be made to fit anything by adjusting the prior and the cost function; the complaint is often dismissed as unfriendly, but there is a theorem behind it.

Ng & Russell (2000) characterised inverse reinforcement learning: given a finite MDP with transitions $P_a$ and a policy $\pi$ that takes action $a^*$ everywhere, the set of reward vectors $R$ for which $\pi$ is optimal is
$$\big\{R : (P_{a^*}-P_a)(I-\gamma P_{a^*})^{-1}R \;\succeq\; 0 \;\;\forall a\big\},$$
a **convex cone** — closed under positive scaling and addition, and containing $R=0$, for which *every* policy is optimal. So observed optimal behaviour never identifies the objective. Extend the argument beyond RL: for any behaviour and any sufficiently rich family of (objective, prior, constraint) triples, a member of the family rationalises it. "The brain is optimal for $X$" is therefore not, by itself, a claim; it becomes one only when the family is restricted in advance.

Three restrictions actually work, and you should insist on at least one:
- **The objective was fixed independently** — by ecology, by physics, by a task the experimenter designed with a known optimum (this is why `../part2-case-studies/C10-evidence-accumulation.md` is strong: Wald's objective is not negotiable).
- **The same objective, with the same parameters, explains behaviour across multiple environments.** Behavioural diversity shrinks the cone. Quantify it: report how much.
- **The objective predicts the *mechanism*, not just the behaviour** — as TD does when it demands a synaptic eligibility trace with a specific sign and time constant.

## 4. What a good level-2 claim looks like: the checklist

Apply this to any claimed circuit→algorithm conversion, including your own. Score each item 0 (absent), 1 (partial), 2 (met). My rule of thumb: below 14/24 it is a story; above 20/24 it is a result.

| # | Criterion | The test | Automatic fail |
|---|---|---|---|
| 1 | **Level-1 statement exists** | Can you write the computational problem — inputs, outputs, objective, constraints — without mentioning neurons? | The "problem" is defined by the circuit's own behaviour |
| 2 | **The algorithm is stated as mathematics** | A reader could implement it from your text alone, with no figures | Only prose and a box-and-arrow diagram |
| 3 | **Parameter accounting** | Number of free parameters, and which dataset each was fitted on, stated explicitly | Refitted per dataset, not reported |
| 4 | **Zero-parameter prediction** | At least one prediction with *no* parameters fitted to the predicted data | Every quantitative claim is a fit |
| 5 | **Prior predictive risk** | What data would the model have excluded a priori? Show the range | "Consistent with" |
| 6 | **Falsification attempt made** | The authors ran an experiment that could have killed it, and report the result | Only confirmatory analyses |
| 7 | **Compression** | Free parameters replaced by a rule that generates them (§3.1); state the ratio | Algorithm has as many knobs as the circuit has synapses |
| 8 | **Degeneracy-robust** | The claim is invariant across the fibre $G_b$ — it does not depend on a particular parameter set | Conclusions change if you re-fit and get an equally good model |
| 9 | **Mechanistic commitment** | The algorithm predicts something about implementation (a plasticity rule sign, a time constant, a connectivity motif) | Purely behavioural |
| 10 | **Causal test** | Perturbation predicted *in advance* by the algorithm, then performed | Correlational only, or perturbation explained post hoc |
| 11 | **Transfer** | The same algorithm, with parameters set independently, applies to another system or task | One system, one task |
| 12 | **Method validated on ground truth** | The analysis pipeline was shown to recover a known algorithm from simulated data (§2) | Pipeline never tested where the answer is known |

A worked application: `../part2-case-studies/C9-adaptive-filters.md` scores near the ceiling — 1,2,4,6,7,9,10 all clearly met (the learning rule's *sign* was predicted and then measured), 11 met via cerebellum and cortical corollary discharge, with 12 the weakest. `../part2-case-studies/C7-dopamine-and-td-learning.md` is comparable, but item 8 is now under pressure from distributional RL and heterogeneity, and item 1 is contested because the state representation is a free choice. `../part2-case-studies/C10-evidence-accumulation.md` aces 1–7 and 11 but fails 10 badly: the causal test was performed and came back negative for LIP. That is not a criticism of the algorithm; it is the checklist doing its job, localising the failure to the circuit assignment rather than the computation.

## 5. Where there really is no clean algorithm

Be prepared for the answer to be "no". Three honest categories:

- **Systems whose behaviour is not a computation.** A homeostat maintaining a set point is well described by control theory but has no interesting level-2 content beyond "it regulates". Do not manufacture an algorithm.
- **Systems where the level-2 description is not shorter.** If the best compression of the circuit's behaviour is the circuit, the circuit *is* the explanation. This is the honest reading of some peripheral transduction and of much of development.
- **Systems that are mid-evolution or mid-learning.** A partly-optimised system implements neither the old algorithm nor the new one, and any single clean description will misfit. This is a real risk for overtrained laboratory animals, and it is the most likely explanation for phenomena like collapsing bounds appearing in monkeys but not humans.

Krakauer, Ghazanfar, Gomez-Marin, MacIver & Poeppel (2017) make the complementary point from the other end: without a precise behavioural/ethological characterisation of what the animal is doing, mechanism-level data cannot be interpreted at all. Level 1 is not optional throat-clearing. It is the load-bearing part, and it is the part most papers skip.

---

## 6. Exercises

**★ 1. Why averaging good models fails.**
Give an explicit two-parameter example where $\theta_1,\theta_2$ both produce the target behaviour but their average does not, and state the general condition.

<details markdown="1"><summary>Solution</summary>

Let the behaviour be "oscillates", and let the model oscillate iff $g_1g_2 \ge 1$ with $g_i>0$ (a caricature of a burster needing sufficient product of an inward and an outward conductance, with an upper bound from another mechanism, say $g_1+g_2\le 5$). Take $\theta_1=(4,1)$ and $\theta_2=(1,4)$: products are 4 and 4, both oscillate. Average is $(2.5,2.5)$, product $6.25$ — fine here, so make the good set an annulus instead: oscillates iff $1\le g_1g_2\le 5$ and $g_1+g_2\le 6$. Then $(5,1)$ and $(1,5)$ have product 5 (good); average $(3,3)$ has product 9 (bad).

General condition: averaging is safe iff $G_b$ is convex. Nothing guarantees this — $G_b$ is a level set of a nonlinear map and is typically a curved, possibly disconnected manifold. Prinz et al.'s pyloric database is visibly non-convex and has multiple components.

Corollary that people miss: the same argument condemns averaging *across animals* whenever the quantity averaged is a parameter rather than an observable. Average firing rates: fine. Average fitted conductances: meaningless. Average connectivity matrices across subjects: meaningless unless you have argued convexity.
</details>

**★★ 2. Sloppiness in the smallest possible model.**
A neuron at rest has $V=\dfrac{g_1E_1+g_2E_2+g_LE_L}{g_1+g_2+g_L}$ and input resistance $R=(g_1+g_2+g_L)^{-1}$. (a) With only $V$ measured, compute the Fisher information matrix in $(g_1,g_2)$ and identify the stiff and sloppy directions. (b) With both $V$ and $R$ measured, compute the Jacobian determinant and the condition number. Interpret.

<details markdown="1"><summary>Solution</summary>

(a) $\partial V/\partial g_i = R\,(E_i - V)$ (differentiate the quotient and use $V=R\sum g_jE_j$). With Gaussian noise of variance $\sigma^2$, $\mathcal I = \sigma^{-2}\nabla V\,\nabla V^\top$, which is **rank 1**. Eigenvalues: $\|\nabla V\|^2/\sigma^2$ and $0$. The stiff direction is $\nabla V \propto (E_1-V,\,E_2-V)$; the sloppy (exactly unidentifiable) direction is orthogonal, $(-(E_2-V),\,E_1-V)$.

Interpretation: the only identifiable combination is $g_1(E_1-V)+g_2(E_2-V)$ — the net ionic current, which at rest must balance the leak. The *split* between the two conductances is invisible, exactly, with infinite data. This is structural non-identifiability in its purest form, and it is the two-parameter cartoon of the entire STG result.

(b) Jacobian $M=\begin{pmatrix} R(E_1-V) & R(E_2-V)\\ -R^2 & -R^2\end{pmatrix}$, so $\det M = -R^3\,(E_1-E_2)$.

Both parameters are now identifiable — but the determinant vanishes as $E_1\to E_2$, and the condition number grows like $|E_1-E_2|^{-1}$. Parameter uncertainty scales as $\sigma/|\det M|$. So two potassium conductances with reversal potentials 5 mV apart are, for practical purposes, one conductance; two with reversals 100 mV apart are separable.

The general lesson worth carrying: **degeneracy is not binary.** It is the eigenvalue spectrum of $\mathcal I$, and the right report is that spectrum, not a claim that parameters "can" or "cannot" be recovered. Report $\lambda_i$ and the associated parameter combinations; that turns hand-waving about degeneracy into a measurement.
</details>

**★★ 3. The reward cone.**
Derive the Ng–Russell characterisation for a finite MDP, verify $R=0$ is in the set, and state the minimal additional assumption that restores identifiability up to scale.

<details markdown="1"><summary>Solution</summary>

Let $\pi(s)=a^*$ for all $s$, and let $P_a$ be the transition matrix under action $a$ taken everywhere. The value of $\pi$ is $V^\pi = (I-\gamma P_{a^*})^{-1}R$. Optimality requires that deviating for one step and then following $\pi$ is not better:
$$R + \gamma P_{a^*}V^\pi \;\succeq\; R+\gamma P_a V^\pi \quad\forall a \;\Longleftrightarrow\; (P_{a^*}-P_a)(I-\gamma P_{a^*})^{-1}R\succeq 0 .$$
(The one-step deviation principle makes this sufficient as well as necessary for the greedy policy w.r.t. $V^\pi$ to be $\pi$.)

The constraint set is defined by homogeneous linear inequalities in $R$, hence a convex cone: closed under $R\mapsto cR$ for $c>0$ and under addition. $R=0$ satisfies every inequality with equality — every policy is optimal when nothing matters — so the cone always contains the degenerate solution and has non-empty interior only when the policy is strictly optimal somewhere.

Restoring identifiability: you need constraints that cut the cone down to a ray. Options that work: (i) observe optimal behaviour in *multiple* MDPs sharing $R$ but with different dynamics — each adds inequalities, and enough of them pin the ray; (ii) observe a *stochastic* policy generated by a known decision rule (e.g. softmax with known temperature), which converts inequalities into equations; (iii) impose a regulariser (max-margin, sparsity) — but note this identifies $R$ by assumption, not by data, and the assumption is now the explanation.

Moral for neuroscience: "the animal is optimising $X$" needs (i) or (ii). Most normative papers rely on (iii) without saying so.
</details>

**★★ 4. Compute a compression ratio.**
For the ELL adaptive filter of `../part2-case-studies/C9-adaptive-filters.md`, estimate both sides of the MDL comparison in §3.1 and say precisely where the compression comes from.

<details markdown="1"><summary>Solution</summary>

*Circuit account.* To predict the medium-ganglion cell's response to a command you need the negative image, i.e. $K$ parallel-fibre weights. Take $K=2\times10^4$ granule cells at ~12 bits of precision (anything finer is below the noise floor): $\mathrm{DL}\approx 2.4\times10^5$ bits, plus connectivity. Residual data cost: small, because the weights predict the response well.

*Algorithm account.* $\mathrm{DL}(\text{algorithm})$: the equations $E=R+w^\top p$, $\Delta w=-\eta E p + \beta(w_0-w)$, plus two real constants — call it $10^2$–$10^3$ bits generously. But the algorithm alone does not predict the specific negative image; it predicts it *given the input statistics* $(C,q)$. Two cases:
- If you must transmit $C$ and $q$ as free parameters, you have $O(K^2)$ numbers and the algorithm is a *loss*.
- If $C$ and $q$ are determined by the animal's environment and body — which you had to describe anyway to state the problem — they cost nothing extra, and the algorithm account costs $\sim10^3$ bits against $\sim10^5$. Ratio $\sim10^{-2}$.

So: **compression is conditional on the environment description being independently required.** This is the precise statement of why level-2 explanations are explanatory. They do not delete information; they relocate it from "arbitrary facts about this animal's synapses" to "facts about the world that any account must include". Any claimed conversion that cannot pass this relocation test — that needs as many free numbers as the circuit has connections — has not explained anything, and you can now say so quantitatively.

Sanity check against the same computation for a lookup-table account ("store the response to each of $n$ commands"): $\mathrm{DL}$ grows linearly in $n$, whereas the algorithm's is constant. Growth rate of description length with dataset size is the cleanest single discriminator between fitting and explaining.
</details>

**★★★ 5. How much is a good fit worth?**
Take a model with $k$ free parameters fitted to $n$ data points with residual s.d. $\hat\sigma$, versus a model making a zero-parameter prediction with the same residual. Quantify the evidential difference, and state the practical rule.

<details markdown="1"><summary>Solution</summary>

Use the BIC/Laplace approximation to the marginal likelihood: $\log p(y\mid M) \approx \log p(y\mid\hat\psi) - \tfrac{k}{2}\log n$. Two models achieving the same maximised likelihood but with $k$ and $0$ parameters differ by $\tfrac{k}{2}\log n$ nats in favour of the parameter-free one. With $n=1000$ and $k=4$, that is $\approx 13.8$ nats, a Bayes factor of $\sim10^{6}$. Fitting four parameters costs you six orders of magnitude of evidence; you should not treat "it fits" as though it were free.

More carefully (and this is the Roberts–Pashler point in Bayesian dress), the penalty is not really about counting parameters but about the *prior predictive spread*: $\log p(y\mid M) = \log p(y\mid\hat\psi) - \log\!\big(\text{Occam factor}\big)$, where the Occam factor is the ratio of prior to posterior parameter volume. A model with four parameters that are tightly constrained a priori (by physiology, say) pays much less than one with four unconstrained parameters. This is why "we constrained $\tau$ to the measured membrane time constant" is a substantive scientific act and not a formality.

Practical rule I would adopt: for every quantitative claim, report the triple $(k_{\text{fitted here}},\,k_{\text{fitted elsewhere}},\,n)$ and compute $\tfrac{k}{2}\log n$. If your model's advantage over a null is smaller than its parameter penalty, you have nothing. In my experience this single computation would remove a large fraction of published "the data are consistent with the model" claims.
</details>

**★★★ 6. Run the microprocessor experiment on something you built.**
Take a circuit model whose algorithm you know exactly — the Wong–Wang network from `../part2-case-studies/C10-evidence-accumulation.md` is ideal — and apply the standard analysis toolkit blind. Report what you recover and what you miss.

<details markdown="1"><summary>Solution</summary>

Protocol: simulate a few hundred units' worth of the spiking or extended rate version. Then, pretending you do not know the mechanism, apply: (a) single-unit tuning curves vs. coherence and choice; (b) lesions (silence random subsets and measure performance); (c) PCA on trial-averaged population activity; (d) pairwise Granger causality between recorded units; (e) a decoder of choice from population activity.

Expected findings, and the honest scoring:
- (a) recovers choice-selective, coherence-graded ramping — *true and useful*. Tuning curves work here because the computation happens to be low-dimensional and aligned with the recorded variables.
- (b) shows graded degradation with no privileged units — *true, and it correctly indicates distributed coding*, but it tells you nothing about the algorithm. This is the 6502 lesion result, in the good case.
- (c) recovers a dominant "decision axis" plus a condition-independent component. This is genuinely informative, but note it does **not** distinguish integration from urgency-gating, nor recover $\lambda$, which is the model's central quantity (see `../part2-case-studies/C10-evidence-accumulation.md` §2.4).
- (d) Granger causality between two units of a network with strong shared input returns dense spurious connectivity, and the recovered "network" bears no resemblance to $J_{11}, J_{12}$. Score: fail.
- (e) works well and tells you nothing mechanistic — decodability is a property of the encoding, not the computation.

The overall lesson to write down for yourself: the analyses that succeeded did so because the network's computation was low-dimensional *and* you looked at variables the task defined. The analyses that failed are exactly the ones people use when they do not know what to look for — which is precisely when they are relied upon. Now repeat with the *leaky* variant ($\lambda\ne0$) and check whether any of (a)–(e) detects the difference. If none does, you have learned that your standard pipeline is blind to the model's key parameter, which is worth more than another paper.
</details>

---

## 7. Reading path

- **Marr (1982)**, *Vision*, Chapter 1 — read it (again) for: the three levels, and Marr's own warning that level-1 analysis is the hard part.
- **Lazebnik (2002)**, *Can a biologist fix a radio?* — read it for: the original version of the argument, short and funny and still true.
- **Golowasch, Goldman, Abbott & Marder (2002)**, *Failure of averaging in the construction of a conductance-based neuron model* — read it for: the two-page proof that a common practice is invalid.
- **Prinz, Bucher & Marder (2004)**, *Similar network activity from disparate circuit parameters* — read it for: the database, and the sheer size of the degenerate set.
- **Marder & Goaillard (2006)**, *Variability, compensation and homeostasis in neuronal and network function* — read it for: the synthesis; the paper that made degeneracy a first-class concept.
- **Gutenkunst et al. (2007)**, *Universally sloppy parameter sensitivities in systems biology models* — read it for: the eigenvalue spectra, and the realisation that this is generic rather than a quirk of neurons.
- **Roberts & Pashler (2000)**, *How persuasive is a good fit?* — read it for: the single most useful corrective in this entire note.
- **Pitt & Myung (2002)**, *When a good fit can be bad* — read it for: why parameter counting understates model flexibility.
- **Ng & Russell (2000)**, *Algorithms for inverse reinforcement learning* — read it for: the theorem that normative explanation is underdetermined, stated as geometry.
- **Bowers & Davis (2012)**, *Bayesian just-so stories in psychology and neuroscience* — read it for: the polemic, then read the replies and decide for yourself.
- **O'Leary, Williams, Franci & Marder (2014)**, *Cell types, network homeostasis, and pathological compensation from a biologically plausible ion channel expression model* — read it for: degeneracy as the output of a controller, not an accident.
- **Jonas & Kording (2017)**, *Could a neuroscientist understand a microprocessor?* — read it for: the benchmark, and then go and use it.
- **Krakauer, Ghazanfar, Gomez-Marin, MacIver & Poeppel (2017)**, *Neuroscience needs behavior: correcting a reductionist bias* — read it for: why level 1 cannot be skipped.
- **Chirimuuta (2018)**, *Explanation in computational neuroscience: causal and non-causal* — read it for: a philosopher taking seriously the idea that some neuroscientific explanations are not mechanistic.
- **Levenstein et al. (2023)**, *On the role of theory and modeling in neuroscience* — read it for: a taxonomy of what models are *for*, which prevents a lot of pointless argument.

---

## 8. Open problems

**A theory of when coarse-graining is exact.** Degeneracy suggests an effective low-dimensional theory, but we have no criterion for when a level-2 description is *derivable* from a circuit rather than merely consistent with it. Renormalisation-group ideas (Transtrum et al.) are the obvious tool and have barely been applied to neural circuits. Settling move: take a system where both levels are known — a trained RNN with a reverse-engineered reduced model — and construct the coarse-graining map explicitly, with error bounds.

**Benchmarks.** There is still no standard, adopted benchmark suite of ground-truth systems for validating computational-neuroscience analysis methods. This is embarrassing and fixable. Anyone reading this could build it.

**Quantifying degeneracy in the wild.** Almost all of the evidence is from models. The direct measurement — how much do the parameters of a given circuit vary across individuals *and* how much of that variation is behaviourally neutral — has been made in the STG and almost nowhere else. Settling move: paired parameter and behaviour measurements in any system where both are accessible, reported as the spectrum of $\mathcal I$ rather than as a variability statistic.

**Is compression the right criterion?** I have leaned on MDL here because it is the only quantitative version of "explanation" I trust. But compression and causal-counterfactual adequacy can diverge — a highly compressed description can get the counterfactuals wrong. Whether they can be reconciled, or whether we need both, is a genuinely open question and a good thing to argue about.

*Next:* `./S3-capstone.md` — put all of this to work on a system where the algorithm is genuinely unsettled, and where you will have to make the checklist above bind on your own hypotheses.
