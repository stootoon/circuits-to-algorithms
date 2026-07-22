---
title: "S1 · RNN as model organism"
parent: Synthesis
nav_order: 1
---

# S1 — Train an RNN, reverse-engineer it, compare to biology

> **The method.** Train a recurrent network on the task your animal solves. Reverse-engineer the trained network completely — you have every weight and every state, so there is no excuse for not understanding it. Then ask, carefully, whether the animal did the same thing.
> **What it is good for.** Generating *specific, mechanistic* hypotheses at level 2 in cases where you cannot guess the algorithm from first principles. It is a hypothesis generator with unusually high resolution.
> **What it is not.** Evidence. A trained RNN is an existence proof that *some* system solves the task a particular way. Turning that into a claim about brains requires work that most papers do not do.
> **Structures thread.** Non-normal linear algebra (left vs. right eigenvectors), slow manifolds and normal forms, invariance classes of similarity measures. Read alongside `../structures/S-07-random-matrices-and-chaos.md` (what the untrained network's spectrum looks like) and `../structures/S-08-low-rank-connectivity.md` (what training adds to it); index in `../structures/README.md`.

---

## 1. Why bother

In `../part2-case-studies/C7-dopamine-and-td-learning.md` and `../part2-case-studies/C10-evidence-accumulation.md` the algorithm came first: someone knew about temporal-difference learning, or about Wald's test, and recognised the neural data as an implementation. That route works when the computational problem has a clean normative solution that a mathematician already found.

It fails for most of what animals do. What is the algorithm for "hold a rule in mind and use it to select which of two noisy sensory streams to integrate"? There is no theorem. But you *can* train a 100-unit RNN to do it in an afternoon, and then you have a complete, transparent, mechanistically specified solution to inspect. Sussillo & Barak (2013) made this a method rather than an anecdote: high-dimensional trained RNNs organise their dynamics around a small number of fixed and slow points, and the linearisations around those points give you a *readable* description of the computation.

The RNN is a model organism in the biologist's sense: a tractable system studied not because it is intrinsically interesting but because you can do experiments on it that you cannot do on the thing you care about. The comparison to *Drosophila* is exact, including the part where you have to argue, separately and with evidence, that the findings transfer.

## 2. The mathematics of reverse-engineering

Take the continuous-time rate RNN

$$\tau\dot x = -x + W\phi(x) + Bu + b + \sqrt{2\tau\sigma^2}\,\xi(t),\qquad z = w_{\text{out}}^\top\phi(x).$$

Write $F(x;u) = \big(-x + W\phi(x) + Bu + b\big)/\tau$.

**Fixed and slow points.** Define $q(x) = \tfrac12\|F(x;u_0)\|^2$ for a *fixed* input $u_0$ (usually the task's mean input in some epoch, or zero). Fixed points are the zeros of $q$; minimise $q$ by gradient descent using the analytic gradient

$$\nabla q = J(x)^\top F(x), \qquad J(x) = \frac{\partial F}{\partial x} = \frac{1}{\tau}\Big(-I + W\,\mathrm{diag}\big(\phi'(x)\big)\Big).$$

Local minima with $q>0$ but small are **slow points**: regions where the flow nearly stalls. These are not artefacts to be discarded. In a network that must hold information for a second, a manifold of slow points is exactly the mechanism, and demanding exact fixed points would make you miss it (Exercise 4).

**Linearisation and the two eigenvectors that matter.** Near a fixed point $x^*$,

$$\delta\dot x = J\,\delta x + \tfrac{1}{\tau}B\,\delta u .$$

$J$ is generically **non-normal**, so its left and right eigenvectors differ. With $J=\sum_a\lambda_a r_a l_a^\top$ and $l_a^\top r_b=\delta_{ab}$,

$$\delta x(t) = \sum_a e^{\lambda_a t}\,r_a\,\big(l_a^\top\delta x_0\big) \;+\; \frac{1}{\tau}\sum_a r_a\!\int_0^t\! e^{\lambda_a(t-s)}\,l_a^\top B\,\delta u(s)\,ds .$$

Read the roles off: **$r_a$ is the direction in state space that mode $a$ moves along; $l_a^\top B$ is the gain with which inputs drive it.** Almost everything interesting in trained RNNs is a statement about one of these two objects.

**Line attractors.** If one eigenvalue $\lambda_1\approx0$ and all others are strongly negative, then over times short compared with $1/|\lambda_1|$,

$$\delta x(t)\approx r_1\,\frac{1}{\tau}\int_0^t l_1^\top B\,\delta u(s)\,ds,$$

a perfect integrator along $r_1$ with input gain $l_1^\top B$. This is the RNN's implementation of the accumulator in `../part2-case-studies/C10-evidence-accumulation.md`, discovered by gradient descent rather than designed.

## 3. Mante et al.: gating is a rotation of the left eigenvector

Mante, Sussillo, Shenoy & Newsome (2013) trained monkeys — and an RNN — on context-dependent integration: two noisy streams (motion, colour), a contextual cue saying which one matters, and a choice based on the integral of the cued stream alone. The intuitive hypothesis is *input gating*: the irrelevant stream is suppressed before it reaches the integrator. In both PFC and the RNN, it is not. Both streams are represented, with similar magnitude, in both contexts.

The mechanism is the geometry above. In each context the network sits near a different region of a line attractor, so the Jacobian differs, so $l_1$ differs — while $r_1$, the integration axis, is essentially unchanged. Context does not change *what* is integrated along; it changes *which input directions project onto the integrating mode*. The irrelevant stream is represented, and pushed in directions that decay.

Two things to appreciate. First, this mechanism is invisible without the left/right distinction: if you only ever look at $r$'s (which is what PCA on population activity gives you), the two contexts look nearly identical and you conclude nothing happens. Second — and this is my favourite fact in this literature — a **linear** time-invariant RNN cannot do this at all. In an LTI system the gain $l_1^\top B$ is a property of the matrices, not of the state, so a constant context input can only add an offset. Context-dependent selection is a strict consequence of nonlinearity, and the specific role the nonlinearity plays is to make the Jacobian state-dependent (Exercise 3). That is a genuine level-2 insight extracted from a network, not from an animal.

## 4. Yang et al.: many tasks, and the question of modularity

Yang, Joglekar, Song, Newsome & Wang (2019) trained a single RNN on twenty cognitive tasks. They defined, for each unit, a task-variance profile, clustered units on it, and found a modest number of clusters that behave like functional modules: lesion a cluster and a specific subset of tasks degrades. They also found compositionality — units and directions reused across tasks that share computational structure (e.g. all the delay tasks share a memory mechanism).

Take the result seriously, but note what it does not establish. Clustering task-variance vectors will return clusters whether or not modules exist; the honest test is the lesion experiment, and the honest analysis is whether the lesion pattern is predicted *in advance* by the clustering. Driscoll, Shenoy & Sussillo (2024) sharpened this considerably by shifting the unit of analysis from neurons to **dynamical motifs** — reusable pieces of the flow field (an attractor, a saddle, a rotation) that are composed differently across tasks. That is the right level: a module defined by dynamics rather than by which neurons are active is invariant to the coordinate changes that plague neuron-level clustering, and it is the level at which one might actually compare a network to a brain.

## 5. Universality, and the reason to be nervous

Do different networks find the same solution? Maheswaranathan, Williams, Golub, Ganguli & Sussillo (2019) trained thousands of RNNs across architectures (vanilla, GRU, LSTM) and seeds, and found the *topology of fixed-point structure* to be strikingly universal — same number and type of fixed points, same qualitative flow — while the representational geometry (which directions, which units) varied wildly. Turner, Dabholkar & Barak (2021) mapped out how the solution space depends on task and regularisation and found genuine multiplicity: different hyperparameters give qualitatively different mechanisms for the same task.

The synthesis I would defend: **topology is often universal; geometry rarely is; and which regime you are in depends on how tightly the task constrains the dynamics.** Tasks with a strong normative solution (integration, delayed comparison) are highly constrained and yield universal topologies. Loosely constrained tasks yield a solution zoo, and any claim that "the RNN found the brain's solution" is then almost meaningless because the RNN would have found something else with a different seed.

This has a practical corollary that almost nobody follows: **train at least 20 seeds and report the distribution of mechanisms.** A single trained network is an anecdote. If 20/20 seeds produce the same fixed-point topology, you have a robust prediction; if 6/20 do, you have found one member of a family and should say so.

## 6. What "the RNN found the same solution as the brain" can mean

Four claims of increasing strength.

1. **Behavioural match.** Model and animal produce the same input–output map, including error patterns. Necessary, extremely weak.
2. **Representational match.** Population geometry is similar under some similarity measure. This is where most papers stop, and it is where the trouble is.
3. **Dynamical match.** The flow fields agree: same fixed-point topology, same slow manifolds, same linearised time constants and input-gain structure. Much stronger, and coordinate-free.
4. **Causal match.** Matched perturbations produce matched deviations. Perturb the model along $r_1$ and along an orthogonal decaying direction; the theory predicts a large behavioural effect for one and none for the other. Do the same experiment optogenetically. This is the only version of the claim that can genuinely fail, and it is therefore the only one worth much.

On (2): the standard measures — RSA (Kriegeskorte, Mur & Bandettini, 2008), CCA/SVCCA (Raghu et al., 2017), CKA (Kornblith et al., 2019), Procrustes-type shape metrics (Williams, Kunz, Kornblith & Linderman, 2021) — differ chiefly in their **invariance class**, and the invariance class is the whole ballgame. CCA is invariant to any invertible linear map, which makes it far too permissive: two representations related by an invertible but wildly ill-conditioned map have CCA similarity 1 while supporting completely different noise-robust decoders (Exercise 5). CKA is invariant only to orthogonal transforms and isotropic scaling, but is dominated by high-variance directions, so a shared high-variance nuisance component drives it to 1 regardless of task-relevant structure. Ding, Denain & Steinhardt (2021) stress-tested these measures against functional behaviour and found the sensitivities do not line up with what we care about. Ostrow, Eisen, Kozachkov & Fiete (2023) proposed comparing *dynamics* rather than geometry (dynamical similarity analysis), which is the right instinct and closer to claim (3).

And read Schaeffer, Khona & Fiete (2022) before you write the sentence "our network developed X-like units, as observed in the brain". They showed that grid-cell-like tuning emerges from trained path-integration networks only under a narrow band of implementation choices that were, in effect, back-fitted to the known answer. A network trained to reproduce a phenomenon reproduces the phenomenon; that is not evidence about mechanism. The discipline this demands: **decide in advance what would count as a mismatch, and report the comparison for a control model that you expect to fail.**

## 7. A recipe you can actually follow

1. **Write the task as the animal experiences it**, including trial timing, stimulus noise statistics, and the response epoch. Get the noise right — the noise level determines whether the network needs to integrate at all.
2. **Choose the smallest architecture that can solve it**, with biological constraints you are willing to defend ($\tau$ matched to membrane/synaptic time constants, private noise injected during training, rate nonlinearity, and — if you care — Dale's law). Regularise: $L^2$ on rates and weights. Unregularised networks find high-dimensional, uninterpretable, non-universal solutions.
3. **Train ≥20 seeds.** Discard nothing.
4. **Find fixed and slow points** for each task-relevant fixed input, from initial conditions sampled from the network's own trajectories (not from random points — you want the fixed points the network actually visits).
5. **Classify.** For each fixed point, compute $J$'s spectrum; identify near-zero eigenvalues (memory/integration), complex pairs (oscillation), unstable directions (decisions/saddles). Compute $l_a^\top B$ for each input channel.
6. **Build the reduced model.** Write down the two- or three-dimensional dynamical system that reproduces the network's behaviour. If you cannot, you have not finished reverse-engineering — the reduced model *is* the extracted algorithm, and by the compression criterion of `./S2-degeneracy-and-limits.md` it had better be much smaller than the network.
7. **Check universality** across seeds and across at least one alternative architecture.
8. **Derive experiment-ready predictions** at the level of claims (3) and (4), then test on neural data — including a pre-registered statement of what would falsify the match.

---

## 8. Model to build

Train a context-dependent integrator with plain numpy BPTT, then find its fixed points.

```python
import numpy as np
rng = np.random.default_rng(0)
N, T, al = 64, 75, 0.2
W  = rng.standard_normal((N, N))/np.sqrt(N); B = rng.standard_normal((N, 4))*0.5
wo = rng.standard_normal(N)/np.sqrt(N);      bb = np.zeros(N)

def batch(bs=64):
    ctx = rng.integers(0, 2, bs)
    coh = rng.choice([-1., -.5, -.25, .25, .5, 1.], (bs, 2))
    u = rng.standard_normal((T, bs, 4))*0.6
    u[:, :, :2] += coh[None]*0.15
    u[:, :, 2:] = 0.0; u[:, np.arange(bs), 2 + ctx] = 1.0
    return u, np.sign(coh[np.arange(bs), ctx])

def fwd(u):
    bs = u.shape[1]; x = np.zeros((T+1, bs, N)); r = np.zeros((T+1, bs, N))
    for t in range(T):
        x[t+1] = (1-al)*x[t] + al*(r[t] @ W.T + u[t] @ B.T + bb)
        r[t+1] = np.tanh(x[t+1])
    return x, r

def grads(u, tgt):
    x, r = fwd(u); bs = u.shape[1]
    e = r[T] @ wo - tgt; L = 0.5*np.mean(e**2)
    gwo = (e[:, None]*r[T]).mean(0)
    dx  = (e[:, None]*wo[None, :])*(1 - r[T]**2)/bs
    gW = np.zeros_like(W); gB = np.zeros_like(B); gb = np.zeros_like(bb)
    for t in range(T, 0, -1):
        gW += al*(dx.T @ r[t-1]); gB += al*(dx.T @ u[t-1]); gb += al*dx.sum(0)
        dx = (1-al)*dx + al*(dx @ W)*(1 - r[t-1]**2)
    return L, gW, gB, gb, gwo

P = [W, B, bb, wo]; M = [np.zeros_like(p) for p in P]; V = [np.zeros_like(p) for p in P]
for it in range(4000):                                   # Adam
    u, tgt = batch(); L, *G = grads(u, tgt)
    for k, (p, g) in enumerate(zip(P, G)):
        M[k] = 0.9*M[k] + 0.1*g; V[k] = 0.999*V[k] + 0.001*g*g
        p -= 2e-3*M[k]/(np.sqrt(V[k]) + 1e-8)
    if it % 500 == 0: print(it, round(L, 4))

def slow_points(u0, n_init=200, steps=4000, lr=0.02):     # minimise q = .5||F||^2
    _, r = fwd(np.repeat(u0[None, None], T, 0).repeat(n_init, 1))
    xs = r[rng.integers(10, T, n_init), np.arange(n_init)].copy()
    xs = np.arctanh(np.clip(xs, -.999, .999))
    for _ in range(steps):
        rr = np.tanh(xs); f = -xs + rr @ W.T + u0 @ B.T + bb
        xs -= lr*(-f + (f @ W)*(1 - rr**2))               # grad q = J^T F
    return xs, 0.5*np.sum((-xs + np.tanh(xs) @ W.T + u0 @ B.T + bb)**2, 1)
```

**What to look for.**
- Loss below ~0.15 within a few thousand iterations; psychometric curves that depend on the *cued* coherence and are flat in the uncued one. If they are not flat, the network is not doing the task, and nothing downstream means anything.
- Run `slow_points` with `u0 = np.array([0,0,1,0])` and `[0,0,0,1]`. Cluster the returned $x$'s: you should find an approximately one-dimensional manifold of slow points (small $q$), not isolated points. That manifold is the line attractor.
- At each slow point compute `J = -np.eye(N) + W*(1-np.tanh(x)**2)` and its eigenvalues. Expect one eigenvalue with $|\lambda|\ll1$ and the rest well into the left half plane.
- Now the payoff: compute the right eigenvector $r_1$ and left eigenvector $l_1$ for the slow mode in each context, and compute $l_1^\top B$ for the two evidence channels. **Prediction:** $r_1$ is nearly the same across contexts (cosine similarity $>0.9$), while $l_1^\top B$ selects the cued channel and suppresses the other. That is Mante et al.'s result, reproduced in 60 lines.
- Control: fit a linear RNN to the same task and confirm it fails to gate (Exercise 3).

---

## 9. Exercises

**★ 1. Two parameterisations, one Jacobian.**
Compare $\tau\dot x=-x+W\phi(x)+Bu$ (current/pre-activation state) with $\tau\dot r = -r + \phi(Wr+Bu)$ (rate/post-activation state). Compute both Jacobians at a fixed point and relate their spectra.

<details markdown="1"><summary>Solution</summary>

Pre-activation: $J_x = \tfrac1\tau(-I + W D)$ with $D=\mathrm{diag}(\phi'(x^*))$.
Post-activation: $J_r = \tfrac1\tau(-I + D'W)$ with $D'=\mathrm{diag}(\phi'(Wr^*+Bu))$.

At corresponding fixed points $D=D'$, and $WD$ and $DW$ are similar whenever $D$ is invertible: $DW = D(WD)D^{-1}$. Hence identical eigenvalues, so **stability, time constants, and topology are parameterisation-independent** — good, since they had better be. The *eigenvectors* are not: right eigenvectors transform as $r\mapsto Dr$. So any claim about "the integration axis in neural space" must state the coordinate convention, and comparisons of eigenvector geometry between papers using different conventions are meaningless unless corrected. If $\phi'$ saturates ($D$ near-singular) the similarity transform degrades and even the spectra can differ in the saturated subspace — worth checking numerically before asserting equivalence.
</details>

**★★ 2. How precisely must a line attractor be tuned?**
The network must hold a memory over a trial of duration $T_{\text{trial}}$. Quantify the tolerance on the slow eigenvalue and hence on the weights.

<details markdown="1"><summary>Solution</summary>

Memory decays as $e^{\lambda_1 t}$; requiring less than a fraction $\epsilon$ of leak over the trial gives $|\lambda_1| < \epsilon/T_{\text{trial}}$. With $\tau=20$ ms, $T_{\text{trial}}=1$ s, $\epsilon=0.05$: $|\lambda_1|<0.05\,\mathrm{s}^{-1}$, i.e. $|\tau\lambda_1| < 10^{-3}$ in dimensionless units, against an $O(1)$ spectral scale. Three orders of magnitude of tuning.

Sensitivity to weights: first-order perturbation theory gives $\delta\lambda_1 = l_1^\top(\delta J)r_1$. For a random perturbation $\delta W$ with i.i.d. entries of s.d. $s$, $\delta J = \delta W D/\tau$, so $\mathrm{sd}(\delta\lambda_1) \approx \tfrac{s}{\tau}\|l_1\|\|D r_1\|$. Because $J$ is non-normal, $\|l_1\|$ can be $\gg1$ (the eigenvector condition number $\kappa_1 = \|l_1\|\|r_1\|$), so a network can be *far* more fragile than the naive estimate. Requiring $\mathrm{sd}(\delta\lambda_1)<\epsilon/T_{\text{trial}}$ gives $s < \epsilon\tau/(T_{\text{trial}}\kappa_1)$.

Two consequences. (i) Biological line attractors need feedback tuning (this is Seung's classic point about the oculomotor integrator), and the condition number tells you how much. (ii) Report $\kappa_1$ whenever you claim a line attractor in a trained network: a solution with $\kappa_1=10^3$ is a numerical curiosity, not a mechanistic hypothesis.
</details>

**★★★ 3. A linear RNN cannot context-gate.**
Prove it, then identify exactly what the nonlinearity supplies.

<details markdown="1"><summary>Solution</summary>

Let $\tau\dot x = -x + Wx + Bu$, and split $u = (s, c)$ into evidence $s$ and a constant context input $c$. Linearity gives superposition: $x(t) = x_s(t) + x_c(t)$, where $x_c$ depends only on $c$ and $x_s$ only on $s$. The readout $z=w^\top x$ is therefore $w^\top x_s(t) + w^\top x_c(t)$: the context contributes an additive, evidence-independent offset and cannot alter the *gain* on $s$. Formally, $\partial z/\partial s$ is independent of $c$ for any LTI system. Multiplying the evidence by a context-dependent factor is a second-order (bilinear) operation, and no LTI system computes one.

What nonlinearity supplies: with $\tau\dot x=-x+W\phi(x)+Bu$, the context input moves the operating point $x^*(c)$, and the Jacobian $J(x^*)=\tfrac1\tau(-I+W\,\mathrm{diag}(\phi'(x^*)))$ therefore *depends on $c$*. The effective input gain onto the slow mode, $l_1(c)^\top B$, is a function of context. So the "multiplication" is implemented as a state-dependent linearisation: the network is locally linear, but *which* linear system it is depends on where it is sitting. This is the general trick — and it is worth naming, because you will meet it again in gain modulation, in attention, and in the synchrony-as-gating hypothesis of `../part2-case-studies/C2-transient-synchrony.md`.

Empirical check: train the linear version. It reaches a loss floor corresponding to integrating a fixed weighted mixture of both streams, and its psychometric curve depends on the *uncued* coherence — a signature you can also look for in animals that have not yet learned the task.
</details>

**★★ 4. Slow points are the point.**
Show $\nabla q = J^\top F$; explain why local minima of $q$ with $q>0$ are functionally meaningful; and give the criterion for deciding whether a slow point is "slow enough".

<details markdown="1"><summary>Solution</summary>

$q=\tfrac12 F^\top F$, so $\nabla q = (\partial F/\partial x)^\top F = J^\top F$. (Note $\nabla q = 0$ also at points where $F\ne0$ but $F\in\ker J^\top$ — these are precisely the interesting non-fixed minima.)

A minimum with small $q>0$ means the flow is slow but not stationary: the state drifts along a low-speed valley. Networks routinely implement memory this way rather than with exact continuous attractors, because gradient descent has no reason to make $\lambda_1$ exactly zero — it only needs the drift to be negligible *on the task's timescale*.

Criterion: the relevant comparison is not $q$ against zero but the drift speed against the task. Convert: the state moves at speed $\|F\|=\sqrt{2q}$, so over the trial it travels $\approx\sqrt{2q}\,T_{\text{trial}}$. Compare that with the extent of the slow manifold (say, the range of states used to encode different evidence levels, $\Delta$). A slow point is "slow enough" iff $\sqrt{2q}\,T_{\text{trial}} \ll \Delta$. Reporting raw $q$ values without this normalisation, as is common, tells the reader nothing.
</details>

**★★ 5. Similarity metrics can be fooled, in both directions.**
Construct (a) two representations with linear-CCA similarity 1 that support very different decoders, and (b) two representations with CKA $\approx1$ that differ completely in task-relevant structure.

<details markdown="1"><summary>Solution</summary>

(a) Let $s,t$ be independent standardised task variables. $X=[\,s,\;t\,]$ and $Y=[\,s,\;10^{-4}t\,]$. $Y=XA$ with $A$ invertible and diagonal, and linear CCA is invariant to invertible linear maps, so CCA similarity is exactly 1. But add readout noise of s.d. $10^{-2}$: $t$ is decodable from $X$ with $d'\approx100$ and from $Y$ with $d'\approx0.01$. Identical by the metric, functionally unrelated. The lesson: any metric invariant to *all* invertible linear maps has thrown away the scale information that determines what a noisy downstream neuron can read.

(b) Let $n$ be a shared nuisance variable with variance 100, and let $X=[\,n,\;s\,]$, $Y=[\,n,\;t\,]$ with $s,t$ independent of each other and of $n$, variance 1. The Gram matrices are dominated by the $n$ term: $XX^\top = nn^\top + ss^\top$, $YY^\top = nn^\top+tt^\top$. Linear CKA is $\frac{\|\tilde X^\top\tilde Y\|_F^2}{\|\tilde X^\top \tilde X\|_F\|\tilde Y^\top\tilde Y\|_F}$; with variance ratio 100:1 the shared component contributes $\sim10^4$ to numerator and denominators and the unshared contributes $\sim1$, giving CKA $\approx 1/(1+O(10^{-2}))\approx0.99$. Yet the task-relevant dimensions are orthogonal and encode different variables.

Practical rule: never report a similarity number alone. Report it alongside (i) a *lower* bound established by a control model that should not match, and (ii) a decoding-based comparison that respects noise. Better still, compare at the level of dynamics (claim 3 in §6) where the objects being compared — fixed-point topology, eigenvalues — are already coordinate-free.
</details>

**★★★ 6. Universality, measured.**
Design and run an experiment on your own trained networks that estimates how universal the solution is, and report it as a number.

<details markdown="1"><summary>Solution</summary>

Protocol: train $S=20$ seeds. For each, extract the fixed/slow-point set for each task epoch and compute a **topological fingerprint**: the number of slow-point clusters, the dimension of each (from the local PCA of the manifold), the number of unstable directions at each, and the sorted real parts of the top-$k$ eigenvalues. Define the fingerprint distance between two networks as the sorted-eigenvalue $L^1$ distance plus a penalty for mismatched counts. Report (i) the fraction of seed pairs with matching discrete structure and (ii) the distribution of eigenvalue distances, against a null obtained by comparing networks trained on *different* tasks.

Expected outcome for the context-dependent integrator: near-universal discrete structure (a slow manifold plus two attractors; ≥18/20 seeds), and eigenvalue distributions overlapping heavily within-task and separated across tasks. Then repeat with the $L^2$ rate penalty removed and watch universality drop — which tells you that "universality" is partly a statement about the regulariser, not about the task. Report *that* too. A paper reporting one seed and asserting universality has done a different experiment from the one it claims.
</details>

---

## 10. Reading path

- **Sussillo & Barak (2013)**, *Opening the black box: low-dimensional dynamics in high-dimensional recurrent neural networks* — read it for: the fixed/slow-point method and the linearisation logic. Everything else builds on this.
- **Mante, Sussillo, Shenoy & Newsome (2013)**, *Context-dependent computation by recurrent dynamics in prefrontal cortex* — read it for: the canonical worked example, and the selection-vector geometry.
- **Sussillo (2014)**, *Neural circuits as computational dynamical systems* — read it for: the method stated as a research programme.
- **Barak (2017)**, *Recurrent neural networks as versatile tools of neuroscience research* — read it for: a clear-eyed account of what this method can and cannot deliver.
- **Song, Yang & Wang (2016)**, *Training excitatory–inhibitory recurrent neural networks for cognitive tasks: a simple and flexible framework* — read it for: how to impose Dale's law and other biological constraints without breaking training.
- **Yang, Joglekar, Song, Newsome & Wang (2019)**, *Task representations in neural networks trained to perform many cognitive tasks* — read it for: multitask training, clustering, and compositionality — and to practise being sceptical about clustering.
- **Maheswaranathan, Williams, Golub, Ganguli & Sussillo (2019)**, *Universality and individuality in neural dynamics across large populations of recurrent networks* — read it for: the topology-universal / geometry-individual distinction.
- **Turner, Dabholkar & Barak (2021)**, *Charting and navigating the space of solutions for recurrent neural networks* — read it for: evidence that solution multiplicity is real and controllable.
- **Driscoll, Shenoy & Sussillo (2024)**, *Flexible multitask computation in recurrent networks utilizes shared dynamical motifs* — read it for: the motif-level unit of analysis, which is where comparisons to brains should be made.
- **Kornblith, Norouzi, Lee & Hinton (2019)**, *Similarity of neural network representations revisited* — read it for: CKA, and the invariance-class argument that organises the whole comparison literature.
- **Williams, Kunz, Kornblith & Linderman (2021)**, *Generalized shape metrics on neural representations* — read it for: how to make similarity a proper metric, which you need if you want to cluster or do statistics.
- **Ding, Denain & Steinhardt (2021)** — read it for: a systematic demonstration that popular similarity measures are not sensitive to the things we care about.
- **Schaeffer, Khona & Fiete (2022)**, *No free lunch from deep learning in neuroscience* — read it for: the strongest available warning against "our network reproduced the phenomenon, therefore mechanism".
- **Ostrow, Eisen, Kozachkov & Fiete (2023)**, *dynamical similarity analysis* — read it for: comparing computations by their dynamics rather than their geometry.

---

## 11. Open problems

**When does training on a task recover the brain's algorithm, and when does it recover an algorithm?** We have no theory of which tasks are constraining enough. A serious attack would characterise the solution space (à la Turner et al.) as a function of task structure and identify a computable constraint index. Until then, universality must be measured case by case, per §9.6.

**Learning rule mismatch.** Animals do not do BPTT. It is an open question — not a rhetorical one — whether biologically plausible learning rules land in the same solution basins. Testable: train matched architectures with node perturbation, feedback alignment, and BPTT, and compare fixed-point topologies. If the topology is invariant to the learning rule, one large worry evaporates; if not, the entire literature inherits a confound.

**Comparison without a ground truth.** Every similarity measure has been validated on synthetic cases chosen by its authors. What the field needs is the Jonas & Kording move from `./S2-degeneracy-and-limits.md`: a system where the algorithm is known exactly, into which we can inject known mechanistic differences, and on which every comparison method is scored. Trained RNNs are the obvious substrate. Someone should just do it.

*Next:* `./S2-degeneracy-and-limits.md`, which asks when there is no clean algorithm to find — and gives you a checklist for deciding whether any claimed circuit→algorithm conversion, including your own, is worth believing.
