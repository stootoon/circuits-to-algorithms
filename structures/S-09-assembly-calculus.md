---
title: S-09 Assembly calculus
parent: Structures
nav_order: 9
---

# S-09 — Assembly calculus: inventing the missing formalism

> **The object.** A formal computational model whose data type is the *assembly* — a set of $k$
> co-firing neurons — with primitive operations **project**, **associate**, **merge**, and
> **pattern-complete**, convergence theorems, and completeness results.
> **Born from.** The absence of any formal language at the level where Hebb's cell assemblies
> live: above single neurons, below symbolic computation.
> **Mathematical home.** Random graph theory, probabilistic analysis of algorithms, theory of
> computation.
> **Situation.** #7, the missing formalism.
> **Novelty.** **Genuinely new.** A computational model with its own semantics and its own
> completeness theorem, constructed to fit a biological substrate.

---

## 1. Where it comes from

Hebb proposed cell assemblies in 1949. For seventy years the idea was central to how
neuroscientists *talk* and absent from how they *compute*. You could model single neurons in
exquisite detail, or you could write down a Turing machine, and there was nothing in between
that respected the biology.

The gap is not cosmetic. When you say "the assembly for *cat* becomes associated with the
assembly for *animal*," you are making a claim that is precise-sounding and formally empty: no
definition of assembly, no definition of associate, no theorem that it works. Papadimitriou,
Vempala and colleagues closed the gap by building the formalism, and the result is a system with
theorems.

The lesson generalizes beyond assemblies, and it is why this note is in the course: **if you
cannot state your hypothesis without hand-waving, the formalism is missing, and writing it down
is a contribution — often a bigger one than the hypothesis.**

---

## 2. The object, precisely

### The substrate (the NEMO model)

- A finite set of **brain areas**, each with $n$ excitatory neurons.
- Within and between areas, directed synapses exist independently with probability $p$
  (Erdős–Rényi). Everything is random; no structure is designed in.
- **$k$-cap:** in each area, at each discrete time step, exactly the $k$ neurons with the highest
  total synaptic input fire. This models local inhibition enforcing a fixed activity level.
- **Hebbian plasticity:** if $j$ fired at $t-1$ and $i$ fires at $t$ and $j \to i$ exists, then
  $w_{ij} \mathrel{\times}= (1+\beta)$.

That is the entire model: random graph, top-$k$, multiplicative Hebb. Three parameters
$(n,p,k)$ plus $\beta$.

Notice immediately that **$k$-cap is exactly the operation of the mushroom body**: the APL/GGN
feedback inhibition in the insect that keeps Kenyon cell activity sparse is a biological $k$-cap.
[C1](../part2-case-studies/C1-expansion-and-sparsening.md) and this note are describing the same
circuit motif at two levels of abstraction — one as a hashing primitive, one as the substrate of
a programming language.

### An assembly

An **assembly** is a set of $k$ neurons in an area that fires stably and recurrently as a unit —
a fixed point of the $k$-cap dynamics with plasticity.

### The operations

**Projection.** Repeatedly fire assembly $x$ in area $A$. Downstream area $B$ initially fires a
somewhat different $k$-set each step, but Hebbian strengthening makes the recruited set
increasingly stable, and it converges to a fixed assembly $y \subset B$.

*This is a theorem, not an observation.* Under conditions on $(n,p,k,\beta)$, the projection
converges with high probability, in $O(\log)$-ish time, and the total number of distinct neurons
ever recruited is bounded. The proof is a probabilistic argument about the tails of the input
distribution under repeated top-$k$ selection — genuinely nontrivial, and the mathematical heart
of the enterprise.

**Association.** Co-firing assemblies $x$ and $y$ in different areas increases their overlap —
a fraction of each assembly's neurons come to be shared, implementing a graded, symmetric
relation.

**Merge.** Given assemblies $x \subset A$ and $y \subset B$, create $z \subset C$ with strong
reciprocal connections to both. Merge is the operation that gives the calculus **recursive
structure** — it is what lets you build a representation of a *pair*, and hence of a tree, and
hence of a parse.

**Pattern completion.** Firing a subset of an assembly recruits the rest.

### The results

- Convergence theorems for projection and the other operations.
- **Completeness:** the calculus can simulate arbitrary computation (under assumptions on the
  parameters and with appropriate space bounds) — so it is a *bona fide* model of computation,
  not a collection of demos.
- **Applications:** a syntactic parser for natural language implemented purely in assembly
  operations (Mitropolsky and colleagues), which is the most striking demonstration that the
  formalism has real expressive reach.

---

## 3. Why the neuroscience forced it

No existing computational model has the right primitives. Turing machines and RAM models assume
addressable memory, which brains do not have. Neural network models operate on real-valued
vectors and provide no notion of a discrete, nameable item. Production systems and symbolic
architectures assume symbols with no account of where they come from.

Assembly calculus is built so that its primitives are (a) implementable by the actual substrate
— random connectivity, inhibition, Hebbian plasticity, all uncontroversial — and (b) rich enough
to compose into cognition. Getting both simultaneously is the achievement, and it required
inventing the language rather than borrowing one.

---

## 4. How to recognize the pattern elsewhere

Situation #7 announces itself when **your field's central concepts are nouns without
definitions.** "Cell assembly." "Binding." "Chunk." "Schema." "Attractor state" used loosely.
And, directly relevant here: **"transient synchrony."**

The diagnostic: try to write your hypothesis as a program. If you cannot say what the data types
and operations are, they do not exist yet, and constructing them is the work.

The construction recipe, from this case:

1. **Fix the substrate** — the minimal set of uncontroversial biological mechanisms.
2. **Define the data type** as an object the substrate naturally produces and stabilizes.
3. **Define operations** and *prove they terminate/converge* on the substrate. This step is where
   most attempts die, and passing it is what separates a formalism from a notation.
4. **Establish expressive power** — what can be computed?
5. **Demonstrate on a real cognitive task.**

**The direct challenge for your work:** there is no assembly-calculus analogue for *temporal*
codes. Assembly calculus is fundamentally about which neurons fire, with time as discrete
steps and no role for phase, synchrony, or oscillation. A calculus whose data type was a
*synchronous sub-assembly bound to an oscillation phase*, with operations for binding,
sequencing, and multiplexing across phases, does not exist. The locust's 20–30 Hz cycle with
odour-specific subsets synchronizing on particular cycles is exactly a substrate that wants such
a formalism. Lisman & Idiart's theta-gamma multiplexing is the closest existing thing and it is
a mechanism, not a calculus. Writing that formalism down would be a genuine contribution and it
is precisely the kind of thing this course exists to point at.

---

## 5. Exercises

**Ex S9.1 ★★** — Implement projection and verify convergence empirically.

<details markdown="1"><summary>Solution</summary>

```python
import numpy as np
rng = np.random.default_rng(0)
n, p, k, beta, T = 10000, 0.01, 100, 0.10, 30

W = (rng.random((n,n)) < p).astype(float)      # A -> B synapses
x = rng.choice(n, k, replace=False)            # assembly in A (fixed, fired every step)

prev = None; overlaps = []
for t in range(T):
    inp = W[:, x].sum(1)
    if prev is not None:                       # recurrent support from B's own assembly
        inp += W[:, prev].sum(1)
    y = np.argpartition(-inp, k)[:k]
    if prev is not None:
        overlaps.append(len(set(y) & set(prev)) / k)
        W[np.ix_(y, x)] *= (1+beta)            # Hebbian: A->B
        W[np.ix_(y, prev)] *= (1+beta)         # Hebbian: B->B
    prev = y
print(np.round(overlaps, 3))
```

The overlap between successive $k$-caps starts well below 1 and climbs toward 1: the assembly
converges. Track the cumulative set of neurons ever recruited — it saturates, which is the
"bounded total recruitment" part of the theorem.

Then set $\beta = 0$ and rerun: convergence fails or is far weaker. Plasticity is doing real
work, and seeing exactly how much is the point of the exercise.
</details>

**Ex S9.2 ★★** — Argue informally why $k$-cap plus Hebbian plasticity produces convergence, in
terms of the input distribution's tail.

<details markdown="1"><summary>Solution</summary>

At each step, neuron $i \in B$ receives input $\sum_{j \in x} w_{ij}$ plus recurrent input from
the previous cap. Across $i$, this is a sum of $\approx pk$ random terms, so approximately
Gaussian with mean $\mu$ and s.d. $\sigma$. The $k$-cap selects the top $k$ of $n$, i.e. the
upper $k/n$ tail — for $k/n$ small, these are neurons several $\sigma$ above the mean.

Hebbian plasticity multiplies the weights *from the winners to the winners* by $(1+\beta)$. On
the next step, a neuron that won previously has systematically larger input — by a factor that
grows geometrically with the number of times it has won. So winning is self-reinforcing: the
distribution develops a separated upper cluster of "incumbents."

Convergence follows once the incumbents' advantage exceeds the typical fluctuation gap at the
cap boundary. The competition is between the geometric growth $(1+\beta)^t$ of incumbent input
and the $O(\sigma)$ spacing of order statistics near the $k$th largest of $n$ Gaussians (which
is $O(\sigma/\sqrt{\log(n/k)})$ locally). Once the former dominates, the cap stops changing.

The full proof handles the dependence between steps and the fact that the recurrent input is
itself changing; that is where the real work is. But the mechanism is exactly this
rich-get-richer dynamic in the tail, and understanding it tells you the parameter regime: you
need $k \ll n$ (a real tail to select from) and $\beta$ large enough to outrun the order-statistic
spacing, but not so large that the first accidental winner locks in immediately.
</details>

**Ex S9.3 ★★★** — *(Situation, not object.)* Sketch the data types and operations for a
*temporal* assembly calculus in which the primitive is a synchronous sub-assembly bound to an
oscillation phase. What must be proven for it to be a formalism rather than a notation?

<details markdown="1"><summary>Solution</summary>

Open, and genuinely unsolved. A sketch of what a serious attempt would need:

*Substrate.* An area with a global oscillation of period $T$ divided into $P$ phase slots;
neurons with a membrane time constant short relative to $T/P$ so that coincidence detection is
selective; local inhibition enforcing a per-slot cap; STDP rather than rate-Hebb.

*Data type.* A **phase-bound assembly**: a pair $(S, \varphi)$ with $S$ a set of $k$ neurons and
$\varphi \in \{1,\dots,P\}$ a phase slot. Two assemblies in different slots coexist without
interference — this is the whole point and is what a rate-based calculus cannot express.

*Operations.* `project` (as before, but preserving phase); `bind` (co-activate $(S_1,\varphi)$
and $(S_2,\varphi)$ to create a conjunctive representation downstream — selective because
downstream coincidence detectors respond only to same-slot input); `sequence` (create an
assembly whose components occupy successive slots, giving ordered structure for free);
`multiplex` (hold $P$ assemblies simultaneously, the Lisman–Idiart capacity story); `rephase`
(move an assembly between slots — needed for any nontrivial manipulation, and probably the
hardest to implement).

*What must be proven.* (i) **Slot stability** — assemblies stay in their slots under plasticity,
with jitter bounded below the slot width; this is the analogue of the projection convergence
theorem and is the make-or-break result. (ii) **Non-interference** — capacity $P$ scales as
$T/(\text{jitter})$ and cross-slot crosstalk is bounded. (iii) **Composition** — `bind` and
`sequence` produce assemblies on which the operations can be applied again, i.e. the type is
closed. (iv) **Expressive power** — at minimum, that `sequence` + `bind` give recursive
structures, so the calculus can represent trees.

*Why it matters for you.* Every one of these is a formalization of something the locust
literature asserts informally. "Synchrony binds features" is operation `bind` without a
convergence proof. "The oscillation provides a temporal frame" is the slot structure without a
capacity theorem. Writing the calculus would convert a body of suggestive experimental claims
into a set of statements with truth values — which is the single highest-leverage thing anyone
could do for that literature, and it is a mathematics problem, not an experimental one.
</details>

---

## 6. Reading

- **Papadimitriou, Vempala, Mitropolsky, Collins & Maass (2020)**, *Brain computation by
  assemblies of neurons* (PNAS) — read it for: the model, the operations, the convergence
  results. This is the paper.
- **Mitropolsky and colleagues**, on a biologically plausible parser implemented in assembly
  operations — read it for: the demonstration that the calculus reaches real linguistic
  structure.
- **Hebb (1949)**, *The Organization of Behavior* — read it for: the original proposal, and to
  appreciate how long a good idea can sit un-formalized.
- **Valiant (1994)**, *Circuits of the Mind* — read it for: the earlier and unjustly neglected
  attempt at exactly this problem (neuroidal model, vicinal algorithms). Reading Valiant then
  Papadimitriou shows you what changed and why the second attempt landed.
- **Lisman & Idiart (1995)**, on storage of multiple items in short-term memory via nested
  theta/gamma oscillations (Science) — read it for: the closest existing thing to a temporal
  calculus, and the starting point for Ex S9.3.
- **Buzsáki (2010)**, on neural syntax and cell assemblies (Neuron) — read it for: the
  experimentalist's case that assemblies with temporal structure are the right unit, which is
  the empirical motivation for the formalism that does not yet exist.
