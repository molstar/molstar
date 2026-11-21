# Fold-back: what the fifth slice taught (draft, collected as it runs, for the author)

The fifth slice's fold-back, started early: lessons are added here as they come, and the full fold-back (what each gate,
lane and use found, the earlier rules checked, the self-check) is written at the close, as the principles' rule "The
fold-back is a step of every close" asks. **Nothing here is applied before the author's line at the end.** Each proposed
change is one line the author decides: "Proposed: … Why: … Where: …".

Where to: **P** principles (the instance's adaptations), **D** a design, **S** a package `_spec`, **C** a check, **L**
the next loop file.

## B. What this slice taught (collected)

1. **A lane's model follows its scope.**
   - Proposed: a lane with a closed item list and named stories (a fix lane over `fixes.md`, a storybook-only lane, a
     pure split, a fast review) runs on a Sonnet subagent; a lane that builds a design's core (language, resolve,
     layouts, runtimes on untrusted input) and its independent strong review keep the strong model; the brief names the
     model, the findings row records it, and the orchestrator's look and checks gate both alike.
   - Why: the author, 2026-11-21: "maybe use this more in general for well-scoped subtasks?", after NS (the Nightingale
     showcase: three stories working first time, light and dark) and FX2 (six of nine items fixed) ran on Sonnet near
     the usage limit; it saves tokens and shows how a weaker model handles the library.
   - Where: L (the loop's Lanes section, with B3's agent count), P (the instance's adaptations: who runs a lane).

2. **A finding that repeats a fixed one is a class: fix it by a rule and a sweep, not by the instance.**
   - Proposed: when a finding matches one already fixed (the same miss on another layer, layout or story), the
     orchestrator doesn't log another per-story item: it writes the rule the instances break (into the design or
     `_spec` it belongs to) and a sweep test over every story and layer that holds it, folds the open instances into
     that one item, and briefs it to a strong-model lane; FX briefs list which items are instances of a class.
   - Why: the author, 2026-11-21, after the hover box missed on a ring's own hover (FX3), the Regions density row, a
     ligand's score box and the ring-at-scale density bands: "we should not be fixing this one-by-one btw. since we
     already fixed a similar issue, we need a general rule for this and global fix." Fixed one at a time, each layer
     type passed its own story's test and the next one missed again.
   - Where: L (the loop's fix-point step: classify before logging), C (a sweep per contract: hover, empty tracks,
     cyclic intervals), P (the instance's adaptations).

## Decided by the author

Not yet.
