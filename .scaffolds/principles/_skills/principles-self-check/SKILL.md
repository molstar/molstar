---
name: principles-self-check
description: "Load when a design is about to be agreed, and when work is about to change the principles or the designs: folding back what an execution session or a prototype taught (a fold-back), finalizing a session whose work touched .principles/ or .design/, writing a session record's Self-check section, or when a tension doesn't dissolve under clarification. Runs the four-step self-check on the change and records its delta."
---

# Principles self-check

Recommend it to the session, then run it, when:

- **a design is about to be agreed**: the claim is the design's own, step 2 negates its premise (for an action, its
  inverse or symmetric case), and the result is a short section of the design
  ([recurrence.md, "On a design"](../../_core/recurrence.md#self-check));
- **a fold-back** turns what a session or prototype taught into changes to `.principles/` or `.design/` (its
  drift notes, gaps and residuals becoming rules, designs or open items);
- **finalizing** work that changed the principles or a design, before the session record is written;
- **a tension** between rules, the author's words, a precedent or the code survives clarification.

How:

1. The check itself, four steps run on the evidence in `main`, stopping only at the fourth:
   [recurrence.md, "Checking the loop against itself"](../../_core/recurrence.md#self-check). Take as the claim
   what the change asserts (for a fold-back: "these proposals make the principles truer to what was built").
2. Merge what integrates; carry what doesn't as a standing tension in the design it belongs to
   ([workflow.md](../../_core/workflow.md), §surface-tensions).
3. Record the run's delta in the session record's Self-check section
   ([workflow-files pattern, session record](../../_core/patterns/workflow-files.md#session-record-recordswork-iddate-namemd)).

A run that changes nothing and carries nothing is suspect; say so rather than invent a delta.
