---
name: principles
description: "Load at the start of any non-trivial work in this repo, and when adopting the principles in a project (bootstrap.md), and before writing a design, a plan, a loop or state file, delegating to subagents, editing AGENTS.md or the spec, or when code and design disagree. Then load principles-router for the area's modules."
---

# Principles

How work is run. Load only what the work needs:

- **Always:** the "In brief" of [workflow.md](../../_core/workflow.md), then the full section for what you are about
  to do (design, close, drift).
- **Running or joining lanes:** its "Lanes" section ([workflow.md](../../_core/workflow.md#lanes)) and the lane brief
  in the [workflow-files pattern](../../_core/patterns/workflow-files.md#lane-brief). A lane reads its brief and the
  rule sections the brief names, not all of `_core`.
- **Changing the principles or a design** (a fold-back, finalizing): [recurrence.md](../../_core/recurrence.md) and
  the `principles-self-check` skill ([../principles-self-check/SKILL.md](../principles-self-check/SKILL.md)).
- **Changing the package itself** (modules, ids, skills): [system.md](../../_core/system.md).
- **Adopting the package in a project:** start at [bootstrap.md](../../bootstrap.md), which runs the
  [adopting](../../_core/patterns/adopting.md) procedure.
- The process patterns (drift triage, workflow files, spec mirror, sources, adopting) are listed in
  [README.md](../../_core/README.md).

Then load `principles-router` ([../principles-router/SKILL.md](../principles-router/SKILL.md)) for the modules of the
area you touch.
