# Core

- **Kind:** rule set (with process patterns)
- **Status:** stable: every project starts with it and rarely changes it; a project may still change it
- **Scope:** how the principles work (modules, kinds, ids, skills, adapting) and how work is run (design, plan, loop,
  state, record; lanes; drift; verification; the gaps; machine-invariant). Not the tools a workspace uses
  (`workspace`) or how code is built (the other modules).
- **Sources:** see the catalogue ([../README.md](../README.md#modules-the-catalogue))
- **Depends on:** none (`_core` stands alone: it cites no other module, so any module can be dropped)

| File | What it is |
|---|---|
| [recurrence.md](recurrence.md) | how work accumulates: a session reads a slice of `main`, works through it, and a merge folds it back; the frame for the rest |
| [system.md](system.md) | how the package is built: modules, kinds, the `_` rule, ids and tombstones, skills, adapting |
| [workflow.md](workflow.md) | how work is run: agent context, design → plan → loop → state → record, iteration, the recursive loop, meeting the author |
| [patterns/drift-triage.md](patterns/drift-triage.md) | drift triage |
| [patterns/workflow-files.md](patterns/workflow-files.md) | the workflow files: designs, plans, loops, state, records |
| [patterns/spec-mirror.md](patterns/spec-mirror.md) | the spec mirror |
| [patterns/adopting.md](patterns/adopting.md) | adopting the package in a project: ask first, create the instance, copy and map, coexist with the host, link skills, wire checks, first lane, fold back ([../bootstrap.md](../bootstrap.md) is the entry point) |
| [patterns/sources.md](patterns/sources.md) | sources: what the principles are drawn from and tested against (a basic guideline) |
