# Architecture

- **Kind:** rule set
- **Status:** stable
- **Scope:** layout (lib, product packages, apps, services), provider runtime packages, ports, API surfaces (sync/async/streaming), composition root and DI, error codes, commands, performance, access enforcement with application policies, identity and sessions, one consent model for agents' unrequested effects. Not the folder tree itself (layout) or the tools (workspace).
- **Sources:** the principles as first written, hardened by a reference prototype in the package's home
- **Depends on:** [_core](../_core/README.md), [representations](../representations/README.md), [testing](../testing/README.md), [ui](../ui/README.md), `workspace` (not adopted here)

The rules are in [rules.md](rules.md). Code patterns: [bootstrap-patterns/](bootstrap-patterns/README.md) (§api-surfaces, §commands, §error-codes, §hot-loop-allocations, §ports-and-adapters).
