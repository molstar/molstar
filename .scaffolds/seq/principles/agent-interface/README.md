# Agent interface

- **Kind:** rule set
- **Status:** stable
- **Scope:** skills (discoverable and searchable, partly generated, small), the tool surface (`<domain>_read`/`_write` taking expressions), agent communication (declared plans, claims, messages) as a queryable domain; outline-first reads, affordances in results, refusals that say how to recover.
- **Sources:** the principles as first written, hardened by the dorothy prototype
- **Depends on:** [_core](../../../principles/_core/README.md), `agent-harness` (the upstream package's, not adopted here), [architecture](../../../principles/architecture/README.md), [queries](../queries/README.md), [representations](../../../principles/representations/README.md), [testing](../../../principles/testing/README.md), [ui](../../../principles/ui/README.md)

The rules are in [rules.md](rules.md). Code patterns: [sketch-patterns/](sketch-patterns/README.md) (§agent-protocol, §skills-and-tools).
