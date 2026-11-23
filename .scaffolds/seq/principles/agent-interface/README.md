# Agent interface

- **Kind:** rule set
- **Status:** stable
- **Scope:** skills (discoverable and searchable, partly generated, small), the tool surface (`<domain>_read`/`_write` taking expressions), agent communication (declared plans, claims, messages) as a queryable domain; outline-first reads, affordances in results, refusals that say how to recover.
- **Sources:** the principles as first written, hardened by the dorothy prototype
- **Depends on:** [_core](../../oz/_core/README.md), `agent-harness` (oz, not adopted here), [architecture](../../oz/architecture/README.md), [queries](../queries/README.md), [representations](../../oz/representations/README.md), [testing](../../oz/testing/README.md), [ui](../../oz/ui/README.md)

The rules are in [rules.md](rules.md). Code patterns: [sketch-patterns/](sketch-patterns/README.md) (§agent-protocol, §skills-and-tools).
