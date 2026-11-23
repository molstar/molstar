# Agent interface: sketch patterns

Directional: the shape is fixed, the details are not; harden them through a design and a prototype, and report drift. What stable and sketch mean and how to use a pattern: [../../_core/system.md](../../../oz/_core/system.md#code-patterns).

§skills-and-tools has a "Hardened by a prototype" section: read it with the sketch.

| Pattern | Name |
|---|---|
| [§agent-protocol](agent-protocol/README.md) | Agent protocol: plans, claims, messages |
| [§skills-and-tools](skills-and-tools/README.md) | Skills and tools: the agent interface |

## Shared names

The names these patterns own (one owner and one name per concept; see `_core/system.md`, "Shared names").

| Concept | Name | Owner |
|---|---|---|
| Agent interface | `tools`, `ToolSession`, `domainRead`/`domainWrite`, `readExpr`, `help` | §skills-and-tools |
| Agent coordination | `agent` domain (`plan`, `step`, `claim`, `msg`), `agentQuery`, `agentErrors` | §agent-protocol |
