# Queries and semantics: sketch patterns

Directional: the shape is fixed, the details are not; harden them through a design and a prototype, and report drift. What stable and sketch mean and how to use a pattern: [../../_core/system.md](../../../../principles/_core/system.md#code-patterns).

§semantic-layer has a "Hardened by a prototype" section: read it with the sketch.

| Pattern | Name |
|---|---|
| [§core-concepts-and-embedding](core-concepts-and-embedding/README.md) | Core concepts and embedded languages |
| [§query-dsls](query-dsls/README.md) | Query DSLs over representations |
| [§semantic-layer](semantic-layer/README.md) | Semantic layer: concepts, spaces, mappings |

## Shared names

The names these patterns own (one owner and one name per concept; see `_core/system.md`, "Shared names").

| Concept | Name | Owner |
|---|---|---|
| Queries | `QueryLang`, `QueryEngine`, `QueryContext`, `QueryCache`, `Selection`/`SelOp`, `Dialect`, `resolveQuery` | §query-dsls |
| Meaning and navigation | `Locator`, `SpaceSpec`, `Mapping`, `SemanticGraph` | §semantic-layer |
| Languages | `LanguageRegistry`, `Language`, `FormSpec`, `Slot`, core concepts | §core-concepts-and-embedding |
