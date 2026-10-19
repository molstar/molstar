# Queries

`src/model/structure/query/` and `src/script/`: selecting parts of a structure.

## 1. The engine

- `StructureQuery = (ctx: QueryContext) => StructureSelection`; `StructureQuery.run(q, structure)` makes the context.
- `QueryContext` (`context.ts`): the input structure, the current `element` Location, the current structure and
  selection, bond info for bond tests, and an optional timeout (`timeoutMs`, 0 for none). The stacks are pushed and
  popped by the queries that set them.
- `StructureSelection` (`selection.ts`): `Singletons` (every element its own selection, one structure holds them all) or
  `Sequence` (a list of structures). `toLociWithSourceUnits` and `unionStructure` turn one into a loci or a structure;
  `LinearBuilder` and `UniqueBuilder` collect results (all-singleton results collapse to `Singletons`).
- Queries (`Queries` in `query.ts`):
  - generators: `all`, `none`, `atoms` (atom, residue, chain, entity tests, optional grouping), `rings`,
    `querySelection`, `bondedAtomicPairs`, `elementSchema`;
  - filters: `pick`, `first`, `withSameAtomProperties`, `areIntersectedBy`, `within`, `isConnectedTo`;
  - modifiers: `wholeResidues`, `includeSurroundings`, `includeConnected`, `surroundingLigands`, `intersectBy`,
    `exceptBy`, `union`, `expandProperty`, `querySelection`;
  - combinators: `merge`, `intersect`; plus `internal` (polymer, water, het, spheres, bundles) and predicates.

## 2. MolScript and languages

- Expressions are `@molstar/query-language`'s MolScript. `script/runtime/query/base.ts` compiles one to a `QueryFn`:
  literals and all-constant applications of constant symbols are evaluated at compile time; an unknown head throws.
- Symbol implementations register in `QueryRuntimeTable` as an import side effect of `table.ts`
  (`runtime/query/compiler.ts` imports it); registrations are counted, so the shared `DefaultQueryRuntimeTable` stays
  valid across plugins. Code that imports `compile` from `base` relies on someone having imported the table.
- `Script` (`script/script.ts`): `{ expression, language }` for `mol-script`, `pymol`, `vmd`, `jmol` (transpiled by
  `@molstar/query-language`); `toQuery`, `toLoci`, `getStructureSelection`. `script/util/id-list.ts` compiles id lists
  (auth, label, atom id, element symbol).
