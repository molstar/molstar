# Targets

How operations name the atoms, residues, chains or structures they act on: permissive in, MolQL underneath. The code is
the truth for how (`src/targets/`, `sources/targets.ts`, `api/target-param.ts`).

## 1. The type

- **`Target = string | SelectorT`**, `Scoped = { target, structure? }` (`targets/types.ts`). A string reads as a keyword
  (`targets/keywords.ts`; `selection`, `focus` are the person's), else an address list (§5, author numbering,
  permissive: `A/120`, `A 120-130`, a trailing comma), else query text; else fails with the part that did not read and a
  hint. `SelectorT` is MVS's: schema items or `{ molql }`; code in the page may pass `{ loci }`, never stored.
- **Query languages are the spec's `languages`** (viewer-spec.md §2.7): only those transpilers are bundled, the first is
  the default (PyMOL in the default spec), `language` picks another and is hidden when none is registered.
- **The structure, one rule** (`workingStructure` in `sources/naming.ts`, `api/scope.ts`): a named one (state ref, else
  entry or label); else the focused one, else the only one holding the selection; else a target applies to every loaded
  one (`selection`: every one holding some), a one-structure operation takes the only one or fails `ambiguous`,
  `EveryStructureParam` takes all. A destructive operation (`structure.remove`) skips focus and selection: with several
  loaded it needs a name. The result says which rule chose (`ResolvedTarget.scope`).

## 2. Underneath: MolQL

- **Everything compiles to MolQL**: keywords to fixed expressions, addresses and schema items through
  `StructureElement.Schema.toExpression`, text through `@molstar/query-language`. One evaluator (`targets/evaluate.ts`).
- **Kept only as MolQL JSON** (`toMolqlJson`): components, snapshots, the history and run arguments hold `{ molql }`, so
  a session replays without a parser.
- **What it selected** (`Selected`: atoms, residues, chains, addresses, one line) comes with every resolved target; no
  second evaluation. A coarse element counts by its label chain and the residues it covers; a molecule is an entity.

## 5. Atoms

Address forms (`targets/addresses.ts`): `C:122@OG1` or `atom:OG1` (`auth_atom_id`), `id:1234` or `@#1234` (`atom_id`),
`index:1233` (`atom_index`, 0-based), `element:FE` (`type_symbol`), `A[ASM-X0-1]:120` (`instance_id`). Alternate
locations: query text only.
