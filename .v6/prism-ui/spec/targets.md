# Targets: how the viewer names what it acts on

One way for the API, commands, tools, the sequence view and URL actions to say which atoms, residues, chains or
structures they mean. Decided with David (2026-10-15): **permissive on the way in, MolQL underneath**, starting from
what Stories Untold already has (its `any-molql` package; in that repository `designs/architecture/selections.md`).
How the agent writes targets is its selections skill's guidance; this file is the contract.

## 1. The type

```ts
type Target = string | SelectorT;
//  string:    a keyword          all, polymer, protein, nucleic, branched, ligand, ion, water, coarse, selection, focus
//             an address list    A, B:120, B:200-210, C:122@OG1   (author numbering; insertion codes as B:120A)
//             a query text       "resn HIS and chain A", "within 5 of resn STI"   (PyMOL by default)
//  SelectorT: schema items       { auth_asym_id: 'A', beg_auth_seq_id: 248, end_auth_seq_id: 256 } or a list
//             { molql }          a MolScript expression as JSON

interface Scoped { target: Target; structure?: StructureRef }   // StructureRef: a state ref, or a label ("1TQN")
```

- **Reading a string**, in order: a keyword, then an address list, then a query text in the given language (the first
  the viewer reads, PyMOL in the default spec, unless the call says `language: 'vmd' | 'jmol' | 'molscript'`), when the
  viewer reads it: its spec lists the languages it bundles (viewer-spec.md §2.7; the default spec all four). The
  `language` param offers only those, and is hidden where there are none. Permissive as Stories reads small models'
  input: spaces, case, `A/120`, `A 120-130`, a trailing comma are all read; what still does not parse fails with the part
  that did not and a hint.
- **`selection` and `focus`** are the person's: Mol\*'s current selection and focus, so "color the selection red" and
  the agent's "what the person points at" are one target.
- **The structure** (one rule for every operation, 2026-10-16; `workingStructure` in `sources/naming.ts`): a
  `structure` given is that one (a state ref, else a label). Absent, it is the one the person works with: the focused
  one, else the one holding the selection when only one does. Past that it depends on what the operation acts on: a
  target applies to every loaded structure (`selection` to every one holding some); an operation on one structure takes
  the only one loaded, else fails `ambiguous` listing them; an operation on many (`EveryStructureParam`) takes every
  one. A destructive operation (`structure.remove`) skips the one worked with: with several loaded, it needs a name.
  The result says which it was.

## 2. Underneath: MolQL

- **Everything compiles to MolQL** (a MolScript expression) before it reaches the engine: keywords to fixed
  expressions (Stories' `KEYWORD_EXPRESSIONS`), address lists to schema items, query texts through Mol\*'s transpilers
  (`@molstar/query-language`: PyMOL, VMD, Jmol), schema items through MVS's element-schema path, which resolves rows
  through the model's indices (`substructureFromRows` in `@molstar/model`'s schema selection, #1972) rather than per atom. One evaluator, one set of
  semantics, whoever asked.
- **Stored as MolQL JSON only** where a target is kept (a component, a snapshot, the history, a run's arguments): no
  PyMOL text in saved state, so a session replays without a parser. A stored query is shown with
  `formatMolScript`, and with what it selected.
- **What it selected** comes back with every resolved target: residues by author number per chain, atom count, chains,
  the structure. Results, the history and the agent's account read this, so "selected 42 residues in chains A and B"
  needs no second evaluation.

## 3. The sequence view and SeqQL

- A Seq\* selection becomes a target through the link (`coordinatesOf`: positions on linked chains to residue
  addresses), so a sequence action names residues as a command does (sequence.md §5).
- SeqQL's `molstar.select` runs a MolQL query on the linked structure and returns positions on Seq\* units: the
  opposite direction, through the same compiler.

## 4. Where it lives

- **In the viewer's package for now** (David, 2026-10-15): `packages/prism/viewer/src/targets/`, ported from Stories'
  `any-molql` (addresses, keywords, selector handling, compile, evaluate, project), DOM-free. It imports
  `@molstar/model` and `@molstar/query-language`, not the plugin. Once it settles, it moves next to
  `@molstar/query-language`, so MVS and Stories share it.
- The viewer's API takes `Target` (or `Scoped`) wherever an operation acts on part of a structure; its tools' fields
  are the same strings, described once in the selections skill.

## 5. Atoms

Atom-level targets use what MVS's component expressions already have:

| Means                                | Schema item                 | Address form (permissive)   |
| ------------------------------------ | --------------------------- | --------------------------- |
| atom name                            | `label_atom_id`, `auth_atom_id` | `C:122@OG1`             |
| atom id (`_atom_site.id`)            | `atom_id`                   | `id:1234`, `@#1234`         |
| source index (0-based, in the file)  | `atom_index`                | `index:1233`                |
| element                              | `type_symbol`               | `element:FE`, `C:122 element:S` |
| instance (assembly or crystal copy)  | `instance_id`               | `A[ASM-X0-1]:120`           |

Alternate locations have no schema field: a MolQL query (`altloc A` in PyMOL text) until a field is asked for.
