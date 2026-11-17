# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.0.0` until the first slice runs). Why: one language, one
  description, one number for readers to check.

## Language

`lang/`: `SeqLanguage`, the description of SeqViewSpec on `@molstar/core/repr` (kinds, parents, params, docs), and the
checks it cannot express.

- **Kinds** (first slice): `spec` (root: `version`, `title`); the containers `data` and `views` (reserved ids of the
  same name); under `data`: `source`, `sequence`, `ranges`, `values`, `map`; under `views`: `view`, then `section`,
  `track`, `layer` (ids optional on a layer only). Reference params (`target`, `of`, `from`, `to`, `follow`, `source`,
  `data`) are `ref` params naming the kinds they take; a reference is an address.
- **Positions** are 0-based, ranges half-open. A range with `start > end` wraps the origin where its sequence is
  `cyclic` and is an error elsewhere; a section's `range` follows the same rule. Strands are `+` or `-`. Exons of a
  codon map are in genomic order and do not overlap (`strand` says which way the protein runs). Why: one convention for
  every coordinate the language writes.
- **Gate order of `checkSeqSpec`**: (1) the version: a newer major is `seq.lang.version.unsupported` alone, since this
  version's description says nothing of a newer spec; (2) the description's validation with normalization (`repr.*`,
  columns typed, so a spec read from JSON is checked as one built); (3) the language's checks (`seq.lang.*`), only when
  step 2 found no error, so they read a valid, normalized tree. Why: each gate assumes what the one before settled.
- **What the checks cover**: one of two params (`letters` or `source`, `rows` or `source`, `values` or `source`, `uri`
  or `value`), an `index` only for a file, letters of the alphabet, row intervals, bounds and strands (the issue's path
  carries the row), what a map method needs, exons, a codon map's alphabets, section ranges, that a reference names a
  node of a kind its param takes, that a layer's data is what its type draws (reported once: a kind the param does not
  take is the reference's issue alone), that `follow` chains end. Extents are checked where the letters are inline.
- **Open**: the theme slot of a layer's `color` and a query in its `data` wait for the sublanguages (TODO).

## Language builder

`createSeqSpec(options)`: the way most specs are written; `build()` is the normalized tree, `check()` is `checkSeqSpec`
of it.

- **Handles are addresses, branded by kind**: `Address<'sequence'>` is the string `/data/p53` that says what it names.
  Every call that adds a node with an id returns its address; a `layer` without an id has none and returns `undefined`.
  Section, track and view builders carry their `address` and a `toString`, and stand where an address of their kind
  does.
- **A reference param takes** a handle (or builder) of a kind its description allows, or an address written out (any
  string); a handle of another kind is a type error. Which params are references is read from the description, in the
  types and when the builder converts a handle to its address. Why: no list of param names to keep in step with the
  language.
- **Misuse** that types do not reach (an id twice under one parent, a malformed id) throws where it is written; the rest
  is found by `check`. Calls are recorded and replayed on the engine of `@molstar/core/repr/builder` at `build`.
- **Type tests** (`_test/types.test.ts`) mark wrong kinds and params with `@ts-expect-error`; they are checked by
  `tsc -p tsconfig.test.json`, which the package build does not run.

## Issues

`seq.lang.*` codes (`LangIssues`), each with a severity, a description and a hint; a test fails on a code a source uses
and the catalog lacks, and on one cataloged and never used, and a misuse for each code is in `_test/lang.test.ts`. The
namespaces are `repr.*` (the tree core's), `seq.lang.*`, and the later `seq.resolve.*`, `seq.source.*`, `seq.store.*`,
`seq.view.*`, `seq.render.*`, `seq.address.*`. An issue's `path` is the step path and its `address` the node's, where it
has one.
