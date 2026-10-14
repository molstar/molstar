# SeqQL: what a query means

SeqQL is a query language over the model, after Mol\*'s MolQL. A query is an expression; compiled, it is a
function from a context to a **sequence selection**. The language (table, reader, binder, printer,
checker) is `@molstar/seq-core/query`, a leaf; the runtime is `@molstar/seq-core/model/query`. Every symbol, its signature,
arguments, defaults, examples and hints is in [reference/seqql.md](reference/seqql.md), generated from
`core/src/query/table.ts`; the agent's skill is [reference/seqql-skill.md](reference/seqql-skill.md).

```lisp
(pick (items :of domains) :where (matches (label) kinase))   ; kinase domains
(onto (sel "1brs/A[10-20]") msa)                              ; their alignment columns
(and 1brs/A (partners 1brs/D :via contacts))                  ; residues of A in contact with D
```

## Overview

- **Generators** make selections from the data (a ref, `all`, `items`, `values`, `letters`, `search`, an
  input); **modifiers** map selections to selections (`range`, `slice`, `expand`, `onto`, `group`,
  `partners`, `each`); **filters** keep groups (`pick`, `first`, `ids`); `and`, `or`, `not` combine
  selections as they combine tests; **properties** read an element inside a test.
- One tree, three spellings: the S-expression (people, agents, links, tests), JSON with full names
  (stored), and the builder `Q` (code).
- A query is deterministic: the same expression on the same set gives the same selection in the same order.

## The context

A query reads a **set** (`SeqSet`), never the document or the view, so it gives the same answer in a
page, a CLI, a server and Nightingale's load. The host adds named **input selections** (`selection`,
`focus`, `structure`), limits, an abort signal and an optional `resolveRef` hook for its own names (a
Context entity's id). Refs are the set's keys. What depends on the view (closure, granularity,
projection) applies after the query, to its loci.

## Elements and groups

- **Elements** are positions in a coordinate system of kind `unit`, `molecule`, `row` or `columns`, and
  annotation **items**, selected with their extent. A complex names its units; a molecule names every
  alive unit over it (itself when it has none). Boundaries and gap cells are not elements: a gap is
  selected as its column.
- A **selection** is a sequence of **groups**, each a set of runs in coordinate systems. A group carries
  its item while its runs are exactly that item's extent: modifiers (`range`, `slice`, `expand`, `onto`,
  `group`, `partners`) and `and`, `not` return untagged groups; filters and `or` keep the tag.
- Groups are what filters test, what `each` binds, what results step through and what becomes one item
  when a selection is saved. A ref to a system gives one group per system; `items` one per item; scans
  one per maximal run; `search` one per hit.

## Slots

A **slot** is the element a test reads, bound by the `:where` it sits in. Each binder makes some slots
readable, in order (its **slot chain**):

| Binder | Its element | Slot chain |
|---|---|---|
| `positions :where` | a position | position, then its system |
| `systems :where`, `rows :where` | a system | system |
| `items :where`, `partners :where` | an item (a pair's, for `partners`) | item |
| `pick :where` | a group | group, its item, its system (when it has exactly one) |

A selection inside a test (`(overlaps Q)`'s `Q`) reads no slot: it is evaluated once per evaluation of the
filter, never per element.

## Types

Checked before evaluation: `query` (a selection), `test` (a `:where`), `bool`, `num`, `int`, `str`,
`json`, `ref`, `label` and enumerations. A word or string where a selection goes is a ref; a number there is an error that says to
quote it (`"123"`). A `:where` takes a `bool` expression read per element, and reads only its binder's
slot chain. Missing values are two-valued: a comparison with a missing value is false, and `not` of it is
true.

## Names

Every symbol has a **full name** (`seq.item.label`), used in the table, the JSON and the builder, and
usually a **short name** (`label`), used in the S-expression. A full name is always accepted. A short
name is resolved where it is written:

1. Every symbol with that short name is a candidate; none is `seqql.symbol.unknown` with suggestions.
2. For `and`, `or`, `not`: all operands selections, the combinators; all tests, the logic; mixed,
   `seqql.type` ("filter with `(pick Q :where T)`").
3. For a property: the first slot along the innermost binder's chain with a candidate whose signature
   fits the arguments given; none, `seqql.slot`, naming where the name works.
4. Otherwise the one candidate.

So in `pick`, `label` is the group's item's label. Resolution is unambiguous by construction:
registration refuses two symbols the steps cannot tell apart.

## Labels and indices

Exactly two ways to point at positions, never mixed in one form:

- **Labels**, for what a person says: the system's default numbering or a named one, inclusive:
  `(range 1brs/A 10-20)`, `1brs/A[10-20]`. A label range covering nothing is a warning. A namespace alone
  names the one numbering under it: `uniprot` is `uniprot:P04637`, as loaders name it; with several, the
  query is an error that lists them (`seqql.numbering.ambiguous`), and a full name picks one,
  `R[uniprot:P04637-2: 120-140]`.
- **Indices**, for what a tool returned: 0-based, half-open, in the system's own index space:
  `(slice 1brs/A 9 20)`, `1brs/A@9:20`. A slice past the domain is clipped with a warning
  (`seqql.slice.clipped`); wholly outside, an error (`seqql.slice.range`).

In the JSON a `range` keeps its labels and a `slice` its indices.

**Cyclic systems.** On a system that is a whole cyclic molecule (the molecule, or a unit over all of it), a
`range` with `beg` after `end`, a `slice` with its start after its stop, and an `expand` past either end run
across the origin, as the format's ranges do (`seqspec:format#positions`): `(range ring 38-3)` is two runs,
38–40 and 1–3. Elsewhere a slice's start after its stop is `seqql.slice.range`, as equal ends are anywhere (they
would select nothing); `sel`'s `R@37:3` is the same slice.

## Combining

- `(and Q Q…)`: each group of the first, intersected with each other operand in turn (each operand's
  groups unioned). `(and A (not B))` removes `B` from each of `A`'s groups. Set algebra is how
  annotations are joined: "positions of 1brs/A inside any domain" is `(and 1brs/A domains)`.
- `(or Q Q…)`: every group of every operand, in order, a group equal to an earlier one dropped.
- `(not Q :in U)`: one group, `U`'s positions (default everything) not in `Q`.
- A scan that is a direct operand of `and` after the first, with no `:in` of its own
  (`(and 1brs/A (letters C))`), scans within the earlier operands' runs only.

## Onto

`(onto Q R)` is each group's image in other coordinate systems: `R` names them as a ref does (a unit, a
row, a complex's units, a molecule's units, an alignment's columns, or its rows with `:kind row`). Images
go through the model's projection rules (`seqspec:model#projection`), including codon steps and partial
maps. A group with no image, or one reached ambiguously, is skipped with a warning.

## Partners

`(partners Q :via A :side S :where T)` returns one side of a pair by its partner: for every pair of `A`
passing `T`, the end whose other end overlaps `Q`. `:side` limits which ends may be returned: `any`
(the default), `a` or `b`, the pair's two column sides. It is MolQL's interface query in domain-free terms: contacts, hydrogen bonds,
base pairs and crosslinks are `pairs` annotations told apart by category.

## Each

`(each Q E)` runs `E` once per group of `Q` with `(current)` bound to it, concatenating the results: the
query monad's bind. Any per-group behavior ("the first ten residues of every domain") is a symbol over
`(current)`, not a new runtime concept. `(current)` outside `each` is an error.

## Order

Groups keep generation order: systems in canonical order (kind, then id, then row), then by start; items
in annotation, then item order; hits in search order; `or` in operand order. Runs within a group are
in canonical system order, then by start, disjoint, merged where they touch.

## Lazy data

A query reads what is loaded, as hover and selection do. Positions of a lazy system exist (its length is
known); its letters and values do not, and each scan skips it with one `seqql.lazy` warning naming how many it skipped (a streamed values
annotation, the same).

## The selection

A `SeqSelection` is groups of runs in typed arrays (system ids within the set's lineage, starts, stops),
each group optionally tagged with its item and a score. It is immutable once built.

- `summary` is the line a person or the agent reads: groups, positions and runs in each system's default
  labels, written so a label range copies into `range`; with `indices: true`, each run's slice too.
- `toLoci` gives the model's loci (unit runs, molecule runs as every alive unit's, rows, columns, items),
  needing no view. In the viewer, loci go to the channels, which close, widen and project them as any
  loci (`seqspec:interaction#channels`).
- `toExtents` feeds focus and the Mol\* link; `toRangesInput` gives a `ranges` annotation, one item per
  group: saving a selection as a location set.

## The reader

The input is one form. Lists `(head arg…)`, whose head names a symbol; JSON strings; JSON numbers; `true`
and `false`; keywords `:name`; and words (any other run of characters, a string without quotes). `;`
starts a comment. `[`, `]`, `,` outside a string are errors pointing at `sel`; a word starting with `$`
points at `(input name)`. Refs are words: `1brs/A`, `1brs/entity:1`, `msa/3`, `ctx.site.hinge`;
quoting such a word changes nothing, but an unquoted `123`, `true` or `false` is a number or boolean,
and a quoted head is an error. Positionals come before keywords. Depth 64 and 64 KB at most; the reader is iterative.
A stored expression, which never passed the reader, is held to the same depth before it is checked
(`seqql.read.depth`).

## Binding and checking

Reading gives a tree; **binding** resolves names (Names) and names positionals, giving the JSON;
checking expands macros and applies the types, slots and arities. Every issue past the reader is at a symbol, with its signature and
the symbol's hint. Macros (`sel`, an application's words) are expanded for checking and compiling and
kept as written in the JSON.

## The printer

`printQuery` gives the canonical text stored beside the JSON: short names wherever they resolve back to
the same symbol, positionals then keywords in signature order, defaults omitted, a string as a word when
it reads back as the same word, one line up to 80 columns, else the head and the leading positionals that fit, then one argument a
line indented two spaces.
Binding a print gives the expression back; printing is a fixed point. In a document, a stored `text` that differs
from `print(expr)` is `svs.query.stale` (a warning), and the `expr` wins.

## The sel shorthand

`(sel "…")` is a selector shorthand for the common case, usable wherever a selection is, expanded when the
query is checked, with no data:

| Shorthand | Expands to |
|---|---|
| `R[10-20,45]` | `(range R 10-20 45)` |
| `R[uniprot: 120-140]` | `(range R 120-140 :numbering uniprot)` |
| `R@9:20`, `R@9` | `(slice R 9 20)`, `(slice R 9 10)` |
| `R[#PF00069,#PF07714]` | `(ids R PF00069 PF07714)` |
| `$selection[10-20]` | `(range (input selection) 10-20)` |
| `t1, t2` | `(or t1 t2)` |

It has no `and`, `not` or tests. The JSON keeps `sel` with its normalized string, so the shorthand round
trips; its errors point into the string.

## Queries in text

In free text (a chat reply), each query is wrapped `[SeqQL]…[/SeqQL]`; several blocks are several queries
and text around them is fine. Without tags, a fenced block with the info string `seqql` is read; without
either, the whole text when it is one form. Where a query is already a value (a tool's argument, a
query box, a document) it is the bare S-expression. Command and link targets that take a query are not
built. Issues point into the whole text.

## Extension

A host assembles its language from entries, as the viewer assembles its registry; there is no module
state.

- An extension's full names are `<ns>.<name>` (`molstar.select`) outside `core`, `seq`, `svs` and
  `seqv` (`seqql.symbol.namespace`); a symbol defined twice is `seqql.symbol.duplicate`.
- Short names are its choice; one that is a base short name or that two symbols cannot be told apart by
  is refused, and the symbols keep their full names.
- A definition is joined to its runtime by the definition object: a definition without one checks but
  does not compile (`seqql.symbol.unavailable`).
- **Domain meaning stays in applications**: every base symbol is stated without a biological word. Meaning
  enters as data (semantic types, categories, fields, letter-class tables), macros (Stories' `chain`) and
  refs the host resolves.

## Errors

Every issue is `{code, severity, message, path}`, placed at a line and column of the text (runtime issues
too), rendered as `line:column code: message` and the symbol's signature, with did-you-mean suggestions
(edit distance at most 2, or 3 for long names) over what is in scope. Codes are stable: `seqql.read.*`
(the reader), `seqql.symbol.*`, `seqql.arg.*`, `seqql.type`, `seqql.slot`, `seqql.sel.*` (binding and
checking), `seqql.ref.*`, `seqql.limit.*`, `seqql.onto.*` and the rest (running), `seqql.reply.*`
(extraction). Each code has one severity, listed in the reference's Issue codes: an error is a query written
wrong or one that cannot run, a warning a part of what was written that selects nothing or a result that may
be incomplete, and `info` a note on how a valid query ran (why its result is empty). Messages show positions as
labels and add the slice where indices matter.

## Empty results

An empty result is not an error, but "selects nothing" alone teaches nothing. When a whole query selects
nothing, the runtime names the first node, in evaluation order, that produced nothing from a non-empty
input, in one `seqql.empty` note (`info`) saying what it knows (a value range, the positions read, how many
items a `:where` tested and which value was missing). A node whose own issue already says why gets no
second one.

## In documents

A query lives in a document as a `query` computed node (`seqspec:format#queries`): its result is a
`ranges` annotation, one item per group (its item's id, else `q<k>`), so layers, highlights, focus and
filters take it.

- **Data only**: `input` and `current` outside `each` are `svs.query.input`; `search` is
  `svs.query.search` (save the hits as an annotation, then name it).
- Its refs (after expanding `sel` and macros) are the node's dependencies: an unknown one is
  `svs.ref.unresolved`, a cycle `svs.ref.cycle`. Query nodes run after the other providers.
- The checker's errors are `svs.query.invalid`, its warnings `svs.query.warning` and its notes `svs.query.note`,
  each message leading with the SeqQL code; a stale `text` is `svs.query.stale`. A run's issues are the same, but
  an error of the environment (a limit, a symbol with no runtime) leaves the node unavailable
  (`svs.provider.unavailable`), and a cancel is none.

## Versioning

SeqQL has its own semver, written in a stored query. Minor versions add symbols, short names, optional
keywords and macros; anything else, a change to `sel`'s expansion, to a short name's resolution or to the
canonical print included, is major. A reader of `1.x` reads any `1.y`, but an unknown symbol or keyword is
an error: a query that silently loses a generator selects something else. Stored JSON uses full names, so
a new short name never changes a stored query.

## Performance

Queries run at edge rate (a click, a command, a load), never per frame. The rules:

- Runs, not positions: every generator emits runs and every set operation works on runs; positions are
  visited only by scans, the general `positions` test and `partners`.
- Generators emit in canonical order, so `and`, `or`, `not` are linear merges per system.
- No allocation per position or per run; a run allocates its result, its scratch and one closure tree.
  A text property compared with constants (`(= (label) "27A")`, `in`, `!=`) is compared as the numbering
  holds it, never as text. `(label)` compared with anything but a constant (`(= (label) (label :numbering
  auth))`), or under `matches` (a regular expression reads only strings), makes its text, a string per
  position.
- **Limits** (context options): positions visited by scans and tests (50,000,000), runs and groups in
  any intermediate selection (1,000,000 each), each reported with a narrowing hint; and a time limit
  (none by default; checked every 1,024 runs and per system). The time limit and an abort signal say
  only that the run stopped.
