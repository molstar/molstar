# Representations {#representations}

- **Implementations:** [TypeScript](typescript.md)

Rules §declarative-representations, §parser-repairs-slips: state is one declarative, recursive, streamable tree with a strict parser and printer in `repr-core`; lenient repair only at agent-facing edges; executors are bound to state in the composition root.

## Problem

When state is a web of objects with behavior, nothing can print, diff, stream or store it, or show it to an agent,
without per-type code, and whatever renders or runs it drifts into being the truth. When the text spelling has no
single strict grammar, every consumer guesses, and an agent's small slip (a missing paren, a code fence) becomes an
ambiguous parse instead of a reported repair.

## Shape

**The tree.** A `Node` is an atom, a symbol or a list of nodes. By convention a list's head symbol is its tag.

| Kind | Holds | Spelling |
|---|---|---|
| string | text data | JSON string syntax: `"a\nb"` |
| number | a float64 | shortest round-trip digits; `-0`, `#nan`, `#inf`, `#-inf` |
| bigint | an arbitrary-precision integer (ids, counters, 64-bit values) | `123n` |
| boolean, null | | `true`, `false`, `nil` |
| symbol (`sym`) | a tag or name, never data | bare: `board`, `x-y`, `<=` |
| list | children in order | `(board (view (at 0 2 5)))` |

In the examples, `(view (at 0 2 5))` is a board's viewport: `at` holds x, y and zoom.

- A symbol name must print and re-parse as itself: no whitespace or delimiters `( ) [ ] { } " ; ,`, no leading digit
  or sign-digit, no leading `#`, not `nil`/`true`/`false`, no smart quotes. `sym(name)` refuses anything else with
  `repr.bad_symbol`.
- A number prints with shortest round-trip digits and, where an exponent is needed, a signed `e` exponent (`1e+21`,
  `5e-324`); that is the canonical printed form (the codec's and the reference prototype's), so other runtimes print
  the same.
- Exact decimals are a tagged form holding a string, `(decimal "12.50")`; arithmetic belongs to a library, not the tree.
- Large numeric arrays are data references (§bulk-data-not-nodes): the tree holds a reference, the bytes live elsewhere;
  lists of homogeneous records are tables (a vocabulary the project defines). Other encodings (JSON, binary) are codecs
  of the same tree, defined by the project's wire and storage layer.

**Addressing.** `Path` is a list of child indices from the root; `nodeAt(root, path)` returns the node or nothing.
`attrOf(n, key)` is the first value of the child list tagged `key`, `attrsOf(n, key)` all its values: on
`(view (at 0 2 5))`, `attrOf(_, at)` is `0` and `attrsOf(_, at)` is `0 2 5`. `walk(n, path, visit)` is a
pre-order walk; `visit` returning false skips that node's children.

**Print and parse.** `print(node) -> text` and `parse(text) -> node | coded error` are the one text spelling. This pattern spells the tree as S-expressions, the example; a project that picks
another spelling keeps the same tree, nodes and invariants, and writes its own `print` and `parse`.
`( ) "` delimit; `[ ] { } ; ,` are reserved and rejected by the strict parser, so repair can treat `;` and `,` as
slips. Exactly one top-level form is allowed. Errors carry a position (offset, line, column) and a hint:

| Code | Example input |
|---|---|
| `repr.parse.unexpected_end` | `(a (b` |
| `repr.parse.bad_string` | `("abc` |
| `repr.parse.bad_char` | `(a ]`, `(a, b)` |
| `repr.parse.bad_token` | `(1abc)` |
| `repr.parse.trailing_input` | `(a) (b)` |

**Lenient parsing**, for agent and user input only. `repair(text) -> { text, repairs }` fixes mechanical, local slips
and reports each one; `parseLenient(text) -> { node, repairs, normalized }` repairs, then runs the strict parser, so
the strict grammar stays the one definition. Repair kinds: `stripped-fence`, `smart-quotes`, `dropped-comment`,
`commas-as-space`, `escaped-newline`, `closed-string`, `closed-parens` (the riskiest: reported with a count, and the
normalized form is echoed back) and `dropped-extra-close`. Not repaired: unknown tags, wrong arity, several top-level
forms, mismatched structure in the middle; those stay coded errors. Commands echo `{ normalized, repairs }` (§commands).

**Streaming.** `parseStream(chunks) -> stream of Node` emits each complete top-level form (list, string or atom) as soon
as it ends; a top-level atom ends at the next delimiter or at the end of input.

**State, snapshot and patches.** `TreeState` is the observable current root (a behavior, §api-surfaces); unchanged
subtrees stay the same node (structural sharing). A positional `Patch` (`replace path value`, `insert path value`,
`remove path`) exists only to stream and mirror a tree over the wire (§streamable): `diff(a, b) -> patches`,
`apply(root, patches) -> root | repr.mirror.bad_path`; a sender emits the patches between consecutive roots, a receiver
starts from a snapshot and applies them. Edits are never patches: they are id-addressed `Op`s
(§committed-transient-identity), which survive concurrent edits; positions do not.

**Executors.** An `Executor` port interprets a description: `reconcile(state, hint{prev, changed paths})` and
`dispose()`. `bindExecutor(state, executor, report)` binds one in the composition root (§executors-by-composition): it
reconciles the current state, then one reconcile at a time, latest wins (intermediate states skipped), and reports
errors instead of throwing. Disposing the executor is its creator's job. State never knows its executors; hit-testing is
a screen-to-doc mapping, not an executor method. Agents need no executor: they inspect state through
commands, which take and return representations (§commands-on-representations): read commands return printed nodes.

**Domain modules.** One module per domain representation (vocabulary, migrations, alternative encodings, query
functions) depends only on the core; executors and apps depend on it, never the reverse (§representation-modules).
Migrations are pure `Node -> Node` functions applied in order from a version (§versioned-representations).

## Invariants

- `parse(print(x))` equals `x` for every node, including `-0`, NaN, both infinities and bigints.
- A symbol that would not re-parse as itself cannot be constructed, so `print` only produces text `parse` accepts.
- The strict parser accepts exactly one form and rejects reserved characters; every rejection is a coded error with a
  position.
- `repair(print(x))` returns `print(x)` with no repairs; every change `repair` makes is reported.
- Streaming parse equals whole parse for any chunking of the same text.
- `apply(a, diff(a, b))` equals `b`; `diff(a, a)` is empty.
- Committed changes never travel as positional patches.
- An executor never runs two reconciles at once, and a failing reconcile does not stop later ones.
- State holds no behavior and references no executor.

## Conformance cases

- Numbers: `-0` parses to negative zero; `#nan`, `#inf`, `#-inf` to NaN and the infinities; `123n` to a bigint and
  `1n` prints as `1n`; `1e21` prints as `1e+21`; `5e-324` round-trips; `1` parses to the number 1.
- `("a\nb" nil true false ())` parses to a string with a newline, null, true, false and an empty list.
- Symbols: `a`, `x-y`, `k.v`, `+`, `<=`, `a1` are valid; the empty name, `1a`, `-1`, `+.5`, `nil`, `true`, `false`,
  `a b`, `a(b`, `#x`, `a;b`, `“q”` are not; `sym("1a")` fails with `repr.bad_symbol`.
- On `(board (view (at 0 2 5)) (shape "t"))`: the tag is `board`; the tag of `"x"` and of `()` is nothing;
  `nodeAt [2 1]` is `"t"`, `nodeAt [9 1]` is nothing; `attrOf(view, nope)` is nothing and `attrsOf(view, nope)` is
  empty; a walk that skips `view`'s children visits paths `()`, `0`, `1`, `2`, `2 0`, `2 1`.
- Each row of the parse table yields its code; `(a` newline ` (b ]` reports offset 7, line 2, column 5.
- One case per repair kind: a ```` ```lisp ```` fence around `(a)` → `stripped-fence`; `(a “b c”)` → `smart-quotes`;
  `(a) ; note` → `dropped-comment`; `(a, b)` → `commas-as-space`; a raw newline inside a string → `escaped-newline`;
  `(a "x` → `closed-string`; `(a (b` → `closed-parens` with detail `2 at end`; `(a))` → `dropped-extra-close` with
  detail `1 unmatched`.
- `parseLenient` on a fenced `(a, b (c` gives normalized `(a b (c))` and repairs `stripped-fence`, `commas-as-space`,
  `closed-parens`; `(a) (b)` still fails with `repr.parse.trailing_input` and `(1abc)` with `repr.parse.bad_token`.
- Streaming: the chunk `(a b) (c` emits `(a b)` at once; ` d)` then emits `(c d)`; `(a (b` at the end fails with
  `repr.parse.unexpected_end`; `a)` fails with `repr.parse.bad_char`.
- Property tests over generated trees of every atom kind: round trip, repair identity, streaming equals whole parse
  for chunk sizes 1 to 7, `apply(a, diff(a, b))` equals `b`.
- A patch whose path does not exist fails with `repr.mirror.bad_path`; a mirror fed a snapshot and the sender's
  patches stays equal to the sender.
- `bindExecutor` reconciles the current state first, then the latest of several quick changes, one at a time, and
  reports a reconcile that throws.

## Open

- The JSON codec's escape spelling (`{"$big":"123"}`, `{"$num":"#nan"}`) is not fixed here: a project decides it in its wire
  and storage layer.
