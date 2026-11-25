# Errors with stable codes {#error-codes}

- **Implementations:** [TypeScript](typescript.md)

Rules §stable-error-codes: boundary errors are structured; callers branch on a stable code, never on message text; one shared shape across parsing, repair, validation, queries, commands and ports; codes are cataloged per module.

## Problem

A caller that branches on message text breaks the day someone improves the wording. When each layer invents its own
error shape, a parse error, a query error and a port error cannot travel through the same streams, logs and command
results. Without an actionable hint an agent cannot correct itself, and without a catalog nobody knows which codes a
module can raise.

## Shape

**`ErrorInfo`**, the one shape for every layer:

| Field | Meaning |
|---|---|
| `code` | stable, namespaced: `repr.parse.unexpected_end`, `store.conflict` |
| `message` | human text; free to improve |
| `hint`? | the actionable next step (known names, the command that lists candidates, the expected arity) |
| `at`? | text position: `offset`, `line`, `col` (1-based), optional `length` |
| `loc`? | semantic position, a locator (a reference, as data, to an entity or range in the domain) |
| `path`? | tree position in a representation |
| `details`? | extra structured data, as a representation node |
| `retryable`? | whether retrying may succeed |

- **`CodedError`**: the raised form, carrying an `ErrorInfo`.
- **`coded(ns)(code, message, rest?)`**: lightweight constructor bound to a namespace.
- **`defineErrors(ns, spec) -> { spec, error }`**: a module's catalog. Each entry is `code -> { description, retryable?,
  example?, retired? }`; `error(code, message?, rest?)` builds the error from the entry.
- **`isCoded(e)`**, **`isCancelled(e)`**: classification for callers.
- **`textPos(src, offset) -> at`**: text position from an offset.
- **`Result<T>`**: `{ ok: true, value }` or `{ ok: false, errors: [ErrorInfo, ...] }`, for agent-facing paths that should
  not raise; one entry per problem, each with its `path` or `loc`.
- **Serialized form**, a representation node:
  `(error (code "repr.parse.unexpected_end") (message "Input ended inside a form.") (hint "close open lists with )") (pos 1 6))`,
  followed by the locator, `(path ...)`, `(details ...)` and `(retryable ...)` when present.

**Namespaces follow modules**: `repr.*`, `query.*`, `command.*`, `store.*`, `doc.*`, `data.*`, `locator.*`, `derive.*`;
domain modules add theirs. `command.*` is about commands themselves (`command.not_found` is an unknown command name); a
missing entity uses its own namespace (`doc.not_found`). A namespace may be split across catalogs when the
sub-prefixes do not overlap. The catalogs render into the docs hub as an "Error codes" table per module (§design-system-storybook).

## Invariants

- Every code starts with its catalog's namespace and a dot; declaring or raising one outside it is a programmer error.
- The code given explicitly always wins; `rest` cannot override it.
- From a catalog, the message defaults to the entry's description and `retryable` to the entry's flag.
- A code is never reused or renamed; a retired code stays in its catalog marked `retired`.
- Every raised code is declared in a catalog (the one required check, §cataloged-codes).
- Callers branch on `code` or `retryable`, never on `message`.
- `isCancelled` is true for a coded error whose code ends in `.cancelled` and for the runtime's abort error;
  cancellations are not reported as failures (§cancellation-codes).
- Programmer errors stay plain exceptions; public module boundaries raise only coded errors (§programmer-errors).
- A `Result` turns only coded errors into values; any other exception propagates.
- The serialized form writes a text position as `(pos line col)`, never `(at ...)`, which is a locator's form; optional
  fields appear only when set; `retryable` appears whenever it is set, true or false.
- Messages and details carry no secrets (§no-secrets-in-errors).

## Conformance cases

- `textPos("abc", 0)` → `{offset: 0, line: 1, col: 1}`; `textPos("ab\ncd", 4)` → `{offset: 4, line: 2, col: 2}`.
- `coded("store")("doc.not_found", ...)` → programmer error (outside the namespace).
- `coded("store")("store.conflict", "m", {code: "x.y", retryable: true})` → code `store.conflict`, retryable true.
- `defineErrors("doc", {"store.x": ...})` → programmer error.
- Catalog `doc.not_found: {description: "No document with this id.", retryable: false}`; `error("doc.not_found")` →
  message "No document with this id.", retryable false.
- `{code: "repr.parse.unexpected_end", message: "Input ended inside a form.", hint: "close open lists with )", at: {line: 1, col: 6}}`
  serialized → `(error (code "repr.parse.unexpected_end") (message "Input ended inside a form.") (hint "close open lists with )") (pos 1 6))`.
- A path that raises a coded error, wrapped as a `Result` → `{ok: false, errors: [that error]}`; one that raises an
  uncoded error → the error propagates.
- `isCancelled`: code `x.cancelled` → true; the runtime's abort error → true; code `x.failed` → false; uncoded error → false.
- Strict parse of `(a ]` → `repr.parse.bad_char` at line 1, col 4; strict parse of `(a` → an error with a hint.
- Each catalog entry with an `example` (not retired) → that entry's code: a parse example through the strict parser
  (lenient repair would fix most of them), a validation example through `validate(parse(x))`.

## Open

- §cancellation-codes names `cancelled`, `deadline_exceeded` and `limit_exceeded`; `isCancelled` recognizes only codes
  ending in `.cancelled`, and nothing says whether a deadline or limit counts as a cancellation (not a failure).
- A catalog's `error` accepts a code marked `retired`; whether raising a retired code should be refused is not said.
