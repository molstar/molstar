# Errors with stable codes: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## The shape (common library)

```ts
export interface ErrorInfo {
  code: string;                       // stable, namespaced: "repr.parse.unexpected_end"
  message: string;                    // human text; free to improve
  hint?: string;                      // actionable next step
  at?: { offset: number; line: number; col: number; length?: number };   // text position (parse errors)
  loc?: Locator;                      // semantic position, a locator: which entity, item, text range
  path?: Path;                        // tree position (§representations)
  details?: Node;                     // extra structured data, as a representation
  retryable?: boolean;
}

export class CodedError extends Error implements ErrorInfo {
  readonly code!: string; readonly hint?: string; readonly at?: ErrorInfo["at"]; readonly loc?: Locator;
  readonly path?: Path; readonly details?: Node; readonly retryable?: boolean;
  constructor(info: ErrorInfo) {
    super(info.message);
    this.name = "CodedError";
    Object.assign(this, info);
  }
}

/** Lightweight form: coded("store")("store.conflict", msg, { retryable: true }). The code cannot be overridden by `rest`. */
export const coded = (ns: string) => (code: string, message: string, rest: Partial<ErrorInfo> = {}) => {
  if (!code.startsWith(`${ns}.`)) throw new Error(`code ${code} outside namespace ${ns}`);   // programmer error
  return new CodedError({ ...rest, code, message });
};

export const isCoded = (e: unknown): e is CodedError => e instanceof CodedError;
export const isCancelled = (e: unknown): boolean =>
  (isCoded(e) && e.code.endsWith(".cancelled")) || (e instanceof DOMException && e.name === "AbortError");

/** Text position for parse errors (used by §representations' parseError). */
export function textPos(src: string, offset: number): NonNullable<ErrorInfo["at"]> {
  const before = src.slice(0, offset);
  return { offset, line: before.split("\n").length, col: offset - (before.lastIndexOf("\n") + 1) + 1 };
}
```

## Catalog: codes are declared once, per module

```ts
export interface ErrorSpec { description: string; retryable?: boolean; example?: string; retired?: boolean }

/** Declares a module's codes and returns a typed thrower: a misspelled code is a compile error, a missing message uses the description. */
export function defineErrors<C extends Record<string, ErrorSpec>>(ns: string, spec: C) {
  for (const k of Object.keys(spec)) if (!k.startsWith(`${ns}.`)) throw new Error(`${k} outside ${ns}`);
  const error = (code: keyof C & string, message?: string, rest: Partial<ErrorInfo> = {}) =>
    new CodedError({ retryable: spec[code]!.retryable, ...rest, code, message: message ?? spec[code]!.description });
  return { spec, error };
}
```

Each module owns its catalog: `reprErrors` in [§representations](../../../representations/bootstrap-patterns/representations/typescript.md#print-and-parse-s-expression-spelling), a domain module's catalog in that module, and so on.

Namespaces, retiring codes and invariants: see [README.md](README.md#shape). Retire a code with `retired: true`.

- A namespace may be split across catalogs when the sub-prefixes do not overlap: `repr.parse.*`, `repr.bad_symbol` and `repr.mirror.*` in §representations, `repr.patch.*` in §committed-transient-identity, and further `repr.*` sub-prefixes in the project's own representation modules.

## Same shape at every layer; errors are data

```ts
throw coded("doc")("doc.not_found", `no doc ${id}`, { hint: "use doc.list" });
throw coded("store")("store.conflict", "version mismatch", { retryable: true, hint: "reload and retry" });
```

Errors serialize as representations, so they travel through the same streams, logs and command results as everything else:

```ts
export const errorToNode = (e: ErrorInfo): Node =>
  [sym("error"), [sym("code"), e.code], [sym("message"), e.message],
   ...(e.hint ? [[sym("hint"), e.hint]] : []),
   ...(e.at ? [[sym("pos"), e.at.line, e.at.col]] : []),   // not (at ...): that form is a locator
   ...(e.loc ? [locatorToNode(e.loc)] : []),                     // locatorToNode: 17
   ...(e.path ? [[sym("path"), ...e.path]] : []),
   ...(e.details !== undefined ? [[sym("details"), e.details]] : []),
   ...(e.retryable !== undefined ? [[sym("retryable"), e.retryable]] : [])];   // remote callers branch on it too
// (error (code "repr.parse.unexpected_end") (message "Input ended inside a form.") (hint "close open lists with )") (pos 1 6))
```

Agent-facing paths that should not throw return a result. One shape covers single errors and validation with many problems (one entry per problem, each with its `path` or `loc`):

```ts
export type Result<T> = { ok: true; value: T } | { ok: false; errors: readonly ErrorInfo[] };

export function tryParseLenient(src: string): Result<Lenient> {
  try { return { ok: true, value: parseLenient(src) }; }
  catch (e) { if (isCoded(e)) return { ok: false, errors: [e] }; throw e; }    // only coded errors are values
}
```

## Callers branch on codes

```ts
const r = tryParseLenient(raw);
if (!r.ok) {
  const [e] = r.errors;
  if (e!.code === "repr.parse.unexpected_end") return askToComplete(e!);      // agent: complete the form
  return report(r.errors);
}
// retries branch on the flag, not the text:  catch (e) { if (isCoded(e) && e.retryable) return retry(); throw e; }
// helpers assumed: askToComplete, report, retry
```

## Checks

The one required check (§cataloged-codes): every thrown code is declared. With `defineErrors(...).error` a typo does not compile; for the lightweight `coded(ns)` form, a gate script greps code literals at throw sites and fails on any that no catalog declares. Example tests are welcome where cheap:

```ts
// parse examples go through the STRICT parser (lenient repair would fix most of them); validation examples through validate(parse(x))
export const reprErrorsSuite = suite<void>("ReprErrors",
  kase("each parse example produces its code", "quick", () => {
    for (const [code, s] of Object.entries(reprErrors.spec))
      if (s.example && !s.retired && code.startsWith("repr.parse.")) expect(codeOf(() => parse(s.example!))).toBe(code);
  }),
  kase("positions point at the offending character", "quick", () => {
    const e = catchCoded(() => parse("(a ]"));
    expect(e.code).toBe("repr.parse.bad_char"); expect(e.at).toMatchObject({ line: 1, col: 4 });
  }),
  kase("agent-facing parse errors carry a hint", "quick", () => expect(catchCoded(() => parse("(a")).hint).toBeTruthy()),
);
// helpers assumed: codeOf(f) returns the thrown code; catchCoded(f) returns the thrown CodedError
```

A lint rejects `throw new Error(` at public module boundaries; internal invariant failures may stay plain (§programmer-errors). The catalogs render into the docs hub ([§design-system-storybook](../../../ui/bootstrap-patterns/design-system-storybook/typescript.md#documentation-hub-including-non-ui-modules)) as an "Error codes" table per module.
