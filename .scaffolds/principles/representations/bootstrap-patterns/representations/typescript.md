# Representations: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Module layout

```
lib/repr/                    # family folder; packages @scope/repr-<name>
  core/                      # tree, print/parse, repair, parseStream, Patch, walk, hashing/ops/doc (§committed-transient-identity)
  schema/                    # one domain representation (same shape for doc, board, ...)
    src/vocabulary.ts        # node tags + shapes the schema accepts
    src/migrations.ts        # v1 -> v2 ... as pure functions Node -> Node
    src/encodings/           # alternative encodings (JSON, binary), each passing encodingSuite (§conformance)
    src/query.ts             # domain query functions
    _spec/module.md
```

Domain modules depend only on `repr-core` and the common base (`lib/core`); executors and apps depend on them, never the reverse (§representation-modules).

```ts
// lib/repr/schema/src/migrations.ts: versioned, pure, testable on samples (§versioned-representations)
export const migrations: ReadonlyArray<(n: Node) => Node> = [
  /* v0 -> v1 */ (n) => mapTag(n, "shape", (s) => [...s, [sym("locked"), false]]),
];
export const migrate = (n: Node, from: number): Node => migrations.slice(from).reduce((x, m) => m(x), n);
declare function mapTag(n: Node, tag: string, f: (list: readonly Node[]) => Node): Node;   // rewrite every list with this tag
```

## The tree

```ts
export type Atom = string | number | bigint | boolean | null;   // strings are data; symbols are tags
export type Sym = { readonly sym: string };
export type Node = Atom | Sym | readonly Node[];                // lists; by convention head = tag

// A symbol must print and re-parse as itself: no delimiters, no leading digit/sign-digit, no '#', not nil/true/false,
// and no smart quotes (repair would read them as string quotes).
const SYM_RE = /^(?![-+]?\.?\d)[^\s()\[\]{}";,#“”][^\s()\[\]{}";,“”]*$/;
const RESERVED = new Set(["nil", "true", "false"]);
export const isSymName = (s: string): boolean => SYM_RE.test(s) && !RESERVED.has(s);
export const sym = (s: string): Sym => {
  if (!isSymName(s)) throw reprErrors.error("repr.bad_symbol", `bad symbol ${JSON.stringify(s)}`, { hint: "use letters and - _ . / : ! ? * + < > =" });
  return { sym: s };
};
export const isSym = (n: Node): n is Sym => typeof n === "object" && n !== null && !Array.isArray(n);
export const tag = (n: Node): string | undefined =>
  Array.isArray(n) && isSym((n as readonly Node[])[0] as Node) ? ((n as readonly Node[])[0] as Sym).sym : undefined;

export type Path = readonly number[];
export const nodeAt = (root: Node, path: Path): Node | undefined => {
  let n: Node | undefined = root;
  for (const i of path) n = Array.isArray(n) ? (n as readonly Node[])[i] : undefined;
  return n;
};
// (view (at 0 2 5)), a viewport at x, y, zoom: attrOf(view, "at") = 0, attrsOf(view, "at") = [0, 2, 5]
const attrList = (n: Node, key: string) =>
  Array.isArray(n) ? ((n as readonly Node[]).find((c) => tag(c) === key) as readonly Node[] | undefined) : undefined;
export const attrOf = (n: Node, key: string): Node | undefined => attrList(n, key)?.[1];
export const attrsOf = (n: Node, key: string): readonly Node[] => attrList(n, key)?.slice(1) ?? [];

/** Pre-order walk; `visit` returning false skips that node's children. Allocates a path per node (fine for scans, §performance). */
export function walk(n: Node, path: Path, visit: (n: Node, p: Path) => boolean | void): void {
  if (visit(n, path) === false || !Array.isArray(n)) return;
  const list = n as readonly Node[];
  for (let i = 0; i < list.length; i++) walk(list[i]!, [...path, i], visit);
}

// Example state: pure data, no behavior
const board: Node = [sym("board"),
  [sym("view"), [sym("at"), 0, 2, 5]],
  [sym("shape"), "box", [sym("fill"), "gray"]]];
```

## Print and parse (S-expression spelling)

`(` `)` and `"` delimit; `[ ] { } ; ,` are reserved and rejected by the strict parser (so `repair` can treat `;` and `,` as slips). Every token is validated in full.

```ts
export function print(n: Node): string {
  if (n === null) return "nil";
  if (typeof n === "string") return JSON.stringify(n);
  if (typeof n === "bigint") return `${n}n`;                       // 123n: never confused with a float
  if (typeof n === "number") {                                      // lossless: String() would lose -0 and print NaN/Infinity as symbols
    if (Object.is(n, -0)) return "-0";
    if (Number.isNaN(n)) return "#nan";
    if (!Number.isFinite(n)) return n > 0 ? "#inf" : "#-inf";
    return String(n);                                               // shortest string that round-trips a float64
  }
  if (typeof n === "boolean") return String(n);
  if (isSym(n)) return n.sym;                                       // safe: sym() validated the name
  return `(${n.map(print).join(" ")})`;
}

export function parse(src: string): Node {
  const p = new Parser(src); const n = p.node(); p.end(); return n;
}

const STR_RE = /"(?:[^"\\]|\\.)*"/y, TOKEN_RE = /[^\s()\[\]{}";,]+/y;      // sticky: no slicing
const NUM_RE = /^-?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i, BIGINT_RE = /^-?\d+n$/, LOOKS_NUMERIC = /^[-+]?\.?\d/;

class Parser {
  private i = 0;
  constructor(private readonly s: string) {}
  private ws() { while (/\s/.test(this.s[this.i] ?? "")) this.i++; }
  private match(re: RegExp): string | undefined { re.lastIndex = this.i; const m = re.exec(this.s)?.[0]; if (m) this.i += m.length; return m; }
  node(): Node {
    this.ws(); const at = this.i, c = this.s[at];
    if (c === undefined) throw parseError("repr.parse.unexpected_end", this.s, at, "close open lists with )");
    if (c === "(") { this.i++; const out: Node[] = [];
      for (;;) { this.ws(); if (this.s[this.i] === ")") { this.i++; return out; } out.push(this.node()); } }
    if (c === '"') {
      const m = this.match(STR_RE);
      try { if (m) return JSON.parse(m) as string; } catch { /* raw control character or bad escape */ }
      throw parseError("repr.parse.bad_string", this.s, at, 'close the string; use JSON escapes (\\n \\" \\\\ \\uXXXX)');
    }
    const t = this.match(TOKEN_RE);
    if (!t) throw parseError("repr.parse.bad_char", this.s, at, `unexpected ${JSON.stringify(c)}; [ ] { } ; , are reserved`);
    if (t === "nil") return null; if (t === "true") return true; if (t === "false") return false;
    if (t === "#nan") return NaN; if (t === "#inf") return Infinity; if (t === "#-inf") return -Infinity;
    if (BIGINT_RE.test(t)) return BigInt(t.slice(0, -1));
    if (NUM_RE.test(t)) return Number(t);
    if (LOOKS_NUMERIC.test(t) || !isSymName(t)) throw parseError("repr.parse.bad_token", this.s, at, "numbers: 1 -2.5 1e-9 12n #nan #inf #-inf");
    return { sym: t };
  }
  end() { this.ws(); if (this.i < this.s.length) throw parseError("repr.parse.trailing_input", this.s, this.i, "exactly one top-level form is allowed"); }
}

// Error catalog for this module (§error-codes). 15 (opErrors) and 18 (embedErrors) declare their own repr.* codes beside their code.
export const reprErrors = defineErrors("repr", {
  "repr.bad_symbol":            { description: "Not a valid symbol name." },
  "repr.parse.unexpected_end":  { description: "Input ended inside a form.",            example: "(a (b" },
  "repr.parse.bad_string":      { description: "Malformed string literal.",             example: '("abc' },
  "repr.parse.bad_char":        { description: "Character not valid at this position.", example: "(a ]" },
  "repr.parse.bad_token":       { description: "Malformed number or symbol.",          example: "(1abc)" },
  "repr.parse.trailing_input":  { description: "Extra input after the first form.",    example: "(a) (b)" },
  "repr.mirror.bad_path":       { description: "Patch path does not exist in the tree." },   // repr.patch.* is §committed-transient-identity's (ops)
});
export const parseError = (code: keyof typeof reprErrors.spec & string, src: string, offset: number, hint?: string) =>
  reprErrors.error(code, undefined, { hint, at: textPos(src, offset) });          // textPos: 14
```

## Lenient parsing and trivial repairs

Agents slip: an unbalanced paren at the end, a fenced code block, a smart quote. `repair` fixes only these mechanical slips, reports each one, and hands the text to the *strict* parser, so the strict grammar stays the single definition.

```ts
export type RepairKind =
  | "stripped-fence" | "smart-quotes" | "dropped-comment" | "commas-as-space"
  | "escaped-newline" | "closed-string" | "closed-parens" | "dropped-extra-close";
export interface Repair { kind: RepairKind; detail?: string }

export function repair(src: string): { text: string; repairs: Repair[] } {
  const repairs: Repair[] = [];
  const note = (r: Repair) => { if (!repairs.some((x) => x.kind === r.kind && x.detail === r.detail)) repairs.push(r); };

  let s = src.trim();
  const fence = /^```[\w-]*\r?\n([\s\S]*?)\r?\n?```$/.exec(s);
  if (fence) { s = fence[1]!; note({ kind: "stripped-fence" }); }

  let out = "", depth = 0, opener = "", esc = false, extra = 0;      // opener: the quote that opened the current string
  for (let i = 0; i < s.length; i++) {
    const c = s[i]!;
    if (opener) {
      if (esc) { out += c; esc = false; }
      else if (c === "\\") { out += c; esc = true; }
      else if (c === '"' && opener === '"') { out += c; opener = ""; }
      else if (c === "”" && opener === "“") { out += '"'; opener = ""; note({ kind: "smart-quotes" }); }
      else if (c === '"') out += '\\"';                                  // plain quote inside a smart-quoted string
      else if (c === "\n" || c === "\r") { out += c === "\n" ? "\\n" : "\\r"; note({ kind: "escaped-newline" }); }
      else out += c;                                                     // includes a closing smart quote inside a "..." string
      continue;
    }
    if (c === '"' || c === "“") { if (c !== '"') note({ kind: "smart-quotes" }); out += '"'; opener = c; }
    else if (c === ";") { while (i + 1 < s.length && s[i + 1] !== "\n") i++; note({ kind: "dropped-comment" }); }
    else if (c === ",") { out += " "; note({ kind: "commas-as-space" }); }
    else if (c === "(") { depth++; out += c; }
    else if (c === ")") { if (depth === 0) extra++; else { depth--; out += c; } }
    else out += c;
  }
  if (opener) { out += '"'; note({ kind: "closed-string" }); }
  if (depth > 0) { out += ")".repeat(depth); note({ kind: "closed-parens", detail: `${depth} at end` }); }
  if (extra > 0) note({ kind: "dropped-extra-close", detail: `${extra} unmatched` });
  return { text: out, repairs };
}

export interface Lenient { node: Node; repairs: Repair[]; normalized: string }

/** Only for agent/user-facing input. Anything the repair cannot fix mechanically still throws from strict parse(). */
export function parseLenient(src: string): Lenient {
  const { text, repairs } = repair(src);
  const node = parse(text);
  return { node, repairs, normalized: print(node) };
}
```

What is and is not repaired, and the cases: see [README.md](README.md#shape). Commands use it through `reprArg` ([§commands](../../../architecture/bootstrap-patterns/commands/README.md)) for representation-valued string fields and echo `{ normalized, repairs }` in their result, e.g. `[{"kind":"closed-parens","detail":"2 at end"}]`.

## Incremental text streaming

Feed chunks; each complete top-level form (list, string or atom) is emitted as soon as it ends. An atom at top level ends at the next delimiter or at completion.

```ts
export function parseStream(chunks: Observable<string>): Observable<Node> {
  return observable((emit, fail, done) => {
    let buf = "", i = 0, start = -1, depth = 0, inStr = false, esc = false;
    const form = (end: number) => { const text = buf.slice(start, end); start = -1; emit(parse(text)); };
    const scan = () => {
      for (; i < buf.length; i++) {
        const c = buf[i]!;
        if (inStr) { if (esc) esc = false; else if (c === "\\") esc = true; else if (c === '"') { inStr = false; if (depth === 0) form(i + 1); } continue; }
        if (depth === 0 && start >= 0 && /[\s()"]/.test(c)) form(i);                  // top-level atom ended; c handled below
        if (c === '"') { inStr = true; if (depth === 0) start = i; }
        else if (c === "(") { if (depth++ === 0) start = i; }
        else if (c === ")") { if (depth === 0) throw parseError("repr.parse.bad_char", buf, i, "unbalanced )"); if (--depth === 0) form(i + 1); }
        else if (depth === 0 && start < 0 && !/\s/.test(c)) start = i;
      }
      const keep = start >= 0 ? start : i;              // drop consumed text once per chunk, not once per form (§performance)
      if (keep > 0) { buf = buf.slice(keep); i -= keep; if (start >= 0) start = 0; }
    };
    return chunks.subscribe({
      next(chunk) { try { buf += chunk; scan(); } catch (e) { fail(e); } },
      error: fail,
      complete() {
        try {
          if (inStr || depth > 0) throw parseError("repr.parse.unexpected_end", buf, buf.length, "stream ended inside a form");
          if (start >= 0) form(buf.length);
          done();
        } catch (e) { fail(e); }
      },
    });
  });
}
```

## Snapshot + patches (wire streaming and mirroring)

State is a behavior of the root: `TreeState = ReadonlyBehavior<Node>` ([§api-surfaces](../../../architecture/bootstrap-patterns/api-surfaces/README.md)). Read it with `state.value`, follow it with `changes(state)`, and subscribe to a subtree with `select(state, (r) => nodeAt(r, path))` (structural sharing keeps unchanged subtrees `===`).

Positional `Patch` exists only to stream and mirror a tree over the wire (§streamable). It is **not** the edit channel: committed changes are id-addressed `Op`s in an `OpEnvelope` ([§committed-transient-identity](../committed-transient-identity/README.md)), which survive concurrent edits; positions do not.

```ts
export type TreeState = ReadonlyBehavior<Node>;

export type Patch =
  | { op: "replace"; path: Path; value: Node }
  | { op: "insert"; path: Path; value: Node }                    // insert at index path[last] of the parent
  | { op: "remove"; path: Path };

export function diff(a: Node, b: Node, path: Path = []): Patch[] {
  if (a === b || hashNode(a) === hashNode(b)) return [];       // hashNode: 15, cached per immutable node
  if (Array.isArray(a) && Array.isArray(b) && tag(a) === tag(b)) {
    const x = a as readonly Node[], y = b as readonly Node[], n = Math.min(x.length, y.length), out: Patch[] = [];
    for (let i = 0; i < n; i++) out.push(...diff(x[i]!, y[i]!, [...path, i]));
    for (let i = x.length - 1; i >= n; i--) out.push({ op: "remove", path: [...path, i] });     // from the end: indices stay valid
    for (let i = n; i < y.length; i++) out.push({ op: "insert", path: [...path, i], value: y[i]! });
    return out;
  }
  return [{ op: "replace", path, value: b }];                  // positional: a reorder is replaces, which is fine for mirroring
}

export const apply = (root: Node, patches: readonly Patch[]): Node => patches.reduce(applyOne, root);

function applyOne(root: Node, p: Patch): Node {
  const bad = () => reprErrors.error("repr.mirror.bad_path", `no ${p.op} target at (${p.path.join(" ")})`, { path: p.path });
  if (p.path.length === 0) { if (p.op !== "replace") throw bad(); return p.value; }
  if (!Array.isArray(root)) throw bad();
  const list = root as readonly Node[], [i, ...rest] = p.path as [number, ...number[]];
  const last = rest.length === 0 && p.op === "insert" ? list.length : list.length - 1;
  if (!Number.isInteger(i) || i < 0 || i > last) throw bad();
  const out = list.slice();
  if (rest.length) out[i] = applyOne(list[i]!, { ...p, path: rest });
  else if (p.op === "replace") out[i] = p.value;
  else if (p.op === "insert") out.splice(i, 0, p.value);
  else out.splice(i, 1);
  return out;
}

/** Sender: patches between consecutive roots. */
export const patchStream = (state: TreeState): Observable<Patch[]> => observable((emit) => {
  let prev = state.value;
  return changes(state).subscribe((next) => { const ps = diff(prev, next); prev = next; if (ps.length) emit(ps); });
});

/** Receiver: snapshot first, then patches. */
export function mirror(snapshot: Node, patches: Observable<Patch[]>): { state: TreeState; stop: Unsubscribe } {
  const s = createBehaviorSubject(snapshot);
  return { state: s, stop: patches.subscribe((ps) => s.next(apply(s.value, ps))) };
}
```

## Executors are bound in composition, not in state

State holds descriptions; an executor interprets them (§executors-by-composition). State never imports an executor. Picking/hit-testing is not an executor method; it is a screen-to-doc mapping.

```ts
// port (repr-core): interprets a declarative description
export interface Executor {
  /** Reconcile the running world with this state. `hint` lets it skip unchanged subtrees (by hash, §committed-transient-identity). */
  reconcile(state: Node, hint?: { prev?: Node; changed?: readonly Path[] }): void | Promise<void>;
  dispose(): void;
}

/** Serialized and latest-wins: one reconcile at a time; intermediate states are skipped; errors are reported, not thrown. */
export function bindExecutor(state: TreeState, ex: Executor, report: (e: unknown) => void = console.error): Unsubscribe {
  let running = false, stopped = false, pending: Node | undefined, prev: Node | undefined;  // undefined = none (null is a Node)
  const run = async (first: Node) => {
    running = true;
    for (let next: Node | undefined = first; next !== undefined && !stopped; next = pending, pending = undefined) {
      try { await ex.reconcile(next, { prev }); prev = next; } catch (e) { report(e); }
    }
    running = false;
  };
  const push = (n: Node) => { if (running) pending = n; else void run(n); };
  push(state.value);
  const off = changes(state).subscribe(push);
  return () => { stopped = true; off(); };                       // disposing the executor is its creator's job (§explicit-lifetimes)
}

// production: bindExecutor(doc.effective, canvasRenderer(canvas))
// tests:      bindExecutor(doc.effective, recordingExecutor())  -> assert on the recorded interpretation
// agents:     no executor at all; inspect state through commands (§commands) and queries
```

Invariants: see [README.md](README.md#invariants).

## Property tests (quick tier)

```ts
// arbNode generates only valid symbols (isSymName) and every atom kind, including -0, NaN, ±Infinity and bigint.
test.prop([arbNode, arbNode])("apply(a, diff(a, b)) == b", (a, b) => {
  expect(print(apply(a, diff(a, b)))).toBe(print(b));
});
test.prop([arbNode])("parse(print(x)) == x", (x) => expect(print(parse(print(x)))).toBe(print(x)));
test.prop([arbNode])("repair(print(x)) is the identity", (x) => expect(repair(print(x))).toEqual({ text: print(x), repairs: [] }));
test.prop([fc.array(arbNode, { maxLength: 4 }), fc.integer({ min: 1, max: 7 })])("stream parse == whole parse (any chunking)", async (xs, n) => {
  const chunks = chunkEvery(xs.map(print).join(" "), n);
  expect((await collect(parseStream(of(...chunks)))).map(print)).toEqual(xs.map(print));
});
// helpers assumed: arbNode, chunkEvery, collect
```
