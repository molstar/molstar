/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's language types (seqspec:seqql#types, seqspec:seqql#names, seqspec:seqql#errors): slots, argument and result
// types, symbol definitions with the metadata the checker, the skill and the reference read, and expressions.

/** The element a test reads, bound by the `:where` it sits in (seqspec:seqql#elements-and-groups). */
export type Slot = "position" | "system" | "item" | "group";

/** Argument and result types, checked before evaluation. */
export type Type =
  | "num"
  | "str"
  | "bool"
  /** A non-negative integer: an index. */
  | "int"
  /** A field's value: a number, string, boolean or null. */
  | "json"
  /** A set ref: a word or string ("1brs/A", "msa/3", "domains"). */
  | "ref"
  /** A label or label range: a word or number ("10-20", "45A", 45). */
  | "label"
  /** A selection: the monadic type; a word or string here is a ref. */
  | "query"
  /** A bool read per element of a slot. */
  | "test"
  | { readonly enum: readonly string[] };

export interface ArgDef {
  /** The keyword, or the positional's name in the signature. */
  readonly name: string;
  readonly type: Type;
  /** Written without its keyword, in order. */
  readonly positional?: boolean;
  /** The last positional takes every remaining one. */
  readonly rest?: boolean;
  readonly optional?: boolean;
  readonly default?: null | boolean | number | string;
  /** A `:where`: the slot its test reads. */
  readonly binds?: Slot;
  readonly doc: string;
}

export type SymbolKind = "core" | "generator" | "modifier" | "combinator" | "filter" | "property" | "macro";

export interface SymbolExample {
  /** An S-expression: "(pick (items :of domains) :where (matches (label) kinase))". */
  readonly text: string;
  /** In English: "kinase domains". */
  readonly says: string;
  /** Its place in the skill's EXAMPLES, from 1 (absent: the reference only). */
  readonly skill?: number;
  /** A committed case it evaluates on ("case7"), for the table's tests. */
  readonly fixture?: string;
}

export interface SymbolHint {
  readonly code: string;
  readonly text: string;
  /** A `seqql.type` hint for one rejected enum value (any case, plural or not): the message when it is written. */
  readonly value?: string;
}

/** A symbol's definition: what `defineSymbol` takes (seqspec:seqql#extension). */
export interface SymbolSpec {
  /** Full: "seq.item.label", "molstar.select". */
  readonly name: string;
  /** Short: "label"; null: the full name only. */
  readonly alias: string | null;
  readonly kind: SymbolKind;
  readonly args: readonly ArgDef[];
  readonly result: Type;
  /** A property's slot. */
  readonly reads?: Slot;
  /** "(items :of A :where T)", in short names: errors and the reference. */
  readonly signature: string;
  /** One line: "items of annotation A (all annotations if no :of), one group each". */
  readonly doc: string;
  /** One or two. */
  readonly examples: readonly SymbolExample[];
  /** Error hints the checker adds. */
  readonly hints: readonly SymbolHint[];
  /**
   * A macro (sel, an app's word): its expansion, for checking and compiling (the JSON keeps the call as written).
   * Takes the bound call; reports issues inside its string arguments by index (seqspec:seqql#the-sel-shorthand).
   */
  readonly expand?: (call: Apply, report: (r: MacroReport) => void) => Expression;
  /** A macro's stored form of its call (sel: its string without spaces), applied by the binder and the builder. */
  readonly normalize?: (call: Apply) => Apply;
}

/**
 * An issue a macro finds in its call: `index` and `open` count UTF-16 units into string argument `arg` ("0"). With a
 * source, the issue is placed there, and the message gains the position of `open` (" at 1:16").
 */
export interface MacroReport {
  readonly code: string;
  readonly message: string;
  readonly arg?: string;
  readonly index?: number;
  readonly open?: number;
}

/** Keyword arguments, the builder's trailing object: `{ where: T }`. */
export type Keywords = Readonly<Record<string, Expression | readonly Expression[] | undefined>>;

/** A defined symbol: its spec, callable as the builder (positionals, then an optional keywords object). */
export interface SymbolDef extends SymbolSpec {
  (...args: readonly (Expression | Keywords)[]): Apply;
}

export type Expression = string | number | boolean | Apply;

export interface Apply {
  /** A full name. */
  readonly head: string;
  /** Positionals under "0", "1", …; a rest positional is an array at its index; keywords by name. */
  readonly args?: Readonly<Record<string, Expression | readonly Expression[]>>;
}

/** What a document or a story stores. */
export interface SeqQLQuery {
  /** The SeqQL version, "1.0" (seqspec:seqql#in-documents). */
  readonly seqql: string;
  /** Full names: canonical. */
  readonly expr: Expression;
  /** print(expr): the canonical S-expression, short names. */
  readonly text: string;
}

export const SEQQL_VERSION = "1.0";

/** Structurally a model `Issue`, plus the place and help. */
export interface QueryIssue {
  /** "seqql.symbol.unknown". */
  readonly code: string;
  /** seqspec:seqql#errors: one per code (`seqqlSeverity`); `info` a note on how a valid query ran (`seqql.empty`). */
  readonly severity: "error" | "warning" | "info";
  /** A JSON pointer into the expression ("" for reader and registration errors). */
  readonly path: string;
  /** What is wrong and what to write instead. */
  readonly message: string;
  /** 1-based, in the text. */
  readonly at: { readonly line: number; readonly column: number } | null;
  /** The unclosed "(" or "[". */
  readonly open: { readonly line: number; readonly column: number } | null;
  /** The symbol's, in short names, for checker errors. */
  readonly signature: string | null;
  /** The names it could have been, nearest first, at most 8. */
  readonly known?: readonly string[];
}

export function isApply(e: unknown): e is Apply {
  return typeof e === "object" && e !== null && !Array.isArray(e) && typeof (e as Apply).head === "string";
}

/**
 * How deep a stored expression may nest: as deep as the reader lets text nest (`read.ts`' `MAX_DEPTH`), so that
 * checking, binding and compiling, which recurse, never meet a deeper one (seqspec:seqql#the-reader).
 */
export const MAX_EXPRESSION_DEPTH = 64;

/**
 * Whether a value has the JSON shape of an expression: a base type, or `{ head, args? }` with expressions inside,
 * nested at most `MAX_EXPRESSION_DEPTH` deep (`expressionTooDeep` says when that is why not).
 */
export function isExpression(value: unknown, depth = 0): value is Expression {
  if (typeof value === "string" || typeof value === "boolean") return true;
  if (typeof value === "number") return Number.isFinite(value);
  if (!isApply(value) || depth >= MAX_EXPRESSION_DEPTH) return false;
  for (const key of Object.keys(value)) if (key !== "head" && key !== "args") return false;
  const args: unknown = value.args;
  if (args === undefined) return true;
  if (typeof args !== "object" || args === null || Array.isArray(args)) return false;
  for (const v of Object.values(args)) {
    if (Array.isArray(v) ? !v.every((x) => isExpression(x, depth + 1)) : !isExpression(v, depth + 1)) return false;
  }
  return true;
}

/** Whether `{ head, args }` values nest deeper than `MAX_EXPRESSION_DEPTH`, found without recursing (any input). */
export function expressionTooDeep(value: unknown): boolean {
  let level: unknown[] = [value];
  for (let depth = 0; level.length > 0; depth++) {
    if (depth >= MAX_EXPRESSION_DEPTH) return true;
    const next: unknown[] = [];
    for (const e of level) {
      const args = isApply(e) ? (e as { args?: unknown }).args : undefined;
      if (typeof args !== "object" || args === null) continue;
      for (const v of Object.values(args)) {
        if (Array.isArray(v)) for (const x of v) next.push(x);
        else next.push(v);
      }
    }
    level = next.filter(isApply);
  }
  return false;
}
