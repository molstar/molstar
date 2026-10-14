/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `query` computed provider (seqspec:seqql#in-documents; seqspec:format#computed-nodes): a SeqQL query stored in a
// document as `{ seqql, expr, text }`, checked statically (the query language's checker, its codes as `svs.query.*`; `search` and the
// host's state refused, a document being data only; the refs it names as the node's dependencies, `sel` expanded), and
// run at load into a `ranges` annotation, one item per group.

import { type AnnotationInput, Issue } from "@molstar/seq-core/model";
import { compileQuery, createQueryEnv, SeqSelection } from "@molstar/seq-core/model/query";
import {
  type Apply,
  checkQuery,
  type Expression,
  expressionTooDeep,
  isApply,
  isExpression,
  MAX_EXPRESSION_DEPTH,
  printQuery,
  type QueryIssue,
  type SeqQLQuery,
} from "@molstar/seq-core/query";
import type { StageContext } from "./contract.ts";
import { report as reportAt } from "./load/issues.ts";
import { type CheckContext, type FieldType, isPlainObject, ptr, report } from "./schema/field-schema.ts";

/** Where a ref in a query may point: anything a selection or an annotation argument names. */
const KINDS = ["sequence", "ligand", "unit", "complex", "alignment", "row", "annotation", "computed"];
const SEARCH =
  "a document query cannot search; run the search, then save its hits as a ranges annotation (results.save) and name that annotation in the query";
const INPUT =
  'a document\'s query is data only: it cannot read what the host passes ((input …)) or (current) outside each; name the data instead, as a ref or (sel "1brs/A[10-20]")';
const ENV = createQueryEnv();

/** The path of `e`'s part at `path` in the stored expression, or its nearest enclosing part (a `sel` expanded). */
function placed(stored: Expression, path: string): string {
  let [x, at, out] = [stored as unknown, "", ""];
  for (const k of path.split("/").slice(1)) {
    x = (x as Record<string, unknown> | undefined)?.[k];
    if (x === undefined) return out;
    at += `/${k}`;
    if (isApply(x)) out = at;
  }
  return at;
}

/** Refs, `search` and `input` in an expanded expression: refs where a selection or a ref goes. */
function walk(e: Expression, path: string, visit: (kind: "ref" | "search" | "input", path: string, v: string) => void) {
  if (!isApply(e)) return;
  if (e.head === "seq.generator.search") visit("search", path, "");
  if (e.head === "seq.generator.input") visit("input", path, "");
  const defs = ENV.language.symbols.get(e.head)?.args ?? [];
  const positional = defs.filter((d) => d.positional);
  for (const [k, v] of Object.entries((e as Apply).args ?? {})) {
    const def = /^\d+$/.test(k) ? positional[Number(k)] : defs.find((d) => d.name === k && !d.positional);
    const refs = def?.type === "query" || def?.type === "ref";
    const each = Array.isArray(v)
      ? v.map((x, i) => [x, `${path}/args/${k}/${i}`] as const)
      : [[v, `${path}/args/${k}`] as const];
    for (const [x, p] of each) {
      if (isApply(x)) walk(x, p, visit);
      else if (refs && typeof x === "string") visit("ref", p, x);
    }
  }
}

/** The heads of the refs a valid query names (`sel` expanded): a query node's dependencies; none when it does not check. */
export function queryRefs(q: SeqQLQuery | undefined): string[] {
  const checked = q && isExpression(q.expr) ? checkQuery(q.expr) : null;
  const out: string[] = [];
  if (checked?.expr) walk(checked.expr, "", (kind, _, ref) => kind === "ref" && out.push(ref.split("/")[0]));
  return out;
}

/**
 * `params.query`: the stored form, checked; its refs pushed as the node's (cycles, unresolved names:
 * seqspec:format#validation).
 */
export const seqqlQuery: FieldType<SeqQLQuery> = {
  name: "SeqQLQuery",
  check(v, path, ctx: CheckContext) {
    if (isPlainObject(v) && expressionTooDeep(v.expr)) {
      const message = `seqql.read.depth: the expression nests deeper than ${MAX_EXPRESSION_DEPTH}; write it flatter`;
      return report(ctx, "svs.query.invalid", ptr(path, "expr"), message);
    }
    if (!isPlainObject(v) || typeof v.seqql !== "string" || typeof v.text !== "string" || !isExpression(v.expr)) {
      return report(ctx, "svs.params.type", path, 'a query is { "seqql": "1.0", "expr": …, "text": "(…)" }');
    }
    if (!/^1\.\d+$/.test(v.seqql))
      return report(
        ctx,
        "svs.query.version",
        ptr(path, "seqql"),
        `this query is SeqQL ${v.seqql}; this reader reads 1.x`,
      );
    const expr = ptr(path, "expr");
    const checked = checkQuery(v.expr as Expression);
    for (const i of checked.issues) {
      const input = i.code === "seqql.current";
      const message = `${i.code}: ${input ? INPUT : i.message}`;
      if (i.severity === "error") report(ctx, input ? "svs.query.input" : "svs.query.invalid", expr + i.path, message);
      else report(ctx, i.severity === "info" ? "svs.query.note" : "svs.query.warning", expr + i.path, message);
    }
    const stored = v.expr as Expression;
    if (!checked.expr) return;
    walk(checked.expr, "", (kind, p, ref) => {
      const at = expr + placed(stored, p);
      if (kind === "search") report(ctx, "svs.query.search", at, SEARCH);
      else if (kind === "input") report(ctx, "svs.query.input", at, INPUT);
      else ctx.refs.push({ ref, kinds: KINDS, path: at, from: ctx.from });
    });
    const printed = printQuery(stored);
    if (printed !== v.text) {
      const message = `text is not the printed expr, which wins: ${printed}`;
      ctx.issues.push(Issue.warning("svs.query.stale", ptr(path, "text"), message));
    }
  },
};

/** Run errors that say the environment cannot run the query now, not that it is written wrong. */
const UNAVAILABLE = /^seqql\.(limit\.|symbol\.unavailable$)/;

/**
 * A run's issue as the node's (seqspec:seqql#in-documents): a warning is `svs.query.warning`, a note
 * `svs.query.note` (an empty result explained); an error the checker
 * cannot see (a regex that does not compile, an unknown field) is `svs.query.invalid`, and one of the environment (a
 * limit, a symbol with no runtime) leaves the node unavailable. A cancel is no issue.
 */
function said(ctx: StageContext, i: QueryIssue): void {
  if (i.code === "seqql.aborted") return;
  const code =
    i.severity === "info"
      ? "svs.query.note"
      : i.severity === "warning"
        ? "svs.query.warning"
        : UNAVAILABLE.test(i.code)
          ? "svs.provider.unavailable"
          : "svs.query.invalid";
  reportAt(ctx, code, `/params/query/expr${i.path}`, `${i.code}: ${i.message}`);
}

/** The query run on the set: a ranges annotation, one item per group (its item's id, else `q<k>`), label and score. */
export async function computeQuery(
  ctx: StageContext,
  o: { query: SeqQLQuery; label: string | null },
): Promise<AnnotationInput | null> {
  const set = ctx.set;
  if (!set) return null;
  const r = compileQuery(o.query.expr, ENV).run({ set, signal: ctx.signal });
  for (const i of r.issues) said(ctx, i);
  if (r.issues.some((i) => i.severity === "error")) return null;
  return SeqSelection.toRangesInput(r.selection, "", { label: o.label ?? undefined });
}
