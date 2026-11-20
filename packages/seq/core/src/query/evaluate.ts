/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A nested value of a form (_spec/README.md "Queries", Forms): a form evaluated against the same prefix (its resolver,
// else child selection), an address from the root; its issues and pending regions carried into the outer form's.
// Internal to the query entry.

import {
  type Issue,
  isNode,
  matchesForm,
  type Named,
  own,
  type Selected,
  type SelectOptions,
  selectAddress,
  type Value,
} from '@molstar/core/repr/tree';
import { checkLimit, type Ctx, Failure, fail, prefixNodes, tick } from './form-kit.js';
import { QueryIssues } from './issues.js';
import { type Location, unionLocations } from './location.js';
import { fromNamed, namedType, type QueryResult } from './result.js';

/** A nested value's result: a form against the prefix (its resolver, else child selection), an address from the root. */
export function evaluate(ctx: Ctx, value: Value, prefix: Named): QueryResult {
  if (isNode(value)) {
    const resolver = own(ctx.resolvers, value.kind);
    if (!resolver) {
      const out: Selected[] = [];
      for (const m of prefixNodes(prefix, value.kind))
        for (const c of m.node.children ?? []) {
          tick(ctx);
          if (c.id !== undefined && matchesForm(value, c)) out.push({ node: c, address: [...m.address, c.id] });
        }
      checkLimit(ctx, out.length);
      return { type: 'nodes', nodes: out };
    }
    const r = resolver({ ...ctx.input, form: value, prefix });
    return carry(ctx, r, r.result);
  }
  if (typeof value === 'string' && value.startsWith('/')) {
    const { input } = ctx;
    const options: SelectOptions = {
      limit: input.limit,
      maxDepth: input.maxDepth,
      resolvers: ctx.resolvers,
      signal: input.signal,
      deadline: input.deadline,
    };
    const sel = selectAddress(input.root, value, options);
    // refs are not resolved yet (a later pass's): an operand naming one is refused, never read as the ref node
    if (sel.matches.some((m) => m.node.kind === 'ref'))
      fail('seq.query.ref.unresolved', 'An operand names a ref', {
        hint: QueryIssues.codes['seq.query.ref.unresolved'].hint,
      });
    return carry(ctx, sel, sel.result ?? { kind: 'nodes', matches: sel.matches });
  }
  return fail('seq.query.form.param', 'A form takes a nested form or an address here');
}

/** Carries a nested result's issues and pending regions; its errors fail the outer form. */
export function carry(
  ctx: Ctx,
  r: { readonly issues?: readonly Issue[]; readonly complete?: boolean; readonly pending?: unknown },
  named: Named | undefined,
): QueryResult {
  const issues = r.issues ?? [];
  if (issues.some((i) => i.severity === 'error') || !named) throw new Failure(issues);
  ctx.issues.push(...issues);
  if (r.complete === false) {
    ctx.complete = false;
    if (r.pending) ctx.pending = unionLocations(ctx.pending, r.pending as Location);
  }
  const result = fromNamed(named);
  return result ?? fail('seq.query.form.type', `A nested form gives a result of type ${namedType(named)}`);
}

export function asType<T extends QueryResult['type']>(
  r: QueryResult,
  type: T,
  kind: string,
): Extract<QueryResult, { type: T }> {
  if (r.type !== type) fail('seq.query.form.type', `A form of kind ${kind} takes ${type}, and is given ${r.type}`);
  return r as Extract<QueryResult, { type: T }>;
}
