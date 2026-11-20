/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What every query form shares (_spec/README.md "Queries", Forms; Staged and bounded): a failure thrown inside a
// resolver and caught at its edge, the shared visit counter read against the signal and the deadline, the caller's
// limit, what is pending, a form's params, and what its prefix names. Internal to the query entry.

import {
  CancelEvery,
  type FormResolver,
  type FormResolverInput,
  type Issue,
  type IssuePlace,
  issueOf,
  type Named,
  type Node,
  own,
  printAddress,
  type Selected,
  TreeIssues,
  type Value,
} from '@molstar/core/repr/tree';
import { type QueryIssueCode, QueryIssues } from './issues.js';
import { isEmptyLocation, type Location, unionLocations } from './location.js';
import type { QueryAnnotation, QueryModel, QuerySystem } from './model.js';
import type { CompileContext } from './expressions.js';
import type { QueryResult } from './result.js';

/** An evaluation failed: its issues (at least one error). Thrown inside a resolver, caught at its edge. */
export class Failure {
  constructor(readonly issues: readonly Issue[]) {}
}

export const fail = (code: QueryIssueCode, message: string, place?: IssuePlace): never => {
  throw new Failure([issueOf(QueryIssues, code, message, place)]);
};

/** What one evaluation shares across its nested forms, nested addresses included. */
export interface Shared {
  /** Visits since the signal and the deadline were last read. */
  visits: number;
  /** Inside a rect's side: a range may not wrap (EC-10). */
  noWrap: boolean;
}

/** One resolver call's state: its input, the issues it carries, what is pending, the shared counter. */
export interface Ctx {
  readonly model: QueryModel;
  readonly input: FormResolverInput;
  readonly resolvers: Readonly<Record<string, FormResolver>>;
  readonly issues: Issue[];
  readonly shared: Shared;
  /** The properties and operators predicates read: the built-ins and the registry's (`query-extensions.ts`). */
  readonly ext: CompileContext;
  complete: boolean;
  pending: Location;
}

/** Counts a visit; every CancelEvery reads the signal and the deadline. */
export function tick(ctx: Ctx) {
  if (++ctx.shared.visits < CancelEvery) return;
  ctx.shared.visits = 0;
  const { signal, deadline } = ctx.input;
  if (signal?.aborted)
    throw new Failure([issueOf(TreeIssues, 'repr.tree.address.cancelled', 'The selection was cancelled')]);
  if (deadline !== undefined && performance.now() > deadline)
    throw new Failure([
      issueOf(TreeIssues, 'repr.tree.address.deadline_exceeded', 'The selection ran past its deadline'),
    ]);
}

export function checkLimit(ctx: Ctx, count: number) {
  const { limit } = ctx.input;
  if (count <= limit) return;
  throw new Failure([
    issueOf(TreeIssues, 'repr.tree.address.limit_exceeded', `More than ${limit} matches`, {
      hint: 'narrow the address, or raise the limit',
    }),
  ]);
}

export function stage(ctx: Ctx, pending: Location | undefined) {
  if (!pending || isEmptyLocation(pending)) return;
  ctx.complete = false;
  ctx.pending = unionLocations(ctx.pending, pending);
}

// ---- params ----

export function checkParams(form: Node, allowed: readonly string[]) {
  for (const key of Object.keys(form.params ?? {}))
    if (!allowed.includes(key))
      fail('seq.query.form.param', `A form of kind ${form.kind} names a param it does not take`, {
        hint: allowed.length ? `params: ${allowed.join(', ')}` : 'no params',
      });
}

export const param = (form: Node, key: string): Value | undefined => (form.params ? own(form.params, key) : undefined);

export function required(form: Node, key: string): Value {
  const v = param(form, key);
  if (v === undefined) fail('seq.query.form.param', `A form of kind ${form.kind} needs :${key}`);
  return v as Value;
}

// ---- what the prefix names ----

export function prefixNodes(prefix: Named, kind: string): readonly Selected[] {
  if (prefix.kind === 'nodes') return prefix.matches;
  return fail('seq.query.form.prefix', `A form of kind ${kind} applies to nodes, and its prefix names ${prefix.type}`);
}

export function systemsOf(ctx: Ctx, prefix: Named, kind: string): QuerySystem[] {
  return prefixNodes(prefix, kind).map((m) => {
    const address = printAddress(m.address);
    const s = ctx.model.system(address);
    return s ?? fail('seq.query.form.prefix', `A form of kind ${kind} applies to systems`, { address });
  });
}

export function annotationsOf(ctx: Ctx, prefix: Named, kind: string): QueryAnnotation[] {
  return prefixNodes(prefix, kind).map((m) => {
    const address = printAddress(m.address);
    const a = ctx.model.annotation(address);
    return a ?? fail('seq.query.form.prefix', `A form of kind ${kind} applies to annotations`, { address });
  });
}

/** A form's evaluation over its prefix. */
export type Resolve = (ctx: Ctx, form: Node, prefix: Named) => QueryResult;
