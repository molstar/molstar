/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The first forms as `repr` form resolvers (_spec/README.md "Queries", Forms): each receives what its prefix names
// and returns one of the three result types, or issues. Location forms read systems (`range`, `at`, `all`), element
// forms annotations (`item`, `items`), `node` the model's dynamic nodes; `extent` and `items :overlapping` convert
// between types explicitly; `and`, `or`, `not` combine results of one type. A nested form is evaluated against the
// same prefix (a resolver's, else child selection), a nested address from the root. Every evaluation is bounded (the
// caller's limit, signal and deadline) and staged: over data read by region it says what is pending.

import {
  type Atom,
  CancelEvery,
  type FormResolver,
  type FormResolverInput,
  type FormResolverResult,
  type Issue,
  type IssuePlace,
  issueOf,
  isNode,
  matchesForm,
  type Named,
  type Node,
  own,
  parseAddress,
  printAddress,
  resolveAddress,
  type Selected,
  type SelectOptions,
  selectAddress,
  TreeIssues,
  type Value,
} from '@molstar/core/repr/tree';
import {
  type ElementKey,
  type Elements,
  EmptyElements,
  elements,
  intersectElements,
  subtractElements,
  unionElements,
} from './elements.js';
import { type QueryIssueCode, QueryIssues } from './issues.js';
import { compareLabels, hasLabelSuffix, isNoInsertionCode, readLabel } from './labels.js';
import {
  EmptyLocation,
  intersectLocations,
  isEmptyLocation,
  type Location,
  location,
  overlapsRanges,
  partOf,
  rangesLocation,
  subtractLocations,
  unionLocations,
} from './location.js';
import type { QueryAnnotation, QueryModel, QuerySystem } from './model.js';
import {
  DefaultNumbering,
  findNumbering,
  labelIndex,
  type NumberingOf,
  nearestLabels,
  numberingIds,
  snapLabel,
} from './numbering.js';
import { fromNamed, namedType, type QueryResult, toNamed } from './result.js';

/** The forms' kinds, in the order the docs list them. */
export const QueryFormKinds = ['range', 'at', 'all', 'item', 'items', 'node', 'extent', 'and', 'or', 'not'] as const;

/** An evaluation failed: its issues (at least one error). Thrown inside a resolver, caught at its edge. */
class Failure {
  constructor(readonly issues: readonly Issue[]) {}
}

const fail = (code: QueryIssueCode, message: string, place?: IssuePlace): never => {
  throw new Failure([issueOf(QueryIssues, code, message, place)]);
};

/** One resolver call's state: its input, the issues it carries, what is pending, a visit counter. */
interface Ctx {
  readonly model: QueryModel;
  readonly input: FormResolverInput;
  readonly resolvers: Readonly<Record<string, FormResolver>>;
  readonly issues: Issue[];
  complete: boolean;
  pending: Location;
  visits: number;
}

/** Reads the signal and the deadline every CancelEvery visits. */
function tick(ctx: Ctx) {
  if (++ctx.visits < CancelEvery) return;
  ctx.visits = 0;
  const { signal, deadline } = ctx.input;
  if (signal?.aborted)
    throw new Failure([issueOf(TreeIssues, 'repr.tree.address.cancelled', 'The selection was cancelled')]);
  if (deadline !== undefined && performance.now() > deadline)
    throw new Failure([
      issueOf(TreeIssues, 'repr.tree.address.deadline_exceeded', 'The selection ran past its deadline'),
    ]);
}

function checkLimit(ctx: Ctx, count: number) {
  const { limit } = ctx.input;
  if (count <= limit) return;
  throw new Failure([
    issueOf(TreeIssues, 'repr.tree.address.limit_exceeded', `More than ${limit} matches`, {
      hint: 'narrow the address, or raise the limit',
    }),
  ]);
}

function stage(ctx: Ctx, pending: Location | undefined) {
  if (!pending || isEmptyLocation(pending)) return;
  ctx.complete = false;
  ctx.pending = unionLocations(ctx.pending, pending);
}

// ---- params ----

function checkParams(form: Node, allowed: readonly string[]) {
  for (const key of Object.keys(form.params ?? {}))
    if (!allowed.includes(key))
      fail('seq.query.form.param', `A form of kind ${form.kind} names a param it does not take`, {
        hint: allowed.length ? `params: ${allowed.join(', ')}` : 'no params',
      });
}

const param = (form: Node, key: string): Value | undefined => (form.params ? own(form.params, key) : undefined);

function required(form: Node, key: string): Value {
  const v = param(form, key);
  if (v === undefined) fail('seq.query.form.param', `A form of kind ${form.kind} needs :${key}`);
  return v as Value;
}

/** A numbering's id: a word (a number is spelled as written). */
function numberingId(form: Node): string {
  const v = param(form, 'numbering');
  if (v === undefined) return DefaultNumbering;
  if (typeof v === 'string' || typeof v === 'number') return `${v}`;
  return fail('seq.query.form.param', 'A numbering is named by a word');
}

function label(value: Value, address: string): string {
  const r = readLabel(value);
  if ('label' in r) return r.label;
  return r.error === 'int32'
    ? fail('seq.query.label.int32', 'A label’s integer is outside the 32-bit range', { address })
    : fail('seq.query.label.syntax', 'A label is not an integer or a word', {
        address,
        hint: QueryIssues.codes['seq.query.label.syntax'].hint,
      });
}

// ---- what the prefix names ----

function prefixNodes(prefix: Named, kind: string): readonly Selected[] {
  if (prefix.kind === 'nodes') return prefix.matches;
  return fail('seq.query.form.prefix', `A form of kind ${kind} applies to nodes, and its prefix names ${prefix.type}`);
}

function systemsOf(ctx: Ctx, prefix: Named, kind: string): QuerySystem[] {
  return prefixNodes(prefix, kind).map((m) => {
    const address = printAddress(m.address);
    const s = ctx.model.system(address);
    return s ?? fail('seq.query.form.prefix', `A form of kind ${kind} applies to systems`, { address });
  });
}

function annotationsOf(ctx: Ctx, prefix: Named, kind: string): QueryAnnotation[] {
  return prefixNodes(prefix, kind).map((m) => {
    const address = printAddress(m.address);
    const a = ctx.model.annotation(address);
    return a ?? fail('seq.query.form.prefix', `A form of kind ${kind} applies to annotations`, { address });
  });
}

function numberingOf(system: QuerySystem, id: string): NumberingOf {
  return (
    findNumbering(system, id) ??
    fail('seq.query.numbering.unknown', 'The system has no numbering of this id', {
      address: system.address,
      hint: `numberings: ${numberingIds(system).join(', ')}`,
    })
  );
}

// ---- nested values ----

/** A nested value's result: a form against the prefix (its resolver, else child selection), an address from the root. */
function evaluate(ctx: Ctx, value: Value, prefix: Named): QueryResult {
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
    return carry(ctx, sel, sel.result ?? { kind: 'nodes', matches: sel.matches });
  }
  return fail('seq.query.form.param', 'A form takes a nested form or an address here');
}

/** Carries a nested result's issues and pending regions; its errors fail the outer form. */
function carry(
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

function asType<T extends QueryResult['type']>(
  r: QueryResult,
  type: T,
  kind: string,
): Extract<QueryResult, { type: T }> {
  if (r.type !== type) fail('seq.query.form.type', `A form of kind ${kind} takes ${type}, and is given ${r.type}`);
  return r as Extract<QueryResult, { type: T }>;
}

// ---- location forms ----

function rangeOn(ctx: Ctx, form: Node, system: QuerySystem): Location {
  const address = system.address;
  if (param(form, 'ins') !== undefined)
    fail('seq.query.label.insertion', 'A range takes no insertion code', {
      address,
      hint: QueryIssues.codes['seq.query.label.insertion'].hint,
    });
  checkParams(form, ['from', 'to', 'numbering']);
  const n = numberingOf(system, numberingId(form));
  const fromV = param(form, 'from');
  const toV = param(form, 'to');
  const from = fromV === undefined ? undefined : label(fromV, address);
  const to = toV === undefined ? undefined : label(toV, address);
  if (!system.length) return EmptyLocation;
  let f = from === undefined ? 0 : labelIndex(n, system, from);
  let t = to === undefined ? system.length - 1 : labelIndex(n, system, to);
  const snapped = f < 0 || t < 0;
  if (snapped) {
    if (from !== undefined && to !== undefined && compareLabels(from, to) > 0)
      fail('seq.query.range.backwards', 'A range ends before it starts', {
        address,
        hint: QueryIssues.codes['seq.query.range.backwards'].hint,
      });
    ctx.issues.push(
      issueOf(QueryIssues, 'seq.query.label.snapped', 'A range end the numbering lacks snapped inward', { address }),
    );
    if (f < 0) f = snapLabel(n, system, from!, 'start');
    if (t < 0) t = snapLabel(n, system, to!, 'end');
    if (f < 0 || t < 0 || f > t) {
      ctx.issues.push(issueOf(QueryIssues, 'seq.query.range.empty', 'A range holds no position', { address }));
      return EmptyLocation;
    }
    return rangesLocation(address, [f, t + 1]);
  }
  if (f <= t) return rangesLocation(address, [f, t + 1]);
  if (!system.cyclic)
    return fail('seq.query.range.backwards', 'A range ends before it starts', {
      address,
      hint: QueryIssues.codes['seq.query.range.backwards'].hint,
    });
  if (t + 1 >= f)
    return fail('seq.query.range.lap', 'A wrapping range covers a whole lap', {
      address,
      hint: QueryIssues.codes['seq.query.range.lap'].hint,
    });
  // a wrapped range is two ranges: to the end, and from the origin
  return rangesLocation(address, [f, system.length, 0, t + 1]);
}

function atOn(form: Node, system: QuerySystem): Location {
  const address = system.address;
  checkParams(form, ['label', 'numbering', 'ins']);
  const n = numberingOf(system, numberingId(form));
  let l = label(required(form, 'label'), address);
  const ins = param(form, 'ins');
  if (ins !== undefined) {
    if (typeof ins !== 'string' || !(isNoInsertionCode(ins) || /^[A-Za-z]+$/.test(ins)))
      fail('seq.query.label.syntax', 'An insertion code is not letters', { address });
    if (!isNoInsertionCode(ins as string)) {
      if (hasLabelSuffix(l))
        fail('seq.query.label.insertion', 'A label with a suffix takes no insertion code', {
          address,
          hint: QueryIssues.codes['seq.query.label.insertion'].hint,
        });
      l += ins as string;
    }
  }
  const i = labelIndex(n, system, l);
  if (i < 0)
    fail('seq.query.label.unknown', 'The numbering has no position of this label', {
      address,
      hint: `labels near it: ${nearestLabels(n, system, l).join(' ')}`,
    });
  return rangesLocation(address, [i, i + 1]);
}

function locationForm(on: (ctx: Ctx, form: Node, system: QuerySystem) => Location): Resolve {
  return (ctx, form, prefix) => {
    let out = EmptyLocation;
    for (const s of systemsOf(ctx, prefix, form.kind)) out = unionLocations(out, on(ctx, form, s));
    return { type: 'location', location: out };
  };
}

// ---- element forms ----

/** Whether an item's value fits a param: the same atom, or an atom held by a list param. */
function fits(v: Value, x: Atom | undefined): boolean {
  if (x === undefined) return false;
  if (Array.isArray(v)) return v.some((y) => Object.is(y, x));
  return Object.is(v, x);
}

const item: Resolve = (ctx, form, prefix) => {
  const params = form.params ?? {};
  const keys = Object.keys(params);
  if (!keys.length) fail('seq.query.form.param', 'A form of kind item needs :id or a field', { hint: '(item :id k)' });
  const out: ElementKey[] = [];
  for (const a of annotationsOf(ctx, prefix, 'item')) {
    const read = a.items();
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      let ok = true;
      for (const k of keys) if (!(ok = fits(params[k], k === 'id' ? it.key : it.fields?.[k]))) break;
      if (ok) out.push({ owner: a.address, key: it.key });
    }
  }
  const e = elements(out);
  checkLimit(ctx, e.keys.length);
  return { type: 'elements', elements: e };
};

/** The nodes of systems by their addresses, so a nested location form reads them as its prefix. */
function systemNodes(ctx: Ctx, addresses: readonly string[]): Named {
  const matches: Selected[] = [];
  for (const a of new Set(addresses)) {
    const p = parseAddress(a).pattern;
    const ids = p?.segments.every((s) => typeof s === 'string') ? (p.segments as string[]) : undefined;
    const node = ids && resolveAddress(ctx.input.root, ids);
    if (!node || !ids) return fail('seq.query.form.prefix', 'An annotation’s target is not a node', { address: a });
    matches.push({ node, address: ids });
  }
  return { kind: 'nodes', matches };
}

const items: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['overlapping']);
  const anns = annotationsOf(ctx, prefix, 'items');
  const over = param(form, 'overlapping');
  const loc =
    over === undefined
      ? undefined
      : asType(
          evaluate(
            ctx,
            over,
            systemNodes(
              ctx,
              anns.map((a) => a.target),
            ),
          ),
          'location',
          'items',
        ).location;
  const out: ElementKey[] = [];
  for (const a of anns) {
    const ranges = loc && partOf(loc, a.target)?.ranges;
    if (loc && !ranges) continue;
    const read = a.items(ranges && rangesLocation(a.target, ranges));
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      if (!ranges || overlapsRanges(ranges, it.start, it.end)) out.push({ owner: a.address, key: it.key });
    }
  }
  const e = elements(out);
  checkLimit(ctx, e.keys.length);
  return { type: 'elements', elements: e };
};

// ---- nodes, conversions, combinators ----

const node: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['key']);
  const key = required(form, 'key');
  if (key !== null && typeof key === 'object') return fail('seq.query.form.param', 'A node key is an atom');
  const out: Selected[] = [];
  for (const m of prefixNodes(prefix, 'node')) {
    const address = printAddress(m.address);
    const n = ctx.model.node?.(address, key as Atom);
    if (n) out.push({ node: n, address: [...m.address, form] });
    else ctx.issues.push(issueOf(QueryIssues, 'seq.query.node.missing', 'No node of this key', { address }));
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

/** The extent of elements: the union of their items' ranges on each annotation's target. */
function extentOf(ctx: Ctx, e: Elements): Location {
  const byOwner = new Map<string, Set<string | number>>();
  for (const k of e.keys) {
    let s = byOwner.get(k.owner);
    if (!s) byOwner.set(k.owner, (s = new Set()));
    s.add(k.key);
  }
  const parts: { system: string; ranges: number[] }[] = [];
  for (const [owner, keys] of byOwner) {
    const a = ctx.model.annotation(owner);
    if (!a) return fail('seq.query.form.prefix', 'An element’s owner is not an annotation', { address: owner });
    const read = a.items();
    stage(ctx, read.pending);
    const ranges: number[] = [];
    for (const it of read.items) {
      tick(ctx);
      if (keys.has(it.key)) ranges.push(it.start, it.end);
    }
    parts.push({ system: a.target, ranges });
  }
  return location({ parts });
}

const extent: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['of']);
  const e = asType(evaluate(ctx, required(form, 'of'), prefix), 'elements', 'extent').elements;
  return { type: 'location', location: extentOf(ctx, e) };
};

const nodeKey = (s: Selected) => printAddress(s.address);

function combineNodes(a: readonly Selected[], b: readonly Selected[], keep: (inA: boolean, inB: boolean) => boolean) {
  const inA = new Set(a.map(nodeKey));
  const inB = new Set(b.map(nodeKey));
  const seen = new Set<string>();
  const out: Selected[] = [];
  for (const s of [...a, ...b]) {
    const k = nodeKey(s);
    if (seen.has(k) || !keep(inA.has(k), inB.has(k))) continue;
    seen.add(k);
    out.push(s);
  }
  return out;
}

function combine(ctx: Ctx, kind: 'and' | 'or', rs: readonly QueryResult[]): QueryResult {
  const type = rs[0].type;
  for (const r of rs) if (r.type !== type) fail('seq.query.form.type', `A form of kind ${kind} mixes result types`);
  const or = kind === 'or';
  let acc = rs[0];
  for (let i = 1; i < rs.length; i++) {
    const r = rs[i];
    if (acc.type === 'nodes' && r.type === 'nodes')
      acc = { type: 'nodes', nodes: combineNodes(acc.nodes, r.nodes, or ? (x, y) => x || y : (x, y) => x && y) };
    else if (acc.type === 'elements' && r.type === 'elements')
      acc = { type: 'elements', elements: (or ? unionElements : intersectElements)(acc.elements, r.elements) };
    else if (acc.type === 'location' && r.type === 'location')
      acc = { type: 'location', location: (or ? unionLocations : intersectLocations)(acc.location, r.location) };
  }
  if (acc.type === 'nodes') checkLimit(ctx, acc.nodes.length);
  if (acc.type === 'elements') checkLimit(ctx, acc.elements.keys.length);
  return acc;
}

function operands(ctx: Ctx, form: Node, prefix: Named): QueryResult[] {
  checkParams(form, ['of']);
  const of = required(form, 'of');
  if (!Array.isArray(of) || !of.length)
    return fail('seq.query.form.param', `A form of kind ${form.kind} takes a list of forms`, {
      hint: `(${form.kind} :of [(range :from 1 :to 9) (at :label 20)])`,
    });
  return (of as readonly Value[]).map((v) => evaluate(ctx, v, prefix));
}

const and: Resolve = (ctx, form, prefix) => combine(ctx, 'and', operands(ctx, form, prefix));
const or: Resolve = (ctx, form, prefix) => combine(ctx, 'or', operands(ctx, form, prefix));

/** The complement within what the prefix names: its systems' positions, its annotations' items, its nodes' children. */
const not: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['of']);
  const r = evaluate(ctx, required(form, 'of'), prefix);
  if (r.type === 'location') {
    if (r.location.rects.length)
      return fail('seq.query.form.unsupported', 'A rect has no complement that is a location');
    const all = location({
      parts: systemsOf(ctx, prefix, 'not').map((s) => ({ system: s.address, ranges: [0, s.length] })),
    });
    return { type: 'location', location: subtractLocations(all, r.location) };
  }
  if (r.type === 'elements') {
    let all = EmptyElements;
    for (const a of annotationsOf(ctx, prefix, 'not')) {
      const read = a.items();
      stage(ctx, read.pending);
      all = unionElements(all, elements(read.items.map((it) => ({ owner: a.address, key: it.key }))));
    }
    const e = subtractElements(all, r.elements);
    checkLimit(ctx, e.keys.length);
    return { type: 'elements', elements: e };
  }
  const children: Selected[] = [];
  for (const m of prefixNodes(prefix, 'not'))
    for (const c of m.node.children ?? [])
      if (c.id !== undefined) children.push({ node: c, address: [...m.address, c.id] });
  const out = combineNodes(children, r.nodes, (x, y) => x && !y);
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

// ---- the table ----

type Resolve = (ctx: Ctx, form: Node, prefix: Named) => QueryResult;

const Forms: Readonly<Record<(typeof QueryFormKinds)[number], Resolve>> = {
  range: locationForm(rangeOn),
  at: locationForm((_, form, s) => atOn(form, s)),
  all: locationForm((_, form, s) => {
    checkParams(form, []);
    return rangesLocation(s.address, [0, s.length]);
  }),
  item,
  items,
  node,
  extent,
  and,
  or,
  not,
};

/**
 * The form-resolver table of the query forms over a model, for `selectAddress` (`SelectOptions.resolvers`). A caller
 * may add its own kinds beside them; nested forms see the whole table given as `extra`.
 */
export function queryResolvers(
  model: QueryModel,
  extra?: Readonly<Record<string, FormResolver>>,
): Readonly<Record<string, FormResolver>> {
  const table: Record<string, FormResolver> = { ...extra };
  for (const kind of QueryFormKinds) {
    const resolve = Forms[kind];
    table[kind] = (input): FormResolverResult => {
      const ctx: Ctx = {
        model,
        input,
        resolvers: table,
        issues: [],
        complete: true,
        pending: EmptyLocation,
        visits: 0,
      };
      try {
        const result = toNamed(resolve(ctx, input.form, input.prefix));
        return ctx.complete
          ? { result, issues: ctx.issues }
          : { result, issues: ctx.issues, complete: false, pending: ctx.pending };
      } catch (e) {
        if (e instanceof Failure) return { issues: e.issues };
        throw e;
      }
    };
  }
  return table;
}

/** A query's outcome: its result (none on an error), every issue, and over data read by region what is pending. */
export interface Evaluation {
  readonly result?: QueryResult;
  readonly issues: readonly Issue[];
  /** False when the result covers only the regions loaded; `pending` is then what is not. */
  readonly complete: boolean;
  readonly pending: Location;
}

export interface QueryOptions {
  readonly limit: number;
  readonly maxDepth?: number;
  readonly signal?: AbortSignal;
  readonly deadline?: number;
}

/** Evaluates an address with the query forms over a model: `selectAddress` with `queryResolvers`, typed. */
export function evaluateQuery(root: Node, model: QueryModel, address: string, options: QueryOptions): Evaluation {
  const sel = selectAddress(root, address, { ...options, resolvers: queryResolvers(model) });
  const failed = sel.issues.some((i) => i.severity === 'error');
  const result = failed
    ? undefined
    : sel.result
      ? fromNamed(sel.result)
      : { type: 'nodes' as const, nodes: sel.matches };
  const complete = sel.complete !== false;
  return {
    ...(result ? { result } : {}),
    issues: sel.issues,
    complete,
    pending: complete ? EmptyLocation : ((sel.pending as Location | undefined) ?? EmptyLocation),
  };
}
