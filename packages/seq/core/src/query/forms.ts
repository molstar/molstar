/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The first forms as `repr` form resolvers (_spec/README.md "Queries", Forms): each receives what its prefix names
// and returns one of the three result types, or issues. Location forms read systems (`range`, `at`, `all`, `rect`),
// element forms annotations (`item`, `items`), `node` and `chain` the model's data-made nodes; `extent` and
// `items :overlapping` convert between types explicitly; `and`, `or`, `not` combine results of one type. A nested form
// is evaluated against the same prefix (a resolver's, else child selection), a nested address from the root. Every
// evaluation is bounded (the caller's limit, signal and deadline, one visit counter across nested forms) and staged:
// over data read by region it says what is pending.

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
  printAddress,
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
import { hasLabelSuffix, isNoInsertionCode, readLabel } from './labels.js';
import {
  EmptyLocation,
  intersectLocations,
  isEmptyLocation,
  isPosition,
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
  anyLabelIn,
  labelIndex,
  type NumberingOf,
  nearestLabels,
  numberingIds,
  numberingPending,
  snapLabel,
} from './numbering.js';
import { fromNamed, namedType, type QueryResult, toNamed } from './result.js';

/** The forms' kinds, in the order the docs list them. */
export const QueryFormKinds = [
  'range',
  'at',
  'all',
  'rect',
  'row',
  'item',
  'items',
  'node',
  'chain',
  'extent',
  'and',
  'or',
  'not',
] as const;

/** An evaluation failed: its issues (at least one error). Thrown inside a resolver, caught at its edge. */
class Failure {
  constructor(readonly issues: readonly Issue[]) {}
}

const fail = (code: QueryIssueCode, message: string, place?: IssuePlace): never => {
  throw new Failure([issueOf(QueryIssues, code, message, place)]);
};

/** What one evaluation shares across its nested forms, nested addresses included. */
interface Shared {
  /** Visits since the signal and the deadline were last read. */
  visits: number;
  /** Inside a rect's side: a range may not wrap (EC-10). */
  noWrap: boolean;
}

/**
 * The evaluation running now, if any: a nested address calls `selectAddress`, which calls the table's resolvers again;
 * they join the outer evaluation's counter instead of starting their own. Evaluation is synchronous, so one slot does.
 */
let active: Shared | undefined;

/** One resolver call's state: its input, the issues it carries, what is pending, the shared counter. */
interface Ctx {
  readonly model: QueryModel;
  readonly input: FormResolverInput;
  readonly resolvers: Readonly<Record<string, FormResolver>>;
  readonly issues: Issue[];
  readonly shared: Shared;
  complete: boolean;
  pending: Location;
}

/** Counts a visit; every CancelEvery reads the signal and the deadline. */
function tick(ctx: Ctx) {
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
  if (r.error === 'int32') fail('seq.query.label.int32', 'A label’s integer is outside the 32-bit range', { address });
  return fail('seq.query.label.syntax', 'A label is not an integer or a word', {
    address,
    hint: r.error === 'fraction' ? 'a .n label is quoted: "100.1"' : QueryIssues.codes['seq.query.label.syntax'].hint,
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

/**
 * The system's numbering of an id; undefined when it is declared and its labels have not arrived (a source's): the
 * form then holds nothing and waits, the whole system pending (`seq.query.numbering.pending`), never unknown.
 */
function numberingOf(ctx: Ctx, system: QuerySystem, id: string): NumberingOf | undefined {
  const found = findNumbering(system, id);
  if (found) return found;
  if (numberingPending(system, id)) {
    ctx.issues.push(
      issueOf(QueryIssues, 'seq.query.numbering.pending', 'The numbering’s labels have not arrived', {
        address: system.address,
      }),
    );
    stage(ctx, rangesLocation(system.address, [0, system.length]));
    return undefined;
  }
  return fail('seq.query.numbering.unknown', 'The system has no numbering of this id', {
    address: system.address,
    hint: `numberings: ${numberingIds(system).join(', ')}`,
  });
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
  const n = numberingOf(ctx, system, numberingId(form));
  if (!n) return EmptyLocation;
  const fromV = param(form, 'from');
  const toV = param(form, 'to');
  const from = fromV === undefined ? undefined : label(fromV, address);
  const to = toV === undefined ? undefined : label(toV, address);
  if (!system.length) return EmptyLocation;
  const step = () => tick(ctx);
  // snap first, then one rule in index order for exact and snapped ends alike: labels need not rise with the index
  let f = from === undefined ? 0 : labelIndex(n, system, from);
  let t = to === undefined ? system.length - 1 : labelIndex(n, system, to);
  const snapped = f < 0 || t < 0;
  if (snapped) {
    ctx.issues.push(
      issueOf(QueryIssues, 'seq.query.label.snapped', 'A range end the numbering lacks snapped inward', { address }),
    );
    if (f < 0) f = snapLabel(n, system, from!, 'start', step);
    if (t < 0) t = snapLabel(n, system, to!, 'end', step);
  }
  // nothing to snap to, or ends snapped past each other over a gap with no label in it
  if (f < 0 || t < 0 || (snapped && f > t && !anyLabelIn(n, t + 1, f, step))) {
    ctx.issues.push(issueOf(QueryIssues, 'seq.query.range.empty', 'A range holds no position', { address }));
    return EmptyLocation;
  }
  if (f <= t) return rangesLocation(address, [f, t + 1]);
  if (!system.cyclic)
    return fail('seq.query.range.backwards', 'A range ends before it starts', {
      address,
      hint: QueryIssues.codes['seq.query.range.backwards'].hint,
    });
  if (ctx.shared.noWrap)
    return fail('seq.query.range.wrap', 'A rect’s side wraps past the origin', {
      address,
      hint: QueryIssues.codes['seq.query.range.wrap'].hint,
    });
  if (t + 1 >= f)
    return fail('seq.query.range.lap', 'A wrapping range covers a whole lap', {
      address,
      hint: QueryIssues.codes['seq.query.range.lap'].hint,
    });
  // a wrapped range is two ranges: to the end, and from the origin
  return rangesLocation(address, [f, system.length, 0, t + 1]);
}

function atOn(ctx: Ctx, form: Node, system: QuerySystem): Location {
  const address = system.address;
  checkParams(form, ['label', 'numbering', 'ins']);
  const n = numberingOf(ctx, system, numberingId(form));
  if (!n) return EmptyLocation;
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
      hint: hintOfLabels(nearestLabels(n, system, l, 6, () => tick(ctx))),
    });
  return rangesLocation(address, [i, i + 1]);
}

const hintOfLabels = (labels: readonly string[]) =>
  labels.length ? `labels near it: ${labels.join(' ')}` : 'see the numbering’s labels';

/** Whole systems; a length past int32 is refused (a Location holds int32 positions). */
function allOn(form: Node, s: QuerySystem): Location {
  checkParams(form, []);
  if (!isPosition(s.length) || s.length < 0)
    fail('seq.query.position.invalid', 'A system’s length is not a 32-bit position', { address: s.address });
  return rangesLocation(s.address, [0, s.length]);
}

/** A rect's side: one system's positions, evaluated with wrapping refused (EC-10: rect sides never wrap). */
function rectSide(ctx: Ctx, value: Value, prefix: Named): Location {
  const before = ctx.shared.noWrap;
  ctx.shared.noWrap = true;
  try {
    const l = asType(evaluate(ctx, value, prefix), 'location', 'rect').location;
    if (l.parts.length !== 1 || l.rects.length)
      fail('seq.query.form.type', 'A rect’s side is the positions of one system', {
        hint: '(rect :a (range :from 10 :to 25) :b /data/x/(range :from 40 :to 55))',
      });
    return l;
  } finally {
    ctx.shared.noWrap = before;
  }
}

const rect: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['a', 'b']);
  const a = rectSide(ctx, required(form, 'a'), prefix).parts[0];
  const b = rectSide(ctx, required(form, 'b'), prefix).parts[0];
  return { type: 'location', location: location({ rects: [{ a, b }] }) };
};

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
  const anns = annotationsOf(ctx, prefix, 'item');
  for (const k of keys)
    if (k !== 'id' && !anns.some((a) => a.fields.includes(k)))
      fail('seq.query.form.param', 'A form of kind item names a field its annotation does not have', {
        hint: `params: ${['id', ...new Set(anns.flatMap((a) => a.fields))].join(', ')}`,
      });
  const out: ElementKey[] = [];
  for (const a of anns) {
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

/**
 * The nodes of systems by their addresses, so a nested location form reads them as its prefix. An address may hold
 * forms (`/data/1tup/(chain :id A)`): it is selected, with the table, as the caller's are.
 */
function systemNodes(ctx: Ctx, addresses: readonly string[]): Named {
  const matches: Selected[] = [];
  for (const a of new Set(addresses)) {
    const { input } = ctx;
    const sel = selectAddress(input.root, a, { limit: 1, maxDepth: input.maxDepth, resolvers: ctx.resolvers });
    if (sel.matches.length !== 1)
      return fail('seq.query.form.prefix', 'An annotation’s target is not a node', { address: a });
    matches.push(sel.matches[0]);
  }
  return { kind: 'nodes', matches };
}

const items: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['overlapping']);
  const anns = annotationsOf(ctx, prefix, 'items');
  const over = param(form, 'overlapping');
  const targets = anns.map((a) => a.target);
  const loc =
    over === undefined
      ? undefined
      : asType(evaluate(ctx, over, systemNodes(ctx, targets)), 'location', 'items').location;
  if (loc?.rects.length)
    fail('seq.query.form.unsupported', 'Items overlap positions, not a rect', { hint: 'one of its sides' });
  // TODO: a location on another system reaches the targets by projection through the maps (the static model's paths);
  // until then it is reported, never silently read as overlapping nothing
  if (loc)
    for (const p of loc.parts)
      if (!targets.includes(p.system))
        ctx.issues.push(
          issueOf(QueryIssues, 'seq.query.location.unreachable', 'A location is on a system the items are not on', {
            address: p.system,
          }),
        );
  const out: ElementKey[] = [];
  for (const a of anns) {
    const ranges = loc && partOf(loc, a.target)?.ranges;
    if (loc && !ranges) continue;
    const length = ctx.model.system(a.target)?.length;
    const read = a.items(ranges && rangesLocation(a.target, ranges));
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      if (!ranges || overlapsRanges(ranges, it.start, it.end, length)) out.push({ owner: a.address, key: it.key });
    }
  }
  const e = elements(out);
  checkLimit(ctx, e.keys.length);
  return { type: 'elements', elements: e };
};

// ---- nodes, conversions, combinators ----

/** Nodes data makes under each prefix node, by one param, addressed by the form itself (no id names them). */
function heldNodes(
  key: 'key' | 'id',
  get: (model: QueryModel, owner: string, value: Atom) => Node | undefined,
): Resolve {
  return (ctx, form, prefix) => {
    checkParams(form, [key]);
    const v = required(form, key);
    if (v !== null && typeof v === 'object')
      return fail('seq.query.form.param', `A form of kind ${form.kind} takes an atom`);
    const out: Selected[] = [];
    for (const m of prefixNodes(prefix, form.kind)) {
      tick(ctx);
      const address = printAddress(m.address);
      const n = get(ctx.model, address, v as Atom);
      if (n) out.push({ node: n, address: [...m.address, form] });
      else ctx.issues.push(issueOf(QueryIssues, 'seq.query.node.missing', 'No node of this key', { address }));
    }
    checkLimit(ctx, out.length);
    return { type: 'nodes', nodes: out };
  };
}

const node = heldNodes('key', (model, owner, key) => model.node?.(owner, key));
/** A chain is a system data makes: `/data/1tup/(chain :id A)` is a prefix location forms read (the model's system). */
/**
 * `(row :key k)`: an alignment's rows by their keys (labels; a list for several) as positions of its rows system, so a
 * rect over rows × columns is written `(rect :a <alignment>/(row :key k) :b <alignment>/(range …))`.
 */
const row: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['key']);
  const key = required(form, 'key');
  const keys = (Array.isArray(key) ? key : [key]) as Atom[];
  const parts: { system: string; ranges: number[] }[] = [];
  for (const m of prefixNodes(prefix, 'row')) {
    const at = printAddress(m.address);
    const rows = at.endsWith('/(rows)') ? at : `${at}/(rows)`;
    const ranges: number[] = [];
    // rows read from a source and not held yet: the form waits over the whole rows system (pending, never unknown)
    if (ctx.model.rowsPending?.(rows)) {
      ctx.issues.push(
        issueOf(QueryIssues, 'seq.query.row.pending', 'The alignment’s rows have not arrived', { address: at }),
      );
      const n = ctx.model.system(rows)?.length;
      if (n !== undefined) stage(ctx, rangesLocation(rows, [0, n]));
      continue;
    }
    for (const k of keys) {
      const i = ctx.model.row?.(rows, k);
      if (i === undefined) return fail('seq.query.row.unknown', `No row has the key ${String(k)}`, { address: at });
      ranges.push(i, i + 1);
    }
    parts.push({ system: rows, ranges });
  }
  return { type: 'location', location: location({ parts }) };
};

const chain = heldNodes('id', (model, owner, id) => model.chain?.(owner, id));

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
      if (!keys.has(it.key)) continue;
      if (!isPosition(it.start) || !isPosition(it.end) || it.start > it.end)
        fail('seq.query.position.invalid', 'An item’s extent is not 32-bit positions with start ≤ end', {
          address: owner,
        });
      ranges.push(it.start, it.end);
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

/**
 * `and`: the intersection of its operands; a location less each `(not :of <location>)` among them, subtracted (what a
 * complement within the systems would intersect to, and the one way to take a rect from a rect: rects subtract, a rect
 * has no complement; _spec/README.md "Queries", Location).
 */
const and: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['of']);
  const of = required(form, 'of');
  if (!Array.isArray(of) || !of.length) return combine(ctx, 'and', operands(ctx, form, prefix));
  const kept: QueryResult[] = [];
  const less: Location[] = [];
  for (const v of of as readonly Value[]) {
    const inner = isNode(v) && v.kind === 'not' ? v.params?.of : undefined;
    const r = inner !== undefined ? evaluate(ctx, inner, prefix) : undefined;
    if (r?.type === 'location') less.push(r.location);
    else kept.push(evaluate(ctx, v, prefix));
  }
  if (!less.length) return combine(ctx, 'and', kept);
  if (!kept.length)
    return fail('seq.query.form.unsupported', 'A location less another needs one to take it from', {
      hint: '(and :of [(rect :a … :b …) (not :of (rect :a … :b …))])',
    });
  const acc = combine(ctx, 'and', kept);
  if (acc.type !== 'location') return fail('seq.query.form.type', 'A form of kind and mixes result types');
  let loc = acc.location;
  for (const l of less) loc = subtractLocations(loc, l);
  return { type: 'location', location: loc };
};
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
  at: locationForm(atOn),
  all: locationForm((_, form, s) => allOn(form, s)),
  rect,
  row,
  item,
  items,
  node,
  chain,
  extent,
  and,
  or,
  not,
};

/**
 * The form-resolver table of the query forms over a model, for `selectAddress` (`SelectOptions.resolvers`). A caller
 * may add its own kinds beside them (`extra`), and nested forms see the whole table; an extra kind that is a query
 * form's is a programmer error and throws: a form's meaning is never replaced silently.
 */
export function queryResolvers(
  model: QueryModel,
  extra?: Readonly<Record<string, FormResolver>>,
): Readonly<Record<string, FormResolver>> {
  const table: Record<string, FormResolver> = { ...extra };
  for (const kind of QueryFormKinds) {
    if (extra && Object.hasOwn(extra, kind)) throw new Error(`A resolver of kind ${kind} clashes with a query form`);
    const resolve = Forms[kind];
    table[kind] = (input): FormResolverResult => {
      const outer = active;
      const shared = outer ?? { visits: 0, noWrap: false };
      active = shared;
      const ctx: Ctx = { model, input, resolvers: table, issues: [], shared, complete: true, pending: EmptyLocation };
      try {
        const result = toNamed(resolve(ctx, input.form, input.prefix));
        return ctx.complete
          ? { result, issues: ctx.issues }
          : { result, issues: ctx.issues, complete: false, pending: ctx.pending };
      } catch (e) {
        if (e instanceof Failure) return { issues: e.issues };
        throw e;
      } finally {
        active = outer;
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
