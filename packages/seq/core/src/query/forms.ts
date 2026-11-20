/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The first forms as `repr` form resolvers (_spec/README.md "Queries", Forms): each receives what its prefix names
// and returns one of the three result types, or issues. Location forms read systems (`range`, `at`, `all`, `rect`),
// key forms annotations (`item`, `items`), `node` and `chain` the model's data-made nodes; `extent` and
// `items :overlapping` convert between types explicitly; `merge`, `intersect`, `except` combine results of one type
// (`combinators.ts`); the structure language's forms are `elements.ts`'s and `structure-forms.ts`'s. A nested form
// is evaluated against the same prefix (a resolver's, else child selection), a nested address from the root. Every
// evaluation is bounded (the caller's limit, signal and deadline, one visit counter across nested forms) and staged:
// over data read by region it says what is pending.

import {
  type Atom,
  type FormResolver,
  type FormResolverResult,
  type Issue,
  issueOf,
  isNode,
  type Named,
  type Node,
  parseAddress,
  printAddress,
  type Selected,
  selectAddress,
  type Value,
} from '@molstar/core/repr/tree';
import { type Key, keySet, subtractKeys } from './keys.js';
import {
  annotationsOf,
  checkLimit,
  checkParams,
  type Ctx,
  Failure,
  fail,
  param,
  prefixNodes,
  type Resolve,
  required,
  type Shared,
  stage,
  systemsOf,
  tick,
} from './form-kit.js';
import { Combinators } from './combinators.js';
import { elements } from './elements-form.js';
import { NoQueryExtensions, type QueryExtensions, type QueryFormEntry } from './extensions.js';
import { StructureForms } from './structure-forms.js';
import { asType, evaluate } from './evaluate.js';
import { keysExtent } from './key-extent.js';
import { QueryIssues } from './issues.js';
import { MoleculeForms } from './molecule-forms.js';
import { hasLabelSuffix, isNoInsertionCode, readLabel } from './labels.js';
import {
  EmptyLocation,
  intersectLocations,
  isPosition,
  type Location,
  location,
  locationsEqual,
  overlapsRanges,
  partOf,
  rangesLocation,
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
import { fromNamed, type QueryResult, toNamed } from './result.js';

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
  'merge',
  'intersect',
  'except',
  'elements',
  'groups',
  'filter',
  'expand',
  'linked',
  'molecule',
  'atoms',
] as const;

/**
 * The evaluation running now, if any: a nested address calls `selectAddress`, which calls the table's resolvers again;
 * they join the outer evaluation's counter instead of starting their own. Evaluation is synchronous, so one slot does.
 */
let active: Shared | undefined;

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
  const out: Key[] = [];
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
  const e = keySet(out);
  checkLimit(ctx, e.keys.length);
  return { type: 'keys', keys: e };
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

/**
 * The systems an annotation's rows are on: its `systems` when it has them (a property its interface does not name: the
 * resolver's pairs annotation sets it, its two ends' systems; internal to core), else `[target]`. Why: reach is
 * decided without reading rows (an empty pairs annotation on two systems reaches both), and no public field is added.
 */
function systemsOfAnnotation(a: QueryAnnotation): readonly string[] {
  return (a as { readonly systems?: readonly string[] }).systems ?? [a.target];
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
  // until then it is reported, never silently read as overlapping nothing. Reach is decided per annotation, by its own
  // systems (a pairs annotation's two, `systemsOfAnnotation`): read with the whole location when a part is on one of
  // them; a part no annotation reaches is said unreachable.
  const reach = anns.map(systemsOfAnnotation);
  if (loc)
    for (const p of loc.parts)
      if (!reach.some((r) => r.includes(p.system)))
        ctx.issues.push(
          issueOf(QueryIssues, 'seq.query.location.unreachable', 'A location is on a system the items are not on', {
            address: p.system,
          }),
        );
  const out: Key[] = [];
  for (let i = 0; i < anns.length; i++) {
    const a = anns[i];
    const systems = reach[i];
    if (loc && !loc.parts.some((p) => systems.includes(p.system))) continue;
    // one system (a ranges annotation): its own part alone, as it was always asked; more (a pair's ends): the whole
    const read = a.items(loc && (systems.length === 1 ? rangesLocation(a.target, partOf(loc, a.target)!.ranges) : loc));
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      const system = it.system ?? a.target;
      const ranges = loc && partOf(loc, system)?.ranges;
      if (loc && !ranges) continue;
      if (!ranges || overlapsRanges(ranges, it.start, it.end, ctx.model.system(system)?.length))
        out.push({ owner: a.address, key: it.key });
    }
  }
  const e = keySet(out);
  checkLimit(ctx, e.keys.length);
  return { type: 'keys', keys: e };
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

const extent: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['of']);
  const e = asType(evaluate(ctx, required(form, 'of'), prefix), 'keys', 'extent').keys;
  return { type: 'location', location: keysExtent(ctx, e) };
};

// ---- the table ----

const Forms: Readonly<Record<(typeof QueryFormKinds)[number], Resolve>> = {
  range: locationForm(rangeOn),
  at: locationForm(atOn),
  all: locationForm((_, form, s) => allOn(form, s)),
  rect,
  row,
  item,
  items,
  node,
  extent,
  ...Combinators,
  elements,
  ...StructureForms,
  ...MoleculeForms,
};

/** A value's forms: itself and those nested in its params and args. */
function* formsIn(v: Value): Generator<Node> {
  if (Array.isArray(v)) for (const x of v as readonly Value[]) yield* formsIn(x);
  else if (isNode(v)) {
    yield v;
    for (const x of Object.values(v.params ?? {}) as Value[]) yield* formsIn(x);
    for (const x of v.args ?? []) yield* formsIn(x);
  }
}

/** The first namespaced kind of an address's forms no resolver has: `seq.query.extension.missing`, else undefined. */
function missingExtension(address: string, table: Readonly<Record<string, FormResolver>>, ext: QueryExtensions) {
  for (const s of parseAddress(address).pattern?.segments ?? [])
    if (typeof s !== 'string')
      for (const f of formsIn(s)) {
        const i = f.kind.indexOf(':');
        if (i > 0 && !Object.hasOwn(table, f.kind) && !ext.operators.has(f.kind) && f.kind !== 'fields')
          return issueOf(QueryIssues, 'seq.query.extension.missing', `No extension provides the form ${f.kind}`, {
            hint: `namespace: ${f.kind.slice(0, i)}`,
          });
      }
  return undefined;
}

/** A positional operand outside an expression form (`seq.query.form.param`): only operators, prop and value take one. */
function checkNoArgs(form: Node) {
  if (form.args?.length)
    fail('seq.query.form.param', `A form of kind ${form.kind} takes keyword params only`, {
      hint: 'positional operands are an expression’s: (filter :test (> (prop b-factor) 60))',
    });
}

/** An extension's form, checked by its role: its result its declared type; a filter's a subset of its prefix. */
function extensionForm(e: QueryFormEntry): Resolve {
  return (ctx, form, prefix) => {
    checkParams(form, e.params);
    const before: QueryResult =
      prefix.kind === 'nodes'
        ? { type: 'nodes', nodes: prefix.matches }
        : (fromNamed(prefix) ??
          fail('seq.query.form.prefix', `A form of kind ${form.kind} is given an unknown prefix`));
    const r = e.resolve({
      model: ctx.model,
      form,
      prefix: before,
      evaluate: (v) => evaluate(ctx, v, prefix),
      tick: () => tick(ctx),
      limit: ctx.input.limit,
    });
    if (r.type !== e.result)
      fail('seq.query.form.type', `A form of kind ${form.kind} gave ${r.type}, and declares ${e.result}`);
    if (e.role === 'filter' && !isSubset(r, before))
      fail('seq.query.form.type', `A filter of kind ${form.kind} gave what its prefix does not hold`);
    if (r.type === 'nodes') checkLimit(ctx, r.nodes.length);
    if (r.type === 'keys') checkLimit(ctx, r.keys.keys.length);
    return r;
  };
}

function isSubset(r: QueryResult, of: QueryResult): boolean {
  if (r.type === 'location' && of.type === 'location')
    return locationsEqual(intersectLocations(r.location, of.location), r.location);
  if (r.type === 'keys' && of.type === 'keys') return !subtractKeys(r.keys, of.keys).keys.length;
  if (r.type === 'nodes' && of.type === 'nodes') {
    const held = new Set(of.nodes.map((n) => printAddress(n.address)));
    return r.nodes.every((n) => held.has(printAddress(n.address)));
  }
  return false;
}

/**
 * The form-resolver table of the query forms over a model, for `selectAddress` (`SelectOptions.resolvers`): the
 * built-ins and the registry's (`extensions`, composed by `queryExtensions`). A caller may add its own kinds beside
 * them (`extra`), and nested forms see the whole table; an extra kind that is a query form's is a programmer error and
 * throws: a form's meaning is never replaced silently.
 */
export function queryResolvers(
  model: QueryModel,
  extra?: Readonly<Record<string, FormResolver>>,
  extensions: QueryExtensions = NoQueryExtensions,
): Readonly<Record<string, FormResolver>> {
  const table: Record<string, FormResolver> = { ...extra };
  const ext = { properties: extensions.properties, operators: extensions.operators };
  const kinds: [string, Resolve][] = [
    ...QueryFormKinds.map((k): [string, Resolve] => [k, Forms[k]]),
    ...[...extensions.forms.values()].map((e): [string, Resolve] => [e.name, extensionForm(e)]),
  ];
  for (const [kind, resolve] of kinds) {
    if (extra && Object.hasOwn(extra, kind)) throw new Error(`A resolver of kind ${kind} clashes with a query form`);
    table[kind] = (input): FormResolverResult => {
      const outer = active;
      const shared = outer ?? { visits: 0, noWrap: false };
      active = shared;
      const ctx: Ctx = {
        model,
        input,
        resolvers: table,
        issues: [],
        shared,
        ext,
        complete: true,
        pending: EmptyLocation,
      };
      try {
        checkNoArgs(input.form);
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
  /** The registry's query forms and properties (`queryExtensions`); none by default. */
  readonly extensions?: QueryExtensions;
}

/** Evaluates an address with the query forms over a model: `selectAddress` with `queryResolvers`, typed. */
export function evaluateQuery(root: Node, model: QueryModel, address: string, options: QueryOptions): Evaluation {
  const { extensions, ...select } = options;
  const resolvers = queryResolvers(model, undefined, extensions);
  const missing = missingExtension(address, resolvers, extensions ?? NoQueryExtensions);
  if (missing) return { issues: [missing], complete: true, pending: EmptyLocation };
  const sel = selectAddress(root, address, { ...select, resolvers });
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
