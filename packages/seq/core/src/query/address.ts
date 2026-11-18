/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Addressing (_spec/README.md "Addressing"): an address resolved to what it names (spec nodes, nodes data makes,
// elements, a Location) through the query forms, its user-defined refs followed first; and back, from what is named to
// its shortest stable address (a ref if one names exactly it, else ids and forms). Bounded by a limit and a depth (more
// matches is an issue, never a truncation) and a ref chain's depth; an address that no longer names anything is
// reported with what it named, never re-bound.

import {
  hasErrors,
  type Issue,
  issueOf,
  isNode,
  type Node,
  parseAddress,
  printAddress,
  type Segment,
  type Selected,
  type Value,
} from '@molstar/core/repr/tree';
import { DefaultRefDepth, type FollowedRef, followRef, RefsId, refPrefix } from '@molstar/seqx-core/lang';
import { AddressIssues } from './address-issues.js';
import { type ElementKey, elements, elementsEqual } from './elements.js';
import { evaluateQuery } from './forms.js';
import {
  EmptyLocation,
  isEmptyLocation,
  type Location,
  type LocationPart,
  type LocationRect,
  location,
  locationsEqual,
  unionLocations,
} from './location.js';
import type { QueryModel } from './model.js';
import { DefaultNumbering, findNumbering } from './numbering.js';
import type { QueryResult } from './result.js';

/** The bounds of a resolution (_spec/README.md "Addressing", Configuration). */
export interface AddressBounds {
  /** At most this many nodes or elements: more is `repr.tree.address.limit_exceeded`, never a truncation. */
  readonly limit: number;
  /** How deep below the root a walk may go. */
  readonly maxDepth: number;
  /** How many refs a chain may pass through. */
  readonly refDepth: number;
}

export const DefaultAddressBounds: AddressBounds = { limit: 10_000, maxDepth: 64, refDepth: DefaultRefDepth };

export interface ResolveAddressOptions extends Partial<AddressBounds> {
  readonly signal?: AbortSignal;
  /** A `performance.now()` time the resolution stops at. */
  readonly deadline?: number;
}

/**
 * What became of an address: `resolved` (its result, possibly incomplete over data read by region), `waiting` (it
 * names a node data makes that is not held yet, or data not loaded: it resolves when they arrive), `missing` (it names
 * nothing: reported with what it named), `failed` (an error: its issues).
 */
export type AddressStatus = 'resolved' | 'waiting' | 'missing' | 'failed';

export interface AddressResolution {
  /** The address as asked. */
  readonly address: string;
  readonly status: AddressStatus;
  /** What it names; absent unless resolved. */
  readonly result?: QueryResult;
  /** The refs followed, in order. */
  readonly via: readonly string[];
  readonly issues: readonly Issue[];
  /** False when the result covers only what is loaded; `pending` is then what is not. */
  readonly complete: boolean;
  readonly pending: Location;
}

const refusalCodes = {
  cycle: 'seq.address.ref.cycle',
  depth: 'seq.address.ref.depth',
  set: 'seq.address.ref.set',
  empty: 'seq.address.ref.empty',
} as const;

const refusalMessages = {
  cycle: 'A ref leads back to itself',
  depth: 'A chain of refs is longer than the bound',
  set: 'A ref to a set is followed by more segments',
  empty: 'A ref names no address',
} as const;

function refused(f: Extract<FollowedRef, { refused: unknown }>, address: string): Issue {
  const code = refusalCodes[f.refused];
  return {
    ...issueOf(AddressIssues, code, refusalMessages[f.refused], {
      address,
      ...(AddressIssues.codes[code].hint ? { hint: AddressIssues.codes[code].hint } : {}),
    }),
    places: f.via.map((a) => ({ address: a })),
  };
}

/** A failure inside the rewrite of nested refs, caught at its edge. */
class Refused {
  constructor(readonly issue: Issue) {}
}

/**
 * The keys whose values are operands (a nested form, or an address from the root): only there is a string an
 * address. Elsewhere a string is a value (an id, a label), which a quoted `"/refs/x"` may well be: never rewritten.
 */
const OperandKeys: ReadonlySet<string> = new Set(['of', 'a', 'b', 'overlapping']);

/**
 * A form value with every nested operand address that starts with a ref rewritten to what the ref names: one address,
 * or a set as `(or :of [...])`, so an operand naming a ref reads as its target (the query forms refuse a ref node
 * itself). `operand`: whether the value is an operand's (a string there is an address).
 */
function followNested(root: Node, value: Value, refDepth: number, via: string[], operand: boolean): Value {
  if (isNode(value)) return followForm(root, value, refDepth, via);
  if (Array.isArray(value)) return value.map((v) => followNested(root, v, refDepth, via, operand));
  if (!operand || typeof value !== 'string' || !value.startsWith(`/${RefsId}/`) || !refPrefix(root, value))
    return value;
  const f = followRef(root, value, refDepth);
  if ('refused' in f) throw new Refused(refused(f, value));
  for (const a of f.via) if (!via.includes(a)) via.push(a);
  const out = f.addresses.map((a) => followAddress(root, a, refDepth, via));
  return out.length === 1 ? out[0] : { kind: 'or', params: { of: out } };
}

function followForm(root: Node, form: Node, refDepth: number, via: string[]): Node {
  const params = form.params;
  if (!params) return form;
  let changed = false;
  const next: Record<string, Value> = {};
  for (const key of Object.keys(params)) {
    const v = followNested(root, params[key], refDepth, via, OperandKeys.has(key));
    if (v !== params[key]) changed = true;
    next[key] = v;
  }
  return changed ? { ...form, params: next } : form;
}

/** An address with the refs nested in its forms followed (its own leading ref is followed before). */
function followAddress(root: Node, address: string, refDepth: number, via: string[]): string {
  const p = parseAddress(address).pattern;
  if (!p || p.segments.every((s) => typeof s === 'string')) return address;
  const segments: Segment[] = p.segments.map((s) => (typeof s === 'string' ? s : followForm(root, s, refDepth, via)));
  return printAddress(segments, p.subtree);
}

const nodeKey = (s: Selected) => printAddress(s.address);

/**
 * Results of one type, together; undefined when the types differ. Linear in the members: their parts, keys or nodes
 * gathered, then normalized once (a set of thousands of addresses never merges pair by pair).
 */
function unionResults(rs: readonly QueryResult[], limit: number): QueryResult | Issue | undefined {
  if (!rs.length) return undefined;
  if (rs.length === 1) return rs[0];
  const type = rs[0].type;
  if (rs.some((r) => r.type !== type)) return undefined;
  if (type === 'location') {
    const parts: LocationPart[] = [];
    const rects: LocationRect[] = [];
    for (const r of rs as Extract<QueryResult, { type: 'location' }>[]) {
      for (const p of r.location.parts) parts.push(p);
      for (const x of r.location.rects) rects.push(x);
    }
    return { type, location: location({ parts, rects }) };
  }
  if (type === 'elements') {
    const keys: ElementKey[] = [];
    for (const r of rs as Extract<QueryResult, { type: 'elements' }>[]) for (const k of r.elements.keys) keys.push(k);
    const e = elements(keys);
    return e.keys.length > limit ? limitIssue(limit) : { type, elements: e };
  }
  const seen = new Set<string>();
  const nodes: Selected[] = [];
  for (const r of rs as Extract<QueryResult, { type: 'nodes' }>[])
    for (const n of r.nodes) {
      const k = nodeKey(n);
      if (seen.has(k)) continue;
      seen.add(k);
      nodes.push(n);
    }
  return nodes.length > limit ? limitIssue(limit) : { type: 'nodes', nodes };
}

function limitIssue(limit: number): Issue {
  return {
    code: 'repr.tree.address.limit_exceeded',
    severity: 'error',
    message: `More than ${limit} matches`,
    hint: 'narrow the address, or raise the limit',
  };
}

function isEmptyResult(r: QueryResult): boolean {
  return r.type === 'nodes'
    ? !r.nodes.length
    : r.type === 'elements'
      ? !r.elements.keys.length
      : isEmptyLocation(r.location);
}

/**
 * What an address names (_spec/README.md "Addressing", Resolving): its leading ref followed (a chain, a set member by
 * member), the refs nested in its forms followed, then the query forms over `model`. A set's members must give one
 * type of result, together. An empty result is `waiting` when a node data makes is not held or data is not loaded,
 * else `missing`, with an issue naming the address and, through `places`, the refs it was reached by.
 */
export function resolveAddressIn(
  root: Node,
  model: QueryModel,
  address: string,
  options: ResolveAddressOptions = {},
): AddressResolution {
  const limit = options.limit ?? DefaultAddressBounds.limit;
  const maxDepth = options.maxDepth ?? DefaultAddressBounds.maxDepth;
  const refDepth = options.refDepth ?? DefaultAddressBounds.refDepth;
  const done = (
    status: AddressStatus,
    issues: readonly Issue[],
    via: readonly string[],
    extra: Partial<AddressResolution> = {},
  ): AddressResolution => ({ address, status, via, issues, complete: true, pending: EmptyLocation, ...extra });

  const parsed = parseAddress(address);
  if (!parsed.pattern) return done('failed', parsed.issues, []);
  const followed = followRef(root, address, refDepth);
  if ('refused' in followed) return done('failed', [refused(followed, address)], followed.via);
  const via = [...followed.via];
  let targets: string[];
  try {
    targets = followed.addresses.map((a) => followAddress(root, a, refDepth, via));
  } catch (e) {
    if (e instanceof Refused) return done('failed', [e.issue], via);
    throw e;
  }

  const issues: Issue[] = [];
  const results: QueryResult[] = [];
  let complete = true;
  let pending = EmptyLocation;
  let notHeld = false;
  for (const t of targets) {
    const ev = evaluateQuery(root, model, t, { limit, maxDepth, signal: options.signal, deadline: options.deadline });
    if (hasErrors(ev.issues) || !ev.result) return done('failed', ev.issues, via);
    for (const i of ev.issues) {
      // a node data makes and not held: the address waits for it (reported once, below)
      if (i.code === 'seq.query.node.missing') notHeld = true;
      else issues.push(i);
    }
    if (!ev.complete) {
      complete = false;
      pending = unionLocations(pending, ev.pending);
    }
    results.push(ev.result);
  }
  const result = unionResults(results, limit);
  if (!result)
    return done(
      'failed',
      [issueOf(AddressIssues, 'seq.address.type', 'A set names results of different types', { address })],
      via,
    );
  if (!('type' in result)) return done('failed', [result], via);
  const places = [...via, ...targets.filter((t) => t !== address)].map((a) => ({ address: a }));
  const extra = { complete, pending };
  if (isEmptyResult(result)) {
    if (notHeld || !complete) {
      const waiting = issueOf(AddressIssues, 'seq.address.waiting', 'The address waits for its data', { address });
      return done('waiting', [...issues, places.length ? { ...waiting, places } : waiting], via, extra);
    }
    const missing = issueOf(AddressIssues, 'seq.address.missing', 'The address names nothing', {
      address,
      hint: AddressIssues.codes['seq.address.missing'].hint,
    });
    return done('missing', [...issues, places.length ? { ...missing, places } : missing], via, extra);
  }
  return done('resolved', issues, via, { ...extra, result });
}

// ---- what is named, as positions ----

export interface Located {
  readonly location: Location;
  /** False when some of it is not loaded yet. */
  readonly complete: boolean;
}

/** The section a node is in, by its address (a section, or a track or layer under one). */
function sectionOf(root: Node, address: readonly Segment[]): Node | undefined {
  let node: Node | undefined = root;
  for (const s of address) {
    if (typeof s !== 'string') return undefined;
    node = node.children?.find((c) => c.id === s);
    if (!node) return undefined;
    if (node.kind === 'section') return node;
  }
  return undefined;
}

/**
 * The positions a result covers (what focus and reveal move to): a Location itself; elements by their items'
 * extents; a system's node by the whole system; an annotation's by all its items; a section's (and what is under it)
 * by the region it shows. Other nodes cover nothing.
 */
export function locationOfResult(root: Node, model: QueryModel, result: QueryResult): Located {
  if (result.type === 'location') return { location: result.location, complete: true };
  let complete = true;
  const parts: { system: string; ranges: number[] }[] = [];
  const itemsOf = (owner: string, keys?: ReadonlySet<string | number>) => {
    const a = model.annotation(owner);
    if (!a) return;
    const read = a.items();
    if (read.pending && !isEmptyLocation(read.pending)) complete = false;
    const ranges: number[] = [];
    for (const it of read.items) if (!keys || keys.has(it.key)) ranges.push(it.start, it.end);
    if (ranges.length) parts.push({ system: a.target, ranges });
  };
  const whole = (system: string, range?: Value) => {
    const s = model.system(system);
    if (!s) return;
    if (Array.isArray(range) && range.length === 2) {
      const [from, to] = range as number[];
      parts.push({ system, ranges: from <= to ? [from, to] : [from, s.length, 0, to] });
    } else parts.push({ system, ranges: [0, s.length] });
  };
  if (result.type === 'elements') {
    const byOwner = new Map<string, Set<string | number>>();
    for (const k of result.elements.keys) {
      let keys = byOwner.get(k.owner);
      if (!keys) byOwner.set(k.owner, (keys = new Set()));
      keys.add(k.key);
    }
    for (const [owner, keys] of byOwner) itemsOf(owner, keys);
  } else {
    for (const n of result.nodes) {
      const a = printAddress(n.address);
      if (model.system(a)) whole(a);
      else if (model.annotation(a)) itemsOf(a);
      else {
        const section = sectionOf(root, n.address);
        const of = section?.params?.of;
        if (typeof of === 'string') whole(of, section?.params?.range);
      }
    }
  }
  return { location: location({ parts }), complete };
}

// ---- and back: the address of what is named ----

/** What an address is printed from: the preference order of its forms (_spec/README.md "Addressing", Back). */
export type AddressForm = 'ref' | 'ids';

export interface AddressForOptions extends ResolveAddressOptions {
  /** The numbering labels are printed in (`ordinal` by default); a label that does not read back falls back to it. */
  readonly numbering?: string;
  /** Which forms are tried, in order: a ref naming exactly it, then ids and forms (`['ref', 'ids']`). */
  readonly prefer?: readonly AddressForm[];
}

export const DefaultAddressPreference: readonly AddressForm[] = ['ref', 'ids'];

function resultsEqual(a: QueryResult, b: QueryResult): boolean {
  if (a.type === 'location' && b.type === 'location') return locationsEqual(a.location, b.location);
  if (a.type === 'elements' && b.type === 'elements') return elementsEqual(a.elements, b.elements);
  if (a.type === 'nodes' && b.type === 'nodes') {
    if (a.nodes.length !== b.nodes.length) return false;
    const keys = new Set(a.nodes.map(nodeKey));
    return b.nodes.every((n) => keys.has(nodeKey(n)));
  }
  return false;
}

/** What each ref names, per spec root and model (both immutable per revision), by its address. */
const RefResults = new WeakMap<Node, WeakMap<QueryModel, Map<string, QueryResult | undefined>>>();

/** The refs of a spec, each with what it names now (resolved once per root and model). */
export function refResults(
  root: Node,
  model: QueryModel,
  options: ResolveAddressOptions = {},
): ReadonlyMap<string, QueryResult | undefined> {
  let byModel = RefResults.get(root);
  if (!byModel) RefResults.set(root, (byModel = new WeakMap()));
  let refs = byModel.get(model);
  if (refs) return refs;
  refs = new Map();
  const container = root.children?.find((c) => c.id === RefsId && c.kind === 'refs');
  for (const r of container?.children ?? []) {
    if (r.kind !== 'ref' || r.id === undefined) continue;
    const a = printAddress([RefsId, r.id]);
    refs.set(a, resolveAddressIn(root, model, a, options).result);
  }
  // a model that changes with the data (an overlay's) is a new object each version: its entries go with it
  byModel.set(model, refs);
  return refs;
}

/** A label as a form value: an integer label is a number, any other a string (quoted when it must be). */
function labelValue(label: string): Value {
  return /^(0|-?[1-9]\d*)$/.test(label) ? Number(label) : label;
}

/**
 * The forms that name a part's ranges on its system, labels in `numbering`; undefined when a label is missing there.
 * A cyclic system's range over the origin (`[s, length)` then `[0, e)`) is one wrapping range unless `wrap` is false.
 */
function partForms(model: QueryModel, p: LocationPart, numbering: string, wrap: boolean): Node[] | undefined {
  const system = model.system(p.system);
  if (!system) return undefined;
  const n = findNumbering(system, numbering);
  if (!n) return undefined;
  const label = (i: number) => n.numbering.labelAt(i - n.shift);
  const r = Array.from(p.ranges);
  if (r.length === 2 && r[0] === 0 && r[1] === system.length) return [{ kind: 'all' }];
  const pairs: [number, number][] = [];
  for (let i = 0; i + 1 < r.length; i += 2) pairs.push([r[i], r[i + 1]]);
  // over the origin: the last range runs to the end and the first starts at 0
  if (wrap && system.cyclic && pairs.length > 1 && pairs[0][0] === 0 && pairs[pairs.length - 1][1] === system.length) {
    const first = pairs.shift()!;
    pairs[pairs.length - 1] = [pairs[pairs.length - 1][0], first[1] + system.length];
  }
  const forms: Node[] = [];
  const params = (extra: Record<string, Value>) =>
    numbering === DefaultNumbering ? extra : { ...extra, numbering: labelValue(numbering) };
  for (const [s, e] of pairs) {
    const last = e - 1 >= system.length ? e - 1 - system.length : e - 1;
    const from = label(s);
    const to = label(last);
    if (from === undefined || to === undefined) return undefined;
    forms.push(
      e - s === 1
        ? { kind: 'at', params: params({ label: labelValue(from) }) }
        : { kind: 'range', params: params({ from: labelValue(from), to: labelValue(to) }) },
    );
  }
  return forms;
}

function segmentsOf(address: string): readonly Segment[] {
  return parseAddress(address).pattern?.segments ?? [];
}

/** One system's ranges as an address: the system's, then its forms (several under `or`). */
function partAddress(model: QueryModel, p: LocationPart, numbering: string, wrap: boolean): string | undefined {
  const forms = partForms(model, p, numbering, wrap);
  if (!forms?.length) return undefined;
  const form = forms.length === 1 ? forms[0] : { kind: 'or', params: { of: forms } };
  return printAddress([...segmentsOf(p.system), form]);
}

/** Ids and forms: one item's address, or several under `or` from the root. */
function joined(items: readonly Value[]): string | undefined {
  if (!items.length) return undefined;
  if (items.length === 1) return typeof items[0] === 'string' ? items[0] : printAddress([items[0] as Node]);
  return printAddress([{ kind: 'or', params: { of: items as Value[] } }]);
}

function idsAddress(model: QueryModel, target: QueryResult, numbering: string): string | undefined {
  if (target.type === 'nodes') return joined(target.nodes.map((n) => printAddress(n.address)));
  if (target.type === 'elements') {
    const byOwner = new Map<string, (string | number)[]>();
    for (const k of target.elements.keys) {
      let keys = byOwner.get(k.owner);
      if (!keys) byOwner.set(k.owner, (keys = []));
      keys.push(k.key);
    }
    const items: string[] = [];
    for (const [owner, keys] of byOwner)
      items.push(
        printAddress([...segmentsOf(owner), { kind: 'item', params: { id: keys.length === 1 ? keys[0] : keys } }]),
      );
    return joined(items);
  }
  const items: Value[] = [];
  for (const p of target.location.parts) {
    const a = partAddress(model, p, numbering, true);
    if (!a) return undefined;
    items.push(a);
  }
  for (const r of target.location.rects) {
    const a = partAddress(model, r.a, numbering, false);
    const b = partAddress(model, r.b, numbering, false);
    if (!a || !b) return undefined;
    items.push({ kind: 'rect', params: { a, b } });
  }
  return joined(items);
}

/**
 * The shortest stable address of what is named (_spec/README.md "Addressing", Back): by `prefer`, a ref of the spec
 * that names exactly it (`/refs/dbd`), else its ids and forms (labels in `numbering`; one that does not read back to
 * the same result falls back to `ordinal`, which always does). Undefined when nothing names it (an empty result).
 */
export function addressFor(
  root: Node,
  model: QueryModel,
  target: QueryResult,
  options: AddressForOptions = {},
): string | undefined {
  if (isEmptyResult(target)) return undefined;
  for (const form of options.prefer ?? DefaultAddressPreference) {
    if (form === 'ref') {
      for (const [a, r] of refResults(root, model, options)) if (r && resultsEqual(r, target)) return a;
      continue;
    }
    const numbering = options.numbering ?? DefaultNumbering;
    for (const n of numbering === DefaultNumbering ? [numbering] : [numbering, DefaultNumbering]) {
      const a = idsAddress(model, target, n);
      if (!a) continue;
      // a label repeated or missing in that numbering does not read back: never print an address naming another thing
      const back = resolveAddressIn(root, model, a, options).result;
      if (back && resultsEqual(back, target)) return a;
    }
  }
  return undefined;
}
