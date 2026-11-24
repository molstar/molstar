/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The core's pass over refs (_spec/README.md "Patches", Refs): what a ref names is the node at its id run, the
// leading ids of its address. Walked are the values a language types `ref` (not those with `in`), in params only,
// inside lists, records, maps, unions (by the member validation would choose), tables, data forms' args and slots.
// A stored patch's slot is walked by its own rule: its `at` and `parent`, the refs in the nodes it embeds by the target
// language, and a `set` value whose key the target language types as a ref in any kind that has it (by the language
// alone, never by reading the spec), each resolved over a dry run of the step's earlier ops. The ref index answers by
// address prefix behind `RefIndex`; this one is built on plain maps, kept in chunks by node, and `derive` walks only
// what changed.

import { type Fields, type KindSpec, type Language, type ParamType, type ValidateOptions } from '../schema/index.js';
import {
  type Address,
  type AddressPattern,
  isNode,
  isRecord,
  type Node,
  own,
  parseAddress,
  printAddress,
  resolveAddress,
  type Segment,
  type TreePath,
  type Value,
} from '../tree/index.js';
import { unionMember } from '../schema/validate.js';
import { patchFromNode, patchTarget } from './language.js';
import { type Draft, isUnder, Work, type WorkListener } from './work.js';

/** A ref in a tree: the address of the node that holds it, its step path from the root, and its text. */
export interface RefSite {
  readonly at?: Address;
  readonly path: readonly (string | number)[];
  readonly value: string;
}

export interface RefIndex {
  /** Refs whose id run (or an address nested in their forms) names `prefix` or below. */
  under(prefix: Address): readonly RefSite[];
  /** Refs that name no node. */
  dangling(): readonly RefSite[];
  /** The next revision's index, by node identity and position: what `after` shares with this one is reused. */
  derive(after: Node): RefIndex;
}

/** Where a stored patch's ref resolves: the spec with the step's ops before `op` applied (or through it, `after`). */
export interface Step {
  /** The stored patch's node. */
  readonly slot: Node;
  /** Where the stored patch is in the node that holds it. */
  readonly rel: TreePath;
  readonly op: number;
  readonly after: boolean;
  /** A stored patch inside another's embedded node or `set` value: that op's step, its dry run's starting state. */
  readonly outer?: Step;
}

/** A ref site with what the pass needs: its parsed runs and, in a stored patch, its step. */
export interface Site extends RefSite {
  /** The leading ids of its address. */
  readonly run: Address;
  /** Every id run in it: its own, then those of addresses nested in its forms. */
  readonly runs: readonly Address[];
  readonly step?: Step;
  /** The node that holds it. */
  readonly node: Node;
  /** The address of that node, or of its nearest ancestor with one. */
  readonly anchor: Address;
  /** The child positions from the anchor down to the node (empty when the node has an address). */
  readonly tail: readonly number[];
  /** Where in the node: `['params', key, ...]`. */
  readonly rel: TreePath;
}

/** The leading ids of an address, then those of the addresses nested in its forms' values. */
function runsOf(pattern: AddressPattern): Address[] {
  const run: string[] = [];
  for (const s of pattern.segments) {
    if (typeof s !== 'string') break;
    run.push(s);
  }
  const out: Address[] = [run];
  const nested = (v: Value | undefined) => {
    if (typeof v === 'string' && v.startsWith('/')) {
      const p = parseAddress(v).pattern;
      if (p) out.push(...runsOf(p));
    } else if (Array.isArray(v)) v.forEach(nested);
    else if (isNode(v)) for (const k of Object.keys(v.params ?? {})) nested(v.params![k]);
  };
  for (const s of pattern.segments)
    if (typeof s !== 'string') for (const k of Object.keys(s.params ?? {})) nested(s.params![k]);
  return out;
}

/**
 * A ref with `from` (an address prefix) moved to `to`: its id run if `from` is a prefix of it, and the addresses
 * nested in its forms the same way. The text unchanged when nothing matches, else printed canonically.
 */
export function rewriteRef(value: string, from: Address, to: Address): string {
  if (from.length && !value.includes(from[0])) return value;
  const pattern = parseAddress(value).pattern;
  if (!pattern) return value;
  let changed = false;
  const nested = (v: Value | undefined): Value | undefined => {
    if (typeof v === 'string' && v.startsWith('/')) {
      const w = rewriteRef(v, from, to);
      if (w !== v) changed = true;
      return w;
    }
    if (Array.isArray(v)) {
      const out = v.map(nested) as Value[];
      return out.some((x, i) => x !== v[i]) ? out : v;
    }
    if (isNode(v)) return form(v);
    return v;
  };
  const form = (f: Node): Node => {
    if (!f.params) return f;
    let params: Record<string, Value> | undefined;
    for (const k of Object.keys(f.params)) {
      const w = nested(f.params[k]);
      if (w !== f.params[k]) (params ??= { ...f.params })[k] = w as Value;
    }
    return params ? { ...f, params } : f;
  };
  let segments: Segment[] = [...pattern.segments];
  let run = 0;
  while (run < segments.length && typeof segments[run] === 'string') run++;
  if (run >= from.length && isUnder(segments.slice(0, run) as Address, from)) {
    segments = [...to, ...segments.slice(from.length)];
    changed = true;
  }
  segments = segments.map((s) => (typeof s === 'string' ? s : form(s)));
  return changed ? printAddress(segments, pattern.subtree) : value;
}

/** Whether values of a type may hold a ref or a slot (memoized per type object). */
const Holds = new WeakMap<object, boolean>();
function holds(t: ParamType): boolean {
  let h = Holds.get(t);
  if (h !== undefined) return h;
  switch (t.type) {
    case 'ref':
      h = t.in === undefined;
      break;
    case 'node':
      h = true;
      break;
    case 'list':
    case 'map':
      h = holds(t.of);
      break;
    case 'record':
      h = Object.values(t.fields).some(holds);
      break;
    case 'table':
      h = Object.values(t.columns).some(holds);
      break;
    case 'union':
      h = t.of.some(holds);
      break;
    default:
      h = false;
  }
  Holds.set(t, h);
  return h;
}

type Emit = (rel: TreePath, value: string, step?: Step) => void;

/**
 * The types a language gives a params key across its kinds, those that may hold a ref or a slot (memoized per language
 * object): a stored `set` value is typed by them, without reading the spec (_spec/README.md "Patches", a stored patch).
 */
const KeyTypes = new WeakMap<Language, Map<string, readonly ParamType[]>>();
function keyTypes(language: Language, key: string): readonly ParamType[] {
  let byKey = KeyTypes.get(language);
  if (!byKey) {
    byKey = new Map();
    for (const kind of Object.keys(language.kinds)) {
      const params = own(language.kinds, kind)?.params;
      if (!params) continue;
      for (const k of Object.keys(params)) {
        const t = own(params, k);
        if (!t || !holds(t)) continue;
        const list = byKey.get(k) as ParamType[] | undefined;
        if (!list) byKey.set(k, [t]);
        else if (!list.includes(t)) list.push(t);
      }
    }
    KeyTypes.set(language, byKey);
  }
  return byKey.get(key) ?? [];
}

/** A ref's runs, parsed once per text (most refs of a spec repeat a few texts); undefined when it is no address. */
const Parsed = new Map<string, readonly Address[] | null>();
export function parsedRuns(value: string): readonly Address[] | undefined {
  let runs = Parsed.get(value);
  if (runs === undefined) {
    const pattern = parseAddress(value).pattern;
    runs = pattern ? runsOf(pattern) : null;
    if (Parsed.size > 50_000) Parsed.clear();
    Parsed.set(value, runs);
  }
  return runs ?? undefined;
}

/** The pass over one root: its sites, and the dry runs of its stored patches. */
export class RefPass {
  private readonly resolve: (id: string) => Language | undefined;
  /** The dry runs, by stored patch node and the root it starts from. */
  private readonly dry = new Map<Node, Map<Node, Node[]>>();
  private rootNode?: Node;
  /** While walking an op's embedded node or `set` value: that op's step, the outer step of a stored patch inside. */
  private outer?: Step;

  constructor(
    root: Node | (() => Node),
    readonly language: Language,
    readonly options: ValidateOptions = {},
  ) {
    if (typeof root === 'function') this.rootOf = root;
    else this.rootNode = root;
    this.resolve = (id) => (id === language.id ? language : options.languages?.(id));
  }

  private rootOf?: () => Node;
  get root(): Node {
    return (this.rootNode ??= this.rootOf!());
  }

  /** The refs of a node of the tree: its params by its kind, its data form's args (`rel` from the node). */
  walkNode(n: Node, emit: Emit) {
    this.node(this.language, n, [], emit);
  }

  /** The refs of one param of a node of kind `kind`. */
  walkKey(kind: string, key: string, v: Value, emit: Emit) {
    const spec: KindSpec | undefined = own(this.language.kinds, kind);
    const t = spec?.params ? own(spec.params, key) : undefined;
    if (t && holds(t)) this.value(t, v, ['params', key], emit);
  }

  /** The refs of a subtree of the tree's nodes, each with its node's child positions below the subtree's root. */
  walkSubtree(n: Node, emit: (tail: readonly number[], rel: TreePath, value: string, step?: Step) => void) {
    const go = (x: Node, tail: number[]) => {
      this.walkNode(x, (rel, v, step) => emit(tail, rel, v, step));
      x.children?.forEach((c, i) => go(c, [...tail, i]));
    };
    go(n, []);
  }

  /** Whether a ref's run names a node: in the root, or over its stored patch's dry run. */
  names(run: Address, step?: Step): boolean {
    const root = step ? this.stateOf(step) : this.root;
    return !!root && resolveAddress(root, run) !== undefined;
  }

  /** The spec a stored patch's ref resolves over: its dry run at the op, started from its outer step's state. */
  private stateOf(step: Step): Node | undefined {
    const from = step.outer ? this.stateOf(step.outer) : this.root;
    const states = from && this.states(step.slot, from);
    return states && states[Math.min(step.op + (step.after ? 1 : 0), states.length - 1)];
  }

  /** The dry run of a stored patch over a root: the root after each op (a refused op changes nothing). */
  private states(slot: Node, root: Node): Node[] | undefined {
    let byRoot = this.dry.get(slot);
    let states = byRoot?.get(root);
    if (states) return states;
    const { patch } = patchFromNode(slot);
    if (!patch) return undefined;
    states = [root];
    const work = new Work(root);
    for (let i = 0; i < patch.ops.length; i++) {
      const failed = work.apply(patch.ops[i], i);
      states.push(failed ? states[states.length - 1] : work.finalize());
    }
    if (!byRoot) this.dry.set(slot, (byRoot = new Map()));
    byRoot.set(root, states);
    return states;
  }

  private node(lang: Language, n: Node, path: TreePath, emit: Emit) {
    const spec: KindSpec | undefined = own(lang.kinds, n.kind);
    if (!spec) return;
    if (spec.params && n.params) this.fields(spec.params, n.params, [...path, 'params'], emit);
    if (spec.data && n.args) {
      const specs = spec.data.args;
      n.args.forEach((v, i) => {
        if (specs.length) this.value(specs[Math.min(i, specs.length - 1)], v, [...path, 'args', i], emit);
      });
    }
  }

  private fields(specs: Fields, values: Readonly<Record<string, Value>>, path: TreePath, emit: Emit) {
    for (const key in values) {
      const spec = own(specs, key);
      const v = values[key];
      if (spec && v !== undefined && holds(spec)) this.value(spec, v, [...path, key], emit);
    }
  }

  private value(t: ParamType, v: Value, path: TreePath, emit: Emit) {
    if (!holds(t)) return;
    switch (t.type) {
      case 'ref':
        if (typeof v === 'string') emit(path, v);
        return;
      case 'list':
        if (Array.isArray(v)) v.forEach((x, i) => this.value(t.of, x, [...path, i], emit));
        return;
      case 'map':
        if (isRecord(v)) for (const k of Object.keys(v)) this.value(t.of, v[k], [...path, k], emit);
        return;
      case 'record':
        if (isRecord(v)) this.fields(t.fields, v, path, emit);
        return;
      case 'table':
        if (isRecord(v)) this.fields(t.columns, v, path, emit);
        return;
      case 'union': {
        const m = unionMember(t.of, v, { ...this.options, languages: this.resolve });
        if (m) this.value(m, v, path, emit);
        return;
      }
      case 'node': {
        const lang = isNode(v) ? this.resolve(t.language) : undefined;
        if (!lang || !isNode(v)) return;
        const patch = patchTarget(lang);
        if (patch) this.storedPatch(patch.target, v, path, emit);
        else this.slot(lang, v, path, emit);
        return;
      }
    }
  }

  /** A slot's body: its nodes' params and args, and its children. */
  private slot(lang: Language, n: Node, path: TreePath, emit: Emit) {
    this.node(lang, n, path, emit);
    n.children?.forEach((c, i) => this.slot(lang, c, [...path, 'children', i], emit));
  }

  /**
   * A stored patch, by the stored-patch rule: `at` and `parent` resolved over the spec with the earlier ops applied;
   * an embedded node's refs and a `set` value's through the op. A `set` value is typed by its key across the target's
   * kinds, never by the node the dry run finds, so whether it is a ref never changes under a patch. A stored patch
   * inside an op's node or value carries its own step, its dry run started from the op's state.
   */
  private storedPatch(target: Language | undefined, n: Node, path: TreePath, emit: Emit) {
    const { patch } = patchFromNode(n);
    if (!patch) return;
    const outer = this.outer;
    const stepOf = (op: number, after: boolean): Step => ({
      slot: n,
      rel: path,
      op,
      after,
      ...(outer ? { outer } : {}),
    });
    patch.ops.forEach((op, i) => {
      const at = [...path, 'children', i, 'params'];
      if ('at' in op) emit([...at, 'at'], printAddress(op.at), stepOf(i, false));
      if ('parent' in op) emit([...at, 'parent'], printAddress(op.parent), stepOf(i, false));
      if (!target) return;
      const through = stepOf(i, true);
      // a ref of the op's own content resolves through the op; a stored patch inside carries its own step
      const inner: Emit = (p, v, step) => emit(p, v, step ?? through);
      this.outer = through;
      try {
        if (op.op === 'insert' || op.op === 'replace') this.slot(target, op.node, [...at, 'node'], inner);
        if (op.op === 'set' && (op.where ?? 'params') === 'params') {
          const seen = new Set<string>();
          const once: Emit = (p, v, step) => {
            const k = p.join('\u0000');
            if (seen.has(k)) return;
            seen.add(k);
            inner(p, v, step);
          };
          for (const t of keyTypes(target, op.key)) this.value(t, op.value, [...at, 'value'], once);
        }
      } finally {
        this.outer = outer;
      }
    });
  }
}

/**
 * The sites of a subtree, kept by the node object that holds them: its own, then a chunk for each child that has any.
 * An unchanged node at an unchanged position is the same chunk in the next revision, so a derived index shares what
 * the ops left alone and builds the nodes on the paths they changed. A chunk is absent for a subtree without sites
 * (valid at any position).
 */
interface Chunk {
  readonly node: Node;
  readonly path: readonly (string | number)[];
  readonly at: Address | undefined;
  readonly anchor: Address;
  readonly tail: readonly number[];
  readonly own: readonly Site[];
  readonly kids: readonly (Chunk | undefined)[] | undefined;
}

/** The chunks a language and options have built, by node object, across revisions. */
interface Family {
  readonly chunks: WeakMap<Node, Chunk>;
  readonly empty: WeakSet<Node>;
}

const newFamily = (): Family => ({ chunks: new WeakMap(), empty: new WeakSet() });

const sameRun = (a: readonly (string | number)[] | undefined, b: readonly (string | number)[] | undefined) => {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
};
const samePos = (
  c: Chunk,
  path: readonly (string | number)[],
  at: Address | undefined,
  anchor: Address,
  tail: readonly number[],
) => sameRun(c.path, path) && sameRun(c.at, at) && sameRun(c.anchor, anchor) && sameRun(c.tail, tail);

/**
 * Builds the chunk tree of a root, reusing the family's chunks, and `prev` (the old root's) for the unchanged
 * children of a path copy.
 */
function buildChunks(pass: RefPass, family: Family, prev?: Chunk): Chunk | undefined {
  const go = (
    n: Node,
    path: readonly (string | number)[],
    at: Address | undefined,
    anchor: Address,
    tail: readonly number[],
    old: Chunk | undefined,
  ): Chunk | undefined => {
    const known = family.chunks.get(n);
    if (known && samePos(known, path, at, anchor, tail)) return known;
    if (family.empty.has(n)) return undefined;
    const own: Site[] = [];
    pass.walkNode(n, (rel, value, step) => {
      const runs = parsedRuns(value);
      if (!runs) return;
      const site: { -readonly [K in keyof Site]: Site[K] } = {
        path: [...path, ...rel],
        value,
        run: runs[0],
        runs,
        node: n,
        anchor,
        tail,
        rel,
      };
      if (at) site.at = at;
      if (step) site.step = step;
      own.push(site);
    });
    const children = n.children;
    let kids: (Chunk | undefined)[] | undefined;
    let any = false;
    if (children) {
      kids = new Array(children.length);
      // at the old node's position, its unchanged children are its old chunks, with no position to compare
      const samePlace = !!old && samePos(old, path, at, anchor, tail);
      const oldKids = old?.node.children;
      for (let i = 0; i < children.length; i++) {
        const c = children[i];
        let oc = old?.kids?.[i];
        if (samePlace && oldKids?.[i] === c) {
          kids[i] = oc;
        } else {
          const cAt = at && c.id !== undefined ? [...at, c.id] : undefined;
          if (oc && oc.node.id !== c.id) oc = undefined;
          kids[i] = go(c, [...path, 'children', i], cAt, cAt ?? anchor, cAt ? [] : [...tail, i], oc);
        }
        if (kids[i]) any = true;
      }
    }
    if (!own.length && !any) {
      family.empty.add(n);
      return undefined;
    }
    const chunk: Chunk = { node: n, path, at, anchor, tail, own, kids: any ? kids : undefined };
    family.chunks.set(n, chunk);
    return chunk;
  };
  return go(pass.root, [], [], [], [], prev);
}

function flatten(c: Chunk | undefined, out: Site[]) {
  if (!c) return;
  for (const s of c.own) out.push(s);
  if (c.kids) for (const k of c.kids) flatten(k, out);
}

class MapRefIndex implements RefIndex {
  private allSites?: Site[];
  private byFirst?: Map<string, Site[]>;
  private danglingSites?: Site[];

  constructor(
    readonly pass: RefPass,
    private readonly family: Family,
    private readonly chunk: Chunk | undefined,
    /** The options as the cache is keyed by them (the pass defaults them). */
    private readonly key: ValidateOptions | undefined,
  ) {}

  /** Every site, depth first. */
  get all(): readonly Site[] {
    if (!this.allSites) flatten(this.chunk, (this.allSites = []));
    return this.allSites;
  }

  /** The sites by the first id of each of their runs, built on the first query with a prefix. */
  private get first(): Map<string, Site[]> {
    if (!this.byFirst) {
      const by = (this.byFirst = new Map<string, Site[]>());
      for (const s of this.all) {
        for (const r of s.runs) {
          const f = r.length ? r[0] : '';
          let list = by.get(f);
          if (!list) by.set(f, (list = []));
          if (list[list.length - 1] !== s) list.push(s); // once a site, however many of its runs start so
        }
      }
    }
    return this.byFirst;
  }

  under(prefix: Address): readonly Site[] {
    if (!prefix.length) return this.all;
    return (this.first.get(prefix[0]) ?? []).filter((s) => s.runs.some((r) => isUnder(r, prefix)));
  }

  dangling(): readonly Site[] {
    return (this.danglingSites ??= this.all.filter((s) => !this.pass.names(s.run, s.step)));
  }

  /**
   * The next revision's index, from this one, by node identity and position: a node `after` shares with this index's
   * root, at the position it was, is this index's chunk, and the nodes on the paths that differ are walked again. It
   * may skip revisions. Cached on `after` like a built one.
   */
  derive(after: Node): RefIndex {
    const { language } = this.pass;
    const options = this.key;
    const had = entryOf(after, language, options);
    if (had instanceof MapRefIndex) return had;
    const pass = new RefPass(after, language, options);
    return store(
      after,
      language,
      options,
      new MapRefIndex(pass, this.family, buildChunks(pass, this.family, this.chunk), options),
    );
  }
}

/**
 * A derive not done yet, cached on the new root: the nearest index that exists (never another pending one, so a run of
 * applies leaves each pointing at the same one). It keeps that index alive: its root, its chunk tree and the caches it
 * made, one revision however long the run, until the next apply needs the new root's index and forces it.
 */
class Pending {
  constructor(readonly from: MapRefIndex) {}
}

type Entry = MapRefIndex | Pending;

/** The first ids of a ref's runs: what the overlay keys it by (`''` for a run naming the root). */
function firsts(runs: readonly Address[]): Set<string> {
  const out = new Set<string>();
  for (const r of runs) out.add(r.length ? r[0] : '');
  return out;
}

const Indices = new WeakMap<Node, WeakMap<Language, Map<ValidateOptions | undefined, Entry>>>();

/** The chunks of a language, by the options object callers pass (held weakly), and the ones for none. */
const Families = new WeakMap<Language, { none?: Family; by: WeakMap<object, Family> }>();

function familyOf(language: Language, options: ValidateOptions | undefined): Family {
  let slot = Families.get(language);
  if (!slot) Families.set(language, (slot = { by: new WeakMap() }));
  if (!options) return (slot.none ??= newFamily());
  let family = slot.by.get(options);
  if (!family) slot.by.set(options, (family = newFamily()));
  return family;
}

function entryOf(root: Node, language: Language, options: ValidateOptions | undefined): Entry | undefined {
  return Indices.get(root)?.get(language)?.get(options);
}

function put(root: Node, language: Language, options: ValidateOptions | undefined, entry: Entry) {
  let byLanguage = Indices.get(root);
  if (!byLanguage) Indices.set(root, (byLanguage = new WeakMap()));
  let byOptions = byLanguage.get(language);
  if (!byOptions) byLanguage.set(language, (byOptions = new Map()));
  byOptions.set(options, entry);
}

function store(root: Node, language: Language, options: ValidateOptions | undefined, index: MapRefIndex) {
  put(root, language, options, index);
  return index;
}

/**
 * Leaves the new root's index to be derived from the old root's when it is next read, if the old root has one (built,
 * derived or pending): nothing is walked now. The cache hits only for the same `options` object the caller passes.
 */
export function deferDerive(root: Node, after: Node, language: Language, options?: ValidateOptions) {
  const old = entryOf(root, language, options);
  if (!old || entryOf(after, language, options)) return;
  put(after, language, options, new Pending(old instanceof Pending ? old.from : old));
}

/** Whether the root has an index, derived later or not. */
export function hasIndex(root: Node, language: Language, options?: ValidateOptions): boolean {
  return !!entryOf(root, language, options);
}

/** The ref index of a root if it has one (a pending derive is done now). */
export function peekIndex(root: Node, language: Language, options?: ValidateOptions): RefIndex | undefined {
  return hasIndex(root, language, options) ? mapIndexOf(root, language, options) : undefined;
}

/** The ref index of a revision, cached per root, language and options. */
export function refIndexOf(root: Node, language: Language, options?: ValidateOptions): RefIndex {
  return mapIndexOf(root, language, options);
}

/** An index built from scratch, sharing nothing and cached nowhere: what a derived one must equal. */
export function freshRefIndex(root: Node, language: Language, options?: ValidateOptions): RefIndex {
  const pass = new RefPass(root, language, options);
  const family = newFamily();
  return new MapRefIndex(pass, family, buildChunks(pass, family), options);
}

function mapIndexOf(root: Node, language: Language, options?: ValidateOptions): MapRefIndex {
  const had = entryOf(root, language, options);
  if (had instanceof MapRefIndex) return had;
  if (had) return had.from.derive(root) as MapRefIndex;
  const family = familyOf(language, options);
  const pass = new RefPass(root, language, options);
  return store(root, language, options, new MapRefIndex(pass, family, buildChunks(pass, family), options));
}

/** A ref of the draft: the draft that holds it, where in it, and what it says now. */
export interface LiveSite {
  readonly holder: Draft;
  readonly rel: TreePath;
  readonly value: string;
  readonly run: Address;
  readonly runs: readonly Address[];
  readonly step?: Step;
}

/**
 * The ref index of an apply, derived across its ops (_spec/README.md "Patches", Refs): the old root's index (built
 * once, cached per revision) for what the ops left alone, and an overlay for each key an op wrote or a rewrite changed
 * and each node an op made, kept by the draft that holds it. A base site counts while its node is where the moves so
 * far put it, still the old root's, and its key untouched.
 */
export class LiveRefs implements WorkListener {
  private baseIndex?: MapRefIndex;
  private readonly superseded = new WeakMap<Node, Set<string>>();
  private readonly overlay = new Map<Draft, Map<string, LiveSite[]>>();
  private readonly byFirst = new Map<string, Set<LiveSite>>();
  private overlayPass?: RefPass;

  constructor(
    readonly work: Work,
    readonly language: Language,
    readonly options: ValidateOptions | undefined,
  ) {}

  /** The old root's pass: names in the old root. */
  get oldPass(): RefPass {
    // names alone need no index: an apply that changes no address builds none
    return this.baseIndex?.pass ?? (this.namesPass ??= new RefPass(this.work.old, this.language, this.options));
  }

  private namesPass?: RefPass;

  /** A ref's id run, or undefined when it is no address. */
  runOf(value: string): Address | undefined {
    return parsedRuns(value)?.[0];
  }

  private get base(): MapRefIndex {
    return (this.baseIndex ??= mapIndexOf(this.work.old, this.language, this.options));
  }

  /** A pass for walking new values: typing reads the language alone, never the spec. */
  private get pass(): RefPass {
    return (this.overlayPass ??= new RefPass(() => this.work.finalize(), this.language, this.options));
  }

  key(d: Draft, key: string) {
    if (d.origin && !d.made)
      (this.superseded.get(d.origin) ?? this.superseded.set(d.origin, new Set()).get(d.origin)!).add(key);
    this.drop(d, key);
    const v = d.params ? own(d.params, key) : undefined;
    if (v === undefined) return;
    this.pass.walkKey(d.kind, key, v, (rel, value, step) => this.add(d, rel, value, step));
  }

  made(d: Draft) {
    this.pass.walkSubtree(d.base, (tail, rel, value, step) => {
      const holder = this.work.below(d, tail);
      if (holder) this.add(holder, rel, value, step);
    });
  }

  private add(holder: Draft, rel: TreePath, value: string, step?: Step) {
    const runs = parsedRuns(value);
    if (!runs) return;
    const site: LiveSite = { holder, rel, value, run: runs[0], runs, ...(step ? { step } : {}) };
    const key = String(rel[1]);
    let byKey = this.overlay.get(holder);
    if (!byKey) this.overlay.set(holder, (byKey = new Map()));
    let list = byKey.get(key);
    if (!list) byKey.set(key, (list = []));
    list.push(site);
    for (const f of firsts(runs)) {
      let set = this.byFirst.get(f);
      if (!set) this.byFirst.set(f, (set = new Set()));
      set.add(site);
    }
  }

  private drop(holder: Draft, key: string) {
    const list = this.overlay.get(holder)?.get(key);
    if (!list) return;
    for (const s of list) for (const f of firsts(s.runs)) this.byFirst.get(f)?.delete(s);
    this.overlay.get(holder)!.delete(key);
  }

  /** The draft a base site's node is now, or undefined (removed, made over, or its key changed). */
  private locate(s: Site): Draft | undefined {
    const key = s.rel[1];
    if (typeof key === 'string' && this.superseded.get(s.node)?.has(key)) return undefined;
    let anchor = s.anchor;
    for (const m of this.work.moves) if (isUnder(anchor, m.from)) anchor = [...m.to, ...anchor.slice(m.from.length)];
    const at = this.work.find(anchor);
    const d = at && this.work.below(at, s.tail);
    return d && !d.made && d.origin === s.node && this.work.attached(d) ? d : undefined;
  }

  /** The live sites whose id run (or an address nested in their forms) names `prefix` or below. */
  under(prefix: Address, withBase = true, withOverlay = true): LiveSite[] {
    const out: LiveSite[] = [];
    if (withBase) {
      for (const s of this.base.under(prefix)) {
        const holder = this.locate(s);
        if (holder)
          out.push({
            holder,
            rel: s.rel,
            value: s.value,
            run: s.run,
            runs: s.runs,
            ...(s.step ? { step: s.step } : {}),
          });
      }
    }
    if (!withOverlay) return out;
    const candidates = prefix.length
      ? this.byFirst.get(prefix[0])
      : new Set([...this.byFirst.values()].flatMap((x) => [...x]));
    for (const s of candidates ?? []) {
      if (s.runs.some((r) => isUnder(r, prefix)) && this.work.attached(s.holder)) out.push(s);
    }
    return out;
  }

  /** Whether a ref of the old root names `prefix` or below (by its id run). */
  namedInOld(prefix: Address): boolean {
    return this.base.under(prefix).some((s) => isUnder(s.run, prefix));
  }

  /** Every live site of the overlay: what the ops wrote, made or rewrote. */
  written(): LiveSite[] {
    const out: LiveSite[] = [];
    for (const [holder, byKey] of this.overlay)
      if (this.work.attached(holder)) for (const list of byKey.values()) out.push(...list);
    return out;
  }

  /** `from` moved to `to`: every live ref under it rewritten, its key's sites derived again. */
  move(from: Address, to: Address) {
    const keys = new Map<Draft, Set<string>>();
    for (const s of this.under(from)) {
      const next = rewriteRef(s.value, from, to);
      if (next === s.value) continue;
      this.work.rewrite(s.holder, s.rel, next);
      let set = keys.get(s.holder);
      if (!set) keys.set(s.holder, (set = new Set()));
      set.add(String(s.rel[1]));
    }
    for (const [d, set] of keys) for (const k of set) this.key(d, k);
  }
}
