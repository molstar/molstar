/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The core's pass over refs (_spec/README.md "Patches", Refs): what a ref names is the node at its id run, the
// leading ids of its address. Walked are the values a language types `ref` (not those with `in`), in params only,
// inside lists, records, maps, unions (by the member validation would choose), tables, data forms' args and slots.
// A stored patch's slot is walked by its own rule: its `at` and `parent`, the refs in the nodes it embeds by the target
// language, and a `set` value typed through its target's kind and key, each resolved over a dry run of the step's
// earlier ops. The ref index answers by address prefix behind `RefIndex`; this one is built on plain maps, and
// `derive` rebuilds it.

import {
  type Changes,
  type Fields,
  type KindSpec,
  type Language,
  type ParamType,
  type ValidateOptions,
} from '../schema/index.js';
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
  /** The next revision's index; a rebuild here. */
  derive(after: Node, changes: Changes): RefIndex;
}

/** Where a stored patch's ref resolves: the spec with the step's ops before `op` applied (or through it, `after`). */
export interface Step {
  /** The stored patch's node. */
  readonly slot: Node;
  /** Where the stored patch is in the node that holds it. */
  readonly rel: TreePath;
  readonly op: number;
  readonly after: boolean;
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

/** A ref's runs, parsed once per text (most refs of a spec repeat a few texts); undefined when it is no address. */
const Parsed = new Map<string, readonly Address[] | null>();
function parsedRuns(value: string): readonly Address[] | undefined {
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
  private readonly dry = new Map<Node, Node[]>();
  private rootNode?: Node;

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

  /** Every site of the tree, depth first, each with the node that holds it and where. */
  sites(): Site[] {
    const out: Site[] = [];
    const go = (n: Node, path: TreePath, at: Address | undefined, anchor: Address, tail: number[]) => {
      this.walkNode(n, (rel, value, step) => {
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
        out.push(site);
      });
      const kids = n.children;
      if (!kids) return;
      for (let i = 0; i < kids.length; i++) {
        const c = kids[i];
        const cAt = at && c.id !== undefined ? [...at, c.id] : undefined;
        go(c, [...path, 'children', i], cAt, cAt ?? anchor, cAt ? [] : [...tail, i]);
      }
    };
    go(this.root, [], [], [], []);
    return out;
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
    let root: Node | undefined = this.root;
    if (step) {
      const states = this.states(step.slot);
      root = states && states[Math.min(step.op + (step.after ? 1 : 0), states.length - 1)];
    }
    return !!root && resolveAddress(root, run) !== undefined;
  }

  /** The dry run of a stored patch over the root: the root after each op (a refused op changes nothing). */
  private states(slot: Node): Node[] | undefined {
    let states = this.dry.get(slot);
    if (states) return states;
    const { patch } = patchFromNode(slot);
    if (!patch) return undefined;
    states = [this.root];
    const work = new Work(this.root);
    for (let i = 0; i < patch.ops.length; i++) {
      const failed = work.apply(patch.ops[i], i);
      states.push(failed ? states[states.length - 1] : work.finalize());
    }
    this.dry.set(slot, states);
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
   * an embedded node's refs and a `set` value's (typed through the target's kind and key) through the op.
   */
  private storedPatch(target: Language | undefined, n: Node, path: TreePath, emit: Emit) {
    const { patch } = patchFromNode(n);
    if (!patch) return;
    patch.ops.forEach((op, i) => {
      const at = [...path, 'children', i, 'params'];
      if ('at' in op) emit([...at, 'at'], printAddress(op.at), { slot: n, rel: path, op: i, after: false });
      if ('parent' in op) emit([...at, 'parent'], printAddress(op.parent), { slot: n, rel: path, op: i, after: false });
      if (!target) return;
      const inner: Emit = (p, v) => emit(p, v, { slot: n, rel: path, op: i, after: true });
      if (op.op === 'insert' || op.op === 'replace') this.slot(target, op.node, [...at, 'node'], inner);
      if (op.op === 'set' && (op.where ?? 'params') === 'params') {
        const states = this.states(n);
        const node = states && resolveAddress(states[i], op.at);
        const spec = node ? own(target.kinds, node.kind) : undefined;
        const t = spec?.params ? own(spec.params, op.key) : undefined;
        if (t) this.value(t, op.value, [...at, 'value'], inner);
      }
    });
  }
}

class MapRefIndex implements RefIndex {
  private readonly byFirst = new Map<string, Site[]>();
  private danglingSites?: Site[];

  constructor(
    readonly pass: RefPass,
    readonly all: readonly Site[],
  ) {
    for (const s of all) {
      for (const f of firsts(s.runs)) {
        let list = this.byFirst.get(f);
        if (!list) this.byFirst.set(f, (list = []));
        list.push(s);
      }
    }
  }

  under(prefix: Address): readonly Site[] {
    if (!prefix.length) return this.all;
    return (this.byFirst.get(prefix[0]) ?? []).filter((s) => s.runs.some((r) => isUnder(r, prefix)));
  }

  dangling(): readonly Site[] {
    return (this.danglingSites ??= this.all.filter((s) => !this.pass.names(s.run, s.step)));
  }

  derive(after: Node): RefIndex {
    return refIndexOf(after, this.pass.language, this.pass.options);
  }
}

/** The first ids of a ref's runs: what the index keys it by (`''` for a run naming the root). */
function firsts(runs: readonly Address[]): Set<string> {
  const out = new Set<string>();
  for (const r of runs) out.add(r.length ? r[0] : '');
  return out;
}

const Indices = new WeakMap<Node, WeakMap<Language, Map<ValidateOptions | undefined, MapRefIndex>>>();

/** The ref index of a revision, cached per root, language and options. */
export function refIndexOf(root: Node, language: Language, options?: ValidateOptions): RefIndex {
  return mapIndexOf(root, language, options);
}

function mapIndexOf(root: Node, language: Language, options?: ValidateOptions): MapRefIndex {
  let byLanguage = Indices.get(root);
  if (!byLanguage) Indices.set(root, (byLanguage = new WeakMap()));
  let byOptions = byLanguage.get(language);
  if (!byOptions) byLanguage.set(language, (byOptions = new Map()));
  let index = byOptions.get(options);
  if (!index) {
    const pass = new RefPass(root, language, options);
    byOptions.set(options, (index = new MapRefIndex(pass, pass.sites())));
  }
  return index;
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

  /** A pass over the draft as it is now, for walking new values (a stored patch's dry run reads the spec). */
  private get pass(): RefPass {
    return (this.overlayPass ??= new RefPass(() => this.work.finalize(), this.language, this.options));
  }

  key(d: Draft, key: string) {
    if (d.origin && !d.made)
      (this.superseded.get(d.origin) ?? this.superseded.set(d.origin, new Set()).get(d.origin)!).add(key);
    this.drop(d, key);
    const v = d.params ? own(d.params, key) : undefined;
    if (v === undefined) return;
    // a stored patch in the value is typed over the spec as it is now
    this.overlayPass = undefined;
    this.pass.walkKey(d.kind, key, v, (rel, value, step) => this.add(d, rel, value, step));
  }

  made(d: Draft) {
    this.overlayPass = undefined;
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
