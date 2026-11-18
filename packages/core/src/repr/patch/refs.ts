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
  unionNodeMember,
  type ValidateOptions,
  validateValue,
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
import { patchFromNode, patchTarget } from './language.js';
import { isUnder, Work } from './work.js';

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
interface Step {
  readonly slot: string;
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

type Emit = (path: TreePath, value: string, step?: Step) => void;

/** The pass over one root: its sites, and the dry runs of its stored patches. */
export class RefPass {
  private readonly resolve: (id: string) => Language | undefined;
  private readonly dry = new Map<string, Node[]>();

  constructor(
    readonly root: Node,
    readonly language: Language,
    readonly options: ValidateOptions = {},
  ) {
    this.resolve = (id) => (id === language.id ? language : options.languages?.(id));
  }

  /** Every site of the tree, depth first. */
  sites(): Site[] {
    const out: Site[] = [];
    const go = (n: Node, path: (string | number)[], at: Address | undefined, lang: Language) => {
      this.node(lang, n, path, (p, v, step) => {
        const site = makeSite(v, p, at, step);
        if (site) out.push(site);
      });
      n.children?.forEach((c, i) =>
        go(c, [...path, 'children', i], at && c.id !== undefined ? [...at, c.id] : undefined, lang),
      );
    };
    go(this.root, [], [], this.language);
    return out;
  }

  /** The sites in one node of the tree (its params and its slots' bodies), at `path`, the node's own step path. */
  nodeSites(n: Node, path: TreePath, at: Address | undefined): Site[] {
    const out: Site[] = [];
    this.node(this.language, n, path, (p, v, step) => {
      const site = makeSite(v, p, at, step);
      if (site) out.push(site);
    });
    return out;
  }

  /** Whether a site names a node: in the root, or over its stored patch's dry run. */
  names(site: Site): boolean {
    const root = site.step ? this.stepRoot(site.step) : this.root;
    return !!root && resolveAddress(root, site.run) !== undefined;
  }

  private stepRoot(step: Step): Node | undefined {
    const states = this.dry.get(step.slot);
    if (!states) return undefined;
    return states[Math.min(step.op + (step.after ? 1 : 0), states.length - 1)];
  }

  /** A node's params by its kind, its data form's args (not its children: the caller walks them). */
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
    for (const key of Object.keys(values)) {
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
        const m = isNode(v)
          ? unionNodeMember(t.of, v.kind, this.resolve)
          : t.of.find(
              (x) => validateValue(x, v, { ...this.options, languages: this.resolve, maxIssues: 1 }).length === 0,
            );
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
    const slot = path.join('/');
    const states: Node[] = [this.root];
    const work = new Work(this.root);
    for (let i = 0; i < patch.ops.length; i++) {
      const failed = work.apply(patch.ops[i], i);
      states.push(failed ? states[states.length - 1] : work.finalize());
    }
    this.dry.set(slot, states);
    patch.ops.forEach((op, i) => {
      const at = [...path, 'children', i, 'params'];
      if ('at' in op) emit([...at, 'at'], printAddress(op.at), { slot, op: i, after: false });
      if ('parent' in op) emit([...at, 'parent'], printAddress(op.parent), { slot, op: i, after: false });
      if (!target) return;
      const inner: Emit = (p, v) => emit(p, v, { slot, op: i, after: true });
      if (op.op === 'insert' || op.op === 'replace') this.slot(target, op.node, [...at, 'node'], inner);
      if (op.op === 'set' && (op.where ?? 'params') === 'params') {
        const node = resolveAddress(states[i], op.at);
        const spec = node ? own(target.kinds, node.kind) : undefined;
        const t = spec?.params ? own(spec.params, op.key) : undefined;
        if (t) this.value(t, op.value, [...at, 'value'], inner);
      }
    });
  }
}

function makeSite(value: string, path: TreePath, at: Address | undefined, step: Step | undefined): Site | undefined {
  const pattern = parseAddress(value).pattern;
  if (!pattern) return undefined; // not an address: its spelling is validation's
  const runs = runsOf(pattern);
  const site: { -readonly [K in keyof Site]: Site[K] } = { path, value, run: runs[0], runs };
  if (at) site.at = at;
  if (step) site.step = step;
  return site;
}

class MapRefIndex implements RefIndex {
  private readonly byFirst = new Map<string, Site[]>();
  private readonly all: Site[];
  private danglingSites?: Site[];

  constructor(
    private readonly pass: RefPass,
    sites: Site[],
  ) {
    this.all = sites;
    for (const s of sites) {
      const firsts = new Set(s.runs.map((r) => (r.length ? r[0] : '')));
      for (const f of firsts) {
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
    return (this.danglingSites ??= this.all.filter((s) => !this.pass.names(s)));
  }

  derive(after: Node): RefIndex {
    return refIndexOf(after, this.pass.language, this.pass.options);
  }
}

const Indices = new WeakMap<Node, WeakMap<Language, Map<ValidateOptions | undefined, RefIndex>>>();

/** The ref index of a revision, cached per root, language and options. */
export function refIndexOf(root: Node, language: Language, options?: ValidateOptions): RefIndex {
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
