/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The transform graph (seqspec:viewer#updates-and-status): one node per document node a node kind is registered for,
// plus the derived stages, keyed by ref. Edges come from the node kinds (the format's references); a node reruns only
// when its own content, an input's result or status, or a part its kind adds to the key (a set entity, a registry
// entry, the skin) changed. Runs are async with an AbortSignal; a run whose key changed is aborted and its result
// dropped. When an update settles, `assembled` emits once. VC1: a throw in a kind's `inputs` or `key` is that node's
// error, and an update always publishes; superseded results are disposed; issues are capped per node; a node whose
// params go back to the previous object keeps its token.
//
// Node records are objects in a Map, not struct-of-arrays: the graph runs at edge rate (a document or registry change).

import {
  canonical,
  type DocIndex,
  forEachOf,
  type IndexedNode,
  indexDoc,
  Layout,
  type Loaders,
  type NormalizedParams,
  type ResolvedView,
  type SchemeSubject,
  type SeqViewSpecDoc,
  type SourceResult,
  type StageContext,
  type SvsRegistries,
  sourceParse,
} from '@molstar/seq-core/format';
import type { Issue, SeqSet } from '@molstar/seq-core/model';
import type { SeqIssue } from './context.js';
import type { DynamicSource } from './dynamic.js';
import type { SeqViewerRegistry } from './registry.js';
import { type SeqState, state, type WritableSeqState } from './state.js';
import type { Theme } from './themes.js';

/** An issue's rank by severity, the most severe first: error, warning, then a note. */
export const severityRank = (i: { readonly severity: string }): number =>
  i.severity === 'error' ? 0 : i.severity === 'warning' ? 1 : 2;

export type Stage = 'load' | 'complex' | 'assemble' | 'compute' | 'section' | 'layer' | 'prepare' | 'theme' | 'view';
/** Per-view stages run for the current view only. */
const PER_VIEW: ReadonlySet<string> = new Set(['section', 'layer', 'prepare', 'theme', 'view']);

/** How a node reads another: its result (keyed by identity), its status (keyed by state), or only after it (no key). */
export interface InputCollector {
  result(ref: string): void;
  status(ref: string): void;
  after(ref: string): void;
}
/**
 * seqspec:viewer#extension-points. `inputs` declares edges (a document node's refs from core's `nodeRefs`, with `svs`,
 * the format registries); the key is its own content, its inputs, and what `key` adds.
 */
export interface NodeKind<R = unknown> {
  readonly name: string;
  readonly stage: Stage;
  inputs(ref: string, node: IndexedNode | null, doc: DocIndex, out: InputCollector, svs: SvsRegistries): void;
  key?(rc: RunContext, out: unknown[]): void;
  run(rc: RunContext): R | Promise<R>;
  equals?(a: R, b: R): boolean;
  dispose?(result: R): void;
}
export interface RunContext {
  readonly ref: string;
  /** The document node (a derived stage's: its layer or view); null for assemblies and data entries. */
  readonly node: IndexedNode | null;
  readonly doc: DocIndex;
  /** The run's signal, the format registries, the set this stage reads, sources and unavailable refs; `report`. */
  readonly stage: StageContext;
  readonly registry: SeqViewerRegistry;
  readonly loaders: Loaders;
  readonly skin: Skin;
  /** The node's previous result, for assemblies and the axis identity. */
  readonly prev: unknown;
  /** A data entry given to `setData` (its `#entry` or `#set` node). */
  readonly payload: unknown;
  /** The stream manager's dynamic nodes, for `#assemble.stream` (null: none). */
  readonly dynamic: DynamicSource | null;
  /** The `result` edges, in declared order. */
  readonly inputs: readonly string[];
  result(ref: string): unknown;
  status(ref: string): NodeStatus | null;
  /** True when the node, or the node its unavailability comes from, failed (`error`), not a missing provider. */
  failed(ref: string): boolean;
  /** A document node's own content: changes when its params or `custom` change (compared normalized). */
  token(ref: string): unknown;
  supported(layout: number): ReadonlySet<string>;
  /** The layer's color scheme picked at runtime (`ui.colorSchemes`); undefined: its own. */
  colorScheme(layer: string): string | undefined;
  /** Every color scheme picked at runtime, by layer ref (a new map when one changes): a template's key. */
  readonly colorSchemes: ReadonlyMap<string, string> | undefined;
  /** The layer's color range set at runtime (`ui.colorRanges`); undefined: its own. */
  colorRange(layer: string): readonly [number, number] | undefined;
  /** Every color range set at runtime, by layer ref (a new map when one changes): a template's key. */
  readonly colorRanges: ReadonlyMap<string, readonly [number, number]> | undefined;
}

export interface NodeStatus {
  readonly state: 'pending' | 'running' | 'ok' | 'error' | 'unavailable';
  readonly issues: readonly Issue[];
  /** Unavailable: the ref whose failure or unavailability caused it. */
  readonly cause: string | null;
}
/**
 * The viewer's look (seqspec:interaction#the-skin): the `--seq-*` tokens as canvas fill styles, the font families, and
 * `dark` (the graph's theme nodes read it). Written by the canvas engine's skin reader, or a host's `setSkin`.
 */
export interface Skin {
  readonly dark: boolean;
  readonly surface: string;
  readonly surfaceAlt: string;
  readonly ink: string;
  readonly inkQuiet: string;
  readonly inkFaint: string;
  readonly grid: string;
  readonly accent: string;
  readonly hover: string;
  readonly highlight: string;
  readonly selection: string;
  readonly selectionBorder: string;
  readonly unobserved: string;
  /** A map whose data has not arrived (the vertical overview of lazy rows): never the ground, so it never reads as empty. */
  readonly pending: string;
  /** A window's end grips: at rest, and lit (the pointer on that end, or on the window's body). */
  readonly grip: string;
  readonly gripHover: string;
  /** A window's fill (a navigator's box, a ring's window, CG3): a neutral wash, never a selection's. */
  readonly window: string;
  readonly conflict: string;
  readonly focusRing: string;
  readonly fontMono: string;
  readonly fontUi: string;
}
export interface PreparedLayer {
  readonly renderer: string;
  readonly state: unknown;
  /** Per channel of its type, in order; null where the resolved view's reference is null. */
  readonly themes: readonly (Theme | null)[];
  /** What its color is about, for the color scheme picker (seqspec:format#color-schemes); null: no scheme fits it. */
  readonly subject: SchemeSubject | null;
}
export interface ViewModel {
  readonly ref: string;
  readonly resolved: ResolvedView;
  /** By resolved layer index. */
  readonly prepared: readonly PreparedLayer[];
  readonly set: SeqSet;
}
export interface Assembled {
  readonly set: SeqSet;
  readonly view: ViewModel | null;
  /** Refs whose result changed (or that left) since the last emission. */
  readonly changed: ReadonlySet<string>;
}
export interface SeqGraph {
  readonly assembled: SeqState<Assembled | null>;
  readonly statuses: SeqState<ReadonlyMap<string, NodeStatus>>;
  readonly issues: SeqState<readonly SeqIssue[]>;
  status(ref: string): NodeStatus | null;
  /** Node runs since creation, for tests. */
  readonly runs: number;
  /** Refs run in the last settled update (superseded passes folded in), for tests. */
  readonly ran: ReadonlySet<string>;
  /** Layer types whose renderers draw the layout (0 tracks, 1 text), implemented themes, `char` and `path`. */
  supported(layout: number): ReadonlySet<string>;
  /** Runs an update, folded with any in flight; resolves when one settles. */
  update(): Promise<void>;
}

/** A data entry of `setData`'s build inputs, by its list. */
export interface DataEntry {
  readonly list: 'molecules' | 'units' | 'complexes' | 'alignments' | 'mappings' | 'annotations';
  readonly entry: { readonly key: string };
}
export type DataMode = { readonly entries: ReadonlyMap<string, DataEntry> } | { readonly set: SeqSet } | null;

interface Run {
  readonly key: unknown[];
  readonly controller: AbortController;
  promise: Promise<void>;
}
const PENDING: NodeStatus = { state: 'pending', issues: [], cause: null };
/** A node's record, kept across updates while its ref and kind stay. */
class GNode {
  readonly ref: string;
  readonly kind: NodeKind;
  node: IndexedNode | null = null;
  view: string | null = null;
  payload: unknown;
  results: string[] = [];
  statuses: string[] = [];
  afters: string[] = [];
  key: unknown[] | null = null;
  result: unknown;
  status: NodeStatus = PENDING;
  issues: readonly SeqIssue[] = [];
  run: Run | null = null;
  cyclic = false;
  /** What its kind's `inputs` threw: the node is an error and never runs. */
  thrown: unknown = null;
  constructor(ref: string, kind: NodeKind) {
    this.ref = ref;
    this.kind = kind;
  }
}
interface Token {
  params: unknown;
  custom: unknown;
  version: number;
  token: object;
  /** The one before, kept: params back to that object (a preview cleared) keep its token. */
  prev: Token | null;
}
interface Pass {
  superseded: boolean;
}

const LISTS = ['molecules', 'units', 'complexes', 'alignments', 'mappings', 'annotations'] as const;
const isEmptyPart = (r: Record<string, unknown>) =>
  r.source === undefined && LISTS.every((k) => !(r[k] as unknown[] | undefined)?.length);
const sameKey = (a: readonly unknown[] | null, b: readonly unknown[]) =>
  a !== null && a.length === b.length && a.every((x, i) => x === b[i]);
const mapEquals = <K, V>(a: ReadonlyMap<K, V>, b: ReadonlyMap<K, V>) =>
  a.size === b.size && [...a].every(([k, v]) => b.get(k) === v);
const listEquals = <T>(a: readonly T[], b: readonly T[]) => a.length === b.length && a.every((x, i) => x === b[i]);

export interface GraphHost {
  readonly registry: SeqViewerRegistry;
  readonly loaders: Loaders;
  readonly skin: SeqState<Skin>;
  /** Layer ref → a color scheme picked at runtime (`ui.colorSchemes`); absent: none. */
  readonly colorSchemes?: SeqState<ReadonlyMap<string, string>>;
  /** Layer ref → a color range set at runtime (`ui.colorRanges`); absent: none. */
  readonly colorRanges?: SeqState<ReadonlyMap<string, readonly [number, number]>>;
  /** The document the graph runs: the shown one, or `setData`'s views. */
  doc(): SeqViewSpecDoc | null;
  currentView(): string | null;
  onIssue(issue: SeqIssue): void;
  /** `Issues.PerNode`. */
  issueLimit(): number;
  /** The dynamic nodes the stream manager holds (seqspec:sources#the-viewers-side). */
  dynamic?(): DynamicSource | null;
}

export class Graph implements SeqGraph {
  private readonly host: GraphHost;
  readonly assembled: WritableSeqState<Assembled | null> = state<Assembled | null>(null);
  readonly statuses: WritableSeqState<ReadonlyMap<string, NodeStatus>> = state(new Map(), mapEquals);
  readonly issues: WritableSeqState<readonly SeqIssue[]> = state<readonly SeqIssue[]>([], listEquals);
  private nodes = new Map<string, GNode>();
  private order: GNode[] = [];
  private readonly tokens = new Map<string, Token>();
  private readonly normCache = new WeakMap<object, NormalizedParams>();
  private readonly entryParts = new WeakMap<object, object>();
  private index: DocIndex = { nodes: new Map(), order: [] };
  private indexed: { doc: SeqViewSpecDoc | null; version: number; mode: DataMode } = {
    doc: null,
    version: -1,
    mode: null,
  };
  private readonly supportedSets: { version: number; set: ReadonlySet<string> }[] = [];
  private readonly sources = new Map<string, SourceResult>();
  private readonly gone = new Set<string>();
  private data: DataMode = null;
  private current: Pass | null = null;
  private queued = false;
  private waiters: (() => void)[] = [];
  private changed = new Set<string>();
  private fresh: SeqIssue[] = [];
  private running = new Set<string>();
  private disposed = false;
  ran: ReadonlySet<string> = new Set();
  runs = 0;

  constructor(host: GraphHost) {
    this.host = host;
  }

  status(ref: string) {
    return this.nodes.get(ref)?.status ?? null;
  }

  supported(layout: number): ReadonlySet<string> {
    const reg = this.host.registry;
    const v = reg.version.value;
    const held = this.supportedSets[layout];
    if (held?.version === v) return held.set;
    const geometry = layout === Layout.text ? 'wrapped' : 'linear';
    const next = new Set([
      ...reg.list('renderers').flatMap((r) => (r.geometries.includes(geometry) ? [r.name] : [])),
      ...reg.list('themes').map((t) => t.name),
      'char',
      'path',
    ]);
    // An unchanged set keeps its identity, so it is no key change.
    const set = held && held.set.size === next.size && [...next].every((x) => held.set.has(x)) ? held.set : next;
    this.supportedSets[layout] = { version: v, set };
    return set;
  }

  /**
   * Data without a document (seqspec:viewer#data-mode); null: back to documents. A change of kind starts a new graph.
   */
  setData(mode: DataMode) {
    const kind = (m: DataMode) => (m === null ? 0 : 'set' in m ? 1 : 2);
    if (kind(mode) !== kind(this.data)) this.reset();
    this.data = mode;
  }

  /** Drops every node (a new lineage follows). */
  reset() {
    for (const g of this.nodes.values()) this.drop(g);
    this.nodes = new Map();
    this.order = [];
    this.indexed = { doc: null, version: -1, mode: null };
  }

  /** Runs an update (folded with any in flight); resolves when one settles. */
  update(): Promise<void> {
    if (this.disposed) return Promise.resolve();
    const done = new Promise<void>((resolve) => this.waiters.push(resolve));
    if (this.current) this.current.superseded = true;
    if (!this.queued) {
      this.queued = true;
      void Promise.resolve().then(() => {
        this.queued = false;
        if (this.disposed) return;
        const pass: Pass = { superseded: false };
        this.current = pass;
        void this.pass(pass);
      });
    }
    return done;
  }

  // --- structure -------------------------------------------------------------------------------------------------

  private structure() {
    const doc = this.host.doc();
    const version = this.host.registry.version.value;
    const i = this.indexed;
    if (i.doc === doc && i.version === version && i.mode === this.data) return;
    this.index = doc ? indexDoc(doc, this.host.registry.svs, this.normCache) : { nodes: new Map(), order: [] };
    const next = new Map<string, GNode>();
    const add = (ref: string, name: string, node: IndexedNode | null, view: string | null, payload?: unknown) => {
      const kind = this.host.registry.get('nodeKinds', name);
      if (!kind || next.has(ref)) return;
      let g = this.nodes.get(ref);
      if (!g || g.kind !== kind) g = new GNode(ref, kind); // the old one is dropped below
      Object.assign(g, { node, view: PER_VIEW.has(kind.stage) ? view : null, payload, cyclic: false, thrown: null });
      const r: string[] = [];
      const s: string[] = [];
      const a: string[] = [];
      try {
        kind.inputs(
          ref,
          node,
          this.index,
          {
            result: (x) => r.push(x),
            status: (x) => s.push(x),
            after: (x) => a.push(x),
          },
          this.host.registry.svs,
        );
      } catch (e) {
        g.thrown = e ?? new Error('threw');
      }
      [g.results, g.statuses, g.afters] = [r, s, a];
      next.set(ref, g);
    };
    const views = new Map<string, string | null>();
    for (const n of this.index.order) {
      const view = n.kind === 'view' ? n.ref : (views.get(n.parent ?? '') ?? null);
      views.set(n.ref, view);
      // A template's layers resolve per instance, in its `#template` node (seqspec:format#placement).
      if (n.kind === 'layer' && forEachOf(this.index.nodes.get(n.parent ?? ''))) continue;
      add(n.ref, n.kind, n, view);
      if (n.kind === 'layer') for (const d of ['#prepare', '#theme']) add(`${d}:${n.ref}`, d, n, view);
      if (forEachOf(n)) add(`#template:${n.ref}`, '#template', n, view);
      if (n.kind === 'view') add(`#view:${n.ref}`, '#view', n, n.ref);
    }
    const data = this.data;
    if (data && 'set' in data) add('#set', '#set', null, null, data.set);
    for (const e of data && 'entries' in data ? data.entries.values() : []) {
      // Memoized by the entry's identity: a host that keeps an entry reruns nothing for it.
      const part = this.entryParts.get(e.entry) ?? { [e.list]: [e.entry] };
      this.entryParts.set(e.entry, part);
      add(e.entry.key, '#entry', null, null, part);
    }
    for (const a of ['#assemble.base', '#assemble.mapped', '#assemble', '#assemble.stream']) add(a, a, null, null);
    // Data mode: entries feed the assemblies (a whole set seeds the first).
    const base = next.get('#assemble.base');
    const mapped = next.get('#assemble.mapped');
    if (data && base && mapped) {
      if ('set' in data) base.results.unshift('#set');
      else
        for (const e of data.entries.values())
          (e.list === 'mappings' || e.list === 'annotations' ? mapped : base).results.push(e.entry.key);
    }
    for (const g of this.nodes.values()) if (next.get(g.ref) !== g) this.drop(g);
    this.nodes = next;
    this.sort();
    this.indexed = { doc, version, mode: this.data };
  }

  /** Topological order (Kahn); nodes in a cycle are `svs.ref.cycle` errors and never run, their readers unavailable. */
  private sort() {
    const all = [...this.nodes.values()];
    const ins = (g: GNode) =>
      [...g.results, ...g.statuses, ...g.afters].flatMap((r) => {
        const n = this.nodes.get(r);
        return n && n !== g ? [n] : [];
      });
    const kahn = (skip: ReadonlySet<GNode>) => {
      const indeg = new Map<GNode, number>();
      const out = new Map<GNode, GNode[]>(all.map((g) => [g, []]));
      for (const g of all) {
        const from = skip.has(g) ? [] : ins(g).filter((n) => !skip.has(n));
        for (const n of from) out.get(n)?.push(g);
        indeg.set(g, skip.has(g) ? -1 : from.length);
      }
      const order = all.filter((g) => indeg.get(g) === 0);
      for (let k = 0; k < order.length; k++) {
        for (const d of out.get(order[k]) ?? []) {
          indeg.set(d, (indeg.get(d) ?? 0) - 1);
          if (indeg.get(d) === 0) order.push(d);
        }
      }
      return order;
    };
    let order = kahn(new Set());
    if (order.length < all.length) {
      // Left over: the cycles' members, and what reads them; only members are errors.
      const sorted = new Set(order);
      const reaches = (from: GNode, to: GNode, seen: Set<GNode>): boolean =>
        ins(from).some((x) => x === to || (!sorted.has(x) && !seen.has(x) && reaches(x, to, seen.add(x))));
      const cyclic = new Set(all.filter((g) => !sorted.has(g) && reaches(g, g, new Set())));
      for (const g of cyclic) {
        g.cyclic = true;
        const message = `"${g.ref}" is in a reference cycle`;
        const issue = { code: 'svs.ref.cycle', severity: 'error' as const, path: g.node?.path ?? '', message };
        this.setStatus(g, { state: 'error', issues: [issue], cause: null });
      }
      order = kahn(cyclic);
    }
    this.order = order;
  }

  private drop(g: GNode) {
    g.run?.controller.abort();
    g.run = null;
    this.release(g, g.result);
    this.sources.delete(g.ref);
    this.gone.delete(g.ref);
    this.tokens.delete(g.ref);
    this.changed.add(g.ref);
  }

  // --- running ---------------------------------------------------------------------------------------------------

  /** Always publishes unless superseded, so every waiter resolves. */
  private async pass(pass: Pass) {
    let view: string | null = null;
    try {
      this.structure();
      view = this.host.currentView();
      const evaluated = new Map<GNode, Promise<void>>();
      const evaluate = (g: GNode): Promise<void> => {
        let p = evaluated.get(g);
        if (!p) {
          p = this.evaluate(g, pass, evaluate);
          evaluated.set(g, p);
        }
        return p;
      };
      await Promise.all(this.order.filter((g) => g.view === null || g.view === view).map(evaluate));
    } catch (e) {
      const message = `the graph's update threw: ${e instanceof Error ? e.message : String(e)}`;
      this.fresh.push({ code: 'seqv.state.error', severity: 'error', path: '', message, ref: null });
    } finally {
      if (!pass.superseded && !this.disposed) {
        this.current = null;
        this.publish(view);
      }
    }
  }

  private async evaluate(g: GNode, pass: Pass, evaluate: (g: GNode) => Promise<void>) {
    const inputs = [...g.results, ...g.statuses, ...g.afters].flatMap((r) => {
      const n = this.nodes.get(r);
      return n && !n.cyclic && n !== g ? [evaluate(n)] : [];
    });
    if (inputs.length > 0) await Promise.all(inputs);
    if (pass.superseded || g.cyclic) return;
    const controller = new AbortController();
    const issues: Issue[] = [];
    const rc = this.context(g, controller.signal, issues);
    let thrown = g.thrown;
    const key: unknown[] = [g.kind, g.node ? rc.token(g.node.ref) : g.payload];
    try {
      for (const r of g.results) key.push(this.nodes.get(r)?.result);
      for (const r of g.statuses) key.push(this.nodes.get(r)?.status.state);
      if (thrown === null) g.kind.key?.(rc, key);
    } catch (e) {
      thrown = e ?? new Error('threw');
    }
    if (thrown === null && sameKey(g.key, key)) {
      // Back to the published key (a preview cleared): a run in flight for another key is dropped.
      g.run?.controller.abort();
      g.run = null;
      return;
    }
    if (thrown === null && g.run && sameKey(g.run.key, key)) return g.run.promise;
    g.run?.controller.abort();
    const run: Run = { key, controller, promise: Promise.resolve() };
    g.run = run;
    run.promise = this.start(g, run, rc, issues, thrown);
    return run.promise;
  }

  private async start(g: GNode, run: Run, rc: RunContext, issues: Issue[], thrown: unknown) {
    this.runs++;
    this.running.add(g.ref);
    let result: unknown;
    let threw: unknown = thrown;
    if (threw === null) {
      try {
        result = await g.kind.run(rc);
      } catch (e) {
        threw = e ?? new Error('threw');
      }
    }
    // A superseded or aborted run is dropped: its result never publishes, and is disposed.
    if (g.run !== run || run.controller.signal.aborted || this.disposed) {
      if (result !== g.result) this.release(g, result);
      return;
    }
    g.run = null;
    g.key = run.key;
    if (threw) {
      const message = threw instanceof Error ? threw.message : String(threw);
      const code = PER_VIEW.has(g.kind.stage) ? 'seqv.renderer.error' : 'svs.source.load';
      issues.push({ code, severity: 'error', path: g.node?.path ?? '', message });
    }
    const status = this.classify(g, threw ? undefined : result, issues, threw !== null);
    if (status.state === 'ok' && g.result !== undefined && g.kind.equals?.(g.result, result)) {
      this.release(g, result);
      result = g.result;
    }
    if (result !== g.result) {
      this.release(g, g.result);
      g.result = result;
      this.changed.add(g.ref);
    }
    const source = (result as { source?: SourceResult } | undefined)?.source;
    if (source) this.sources.set(g.ref, source);
    else this.sources.delete(g.ref);
    this.setStatus(g, status);
  }

  /** Status: failures are errors, their readers unavailable with a cause; a missing provider is not a failure. */
  private classify(g: GNode, result: unknown, issues: Issue[], threw: boolean): NodeStatus {
    const stage = g.kind.stage;
    const data = stage === 'load' || stage === 'complex' || stage === 'compute';
    const r = result as Record<string, unknown> | null | undefined;
    const bad = [...g.results, ...g.statuses].find((x) => {
      const s = this.nodes.get(x)?.status.state;
      return s === 'error' || s === 'unavailable';
    });
    let state: NodeStatus['state'] = 'ok';
    const kind = (r?.source as SourceResult | undefined)?.kind;
    if (threw || kind === 'failed') state = 'error';
    else if (data || stage === 'layer') {
      if (r == null || r.unavailable === true || kind === 'unavailable') state = 'unavailable';
      else if (data && isEmptyPart(r) && bad) state = 'unavailable';
      else if (data && isEmptyPart(r) && issues.some((i) => i.severity === 'error')) state = 'error';
    }
    const own = issues.some((i) => i.code === 'svs.provider.unavailable');
    return { state, issues, cause: state === 'unavailable' && !own ? (bad ?? null) : null };
  }

  /** A node's result, disposed by its kind; a throw there is reported, never thrown. */
  private release(g: GNode, result: unknown) {
    try {
      if (result !== undefined) g.kind.dispose?.(result);
    } catch (e) {
      const message = `"${g.ref}"'s dispose threw: ${e instanceof Error ? e.message : String(e)}`;
      this.fresh.push({ code: 'seqv.renderer.error', severity: 'error', path: '', message, ref: g.ref });
    }
  }

  /**
   * Its issues, at most `Issues.PerNode`, the most severe kept first (in order within a severity), the rest counted in
   * one more that says how many of each severity.
   */
  private setStatus(g: GNode, status: NodeStatus) {
    const limit = this.host.issueLimit();
    if (status.issues.length > limit) {
      const all = [...status.issues].sort((a, b) => severityRank(a) - severityRank(b));
      const n = [0, 0, 0];
      for (let k = limit; k < all.length; k++) n[severityRank(all[k])]++;
      const each = [`${n[0]} errors`, `${n[1]} warnings`, `${n[2]} notes`].filter((_, k) => n[k] > 0);
      const message = `and ${all.length - limit} more (${each.join(', ')})`;
      const more: Issue = { code: 'seqv.issues.truncated', severity: 'warning', path: g.node?.path ?? '', message };
      status = { ...status, issues: [...all.slice(0, limit), more] };
    }
    g.status = status;
    g.issues = status.issues.map((i) => ({ ...i, ref: g.ref }));
    for (const i of g.issues) this.fresh.push(i);
    if (status.state === 'ok') this.gone.delete(g.ref);
    else this.gone.add(g.ref);
  }

  private context(g: GNode, signal: AbortSignal, issues: Issue[]): RunContext {
    const nodes = this.nodes;
    const index = this.index;
    const seed = [...g.results, ...g.afters].map((r) => nodes.get(r)?.result).find((x) => (x as { set?: SeqSet })?.set);
    const final = this.final();
    const set = (seed as { set: SeqSet } | undefined)?.set ?? (PER_VIEW.has(g.kind.stage) ? final : null);
    const sources = this.sources;
    const stage: StageContext = {
      signal,
      registries: this.host.registry.svs,
      set,
      sources,
      unavailable: this.gone,
      source: (ref) => sourceParse(index, sources, ref),
      report: (i) => issues.push(i),
    };
    const status = (ref: string) => nodes.get(ref)?.status ?? null;
    const failed = (ref: string): boolean => {
      for (let s = status(ref), k = 0; s && k < 64; s = s.cause ? status(s.cause) : null, k++) {
        if (s.state === 'error') return true;
        if (!s.cause) return false;
      }
      return false;
    };
    return {
      ref: g.ref,
      node: g.node,
      doc: index,
      stage,
      registry: this.host.registry,
      loaders: this.host.loaders,
      skin: this.host.skin.value,
      prev: g.result,
      payload: g.payload,
      dynamic: this.host.dynamic?.() ?? null,
      inputs: g.results,
      result: (ref) => nodes.get(ref)?.result,
      status,
      failed,
      token: (ref) => this.token(ref),
      supported: (layout) => this.supported(layout),
      colorScheme: (layer) => this.host.colorSchemes?.value.get(layer),
      colorSchemes: this.host.colorSchemes?.value,
      colorRange: (layer) => this.host.colorRanges?.value.get(layer),
      colorRanges: this.host.colorRanges?.value,
    };
  }

  /** The final set: the document's with the dynamic nodes held (`#assemble.stream`). */
  private final(): SeqSet | null {
    const g = this.nodes.get('#assemble.stream') ?? this.nodes.get('#assemble');
    return (g?.result as { set?: SeqSet } | undefined)?.set ?? null;
  }

  /** A node's own content identity: new when its normalized params or `custom` change, kept across a registration. */
  private token(ref: string): unknown {
    const n = this.index.nodes.get(ref);
    if (!n) return undefined;
    const version = this.host.registry.svs.version;
    const t = this.tokens.get(ref);
    if (t && t.params === n.params && t.custom === n.node.custom) return t.token;
    const p = t?.prev;
    if (t && p && p.params === n.params && p.custom === n.node.custom) {
      this.tokens.set(ref, { ...p, prev: { ...t, prev: null } });
      return p.token;
    }
    const same = t && t.version !== version && canonical([t.params, t.custom]) === canonical([n.params, n.node.custom]);
    const prev = t ? { ...t, prev: null } : null;
    const next = { params: n.params, custom: n.node.custom, version, token: same && t ? t.token : {}, prev };
    this.tokens.set(ref, next);
    return next.token;
  }

  private publish(view: string | null) {
    const statuses = new Map<string, NodeStatus>();
    const issues: SeqIssue[] = [];
    for (const g of this.nodes.values()) {
      statuses.set(g.ref, g.status);
      for (const i of g.issues) issues.push(i);
    }
    this.statuses.set(statuses, 'graph');
    this.issues.set(issues, 'graph');
    const fresh = this.fresh;
    this.fresh = [];
    for (const i of fresh) this.host.onIssue(i);
    this.ran = this.running;
    this.running = new Set();
    const set = this.final();
    const model = (view ? this.nodes.get(`#view:${view}`)?.result : null) as ViewModel | undefined;
    const prev = this.assembled.value;
    if (!set) this.assembled.set(null, 'graph');
    else if (this.changed.size > 0 || prev?.set !== set || prev?.view !== (model ?? null)) {
      this.assembled.set({ set, view: model ?? null, changed: this.changed }, 'graph');
    }
    this.changed = new Set();
    const waiters = this.waiters;
    this.waiters = [];
    for (const w of waiters) w();
  }

  /** seqspec:viewer#the-context: running work aborted, results disposed, states completed. */
  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    for (const g of this.nodes.values()) this.drop(g);
    this.nodes.clear();
    for (const w of this.waiters) w();
    this.waiters = [];
    this.assembled.complete();
    this.statuses.complete();
    this.issues.complete();
  }
}
