/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { validateAndNormalize } from '@molstar/core/repr/schema';
import {
  CategoryColumn,
  type Issue,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Segment,
  type Value,
} from '@molstar/core/repr/tree';
import type { FormatRuntime, NodeEntry } from '@molstar/seqx-core/data';
import { runtimeFor } from '@molstar/seqx-core/data/channels';
import { checkDataNode, LangIssues, SeqLanguage, SeqLanguageVersion } from '@molstar/seqx-core/lang';
import { ResolveIssues } from './issues.js';
import {
  type AnnotationEntry,
  type BoundSourceEntry,
  DefaultResolveConfig,
  type MapEntry,
  type PathResult,
  pathIssue,
  type ResolveConfig,
  resolveAnnotationNode,
  resolveMapNode,
  type StaticModel,
  type SystemEntry,
} from './model.js';
import { type Edge, findPath, type Step } from './paths.js';

// The overlay (_spec/README.md "Overlay"): tier 2, the dynamic nodes the data store holds, over the static model and
// never part of it. Each node is checked by the language description as a node of its kind under `data`; its refs
// name spec nodes or held dynamic nodes, else it is left out with an issue. Dynamic annotations resolve as the spec's
// do; an item of a ranges node with `map: 'codon'` gets a coordinate system of its residues and a codon map onto the
// target, made when something names it; paths with a dynamic end are searched over the spec's maps and the overlay's.
// Nested sources (source nodes) are listed only when a template needs their depth, never past the configured limit
// (EC-179). One overlay a node version of the store: its version is part of the draw lists' key.

/** What the overlay reads of the store: `DataStore` is one. */
export interface NodeReader {
  /** The node entries held now. */
  readonly nodes: ReadonlyMap<string, NodeEntry>;
  /** Moves when a node entry comes or goes. */
  readonly nodeVersion: { readonly value: number };
}

/** A dynamic node: where it came from, its address, its checked (normalized) form and how deep it is nested. */
export interface OverlayNode {
  /** The address of the source that brought it. */
  readonly source: string;
  readonly key: string;
  /** `<source>/(node :key <key>)`. */
  readonly address: string;
  readonly kind: string;
  readonly node: Node;
  /** 0 for a node of a spec source; one more for each source node between. */
  readonly depth: number;
  /** The spec source at the top of its chain. */
  readonly root: string;
  /** The store's entry it was made from (its identity while held). */
  readonly entry: NodeEntry;
}

/** A source node a template needs opened: its entry for the runtime bound to its format (none: unavailable). */
export interface NestedSource {
  readonly address: string;
  readonly entry: BoundSourceEntry;
  readonly runtime?: FormatRuntime;
  /** Its own depth: 1 for a source node of a spec source. */
  readonly depth: number;
  readonly root: string;
  /** The source that brought its node: the item systems it serves are of that source's nodes. */
  readonly parent: string;
}

export interface Overlay {
  /** The store's node version it was made at (0: no node). */
  readonly version: number;
  /** The valid nodes, by depth, then source, then key. */
  readonly nodes: readonly OverlayNode[];
  /** What was left out and why, each once (grows as item systems and nested sources are asked for). */
  readonly issues: readonly Issue[];
  node(address: string): OverlayNode | undefined;
  /** An item's coordinate system (`<ranges node>/(item :id <id>)`), made on first use. */
  system(address: string): SystemEntry | undefined;
  annotation(address: string): AnnotationEntry | undefined;
  /** A path with a dynamic end over the spec's maps and the overlay's; between spec systems, the static model's. */
  path(from: string, to: string): PathResult;
  /**
   * The source nodes to open: those at a depth some template on their root needs (`needs`: root → depth), within
   * `maxNestedDepth`; deeper ones are an issue and not listed.
   */
  nested(needs: ReadonlyMap<string, number>): readonly NestedSource[];
  /**
   * The system a source names `name` (a nested source's coverage, by name): an item system of a node its source's
   * parent brought (`near`, that parent's address), else a spec system of that name.
   */
  systemNamed(name: string, near?: string): SystemEntry | undefined;
  /**
   * The nodes the given addresses (nodes, or item systems: their node) hold, with every node a held node names
   * (its `target`, its `source`), transitively: what is shown holds what it is drawn through.
   */
  holds(addresses: Iterable<string>): readonly OverlayNode[];
}

/** No dynamic node: what a frame without a store's nodes reads. */
export const EmptyOverlay: Overlay = {
  version: 0,
  nodes: [],
  issues: [],
  node: () => undefined,
  system: () => undefined,
  annotation: () => undefined,
  path: () => ({}),
  nested: () => [],
  systemNamed: () => undefined,
  holds: () => [],
};

function segmentsOf(address: string): readonly Segment[] | undefined {
  return parseAddress(address).pattern?.segments;
}

/** A node's address: its source's, then `(node :key <key>)`. */
export function nodeAddress(source: string, key: string): string {
  return printAddress([...(segmentsOf(source) ?? []), { kind: 'node', params: { key } }]);
}

/** An item's coordinate system: its ranges node's address, then `(item :id <id>)`. */
export function itemAddress(annotation: string, id: string | number): string {
  return printAddress([...(segmentsOf(annotation) ?? []), { kind: 'item', params: { id } }]);
}

/** An address as its canonical spelling (forms too), or undefined. */
function canonical(v: Value | undefined): string | undefined {
  if (typeof v !== 'string') return undefined;
  const p = parseAddress(v);
  if (p.issues.length || !p.pattern || p.pattern.subtree || !p.pattern.segments.length) return undefined;
  return printAddress(p.pattern.segments);
}

const checked = new WeakMap<Node, { node?: Node; why?: string }>();

/** A dynamic node checked as a node of its kind under `data`, normalized; kept per node object. */
function check(node: Node): { node?: Node; why?: string } {
  let c = checked.get(node);
  if (c) return c;
  const spec: Node = {
    kind: 'spec',
    params: { version: SeqLanguageVersion },
    children: [{ kind: 'data', id: 'data', children: [{ ...node, id: 'node' }] }],
  };
  const r = validateAndNormalize(SeqLanguage, spec);
  const n = r.node?.children?.[0]?.children?.[0];
  const { id: _, ...rest } = n ?? node;
  c = n ? { node: rest as Node } : { why: r.issues.find((i) => i.severity === 'error')?.message ?? 'invalid' };
  checked.set(node, c);
  return c;
}

export interface OverlayOptions {
  /** The format runtimes nested sources are bound to (the application's registry). */
  readonly runtimes?: readonly FormatRuntime[];
  readonly config?: Partial<ResolveConfig>;
}

const NoRuntimes: readonly FormatRuntime[] = [];

const SourceKeys = ['uri', 'ref', 'value', 'index', 'fields', 'assembly'] as const;

/** The overlay of the store's nodes now, over a static model. */
export function buildOverlay(model: StaticModel, reader: NodeReader, options: OverlayOptions = {}): Overlay {
  return new DynamicOverlay(model, reader, options);
}

/** Keeps the last overlay while the model and the store's node version hold (one overlay a node version). */
export class OverlayBuilder {
  private last: { model: StaticModel; version: number; overlay: Overlay } | undefined;

  constructor(readonly options: OverlayOptions = {}) {}

  of(model: StaticModel, reader: NodeReader): Overlay {
    const version = reader.nodeVersion.value;
    if (this.last && this.last.model === model && this.last.version === version) return this.last.overlay;
    const overlay = reader.nodes.size ? buildOverlay(model, reader, this.options) : { ...EmptyOverlay, version };
    this.last = { model, version, overlay };
    return overlay;
  }
}

class DynamicOverlay implements Overlay {
  readonly version: number;
  readonly nodes: OverlayNode[] = [];
  readonly issues: Issue[] = [];
  private readonly config: ResolveConfig;
  private readonly byAddress = new Map<string, OverlayNode>();
  private readonly annotations = new Map<string, AnnotationEntry>();
  /** Ranges nodes whose items are coding sequences, and their item systems and maps once made. */
  private readonly coding = new Map<string, OverlayNode>();
  private readonly items = new Map<string, { systems: Map<string, SystemEntry>; maps: MapEntry[] }>();
  private edges: Map<string, Edge[]> | undefined;
  private readonly paths = new Map<string, PathResult>();
  private readonly reported = new Set<string>();

  constructor(
    private readonly model: StaticModel,
    reader: NodeReader,
    private readonly options: OverlayOptions,
  ) {
    this.version = reader.nodeVersion.value;
    this.config = { ...DefaultResolveConfig, ...model.config, ...options.config };
    // parents first: a nested source's nodes are under a longer source address
    const entries = [...reader.nodes.values()].sort(
      (a, b) => a.source.length - b.source.length || cmp(a.source, b.source) || cmp(a.key, b.key),
    );
    for (const entry of entries) {
      const parent = this.byAddress.get(entry.source);
      const place = model.sources.has(entry.source)
        ? { depth: 0, root: entry.source }
        : parent?.kind === 'source'
          ? { depth: parent.depth + 1, root: parent.root }
          : undefined;
      // a node whose source is no longer held waits for it (its parent left first)
      if (!place) continue;
      const address = nodeAddress(entry.source, entry.key);
      const c = check(entry.node);
      if (!c.node) {
        this.report('seq.resolve.node.invalid', `The node is not valid: ${c.why}`, address);
        continue;
      }
      const n: OverlayNode = {
        source: entry.source,
        key: entry.key,
        address,
        kind: c.node.kind,
        node: c.node,
        ...place,
        entry,
      };
      this.byAddress.set(address, n);
      if (n.kind === 'ranges' && n.node.params?.map === 'codon') this.coding.set(address, n);
    }
    // in depth order, coding nodes first: a node may be on an item of a shallower one (or one of its own depth)
    const order = (n: OverlayNode) => n.depth * 2 + (this.coding.has(n.address) ? 0 : 1);
    for (const n of [...this.byAddress.values()].sort((a, b) => order(a) - order(b))) {
      const fails = !this.rules(n) || ((n.kind === 'ranges' || n.kind === 'values') && !this.resolveAnnotation(n));
      if (fails) {
        this.byAddress.delete(n.address);
        this.coding.delete(n.address);
        continue;
      }
      this.nodes.push(n);
    }
    this.nodes.sort((a, b) => a.depth - b.depth || cmp(a.source, b.source) || cmp(a.key, b.key));
  }

  /** What a ref of a node names, as the language's checks read it: a system (its alphabet, cyclic flag) or a source. */
  private targetOf(n: OverlayNode, key: string): { node: Node; length?: number } | undefined {
    const a = canonical(n.node.params?.[key]);
    if (a === undefined) return undefined;
    const sys = this.model.systems.get(a) ?? this.system(a);
    if (sys)
      return { node: { kind: 'sequence', params: { alphabet: sys.alphabet, cyclic: sys.cyclic } }, length: sys.length };
    if (this.model.sources.has(a) || this.byAddress.get(a)?.kind === 'source') return { node: { kind: 'source' } };
    return undefined;
  }

  /** The language's own rules of a data node (`checkDataNode`): an error leaves the node out, with its issue. */
  private rules(n: OverlayNode): boolean {
    let ok = true;
    checkDataNode(
      n.node,
      { target: (key) => this.targetOf(n, key)?.node, length: (key) => this.targetOf(n, key)?.length },
      (code, message, ...tail) => {
        const issue = issueOf(LangIssues, code, message, { address: n.address, path: ['params', ...tail] });
        if (issue.severity === 'error') ok = false;
        this.push(issue);
      },
    );
    return ok;
  }

  private push(issue: Issue): void {
    const id = `${issue.code}|${issue.address}`;
    if (this.reported.has(id)) return;
    this.reported.add(id);
    this.issues.push(issue);
  }

  private report(code: keyof typeof ResolveIssues.codes, message: string, address: string): void {
    const id = `${code}|${address}`;
    if (this.reported.has(id)) return;
    this.reported.add(id);
    const hint = (ResolveIssues.codes[code] as { hint?: string }).hint;
    this.issues.push(issueOf(ResolveIssues, code, message, { address, ...(hint ? { hint } : {}) }));
  }

  private resolveAnnotation(n: OverlayNode): boolean {
    const params = n.node.params ?? {};
    const target = canonical(params.target);
    const system = target === undefined ? undefined : (this.model.systems.get(target) ?? this.system(target));
    if (!system) {
      this.report('seq.resolve.node.ref', 'The node’s target is no system held', n.address);
      return false;
    }
    let source: string | undefined;
    if (params.source !== undefined) {
      source = canonical(params.source);
      if (source === undefined || !(this.model.sources.has(source) || this.byAddress.get(source)?.kind === 'source')) {
        this.report('seq.resolve.node.ref', 'The node’s source is no source held', n.address);
        return false;
      }
    }
    const r = resolveAnnotationNode(n.address, n.node, system, source);
    for (const i of r.issues) this.report(i.code as keyof typeof ResolveIssues.codes, i.message, n.address);
    this.annotations.set(n.address, r.annotation);
    return true;
  }

  node(address: string): OverlayNode | undefined {
    return this.byAddress.get(address);
  }

  annotation(address: string): AnnotationEntry | undefined {
    return this.annotations.get(address);
  }

  system(address: string): SystemEntry | undefined {
    const segments = segmentsOf(address);
    const last = segments?.[segments.length - 1];
    if (!segments || typeof last === 'string' || last?.kind !== 'item') return undefined;
    const owner = printAddress(segments.slice(0, -1));
    const n = this.coding.get(owner);
    return n ? this.itemsOf(n)?.systems.get(address) : undefined;
  }

  /** A coding ranges node's item systems and their codon maps, made on first use (EC-177). */
  private itemsOf(n: OverlayNode) {
    const made = this.items.get(n.address);
    if (made) return made;
    const target = this.annotationTarget(n);
    if (!target) return undefined;
    const rows = n.node.params?.rows as Record<string, Value> | undefined;
    const out = { systems: new Map<string, SystemEntry>(), maps: [] as MapEntry[] };
    this.items.set(n.address, out);
    if (!rows) return out;
    const start = rows.start as Int32Array;
    const end = rows.end as Int32Array;
    const id = rows.id instanceof CategoryColumn ? rows.id : undefined;
    const strand = rows.strand instanceof CategoryColumn ? rows.strand : undefined;
    const pieces = new Map<string | number, { exons: [number, number][]; strand: string }>();
    for (let i = 0; i < start.length; i++) {
      const key = id?.get(i) ?? i;
      let p = pieces.get(key);
      const here = strand?.get(i) === '-' ? '-' : '+';
      if (!p) pieces.set(key, (p = { exons: [], strand: here }));
      else if (p.strand !== here)
        this.report(
          'seq.resolve.item.strand',
          'An item’s rows are on both strands: the first row’s is taken',
          n.address,
        );
      p.exons.push([start[i], end[i]]);
    }
    for (const [key, p] of pieces) {
      p.exons.sort((a, b) => a[0] - b[0]);
      let coding = 0;
      for (const [a, b] of p.exons) coding += b - a;
      const address = itemAddress(n.address, key);
      const system: SystemEntry = {
        address,
        kind: 'sequence',
        alphabet: 'protein',
        length: Math.floor(coding / 3),
        cyclic: false,
        name: String(key),
        numberings: ['index', 'ordinal'],
        // its letters would come from the source that brought it, under the item's name
        access: { kind: 'source', source: n.source, system: String(key), shape: 'letters' },
      };
      const map: Node = {
        kind: 'map',
        params: {
          method: 'codon',
          strand: p.strand,
          exons: { start: Int32Array.from(p.exons, (e) => e[0]), end: Int32Array.from(p.exons, (e) => e[1]) },
        },
      };
      const r = resolveMapNode(address, map, system, target, this.config);
      for (const i of r.issues) this.report(i.code as keyof typeof ResolveIssues.codes, i.message, address);
      out.systems.set(address, system);
      out.maps.push(r.map);
    }
    return out;
  }

  private annotationTarget(n: OverlayNode): SystemEntry | undefined {
    const t = this.annotations.get(n.address)?.target;
    return t === undefined ? undefined : (this.model.systems.get(t) ?? this.system(t));
  }

  /** The overlay's steps: every item map both ways, made when a path is first asked for. */
  private dynamicEdges(): Map<string, Edge[]> {
    if (this.edges) return this.edges;
    const edges = new Map<string, Edge[]>();
    const add = (e: Edge) => {
      let list = edges.get(e.from);
      if (!list) edges.set(e.from, (list = []));
      list.push(e);
    };
    for (const n of this.coding.values())
      for (const m of this.itemsOf(n)?.maps ?? []) {
        add({ map: m.address, from: m.from, to: m.to, backward: false, step: m.step, linear: m.linear });
        add({
          map: m.address,
          from: m.to,
          to: m.from,
          backward: true,
          step: -m.step as Step,
          linear: m.linear.inverse(),
        });
      }
    this.edges = edges;
    return edges;
  }

  path(from: string, to: string): PathResult {
    if (this.model.systems.has(from) && this.model.systems.has(to)) return this.model.path(from, to);
    const k = `${from}\u0000${to}`;
    let r = this.paths.get(k);
    if (r) return r;
    const dynamic = this.dynamicEdges();
    const found = findPath({ edgesFrom: (s) => [...this.model.edgesFrom(s), ...(dynamic.get(s) ?? [])] }, from, to, {
      maxSteps: this.config.maxPathSteps,
      maxCandidates: this.config.maxPathCandidates,
    });
    r = found.ok ? { path: found.path } : { issue: pathIssue(found.why, from, to) };
    this.paths.set(k, r);
    return r;
  }

  systemNamed(name: string, near?: string): SystemEntry | undefined {
    for (const n of this.coding.values()) {
      if (near !== undefined && n.source !== near) continue;
      for (const sys of this.itemsOf(n)?.systems.values() ?? []) if (sys.name === name) return sys;
    }
    for (const sys of this.model.systems.values()) if (sys.name === name) return sys;
    return undefined;
  }

  /** The node an address is or is an item of. */
  private ownerOf(address: string): OverlayNode | undefined {
    const n = this.byAddress.get(address);
    if (n) return n;
    const segments = segmentsOf(address);
    const last = segments?.[segments.length - 1];
    if (!segments || typeof last === 'string' || last?.kind !== 'item') return undefined;
    return this.byAddress.get(printAddress(segments.slice(0, -1)));
  }

  holds(addresses: Iterable<string>): readonly OverlayNode[] {
    const out = new Set<OverlayNode>();
    const todo: OverlayNode[] = [];
    const add = (a: string | undefined) => {
      const n = a === undefined ? undefined : this.ownerOf(a);
      if (n && !out.has(n)) {
        out.add(n);
        todo.push(n);
      }
    };
    for (const a of addresses) add(a);
    while (todo.length) {
      const n = todo.pop()!;
      for (const key of ['target', 'source']) add(canonical(n.node.params?.[key]));
    }
    return [...out];
  }

  nested(needs: ReadonlyMap<string, number>): readonly NestedSource[] {
    const out: NestedSource[] = [];
    for (const n of this.nodes) {
      if (n.kind !== 'source') continue;
      const depth = n.depth + 1;
      if ((needs.get(n.root) ?? -1) < depth) continue;
      if (depth > this.config.maxNestedDepth) {
        this.report('seq.resolve.node.depth', `The source is nested ${depth} levels deep`, n.address);
        continue;
      }
      const params = n.node.params ?? {};
      const format = typeof params.format === 'string' ? params.format : '';
      const entry: Record<string, unknown> = { address: n.address, format };
      for (const key of SourceKeys) if (params[key] !== undefined) entry[key] = params[key];
      // one runtime a format, never the first listed (_spec/README.md "Context", Registry entries)
      const found = runtimeFor(this.options.runtimes ?? NoRuntimes, format);
      const runtime = found === 'conflict' ? undefined : found;
      if (!runtime) this.report('seq.resolve.source.runtime', 'No runtime reads the format of the source', n.address);
      out.push({
        address: n.address,
        entry: entry as unknown as BoundSourceEntry,
        ...(runtime ? { runtime } : {}),
        depth,
        root: n.root,
        parent: n.source,
      });
    }
    return out;
  }
}

function cmp(a: string, b: string): number {
  return a < b ? -1 : a > b ? 1 : 0;
}
