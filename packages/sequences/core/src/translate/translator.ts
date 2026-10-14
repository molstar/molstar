/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The registry of coordinate systems and the routes between them (seqspec:libraries#translate). Platform-free and
// domain-free: systems are nodes, coordinate maps are edges, a route is the composed map of a shortest path. Routes are
// composed on first use per (from, to) and cached until a system or a map changes; nothing on the translating paths
// allocates after that (the visitors' state lives on the instance, as Mol*'s iterators keep theirs).

import {
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  Extents,
  type MapHit,
} from "@molstar/seq-core/model";

export type SystemRefKind = "molecule" | "unit" | "complex" | "alignment" | "annotation" | "item" | "section";
export type SystemStructureField = "structureRef" | "entity" | "labelAsymId" | "authAsymId" | "instanceId";

/** Plain data, no host types: what a host and the viewer say about a system (seqspec:interaction#focus). */
export interface SystemDescriptor {
  /** Stable across patches and reloads: "seqstar:1brs/A", "molstar:3d-1/A@1". */
  readonly key: string;
  readonly kind: "molecule" | "unit" | "row" | "columns" | "item" | "complex" | "layout" | "external";
  /** Its domain. */
  readonly start: number;
  readonly stop: number;
  readonly refs: Partial<Record<SystemRefKind, string>> & { row?: number };
  /** A complex's: its members' keys, in order. A descriptor with members is a composite system. */
  readonly members: readonly string[] | null;
  /** A unit's `mvs` fields, where known, passed on unread. */
  readonly structure: (Partial<Record<SystemStructureField, string>> & { model?: number }) | null;
  /** Names it is labeled or matched by. */
  readonly numberings: readonly string[];
  /** "seqstar", or the namespace of the host that registered it. */
  readonly owner: string;
}

/** An extra index over one system, built on first use and cached until the system is registered again or goes. */
export interface IndexKind<T> {
  readonly name: string;
  build(system: CoordinateSystem, t: Translator): T;
}

/** At most this many maps in a route (the model's bound on a path's counted steps). */
const MaxSteps = 4;

/** An edge from a node to `to`: the map read that way (`fwd`) and the other way (`rev`); null where a span map cannot. */
type Edge = { readonly to: number; readonly fwd: CoordinateMap | null; readonly rev: CoordinateMap | null };
/** A composed route: `map` read forward, or (`reverse`, a span step taken from its to side) its inverse image. */
type Route = { readonly map: CoordinateMap | null; readonly reverse: boolean; readonly target: CoordinateSystem };
/** `forEachAt`'s visitor and state, on the instance. */
type AtState = { visit: TranslatorVisit<unknown>; c: unknown; target: CoordinateSystem | null; count: number };
export type TranslatorVisit<C> = (i: number, c: C, system: CoordinateSystem) => void;

const NoRoutes: readonly Route[] = [];
const NoEdges: readonly Edge[] = [];
const NoSystems: readonly CoordinateSystem[] = [];
const Wide = [-0x80000000, 0x7fffffff] as const;
const inverses = new WeakMap<CoordinateMap, CoordinateMap>();
const inverse = (m: CoordinateMap) => inverses.get(m) ?? inverses.set(m, Coords.invert(m)).get(m);

const addRun = (s: number, t: number, k: { out: Extents; target: CoordinateSystem }) =>
  Extents.add(k.out, k.target, s, t);
function atRun(s: number, t: number, a: AtState) {
  for (let i = s; i < t; i++, a.count++) a.visit(i, a.c, a.target as CoordinateSystem);
}
function atHit(h: MapHit, a: AtState) {
  a.count++;
  a.visit(h.element, a.c, a.target as CoordinateSystem);
}

export class Translator {
  private readonly keys = new Map<string, CoordinateSystem>();
  private readonly owned = new Map<string, readonly CoordinateMap[]>();
  /** By kind: id → row → node (small integer keys: a double key would be boxed per lookup). */
  private readonly nodes: (Map<number, Map<number, number>> | undefined)[] = [];
  private readonly systems: CoordinateSystem[] = [];
  private readonly descriptors: (SystemDescriptor | null)[] = [];
  private readonly indices: (Map<IndexKind<unknown>, unknown> | null)[] = [];
  private maps: CoordinateMap[] = [];
  private adjacency: (Edge[] | undefined)[] | null = null;
  /** By from node: to node → its routes. */
  private routes: (Map<number, readonly Route[]> | undefined)[] = [];
  /** By node: the systems one map away. */
  private near: (readonly CoordinateSystem[] | undefined)[] = [];
  private nextId = 0;
  private readonly sink = { out: Extents.create(1), target: { kind: 0, id: 0, row: -1 } as CoordinateSystem };
  private readonly at: AtState = { visit: () => {}, c: null, target: null, count: 0 };
  private readonly hit: MapHit = { element: 0, phase: 0 };

  /**
   * Registers a system under its key, with the maps it brings (they replace those it brought before; a function of the
   * system when they need it). The key keeps its system once it has one; else `system` (a Seq* system the viewer
   * names), else a new `external` id.
   */
  register(
    d: SystemDescriptor,
    maps: readonly CoordinateMap[] | ((s: CoordinateSystem) => readonly CoordinateMap[]) = [],
    system?: CoordinateSystem,
  ): CoordinateSystem {
    const s = this.keys.get(d.key) ?? system ?? { kind: CoordinateSystemKind.external, id: this.nextId++, row: -1 };
    this.keys.set(d.key, s);
    const n = this.node(s);
    this.descriptors[n] = d;
    this.indices[n] = null;
    const before = this.owned.get(d.key);
    if (before) this.maps = this.maps.filter((m) => !before.includes(m));
    const own = typeof maps === "function" ? maps(s) : maps;
    this.owned.set(d.key, own);
    for (const m of own) this.addMap(m);
    this.changed();
    return s;
  }

  /** Its descriptor, indices and every map touching it go; the key keeps its system for a later `register`. */
  unregister(key: string): void {
    const s = this.keys.get(key);
    if (!s) return;
    const n = this.node(s);
    this.descriptors[n] = null;
    this.indices[n] = null;
    this.owned.delete(key);
    this.maps = this.maps.filter((m) => !Coords.sameSystem(m.from, s) && !Coords.sameSystem(m.to, s));
    this.changed();
  }

  /** A map between any two systems, registered or not; it goes when either end is unregistered. */
  addMap(map: CoordinateMap): void {
    this.node(map.from);
    this.node(map.to);
    this.maps.push(map);
    this.changed();
  }

  /** The key's system, also while unregistered; null: never registered. */
  system(key: string): CoordinateSystem | null {
    return this.keys.get(key) ?? null;
  }
  descriptor(key: string): SystemDescriptor | null {
    const s = this.keys.get(key);
    return s ? this.descriptorOf(s) : null;
  }
  descriptorOf(system: CoordinateSystem): SystemDescriptor | null {
    const n = this.find(system.kind, system.id, system.row);
    return n < 0 ? null : this.descriptors[n];
  }

  /**
   * Runs in `from` appended to `out` as runs in `to`, through every shortest route (a union; copies included). For a
   * composite `from` each entry names its member; for a composite `to` the routes reach its nearest members, and each
   * run names its member. True when `out` gained a run.
   */
  translate(runs: Extents, from: CoordinateSystem, to: CoordinateSystem, out: Extents): boolean {
    const f = this.find(from.kind, from.id, from.row);
    const t = this.find(to.kind, to.id, to.row);
    if (f < 0 || t < 0) return false;
    const composite = this.descriptors[f]?.members != null;
    const before = out.count;
    const sink = this.sink;
    sink.out = out;
    for (let i = 0; i < runs.count; i++) {
      const src = composite ? this.find(runs.kind[i], runs.id[i], runs.row[i]) : f;
      const routes = src < 0 ? NoRoutes : this.routesOf(src, t);
      const s = runs.start[i];
      const e = runs.stop[i];
      for (let r = 0; r < routes.length; r++) {
        const route = routes[r];
        sink.target = route.target;
        if (route.map === null) Extents.add(out, route.target, s, e);
        else if (route.reverse) Coords.forEachFrom(route.map, s, e, addRun, sink);
        else Coords.forEachTo(route.map, s, e, addRun, sink);
      }
    }
    return out.count > before;
  }

  /** Every index in `to` (or its members) one coordinate reaches, copies included: the hover-rate path. The count. */
  forEachAt<C>(from: CoordinateSystem, index: number, to: CoordinateSystem, visit: TranslatorVisit<C>, c: C): number {
    const f = this.find(from.kind, from.id, from.row);
    const t = this.find(to.kind, to.id, to.row);
    if (f < 0 || t < 0) return 0;
    const routes = this.routesOf(f, t);
    const a = this.at;
    // A visitor may translate again: the outer call's state is restored after it (locals, no array: the hot path).
    const pv = a.visit;
    const pc = a.c;
    const pt = a.target;
    const pn = a.count;
    a.visit = visit as TranslatorVisit<unknown>;
    a.c = c;
    a.count = 0;
    try {
      for (let r = 0; r < routes.length; r++) {
        const route = routes[r];
        a.target = route.target;
        if (route.map === null) atRun(index, index + 1, a);
        else if (route.reverse) Coords.forEachFromAt(route.map, index, this.hit, atHit, a);
        else Coords.forEachTo(route.map, index, index + 1, atRun, a);
      }
      return a.count;
    } finally {
      a.visit = pv;
      a.c = pc;
      a.target = pt;
      a.count = pn;
    }
  }

  /** The systems one map from `system`, either way, each once: where the viewer's external routes start. Cached. */
  neighbors(system: CoordinateSystem): readonly CoordinateSystem[] {
    const n = this.find(system.kind, system.id, system.row);
    if (n < 0) return NoSystems;
    const near = this.near[n];
    if (near) return near;
    const to = new Set(((this.adjacency ?? this.link())[n] ?? NoEdges).map((e) => e.to));
    this.near[n] = Array.from(to, (t) => this.systems[t]);
    return this.near[n];
  }

  /** An extra index over a system, built on first use and cached. */
  index<T>(system: CoordinateSystem, kind: IndexKind<T>): T {
    const n = this.node(system);
    const cache = this.indices[n] ?? new Map();
    this.indices[n] = cache;
    if (!cache.has(kind)) cache.set(kind, kind.build(system, this));
    return cache.get(kind) as T;
  }

  private find(kind: number, id: number, row: number): number {
    return this.nodes[kind]?.get(id)?.get(row) ?? -1;
  }
  private node(s: CoordinateSystem): number {
    const n = this.find(s.kind, s.id, s.row);
    if (n >= 0) return n;
    const ids = this.nodes[s.kind] ?? new Map<number, Map<number, number>>();
    const rows = ids.get(s.id) ?? new Map<number, number>();
    this.nodes[s.kind] = ids;
    ids.set(s.id, rows);
    rows.set(s.row, this.systems.length);
    this.descriptors.push(null);
    this.indices.push(null);
    return this.systems.push(s) - 1;
  }
  private changed() {
    this.adjacency = null;
    this.routes = [];
    this.near = [];
  }

  private routesOf(src: number, dst: number): readonly Route[] {
    const cache = this.routes[src];
    const r = cache?.get(dst);
    return r ?? this.cache(src, dst);
  }
  private cache(src: number, dst: number): readonly Route[] {
    const r = this.search(src, dst);
    const cache = this.routes[src] ?? new Map<number, readonly Route[]>();
    this.routes[src] = cache.set(dst, r);
    return r;
  }

  /** Every simple path of the fewest maps (at most one span map) from `src` to `dst` or, for a composite, a member. */
  private search(src: number, dst: number): readonly Route[] {
    const adj = this.adjacency ?? this.link();
    const members = this.descriptors[dst]?.members;
    const targets = new Set<number>();
    for (const key of members ?? []) {
      const s = this.keys.get(key);
      if (s) targets.add(this.node(s));
    }
    if (!members) targets.add(dst);
    if (targets.has(src)) return [{ map: null, reverse: false, target: this.systems[src] }];
    const found: Edge[][] = [];
    const path: Edge[] = [];
    const walk = (n: number, depth: number, spans: number) => {
      if (path.length === depth) {
        if (targets.has(n)) found.push(path.slice());
        return;
      }
      for (const e of adj[n] ?? NoEdges) {
        const span = ((e.fwd ?? e.rev) as CoordinateMap).span > 1 ? 1 : 0;
        if (spans + span > 1 || e.to === src || path.some((p) => p.to === e.to)) continue;
        path.push(e);
        walk(e.to, depth, spans + span);
        path.pop();
      }
    };
    for (let d = 1; d <= MaxSteps && found.length === 0; d++) walk(src, d, 0);
    return found.length ? found.flatMap((p) => this.compose(p)) : NoRoutes;
  }

  /** A path's maps composed: forward, or, with its span map taken backward, the reversed path read as an inverse. */
  private compose(path: readonly Edge[]): Route[] {
    const reverse = path.some((e) => e.fwd === null);
    const seq = (reverse ? path.map((e) => e.rev).reverse() : path.map((e) => e.fwd)) as CoordinateMap[];
    try {
      let m = seq[0];
      for (let i = 1; i < seq.length; i++) m = Coords.compose(m, seq[i], seq[i].to, this.domain(seq[i].to));
      return [{ map: m, reverse, target: this.systems[path[path.length - 1].to] }];
    } catch {
      return []; // a strand - step before a codon step, which Coords.compose refuses: no route
    }
  }

  private domain(s: CoordinateSystem): readonly [number, number] {
    const d = this.descriptorOf(s);
    return d ? [d.start, d.stop] : Wide;
  }

  private link(): (Edge[] | undefined)[] {
    const adj: (Edge[] | undefined)[] = [];
    for (const m of this.maps) {
      const [a, b] = [this.node(m.from), this.node(m.to)];
      const inv = m.span === 1 ? (inverse(m) as CoordinateMap) : null;
      const [ea, eb] = [adj[a] ?? [], adj[b] ?? []];
      ea.push({ to: b, fwd: m, rev: inv });
      eb.push({ to: a, fwd: inv, rev: m });
      [adj[a], adj[b]] = [ea, eb];
    }
    this.adjacency = adj;
    return adj;
  }
}
