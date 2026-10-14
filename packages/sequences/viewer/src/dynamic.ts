/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-a: the dynamic nodes the stream manager holds (seqspec:sources#the-viewers-side; seqspec:sources#dynamic-nodes),
// beside its regions. A patch's nodes are held by (source, key): the same key again with equal params is a no-op, with
// other params a conflict (`svs.stream.conflict`, the held node kept). A node stays while something holds it: a region
// whose patch carried it, a shown placement, or a held node that names it; when the last lets go it is evicted, and the
// nodes it named may go with it. `#assemble.stream` reads the entries and their version; nothing here reaches the
// document. Arrival-rate.

import {
  canonical,
  type DocNode,
  nodeRef,
  nodeRefs,
  type StreamNode,
  type SvsRegistries,
} from "@molstar/seq-core/format";

/** A dynamic node held, and what holds it. */
export interface HeldNode {
  readonly source: string;
  readonly ref: string;
  readonly node: StreamNode;
  /** Its params' canonical form: a later patch's same key compares by it. */
  readonly canon: string;
  /** The held nodes it names, each held by it (resolved once its patch is in). */
  names: readonly string[];
  /** Counted in the cache (`Stream.CacheBytes`), about its JSON's size. */
  readonly bytes: number;
  /** Regions whose patches carried it, held nodes naming it, and a shown placement (the view's data or blocks). */
  carriers: number;
  referrers: number;
  shown: boolean;
  /** The viewer's own (an item's coordinate system and its `.map`, ST7): read only when something names the item. */
  readonly own: boolean;
}

/** What `#assemble.stream` reads (graph.ts). */
export interface DynamicSource {
  readonly version: number;
  readonly entries: ReadonlyMap<string, HeldNode>;
}

/** The refs a held node names. */
export type NodeRefsOf = (ref: string, node: StreamNode) => Iterable<string>;
/** The kinds a provider's patch may hold (`StreamNode`); the viewer's own (an item coordinate system) name nothing. */
const FORMAT_KINDS = new Set<string>(["sequence", "ligand", "unit", "mapping", "annotation", "source"]);
/** Core's `nodeRefs` over the format registries: a held node's typed params, a table's targets in any form. */
export const heldNodeRefs =
  (svs: SvsRegistries): NodeRefsOf =>
  (ref, n) => {
    if (!FORMAT_KINDS.has(n.kind)) return [];
    const node = { kind: n.kind, ref, params: n.params } as DocNode;
    return nodeRefs({ node, ref, kind: n.kind, params: n.params, parent: null, path: "" }, svs).map((u) => u.ref);
  };

export class DynamicNodes implements DynamicSource {
  readonly entries = new Map<string, HeldNode>();
  version = 0;
  bytes = 0;
  private readonly changed: () => void;
  private readonly report: (code: string, message: string) => void;
  private readonly refs: NodeRefsOf;

  constructor(changed: () => void, report: (code: string, message: string) => void, refs: NodeRefsOf) {
    this.changed = changed;
    this.report = report;
    this.refs = refs;
  }

  /**
   * A patch's `nodes` from `source`, carried by a region: held once each (a conflict reported, the held one kept). The
   * refs the region now carries (it lets go of them when it leaves, `release`); null: none.
   */
  add(source: string, nodes: readonly StreamNode[], own = false): string[] | null {
    let carried: string[] | null = null;
    const fresh: HeldNode[] = [];
    for (const node of nodes) {
      if (typeof node?.key !== "string") continue;
      // `#` kinds are the viewer's own (`own`: an item's coordinate system), never a provider's.
      if (!own && String(node.kind).startsWith("#")) {
        this.report("svs.node.kind", `${source}: "${node.kind}" is no dynamic node's kind a patch may add`);
        continue;
      }
      const ref = nodeRef(source, node.key);
      const held = this.entries.get(ref);
      const canon = held?.node === node ? held.canon : canonical([node.kind, node.params ?? {}]);
      if (held && held.canon !== canon) {
        this.report(
          "svs.stream.conflict",
          `${source}: "${node.key}" came again with other params; the held one is kept`,
        );
        continue;
      }
      if (held) held.carriers++;
      else {
        const bytes = 2 * canon.length;
        const e: HeldNode = {
          source,
          ref,
          node,
          canon,
          names: [],
          bytes,
          carriers: 1,
          referrers: 0,
          shown: false,
          own,
        };
        this.entries.set(ref, e);
        this.bytes += e.bytes;
        fresh.push(e);
      }
      if (carried === null) carried = [];
      carried.push(ref);
    }
    // Names resolved once the patch is in: a node may name another of the same patch.
    const names = new Set<string>();
    for (const e of fresh) {
      names.clear();
      for (const r of this.refs(e.ref, e.node)) names.add(r);
      const held: string[] = [];
      for (const n of names) {
        const m = n !== e.ref ? this.entries.get(n) : undefined;
        if (m === undefined) continue;
        m.referrers++;
        held.push(n);
      }
      e.names = held;
    }
    if (fresh.length > 0) this.bump();
    return carried;
  }

  /** A region left (evicted, cancelled, closed): the nodes it carried that nothing else holds are evicted. */
  release(carried: readonly string[]) {
    const v = this.version;
    for (const ref of carried) {
      const e = this.entries.get(ref);
      if (e === undefined) continue;
      e.carriers--;
      this.drop(e);
    }
    if (this.version !== v) this.changed();
  }

  /** The shown view's keys (its layers' data, its blocks' units and molecules): edge-rate. */
  show(keys: ReadonlySet<string>) {
    const v = this.version;
    // An item's `.map` is shown with it (one unit: a shown item keeps its path onto its source).
    for (const e of this.entries.values())
      e.shown = keys.has(e.ref) || (e.own && e.ref.endsWith(".map") && keys.has(e.ref.slice(0, -4)));
    for (const e of [...this.entries.values()]) this.drop(e);
    if (this.version !== v) this.changed();
  }

  /** Every node evicted (the manager closed). */
  clear() {
    if (this.entries.size === 0) return;
    this.entries.clear();
    this.bytes = 0;
    this.bump();
  }

  /** Evicted when nothing holds it; the nodes it named, then, when nothing else holds them. */
  private drop(e: HeldNode) {
    if (e.carriers > 0 || e.referrers > 0 || e.shown || this.entries.get(e.ref) !== e) return;
    this.entries.delete(e.ref);
    this.bytes -= e.bytes;
    this.version++;
    for (const n of e.names) {
      const m = this.entries.get(n);
      if (m === undefined) continue;
      m.referrers--;
      this.drop(m);
    }
  }

  private bump() {
    this.version++;
    this.changed();
  }
}
