/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseAddress, printAddress, type Value } from '@molstar/core/repr/tree';
import type { Overlay, OverlayNode, StaticModel } from '@molstar/seqx-core/resolve';
import { LayerTypes, type ResolvedLayer, type ResolvedSection, type ResolvedTrack, type ResolvedView } from './view.js';

// Track templates (_spec/README.md "Resolved view", Templates): the resolved view holds a template, never its expansion;
// at layout a template is expanded over the overlay into a track per matching node, in order. An instance's layers
// without data draw its node (its annotation on its target), carried to the section along the overlay's path. One
// expansion per template per overlay: an arrival changes the overlay, and the next layout expands again.

const expansions = new WeakMap<Overlay, Map<string, readonly ResolvedTrack[]>>();

/** The depth of nested sources the view's templates need, by their source: what the overlay opens (EC-179). */
export function templateNeeds(view: ResolvedView): Map<string, number> {
  const out = new Map<string, number>();
  for (const s of view.sections)
    for (const t of s.tracks)
      if (t.template) out.set(t.template.source, Math.max(out.get(t.template.source) ?? 0, t.template.depth));
  return out;
}

function atom(v: Value | undefined): string | number | undefined {
  return typeof v === 'string' || typeof v === 'number' ? v : typeof v === 'boolean' ? String(v) : undefined;
}

function compare(a: string | number | undefined, b: string | number | undefined): number {
  if (a === b) return 0;
  if (a === undefined) return 1;
  if (b === undefined) return -1;
  if (typeof a === 'number' && typeof b === 'number') return a - b;
  const x = String(a);
  const y = String(b);
  return x < y ? -1 : x > y ? 1 : 0;
}

/** The node's top ancestor's key (a nested node is named by the node of its spec source it came under). */
function rootKey(overlay: Overlay, n: OverlayNode): string {
  let at = n;
  for (let up = overlay.node(at.source); up; up = overlay.node(up.source)) at = up;
  return at.key;
}

/** The nodes a template selects: of its kind, at its depth under its source, filtered, in its order (then key). */
export function templateNodes(track: ResolvedTrack, overlay: Overlay): OverlayNode[] {
  const t = track.template;
  if (!t) return [];
  const by = t.order.startsWith('-') ? t.order.slice(1) : t.order;
  const sign = t.order.startsWith('-') ? -1 : 1;
  const value = (n: OverlayNode) => (by === 'key' ? n.key : atom(n.node.params?.[by]));
  return overlay.nodes
    .filter(
      (n) =>
        n.root === t.source &&
        n.depth === t.depth &&
        n.kind === t.kind &&
        Object.entries(t.filter).every(([k, v]) => String(atom(n.node.params?.[k])) === v),
    )
    .sort((a, b) => sign * compare(value(a), value(b)) || compare(a.address, b.address));
}

/** A template layer on one node: its annotation on its target, carried to the section; not placed when it cannot be. */
function instanceLayer(
  model: StaticModel,
  overlay: Overlay,
  layer: ResolvedLayer,
  node: OverlayNode,
  track: string,
  section: string,
): ResolvedLayer {
  const address = layer.address ? `${track}/${layer.address.slice(layer.address.lastIndexOf('/') + 1)}` : undefined;
  const { address: _, ...rest } = layer;
  const base = { ...rest, ...(address ? { address } : {}) };
  const a = overlay.annotation(node.address);
  if (!a || LayerTypes[layer.type]?.shape !== a.shape) return { ...base, placed: false };
  const length = (model.systems.get(a.target) ?? overlay.system(a.target))?.length;
  const data = {
    kind: 'address' as const,
    address: node.address,
    system: a.target,
    shape: a.shape,
    access: a.access,
    ...(length !== undefined ? { length } : {}),
  };
  if (a.target === section) return { ...base, data, placed: true };
  const path = overlay.path(a.target, section).path;
  return path ? { ...base, data, path, placed: true } : { ...base, data, placed: false };
}

/**
 * A template track's instances over the overlay: a track per selected node (its address the template's, then
 * `(node :at <node>)`; its label the template's and its node's key), its layers without data on that node. Kept per
 * overlay, so a frame that keeps its overlay keeps its instances.
 */
export function expandTemplate(
  model: StaticModel,
  overlay: Overlay,
  section: ResolvedSection,
  track: ResolvedTrack,
): readonly ResolvedTrack[] {
  if (!track.template) return [];
  let byTrack = expansions.get(overlay);
  if (!byTrack) expansions.set(overlay, (byTrack = new Map()));
  const k = `${section.address}\u0000${track.address}`;
  const hit = byTrack.get(k);
  if (hit) return hit;
  const segments = parseAddress(track.address).pattern?.segments ?? [];
  const out: ResolvedTrack[] = [];
  for (const node of templateNodes(track, overlay)) {
    const address = printAddress([...segments, { kind: 'node', params: { at: node.address } }]);
    const layers = track.layers.map((l) =>
      l.placed || l.data ? l : instanceLayer(model, overlay, l, node, address, section.axis.system),
    );
    const { template: _, ...rest } = track;
    out.push({ ...rest, address, label: `${track.label} ${rootKey(overlay, node)}`, layers });
  }
  byTrack.set(k, out);
  return out;
}
