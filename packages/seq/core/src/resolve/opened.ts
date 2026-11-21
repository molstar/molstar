/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The opened model (_spec/README.md "Overlay", Opened sources): the static model with the molecules, complexes and
// units its opened sources hold read as its own, so a section `of` a source's complex resolves (the static-model rule's
// one exception). One a model and overlay (kept): an overlay is one a node version, so a source's update is one new
// opened model, never one a frame or a hover. The static model itself while the overlay holds no node.

import type { StaticModel } from './model.js';
import { type Molecular, modifiedOf } from './molecule-query.js';
import type { Overlay } from './overlay.js';
import { overlayQueryModel } from './query-model.js';

/** Two maps read as one, the first's entries first (an overlay's addresses are never the spec's). */
class Chained<K, V> implements ReadonlyMap<K, V> {
  private all: Map<K, V> | undefined;
  constructor(
    private readonly a: ReadonlyMap<K, V>,
    private readonly b: ReadonlyMap<K, V>,
  ) {}
  get size() {
    return this.a.size + this.b.size;
  }
  get(k: K) {
    return this.a.get(k) ?? this.b.get(k);
  }
  has(k: K) {
    return this.a.has(k) || this.b.has(k);
  }
  // a walk is rare (an export's plans): the two copied once
  private both() {
    return (this.all ??= new Map([...this.a, ...this.b]));
  }
  entries() {
    return this.both().entries();
  }
  keys() {
    return this.both().keys();
  }
  values() {
    return this.both().values();
  }
  [Symbol.iterator]() {
    return this.both()[Symbol.iterator]();
  }
  forEach(fn: (v: V, k: K, m: ReadonlyMap<K, V>) => void) {
    for (const [k, v] of this.both()) fn(v, k, this);
  }
}

const chain = <K, V>(a: ReadonlyMap<K, V>, b: ReadonlyMap<K, V> | undefined): ReadonlyMap<K, V> =>
  b?.size ? (a.size ? new Chained(a, b) : b) : a;

const Kept = new WeakMap<StaticModel, WeakMap<Overlay, StaticModel>>();

/**
 * The static model with the molecules, complexes, units and their systems the overlay's sources hold as its own
 * (`systems`, `molecular`), paths over the overlay's steps, its query forms the overlay's query model's, and `opened`
 * the sources holding a node; the model itself when the overlay holds none. Kept per model and overlay.
 */
export function openedModel(model: StaticModel, overlay: Overlay): StaticModel {
  if (!overlay.nodes.length) return model;
  let byOverlay = Kept.get(model);
  if (!byOverlay) Kept.set(model, (byOverlay = new WeakMap()));
  let m = byOverlay.get(overlay);
  if (!m) byOverlay.set(overlay, (m = makeOpened(model, overlay)));
  return m;
}

function makeOpened(model: StaticModel, overlay: Overlay): StaticModel {
  const own = model.molecular;
  const dyn = overlay.molecular;
  const molecular: Molecular = {
    molecules: chain(own.molecules, dyn.molecules),
    complexes: chain(own.complexes, dyn.complexes),
    units: chain(own.units, dyn.units),
    atoms: chain(own.atoms, dyn.atoms),
    modified: chain(own.modified, dyn.modified),
    nodes: chain(own.nodes, dyn.nodes),
    ...(dyn.numberings ? { numberings: dyn.numberings } : {}),
    ...(dyn.compounds ? { compounds: dyn.compounds } : {}),
  };
  const both = (from: string, to: string) => model.systems.has(from) && model.systems.has(to);
  return {
    ...overlayQueryModel(model, overlay),
    spec: model.spec,
    hash: `${model.hash}|opened@${overlay.version}`,
    config: model.config,
    systems: chain(model.systems, dyn.systems),
    maps: model.maps,
    sources: model.sources,
    annotations: model.annotations,
    assemblies: model.assemblies,
    ...(model.loading ? { loading: model.loading } : {}),
    numberings: model.numberings,
    hasSourceNumberings: model.hasSourceNumberings,
    molecular,
    compounds: model.compounds,
    ...(model.describedMaps ? { describedMaps: model.describedMaps } : {}),
    opened: new Set(overlay.nodes.map((n) => n.source)),
    modified: (s) => model.modified(s) ?? modifiedOf(dyn, s),
    compound: (code) => model.compound(code) ?? dyn.compounds?.get(code),
    issues: model.issues,
    numberingOf: (s, id) => model.numberingOf(s, id),
    // between two spec systems the model's own path (as the overlay's): a dynamic end searches the overlay's steps too
    path: (from, to, via) => (both(from, to) ? model.path(from, to, via) : overlay.path(from, to, via)),
    edgesFrom: (s) => overlay.edgesFrom?.(s) ?? model.edgesFrom(s),
    // a Location is carried and closed along the spec's maps: what crosses sections goes by the view's links
    project: (loc, to, options) => model.project(loc, to, options),
    close: (loc, among) => model.close(loc, among),
  };
}
