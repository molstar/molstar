/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b: track templates (seqspec:format#placement). A track with `for_each` stands for one track per dynamic node it
// selects: `templateNodes` selects them from the set (the nodes of a kind under a streamed source, within `depth`
// levels of the sources it nests, in key order), each instance's layers resolve through `resolveLayer` with
// the node for `"$node"`, and `instantiate` puts the instances at the template's place in the section, so
// `composeView` composes them as any track. Instance refs are `<template ref>/<node ref>`, a track's and a layer's.

import {
  Annotations,
  CoordinateSystemKind,
  KeyKind,
  type KeyTarget,
  MoleculeKind,
  type NameId,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { DocIndex } from "../load/compile.js";
import { ItemSystemProp } from "../load/nodes.js";
import type { SvsRegistries } from "../registry.js";
import { isPlainObject } from "../schema/field-schema.js";
import type { IndexedNode } from "../schema/tree-schema.js";
import { type KeyedSay, type KeyFilterValue, keyed, type Members, type OrderKeyValue } from "./keyed.js";
import type { ResolvedSection, TrackRow } from "./section.js";

/** A template track's `for_each`, normalized. */
export interface ForEach {
  readonly kind: "sequence" | "ligand" | "unit" | "mapping" | "annotation" | "item_system";
  readonly source: string;
  readonly depth: number | null;
  readonly reaching: boolean;
  /** Key order by default (`null`, `"data"`); else keys over key tables keyed by `node` (keyed sets, keyed.ts). */
  readonly order: OrderKeyValue | readonly OrderKeyValue[] | null;
  readonly filter?: readonly KeyFilterValue[] | null;
}

/** A track's `for_each`; null when it is no template. */
export const forEachOf = (n: IndexedNode | undefined): ForEach | null =>
  n?.kind === "track" && isPlainObject(n.params.for_each) ? (n.params.for_each as unknown as ForEach) : null;

/** An instance's ref: its template's (a track's or a layer's), then its node's. */
export const instanceRef = (template: string, node: string) => `${template}/${node}`;

const NONE: { has(ref: string): boolean } = new Set<string>();

/** How many of `nested` lie between `source` and the node `ref` (its prefixes at a `/`): its depth under `source`. */
function depthUnder(ref: string, source: string, nested: { has(ref: string): boolean }): number {
  let d = 0;
  for (let i = ref.indexOf("/", source.length + 1); i > 0; i = ref.indexOf("/", i + 1))
    if (nested.has(ref.slice(0, i))) d++;
  return d;
}

/**
 * A template track's selected nodes, in key order: the set's nodes of its `kind` under its `source`, at `depth` levels
 * of nested sources below it (`nested`: the held nested sources' refs; null depth: any). `reaching` is the caller's: an
 * instance whose layers place nothing is left out (`instantiate`'s input).
 */
export function templateNodes(
  track: string,
  doc: DocIndex,
  set: SeqSet,
  nested: { has(ref: string): boolean } = NONE,
  reg?: SvsRegistries,
  say?: KeyedSay,
): readonly string[] {
  const fe = forEachOf(doc.nodes.get(track));
  if (!fe) return [];
  const prefix = `${fe.source}/`;
  const out: string[] = [];
  const see = (key: NameId) => {
    const k = set.strings.get(key);
    if (k.startsWith(prefix) && (fe.depth === null || depthUnder(k, fe.source, nested) === fe.depth)) out.push(k);
  };
  // Item coordinate systems (ST7) only by `item_system`: never as the source's sequences, units or mappings.
  const item = (mol: number) => set.molecules[mol]?.custom.props.get(ItemSystemProp) === true;
  const K = CoordinateSystemKind;
  if (fe.kind === "annotation") for (const a of set.annotations) a && see(a.key);
  else if (fe.kind === "mapping")
    for (const m of set.mappings) m && !(m.from.kind === K.molecule && item(m.from.id)) && see(m.key);
  else if (fe.kind === "unit")
    for (let u = 0; u < set.units.count; u++)
      set.units.alive[u] && !item(set.units.molecule[u]) && see(set.units.key[u] as NameId);
  else if (fe.kind === "item_system")
    for (const m of set.molecules) m && m.custom.props.get(ItemSystemProp) === true && see(m.key);
  else {
    const kind = fe.kind === "ligand" ? MoleculeKind.ligand : MoleculeKind.polymer;
    for (const m of set.molecules) if (m?.kind === kind && m.custom.props.get(ItemSystemProp) !== true) see(m.key);
  }
  out.sort();
  const keyedOrder = fe.order != null && fe.order !== "data" && fe.order !== "tree";
  if (!reg || !(keyedOrder || fe.filter?.length)) return out;
  // Key tables keyed by `node` (a `rows` annotation's `node` column) order and filter the instances
  // (seqspec:format#layer-types).
  const at = new Map(out.map((k, i) => [k, i]));
  const members: Members = {
    count: out.length,
    tree: null,
    // A row with no node names none (GV-11: `String(undefined)` named a node called "undefined").
    of: (ann, i) => {
      const v = Annotations.rowValue(ann, set, i, "node");
      return v === undefined ? -1 : (at.get(String(v)) ?? -1);
    },
  };
  const k = keyed(set, reg, members, { order: fe.order, filter: fe.filter }, say ?? (() => {}));
  return Array.from(k.permutation, (m) => out[m]).filter((_, i) => k.kept[k.permutation[i]]);
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
/** An instance's label: its node's (an annotation's, a molecule's, a unit's), else its key under the template's source. */
function labelOf(set: SeqSet, node: string, source: string): number {
  let label: number = -1;
  if (set.byName(node, t)) {
    if (t.kind === KeyKind.annotation) label = set.annotations[t.id]?.label ?? -1;
    else if (t.kind === KeyKind.molecule) label = set.molecules[t.id]?.label ?? -1;
    else if (t.kind === KeyKind.unit) label = set.units.label[t.id];
  }
  // A node given no label is labeled by its key (the model's default): its key under the source reads better.
  return label >= 0 && set.strings.get(label as NameId) !== node
    ? label
    : set.strings.intern(node.slice(source.length + 1));
}

/**
 * The section with each template track (by its index in `section.tracks`) replaced by its instances' tracks, at its
 * place and in the given order, and their layers (`instanceRef` of each template layer) in the section's layer order;
 * a template with no instance shows no track. Split layers are matched again by ref.
 */
export function instantiate(
  section: ResolvedSection,
  instances: ReadonlyMap<number, readonly string[]>,
  doc: DocIndex,
  set: SeqSet,
): ResolvedSection {
  if (instances.size === 0) return section;
  const tracks: TrackRow[] = [];
  const layers: { ref: string; track: number; type: string }[] = [];
  const at = new Map<string, number>();
  const split: string[] = [];
  section.tracks.forEach((tr, k) => {
    const own = section.layers.filter((l) => l.track === k);
    const splitRef = tr.split >= 0 ? section.layers[tr.split]?.ref : undefined;
    const nodes = instances.get(k);
    const ref = set.strings.get(tr.ref);
    const fe = nodes ? forEachOf(doc.nodes.get(ref)) : null;
    for (const node of nodes ?? [null]) {
      const i = tracks.length;
      const label = node === null || tr.params.label != null ? tr.label : labelOf(set, node, fe?.source ?? "");
      tracks.push(node === null ? tr : { ...tr, ref: set.strings.intern(instanceRef(ref, node)), label });
      for (const l of own) {
        const r = node === null ? l.ref : instanceRef(l.ref, node);
        at.set(r, layers.push({ ref: r, track: i, type: l.type }) - 1);
      }
      split.push(splitRef === undefined ? "" : node === null ? splitRef : instanceRef(splitRef, node));
    }
  });
  for (let i = 0; i < tracks.length; i++)
    if (tracks[i].split >= 0) tracks[i] = { ...tracks[i], split: at.get(split[i]) ?? -1 };
  return { ...section, tracks, layers };
}
