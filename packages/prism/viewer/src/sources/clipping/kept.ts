/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { State, StateBuilder, StateObjectCell } from '@molstar/core/state';
import { StateSelection } from '@molstar/core/state';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { lociLabel } from '@molstar/model/model/label';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { ClippingStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/clipping';

// Atoms kept whole by a cut (clipping.md §5): Mol*'s clipping layer on a structure representation (the cell the
// clipping helpers and Prism's `representation.clipping` write, tagged `clipping-controls`), its group bit i meaning
// "not clipped by the representation's clip object i". The layer is the record: a cut's kept atoms are read from it,
// and every change writes it back. An atom's flags are those of the last layer covering it (`Clipping.merge`), so a
// bit is set or cleared on the atoms' current flags, never by a layer over them.

/** The tag Mol*'s clipping helpers give the layer cell (`structure-clipping.ts`, private there) */
export const ClippingLayerTag = 'clipping-controls';

/** Clip objects past this index have no group (Mol*'s groups one to six): their cuts cannot keep atoms whole */
export const MaxKeptIndex = 5;

/** The structure a representation draws (its root keeps the layer valid while the structure changes), or undefined */
export function structureOf(repr: StateObjectCell): Structure | undefined {
  const data = (repr.obj?.data as { sourceData?: unknown } | undefined)?.sourceData;
  return data instanceof Structure ? data : undefined;
}

export function layerCell(state: State, reprRef: string): StateObjectCell | undefined {
  return state.select(
    StateSelection.Generators.ofTransformer(ClippingStructureRepresentation3DFromBundle, reprRef).withTag(
      ClippingLayerTag,
    ),
  )[0];
}

/** The representation's layer, merged (an atom in one layer at most), on its structure's root; empty without one */
export function layerOf(state: State, repr: StateObjectCell): Clipping {
  const structure = structureOf(repr);
  const cell = layerCell(state, repr.transform.ref);
  const layers = (cell?.transform.params as { layers?: Clipping.BundleLayer[] } | undefined)?.layers;
  if (!structure || !layers?.length) return Clipping.Empty;
  return Clipping.merge(Clipping.ofBundle(layers, structure.root));
}

/** The atoms whose flags have `bit` (kept whole by the clip object of that bit), or undefined */
export function keptLoci(layer: Clipping, bit: number): StructureElement.Loci | undefined {
  let out: StructureElement.Loci | undefined;
  for (const l of layer.layers) {
    if (!(l.groups & bit)) continue;
    const loci = l.loci as StructureElement.Loci;
    out = out ? StructureElement.Loci.union(out, loci) : loci;
  }
  return out && !StructureElement.Loci.isEmpty(out) ? out : undefined;
}

/** `loci` (on the layer's structure, its root) with `bit` added to their flags, the others' kept */
export function withKept(layer: Clipping, loci: StructureElement.Loci, bit: number): Clipping {
  const layers: Clipping.Layer[] = [];
  let rest = loci;
  for (const l of layer.layers) {
    const at = l.loci as StructureElement.Loci;
    const both = StructureElement.Loci.intersect(at, loci);
    const only = StructureElement.Loci.subtract(at, loci);
    if (!StructureElement.Loci.isEmpty(only)) layers.push({ loci: only, groups: l.groups });
    if (!StructureElement.Loci.isEmpty(both)) layers.push({ loci: both, groups: (l.groups | bit) as Clipping.Groups });
    rest = StructureElement.Loci.subtract(rest, at);
  }
  if (!StructureElement.Loci.isEmpty(rest)) layers.push({ loci: rest, groups: bit as Clipping.Groups });
  return Clipping.merge({ kind: 'element-loci', layers });
}

/** Every flag's bits mapped by `map`; atoms left with none drop out (clipped by everything again) */
export function mapGroups(layer: Clipping, map: (groups: number) => number): Clipping {
  const layers: Clipping.Layer[] = [];
  for (const l of layer.layers) {
    const groups = map(l.groups);
    if (groups) layers.push({ loci: l.loci, groups: groups as Clipping.Groups });
  }
  return Clipping.merge({ kind: layer.kind, layers } as Clipping);
}

/** `bit` cleared everywhere: the cut keeps nothing whole */
export const withoutKept = (layer: Clipping, bit: number) => mapGroups(layer, (g) => g & ~bit);

/**
 * The flags after clip object `index` is removed: its bit dropped, the bits of the objects after it one lower (they
 * moved up a place). Without this an exclusion meant for one cut would apply to the next.
 */
export const afterRemoval = (layer: Clipping, index: number) =>
  mapGroups(layer, (g) => (g & ((1 << index) - 1)) | ((g >> (index + 1)) << index));

/** The layer written to the representation in `b`: its cell updated, made (tagged), or deleted when empty */
export function writeLayer(state: State, b: StateBuilder.Root, repr: StateObjectCell, layer: Clipping) {
  const structure = structureOf(repr);
  const cell = layerCell(state, repr.transform.ref);
  if (!structure || Clipping.isEmpty(layer)) {
    if (cell) b.delete(cell.transform.ref);
    return;
  }
  const filtered = Clipping.filter(layer, structure) as Clipping<StructureElement.Loci>;
  const params = Clipping.toBundle(filtered) as never;
  if (cell) b.to(cell.transform.ref).update(params);
  else b.to(repr.transform.ref).apply(ClippingStructureRepresentation3DFromBundle, params, { tags: ClippingLayerTag });
}

/** "REA 200", "12 residues": what is kept, short */
export function keptLabel(loci: StructureElement.Loci): string {
  return lociLabel(loci, { htmlStyling: false, condensed: true });
}
