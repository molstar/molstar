/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, type StateTransform, type StateTransformer } from '@molstar/core/state';
import { Color } from '@molstar/core/util/color';
import { Material } from '@molstar/core/util/material';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { ClippingStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/clipping';
import { EmissiveStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/emissive';
import { OverpaintStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { SubstanceStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/substance';
import { TransparencyStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/transparency';
import { summarize } from '../targets/index.js';
import { structureNames } from './naming.js';

// The style layers in words (state.md §1, Mol*'s overpaint, transparency, emissive, substance and clipping cells under
// each representation, the ones the plugin's `state/helpers/structure-*.ts` manage by their tags): what
// `representation.layers` returns, the operations that change one find, and the Selection tab's list of applied
// effects reads. One layer here is one value of a kind on one structure ("red on A:10-17 of 1CBS"): Mol* keeps it per
// representation, merged by value (`Overpaint.merge`), so its atoms are the union over the structure's representations.

export type StyleLayerKind = 'color' | 'transparency' | 'emissive' | 'substance' | 'clipping';

/** Where each kind lives: the transformer of its cells and the tag the plugin's helpers give them. */
export const StyleLayerCells: Readonly<
  Record<StyleLayerKind, { readonly transformer: StateTransformer; readonly tag: string; readonly name: string }>
> = {
  color: { transformer: OverpaintStructureRepresentation3DFromBundle, tag: 'overpaint-controls', name: 'color' },
  transparency: {
    transformer: TransparencyStructureRepresentation3DFromBundle,
    tag: 'transparency-controls',
    name: 'transparency',
  },
  emissive: { transformer: EmissiveStructureRepresentation3DFromBundle, tag: 'emissive-controls', name: 'emissive' },
  substance: { transformer: SubstanceStructureRepresentation3DFromBundle, tag: 'substance-controls', name: 'material' },
  clipping: { transformer: ClippingStructureRepresentation3DFromBundle, tag: 'clipping-controls', name: 'clipping' },
};

export const StyleLayerKinds = Object.keys(StyleLayerCells) as StyleLayerKind[];

/** A layer as Mol* keeps it in a cell's `layers` (the fields of every kind). */
interface BundleLayer {
  readonly bundle: StructureElement.Bundle;
  readonly color?: Color;
  readonly clear?: boolean;
  readonly value?: number;
  readonly material?: Material;
  readonly groups?: Clipping.Groups;
}

/** One applied layer: a value of a kind on (part of) one structure. */
export interface StyleLayerInfo {
  /** `<kind>:<structure ref>:<value>`: stable while its value stays (an edit gives a new id, `styleLayerId`). */
  readonly id: string;
  readonly kind: StyleLayerKind;
  readonly structure: { readonly ref: string; readonly label: string };
  /** The value for people: `#e0523a`, `0.5`, `one, two`, `metalness 0 …`; `reset` for a reset layer. */
  readonly value: string;
  /** A color layer's color (`#rrggbb`); undefined for a reset one. */
  readonly color?: string;
  /** Transparency (0 opaque, 1 hidden) or emissive strength. */
  readonly amount?: number;
  /** A reset layer (color, material): it puts the representation's own back where it is. */
  readonly clear?: boolean;
  /** Where it is, short: its addresses (`A:10-17`, a target) when few, else the counts. */
  readonly where: string;
  /** Where it is, in full: "8 residues (60 atoms) in chain A". */
  readonly line: string;
  readonly atoms: number;
  /** The cells that hold it (one per representation it is on). */
  readonly cells: readonly string[];
}

/** The id of a layer: its kind, structure and value key (`styleLayerKey`). */
export function styleLayerId(kind: StyleLayerKind, structure: string, key: string): string {
  return `${kind}:${structure}:${key}`;
}

/** What tells a layer of a kind apart in its cell: its value (a cell's layer, or the fields of one). */
export function styleLayerKey(kind: StyleLayerKind, fields: object): string {
  const layer = fields as Omit<BundleLayer, 'bundle'>;
  switch (kind) {
    case 'color':
      return layer.clear ? 'reset' : Color.toHexStyle(layer.color ?? Color(0));
    case 'transparency':
    case 'emissive':
      return String(round(layer.value ?? 0));
    case 'substance':
      return layer.clear ? 'reset' : Material.toString(layer.material ?? Material());
    case 'clipping':
      return Clipping.Groups.toNames(layer.groups ?? Clipping.Groups.Flag.None).join(', ') || 'none';
  }
}

const round = (x: number) => Math.round(x * 1000) / 1000;

/** The cells of a kind under a representation, as `state.select` gives them. */
export function styleLayerCellsOf(plugin: PluginContext, kind: StyleLayerKind, representation: string) {
  const { transformer, tag } = StyleLayerCells[kind];
  return plugin.state.data.select(StateSelection.Generators.ofTransformer(transformer, representation).withTag(tag));
}

/** What one layer collects while the cells are read. */
interface Gathering {
  readonly kind: StyleLayerKind;
  readonly structure: StructureRef;
  readonly key: string;
  readonly layer: BundleLayer;
  readonly root: Structure;
  loci: StructureElement.Loci;
  readonly cells: string[];
}

/** The layer cells the last read saw, by ref, and what it made of them: equal cells, the same list. */
interface ReadMemo {
  readonly transforms: ReadonlyMap<string, StateTransform>;
  readonly names: string;
  readonly list: readonly StyleLayerInfo[];
}

/**
 * The applied layers of every loaded structure, in the structures' order, by kind (color, transparency, emissive,
 * material, clipping), then in the order the cells hold them. `memo`: what a previous read made; while no layer cell
 * changed (a camera move, a selection, another part of the tree), it is given back as it was, nothing walked again.
 */
export function readStyleLayers(plugin: PluginContext, memo?: { current?: ReadMemo }): readonly StyleLayerInfo[] {
  const structures = plugin.managers.structure.hierarchy.current.structures;
  const names = structureNames(structures);
  const nameKey = structures.map((s) => names.get(s.cell.transform.ref)).join('\n');
  const transforms = new Map<string, StateTransform>();
  const found: { kind: StyleLayerKind; structure: StructureRef; repr: string; transform: StateTransform }[] = [];
  for (const s of structures) {
    for (const c of s.components) {
      for (const r of c.representations) {
        const repr = r.cell.transform.ref;
        for (const kind of StyleLayerKinds) {
          for (const cell of styleLayerCellsOf(plugin, kind, repr)) {
            transforms.set(cell.transform.ref, cell.transform);
            found.push({ kind, structure: s, repr, transform: cell.transform });
          }
        }
      }
    }
  }
  const last = memo?.current;
  if (last && last.names === nameKey && sameTransforms(last.transforms, transforms)) return last.list;

  const groups = new Map<string, Gathering>();
  for (const { kind, structure, repr, transform } of found) {
    const root = (plugin.state.data.cells.get(repr)?.obj?.data as { sourceData?: Structure } | undefined)?.sourceData
      ?.root;
    if (!root) continue;
    const layers = ((transform.params as { layers?: readonly BundleLayer[] } | undefined)?.layers ??
      []) as BundleLayer[];
    for (const layer of layers) {
      const key = styleLayerKey(kind, layer);
      const id = styleLayerId(kind, structure.cell.transform.ref, key);
      const loci = StructureElement.Bundle.toLoci(layer.bundle, root);
      const g = groups.get(id);
      if (!g) {
        groups.set(id, { kind, structure, key, layer, root, loci, cells: [transform.ref] });
        continue;
      }
      if (loci.structure === g.loci.structure) g.loci = StructureElement.Loci.union(g.loci, loci);
      if (!g.cells.includes(transform.ref)) g.cells.push(transform.ref);
    }
  }

  // Ordered by structure, then kind, then as found.
  const order = new Map(structures.map((s, i) => [s.cell.transform.ref, i]));
  const gathered = [...groups.entries()].sort(
    ([, a], [, b]) =>
      (order.get(a.structure.cell.transform.ref) ?? 0) - (order.get(b.structure.cell.transform.ref) ?? 0) ||
      StyleLayerKinds.indexOf(a.kind) - StyleLayerKinds.indexOf(b.kind),
  );
  const list = gathered.map(([id, g]) => describe(id, g, names.get(g.structure.cell.transform.ref)));
  if (memo) memo.current = { transforms, names: nameKey, list };
  return list;
}

function sameTransforms(a: ReadonlyMap<string, StateTransform>, b: ReadonlyMap<string, StateTransform>): boolean {
  if (a.size !== b.size) return false;
  for (const [ref, t] of b) if (a.get(ref) !== t) return false;
  return true;
}

/** Few addresses read as a target (`A:10-17, B:5`); more as the counts. */
const MaxAddresses = 3;

function describe(id: string, g: Gathering, label: string | undefined): StyleLayerInfo {
  const selected = summarize(g.loci, { maxRangesPerChain: MaxAddresses + 1 });
  const where =
    selected.addresses.length > 0 &&
    selected.addresses.length <= MaxAddresses &&
    selected.byChain.every((c) => c.more === 0)
      ? selected.addresses.join(', ')
      : `${selected.residues.toLocaleString('en-US')} residue${selected.residues === 1 ? '' : 's'}${inChains(selected.chains)}`;
  const { layer, kind } = g;
  const base = {
    id,
    kind,
    structure: { ref: g.structure.cell.transform.ref, label: label ?? g.structure.cell.transform.ref },
    value: g.key,
    where,
    line: selected.line,
    atoms: selected.atoms,
    cells: g.cells,
  };
  switch (kind) {
    case 'color':
      return layer.clear ? { ...base, clear: true } : { ...base, color: g.key };
    case 'transparency':
    case 'emissive':
      return { ...base, amount: layer.value ?? 0 };
    case 'substance':
      return layer.clear ? { ...base, clear: true } : base;
    case 'clipping':
      return base;
  }
}

function inChains(chains: readonly string[]): string {
  if (!chains.length) return '';
  if (chains.length === 1) return ` in chain ${chains[0]}`;
  if (chains.length > 4) return ` in ${chains.length} chains`;
  return ` in chains ${chains.join(', ')}`;
}

/** Equal lists of layers (for a state that emits only on change). */
export function sameStyleLayers(a: readonly StyleLayerInfo[], b: readonly StyleLayerInfo[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (x.id !== y.id || x.where !== y.where || x.atoms !== y.atoms || x.structure.label !== y.structure.label) {
      return false;
    }
  }
  return true;
}
