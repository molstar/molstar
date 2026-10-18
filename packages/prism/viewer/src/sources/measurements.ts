/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { stripTags } from '@molstar/core/util/string';
import {
  angleLabel,
  dihedralLabel,
  distanceLabel,
  lociLabel,
  structureElementLociLabelMany,
} from '@molstar/model/model/label';
import type { Loci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureMeasurementCell } from '@molstar/plugin/state/manager/structure/measurement';
import type { IconName } from '@molstar/prism-ui/core';

// The measurements in words (state.md §1, `managers.structure.measurement`): what the scene tree, the `measurements`
// feature's operations and its Measure section read.

export type MeasureKind = 'distance' | 'angle' | 'dihedral' | 'label' | 'orientation' | 'plane';

/**
 * What each kind is called and drawn as, once: the scene tree's nodes, the Measure section's buttons and the
 * operations' texts read it. `points` is how many picks it takes (0: the selection as a whole).
 */
export const MeasureKinds: Readonly<Record<MeasureKind, { title: string; icon: IconName; points: number }>> = {
  distance: { title: 'Distance', icon: 'measure-distance', points: 2 },
  angle: { title: 'Angle', icon: 'measure-angle', points: 3 },
  dihedral: { title: 'Dihedral', icon: 'measure-angle', points: 4 },
  label: { title: 'Label', icon: 'tag', points: 1 },
  orientation: { title: 'Orientation', icon: 'move-3d', points: 0 },
  plane: { title: 'Plane', icon: 'square', points: 0 },
};

/** One measurement, for `measure.get`, results and the agent's account. */
export interface MeasurementInfo {
  /** The state ref of its representation. */
  readonly ref: string;
  readonly kind: MeasureKind;
  /** What it measures between: "A 10 CA – A 20 CA", "12 Residues". */
  readonly label: string;
  /** What it reads: "12.31 Å", "104.20°", a label's text; empty for an orientation and a plane. */
  readonly value: string;
  /** Each atom or residue measured, in order. */
  readonly atoms: readonly string[];
  readonly hidden: boolean;
}

/** A measurement with what the identity check needs. */
export interface MeasuredCell extends MeasurementInfo {
  readonly cell: StructureMeasurementCell;
  /** The state ref of the selection cell it measures; removing it removes the measurement. */
  readonly selectionRef: string;
  /** The loci it measures, in order. */
  readonly locis: readonly StructureElement.Loci[];
  /** A custom text it was made with. */
  readonly text: string;
}

const Order: readonly [MeasureKind, keyof ReturnType<typeof listsOf>][] = [
  ['label', 'labels'],
  ['distance', 'distances'],
  ['angle', 'angles'],
  ['dihedral', 'dihedrals'],
  ['orientation', 'orientations'],
  ['plane', 'planes'],
];

const listsOf = (plugin: PluginContext) => plugin.managers.structure.measurement.state;

const Separator = ' – ';
const lociName = (l: Loci) => stripTags(lociLabel(l, { condensed: true, htmlStyling: false })).trim();

/** The measurements made, in the order the classic list has them (labels, distances, angles, dihedrals, …). */
export function readMeasurements(plugin: PluginContext): MeasuredCell[] {
  const state = listsOf(plugin);
  const out: MeasuredCell[] = [];
  for (const [kind, key] of Order) {
    for (const cell of state[key] as StructureMeasurementCell[]) {
      const m = describe(plugin, cell, kind);
      if (m) out.push(m);
    }
  }
  return out;
}

/** What a measurement cell's object holds: one of these, by kind. */
interface MeasuredData {
  infos?: { loci: Loci }[];
  pairs?: Loci.Bundle<2>[];
  triples?: Loci.Bundle<3>[];
  quads?: Loci.Bundle<4>[];
  locis?: Loci[];
}

/** The loci a measurement measures, in order: no labels. */
function lociOfData(data: MeasuredData): Loci[] {
  if (data.infos) return [data.infos[0].loci];
  if (data.pairs) return [...data.pairs[0].loci];
  if (data.triples) return [...data.triples[0].loci];
  if (data.quads) return [...data.quads[0].loci];
  return data.locis ? [...data.locis] : [];
}

/**
 * The structure-element loci of one measurement, by its ref, without describing it (the labels of every measurement
 * are for lists, not for a pointer crossing a row).
 */
export function measurementLocis(plugin: PluginContext, ref: string): StructureElement.Loci[] {
  for (const [, key] of Order) {
    for (const cell of listsOf(plugin)[key] as StructureMeasurementCell[]) {
      if (cell.transform.ref !== ref) continue;
      const data = cell.obj?.data.sourceData as MeasuredData | undefined;
      return data ? lociOfData(data).filter(StructureElement.Loci.is) : [];
    }
  }
  return [];
}

/**
 * A measurement as last described, by its cell: kept while the cell's object and transform are the same (its atoms are
 * labelled once, not at every read: each operation and each scene sync reads every measurement, so M added one by one
 * were M² labels; AGENTS.md, bulk operations). Only whether it is hidden is read anew.
 */
interface Described {
  readonly obj: unknown;
  readonly transform: unknown;
  readonly kind: MeasureKind;
  last: MeasuredCell;
}
const described = new WeakMap<StructureMeasurementCell, Described>();

function describe(plugin: PluginContext, cell: StructureMeasurementCell, kind: MeasureKind): MeasuredCell | undefined {
  const data = cell.obj?.data.sourceData as MeasuredData | undefined;
  if (!data) return undefined;
  const parent = plugin.state.data.cells.get(cell.transform.parent);
  const hidden = !!parent?.state.isHidden || !!cell.state.isHidden;
  const known = described.get(cell);
  if (known && known.obj === cell.obj && known.transform === cell.transform && known.kind === kind) {
    if (known.last.hidden !== hidden) known.last = { ...known.last, hidden };
    return known.last;
  }
  const params = (cell.transform.params ?? {}) as { customText?: string; unitLabel?: string };
  const text = (params.customText ?? '').trim();
  const loci = lociOfData(data);
  let value = '';
  if (data.infos) value = text;
  else if (data.pairs) {
    value = text || distanceLabel(data.pairs[0], { measureOnly: true, unitLabel: params.unitLabel ?? 'Å' });
  } else if (data.triples) value = text || angleLabel(data.triples[0], { measureOnly: true });
  else if (data.quads) value = text || dihedralLabel(data.quads[0], { measureOnly: true });
  const elements = loci.filter(StructureElement.Loci.is);
  const atoms =
    kind === 'orientation' || kind === 'plane' ? [summaryOf(elements)] : loci.map((l) => lociName(l) || 'Nothing');
  const last: MeasuredCell = {
    ref: cell.transform.ref,
    kind,
    label: kind === 'orientation' || kind === 'plane' ? atoms[0] : atoms.join(Separator),
    value,
    atoms,
    hidden,
    cell,
    selectionRef: cell.transform.parent,
    locis: elements,
    text,
  };
  described.set(cell, { obj: cell.obj, transform: cell.transform, kind, last });
  return last;
}

const summaryOf = (locis: StructureElement.Loci[]) =>
  stripTags(structureElementLociLabelMany(locis, { countsOnly: true, htmlStyling: false })).trim();

/** The measurement of this kind that measures these loci in this order (with this text), if there is one. */
export function findMeasurement(
  plugin: PluginContext,
  kind: MeasureKind,
  locis: readonly StructureElement.Loci[],
  text: string,
): MeasuredCell | undefined {
  const parentOf = (l: StructureElement.Loci) => plugin.helpers.substructureParent.get(l.structure)?.transform.ref;
  // The points asked for, once; each measurement's, once (`bundleOf`): not per comparison.
  const parents = locis.map(parentOf);
  const bundles = locis.map(bundleOf);
  const same = (a: StructureElement.Loci, i: number) =>
    parentOf(a) === parents[i] && StructureElement.Bundle.areEqual(bundleOf(a), bundles[i]);
  return readMeasurements(plugin).find(
    (m) =>
      m.kind === kind &&
      m.text === text.trim() &&
      m.locis.length === locis.length &&
      m.locis.every((l, i) => same(l, i)),
  );
}

/** A loci's bundle, made once per loci (a measurement's loci are kept with its description). */
const bundles = new WeakMap<StructureElement.Loci, StructureElement.Bundle>();
function bundleOf(l: StructureElement.Loci): StructureElement.Bundle {
  let b = bundles.get(l);
  if (!b) {
    b = StructureElement.Bundle.fromLoci(l);
    bundles.set(l, b);
  }
  return b;
}
