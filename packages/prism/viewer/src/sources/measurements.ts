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

// The measurements in words (state.md §1, `managers.structure.measurement`): what the scene tree, the `measurements`
// feature's operations and its Measure section read.

export type MeasureKind = 'distance' | 'angle' | 'dihedral' | 'label' | 'orientation' | 'plane';

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

function describe(plugin: PluginContext, cell: StructureMeasurementCell, kind: MeasureKind): MeasuredCell | undefined {
  const data = cell.obj?.data.sourceData as
    | {
        infos?: { loci: Loci }[];
        pairs?: Loci.Bundle<2>[];
        triples?: Loci.Bundle<3>[];
        quads?: Loci.Bundle<4>[];
        locis?: Loci[];
      }
    | undefined;
  if (!data) return undefined;
  const params = (cell.transform.params ?? {}) as { customText?: string; unitLabel?: string };
  const text = (params.customText ?? '').trim();
  let loci: Loci[] = [];
  let value = '';
  if (data.infos) {
    loci = [data.infos[0].loci];
    value = text;
  } else if (data.pairs) {
    loci = [...data.pairs[0].loci];
    value = text || distanceLabel(data.pairs[0], { measureOnly: true, unitLabel: params.unitLabel ?? 'Å' });
  } else if (data.triples) {
    loci = [...data.triples[0].loci];
    value = text || angleLabel(data.triples[0], { measureOnly: true });
  } else if (data.quads) {
    loci = [...data.quads[0].loci];
    value = text || dihedralLabel(data.quads[0], { measureOnly: true });
  } else if (data.locis) loci = [...data.locis];
  const elements = loci.filter(StructureElement.Loci.is);
  const atoms =
    kind === 'orientation' || kind === 'plane' ? [summaryOf(elements)] : loci.map((l) => lociName(l) || 'Nothing');
  const parent = plugin.state.data.cells.get(cell.transform.parent);
  return {
    ref: cell.transform.ref,
    kind,
    label: kind === 'orientation' || kind === 'plane' ? atoms[0] : atoms.join(Separator),
    value,
    atoms,
    hidden: !!parent?.state.isHidden || !!cell.state.isHidden,
    cell,
    selectionRef: cell.transform.parent,
    locis: elements,
    text,
  };
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
  const same = (a: StructureElement.Loci, b: StructureElement.Loci) =>
    parentOf(a) === parentOf(b) &&
    StructureElement.Bundle.areEqual(StructureElement.Bundle.fromLoci(a), StructureElement.Bundle.fromLoci(b));
  return readMeasurements(plugin).find(
    (m) =>
      m.kind === kind &&
      m.text === text.trim() &&
      m.locis.length === locis.length &&
      m.locis.every((l, i) => same(l, locis[i])),
  );
}

// --- the Measure section's model -------------------------------------------------------------------------------------
