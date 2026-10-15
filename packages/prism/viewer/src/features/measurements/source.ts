/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
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
import {
  type StructureMeasurementCell,
  StructureMeasurementParams,
} from '@molstar/plugin/state/manager/structure/measurement';
import { type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { merge } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { MeasurementKind, MeasurementPoint, MeasurementsModel } from '../../model/measurements.js';
import type { SceneNode } from '../../model/scene.js';
import { fromChanged, SourceScope, sameItems } from '../../sources/state.js';

// The measurements as the views read them (state.md §1, `managers.structure.measurement` ↔ `MeasurementsModel`): the
// measurement manager's cells described in words, and the Measure section's model over them and over the selection's
// additions history (the ordered picks). Reads the plugin, as `sources/` does; every change of measurements goes
// through `measure.*`. The picks themselves are view state, not steps (state.md §3): ordering and dropping them is done
// here on the selection manager, and nothing records it.

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

const Kinds: readonly MeasurementKind[] = [
  { id: 'distance', title: 'Distance', icon: 'measure-distance', points: 2 },
  { id: 'angle', title: 'Angle', icon: 'measure-angle', points: 3 },
  { id: 'dihedral', title: 'Dihedral', icon: 'measure-angle', points: 4 },
  { id: 'label', title: 'Label', icon: 'tag', points: 1 },
  { id: 'orientation', title: 'Orientation', icon: 'move-3d', points: 0 },
  { id: 'plane', title: 'Plane', icon: 'square', points: 0 },
];

const KindIcons: Readonly<Record<MeasureKind, IconName>> = {
  distance: 'measure-distance',
  angle: 'measure-angle',
  dihedral: 'measure-angle',
  label: 'tag',
  orientation: 'move-3d',
  plane: 'square',
};

/** The model, and what releases its subscriptions. */
export interface MeasurementsSource {
  readonly model: MeasurementsModel;
  dispose(): void;
}

class MeasurementNode implements SceneNode {
  readonly kind = 'measurement' as const;
  readonly #label = observableState('');
  readonly #description = observableState<string | undefined>(undefined);
  readonly #hidden = observableState(false);
  readonly #expanded = observableState(false);
  readonly #children = observableState<readonly SceneNode[]>([]);
  icon: IconName = 'ruler';
  removable = true;

  constructor(
    readonly id: string,
    private readonly actions: { visibility(hide: boolean): void; remove(): void; highlight(on: boolean): void },
  ) {}

  get label(): ReadonlyState<string> {
    return this.#label.asReadonly;
  }
  get description(): ReadonlyState<string | undefined> {
    return this.#description.asReadonly;
  }
  get isHidden(): ReadonlyState<boolean> {
    return this.#hidden.asReadonly;
  }
  get isExpanded(): ReadonlyState<boolean> {
    return this.#expanded.asReadonly;
  }
  get children(): ReadonlyState<readonly SceneNode[]> {
    return this.#children.asReadonly;
  }

  set(m: MeasurementInfo): void {
    this.icon = KindIcons[m.kind];
    this.#label.set(m.label);
    this.#description.set(m.value || undefined);
    this.#hidden.set(m.hidden);
  }
  toggleVisibility(): void {
    this.actions.visibility(!this.#hidden.value);
  }
  toggleExpanded(): void {
    this.#expanded.set(!this.#expanded.value);
  }
  highlight(on: boolean): void {
    this.actions.highlight(on);
  }
  focus(): void {}
  remove(): void {
    this.actions.remove();
  }
}

/** The values of the options form: Mol*'s own (`StructureMeasurementParams`). */
const OptionsParams = {
  distanceUnitLabel: StructureMeasurementParams.distanceUnitLabel,
  textColor: StructureMeasurementParams.textColor,
};

const sameEntries = (a: { id: string; label: string }, b: { id: string; label: string }) =>
  a.id === b.id && a.label === b.label;

/**
 * The model of the Measure section over the plugin. `api` is bound to the person: adding, removing and showing or
 * hiding a measurement are its operations; the points are the selection's additions history (Mol*'s, newest first).
 */
export function createMeasurementsSource(plugin: PluginContext, api: ViewerApi): MeasurementsSource {
  const scope = new SourceScope();
  const selection = plugin.managers.structure.selection;
  const manager = plugin.managers.structure.measurement;

  const points = fromChanged(
    scope,
    (): readonly MeasurementPoint[] =>
      selection.additionsHistory.map((e) => ({ id: e.id, label: stripTags(e.label).trim() })),
    selection.events.additionsHistoryUpdated,
    (a, b) => sameItems(a, b, sameEntries),
  );
  const entry = (id: string) => selection.additionsHistory.find((e) => e.id === id);

  const nodes = new Map<string, MeasurementNode>();
  const added = fromChanged(
    scope,
    (): readonly SceneNode[] => {
      const list = readMeasurements(plugin);
      const live = new Set(list.map((m) => m.ref));
      for (const ref of [...nodes.keys()]) if (!live.has(ref)) nodes.delete(ref);
      return list.map((m) => {
        let node = nodes.get(m.ref);
        if (!node) {
          node = new MeasurementNode(m.ref, {
            visibility: (hide) => void api.measure?.visibility({ measurement: m.ref, visible: !hide }),
            remove: () => void api.measure?.remove({ measurement: m.ref }),
            highlight: (on) => {
              const loci = readMeasurements(plugin).find((x) => x.ref === m.ref)?.locis[0];
              if (on && loci) void api.highlight.set({ target: { loci } });
              else void api.highlight.clear();
            },
          });
          nodes.set(m.ref, node);
        }
        node.set(m);
        return node;
      });
    },
    merge(manager.behaviors.state, plugin.state.data.events.cell.stateUpdated, plugin.state.data.events.changed),
    (a, b) => sameItems(a, b),
  );

  const options = new ParamsFormModel(OptionsParams, manager.state.options, {
    onChange: (values) =>
      void api.measure?.options({
        distanceUnitLabel: values.distanceUnitLabel,
        textColor: Color.toHexStyle(values.textColor),
      }),
  });
  scope.add(manager.behaviors.state.subscribe((s) => options.load(s.options as PD.Values<typeof OptionsParams>)));

  const model: MeasurementsModel = {
    points,
    movePoint: (id, by) => {
      const e = entry(id);
      if (e) selection.modifyHistory(e, by < 0 ? 'up' : 'down', 4);
    },
    movePointTo: (id, index) => {
      // Dragged to a position: the history moves the entry there, in place (Mol*'s list is the state).
      const e = entry(id);
      if (!e) return;
      const history = selection.additionsHistory;
      history.splice(history.indexOf(e), 1);
      history.splice(Math.max(0, Math.min(history.length, index)), 0, e);
      selection.events.additionsHistoryUpdated.next(undefined);
    },
    removePoint: (id) => {
      const e = entry(id);
      if (e) selection.modifyHistory(e, 'remove');
    },
    focusPoint: (id) => {
      const e = entry(id);
      if (e) void api.camera.frame({ target: { loci: e.loci } });
    },
    highlightPoint: (id) => {
      const e = id ? entry(id) : undefined;
      if (e) void api.highlight.set({ target: { loci: e.loci } });
      else void api.highlight.clear();
    },
    kinds: Kinds,
    add: (kind) => {
      const op = api.measure?.[kind];
      if (op) void op({});
    },
    options,
    added,
    styleOf: () => undefined,
  };
  return { model, dispose: () => scope.dispose() };
}
