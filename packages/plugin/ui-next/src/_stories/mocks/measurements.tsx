/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { MeasurementKind, MeasurementPoint, MeasurementsModel } from '../../model/measurements.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';
import { MeasurementsAdd } from '../../views/structure/measurements-add.js';
import type { MockScene } from './scene.js';

const Kinds: MeasurementKind[] = [
  { id: 'distance', title: 'Distance', icon: 'measure-distance', points: 2 },
  { id: 'angle', title: 'Angle', icon: 'measure-angle', points: 3 },
  { id: 'dihedral', title: 'Dihedral', icon: 'measure-angle', points: 4 },
  { id: 'label', title: 'Label', icon: 'tag', points: 1 },
  { id: 'orientation', title: 'Orientation', icon: 'move-3d', points: 0 },
  { id: 'plane', title: 'Plane', icon: 'square', points: 0 },
];

// What each would read for the stories' points (CYS 442 SG, HEM 508 FE, …).
const Values: Record<string, string> = {
  distance: '2.31 Å',
  angle: '104.2°',
  dihedral: '−62.8°',
  label: '',
  orientation: 'principal axes',
  plane: 'best fit',
};

/** Measurements for the stories: two points selected (a distance is ready), adding creates rows under Measurements. */
export function mockMeasurements(scene: MockScene, post: (text: string) => void): MeasurementsModel {
  const points = observableState<readonly MeasurementPoint[]>([
    { id: 'p1', label: 'CYS 442 SG' },
    { id: 'p2', label: 'HEM 508 FE' },
    { id: 'p3', label: 'ARG 105 NH1' },
  ]);
  const label = (p: readonly MeasurementPoint[], n: number) =>
    p
      .slice(0, n)
      .map((x) => x.label)
      .join(' – ');
  return {
    points: points.asReadonly,
    movePoint: (id, by) => {
      const list = [...points.value];
      const i = list.findIndex((p) => p.id === id);
      const j = i + by;
      if (i < 0 || j < 0 || j >= list.length) return;
      [list[i], list[j]] = [list[j], list[i]];
      points.set(list);
    },
    removePoint: (id) => points.set(points.value.filter((p) => p.id !== id)),
    focusPoint: (id) => post(`Focused on ${points.value.find((p) => p.id === id)?.label}`),
    kinds: Kinds,
    add: (kind) => {
      const k = Kinds.find((x) => x.id === kind)!;
      const p = points.value;
      const what = k.points === 0 ? 'Selection' : label(p, k.points);
      scene.addMeasurement(
        k.points === 0 ? `${k.title} · selection` : what,
        k.points === 0 ? Values[kind] : `${k.title}${Values[kind] ? ` · ${Values[kind]}` : ''}`,
        k.icon,
      );
      post(`Added ${k.title.toLowerCase()}: ${what}`);
    },
    options: new ParamsFormModel({
      distanceUnitLabel: PD.Text('Å', { label: 'Distance unit' }),
      textColor: PD.Color(0x000000 as any, { label: 'Text color' }),
      textSize: PD.Numeric(0.5, { min: 0.1, max: 5, step: 0.1 }, { label: 'Text size' }),
    }),
  };
}

/** The Measurements row's inspector: adding measurements (contributed to its slot, as the plugin's UI would). */
export function registerMockMeasurements(shell: ShellModel, model: MeasurementsModel) {
  shell.registerSections([
    {
      id: 'measurements-add',
      slot: 'inspector:measurements',
      title: 'Add',
      icon: 'plus',
      order: 0,
      render: domView(() => <MeasurementsAdd model={model} />),
    },
  ]);
}
