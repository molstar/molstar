/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { StructureElement } from '@molstar/model/model/structure';
import { PluginCommands } from '@molstar/plugin/commands';
import { StructureMeasurementParams } from '@molstar/plugin/state/manager/structure/measurement';
import type { PluginContext } from '@molstar/plugin/context';
import {
  ApiArea,
  type ApiAreaContext,
  type ApiContext,
  type ApiFailure,
  ApiResult,
  type ViewerApi,
} from '../../api/area.js';
import { pickNamed } from '../../api/scope.js';
import { StructureParam } from '../../api/params.js';
import { isLociTarget, LanguageParam, TargetParam, type TargetValue } from '../../api/target-param.js';
import {
  createMeasurementsSource,
  findMeasurement,
  type MeasuredCell,
  type MeasurementInfo,
  type MeasurementsSource,
  type MeasureKind,
  MeasureKinds,
  readMeasurements,
} from './source.js';

// The `measure` area (operations.md §2): Mol*'s measurement manager (`managers.structure.measurement`) as operations.
//
// A point is the center of what its target selects (Mol* measures between centers): an atom for `A:10@CA`, a residue's
// or a ligand's center for `A:10` or `ligand`. A target that selects nothing is `invalid`; one that is on several
// structures is `ambiguous` (name the structure). With no targets at all, a measurement takes the latest picks of the
// selection (its additions history, newest first: distance = the last two), as the classic "Distance from the last
// two picks"; an orientation and a plane take the current selection, as the classic ones do.
//
// Identity (operations.md §1, decision 30): a measurement is the same when its kind, its points in order and its text
// are (so A→B and B→A are two); adding it again is "already so" unless `copy` asks for another.

export type { MeasureKind, MeasurementInfo } from './source.js';

const TextParam = (what: string) => PD.Text('', { description: `${what}; empty: the measured value` });

const Pick = (n: string) => TargetParam(`The ${n} point (A:10@CA, ligand, resn HEM, …); empty: from the latest picks`);

/** Which points the picks give, how many a kind takes. */
const PointCount: Readonly<Record<'distance' | 'angle' | 'dihedral' | 'label', number>> = {
  distance: 2,
  angle: 3,
  dihedral: 4,
  label: 1,
};

const isGiven = (t: unknown) => t !== undefined && t !== null && t !== '';

const word = (n: number) => (n === 1 ? 'point' : 'points');

export class MeasureApi extends ApiArea {
  readonly id = 'measure';

  #source: MeasurementsSource | undefined;
  readonly #unregister: () => void;

  constructor(ctx: ApiAreaContext) {
    super(ctx);
    // The manager's options are not in the plugin's snapshot: registered, a checkpoint holds them, so Undo of
    // `measure.options` (and a session) restores them with the measurements they style.
    this.#unregister = ctx.plugin.state.registerCustomState('prism-viewer.measure', {
      get: () => {
        const o = ctx.plugin.managers.structure.measurement.state.options;
        return { version: 1, distanceUnitLabel: o.distanceUnitLabel, textColor: Color.toHexStyle(o.textColor) };
      },
      set: async (value) => {
        const v = value as { version?: number; distanceUnitLabel?: string; textColor?: string } | undefined;
        const defaults = PD.getDefaultValues(StructureMeasurementParams);
        const options =
          v?.version === 1 && typeof v.distanceUnitLabel === 'string' && /^#[0-9a-f]{6}$/i.test(v.textColor ?? '')
            ? { distanceUnitLabel: v.distanceUnitLabel, textColor: Color.fromHexStyle(v.textColor!) }
            : defaults;
        const manager = ctx.plugin.managers.structure.measurement;
        const now = manager.state.options;
        if (now.distanceUnitLabel === options.distanceUnitLabel && now.textColor === options.textColor) return;
        await manager.setOptions(options);
      },
    });
  }

  /** The model of the Measure section, over the plugin (made once; released with the area). */
  model(api: ViewerApi): import('../../model/measurements.js').MeasurementsModel {
    this.#source ??= createMeasurementsSource(this.ctx.plugin, api);
    return this.#source.model;
  }

  override dispose(): void {
    this.#unregister();
    this.#source?.dispose();
    this.#source = undefined;
  }

  readonly distance = this.operation({
    id: 'measure.distance',
    title: 'Measure a distance',
    params: {
      a: Pick('first'),
      b: Pick('second'),
      text: TextParam('A text to show instead'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measurePoints(ctx, 'distance', [v.a, v.b], v);
    },
  });

  readonly angle = this.operation({
    id: 'measure.angle',
    title: 'Measure an angle',
    params: {
      a: Pick('first'),
      b: Pick('middle (the vertex)'),
      c: Pick('third'),
      text: TextParam('A text to show instead'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measurePoints(ctx, 'angle', [v.a, v.b, v.c], v);
    },
  });

  readonly dihedral = this.operation({
    id: 'measure.dihedral',
    title: 'Measure a dihedral angle',
    params: {
      a: Pick('first'),
      b: Pick('second'),
      c: Pick('third'),
      d: Pick('fourth'),
      text: TextParam('A text to show instead'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measurePoints(ctx, 'dihedral', [v.a, v.b, v.c, v.d], v);
    },
  });

  readonly label = this.operation({
    id: 'measure.label',
    title: 'Add a label',
    params: {
      a: Pick('labelled'),
      text: TextParam('The label text'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measurePoints(ctx, 'label', [v.a], v);
    },
  });

  readonly orientation = this.operation({
    id: 'measure.orientation',
    title: 'Show the principal axes',
    params: {
      target: TargetParam('What to orient (a ligand, A:10-40, …); empty: the current selection'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measureWhole(ctx, 'orientation', v);
    },
  });

  readonly plane = this.operation({
    id: 'measure.plane',
    title: 'Show the best-fit plane',
    params: {
      target: TargetParam('What to fit a plane to (a ligand, A:10-40, …); empty: the current selection'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run(v, ctx) {
      return measureWhole(ctx, 'plane', v);
    },
  });

  readonly remove = this.operation({
    id: 'measure.remove',
    title: 'Remove a measurement',
    params: {
      measurement: PD.Text('', {
        description: 'A measurement by state ref or by its label ("A 10 CA – A 20 CA"); empty: the only one',
      }),
    },
    approval: 'destructive',
    async run({ measurement }, ctx) {
      const found = oneMeasurement(ctx, measurement);
      if (ApiResult.isFailure(found)) return found;
      const state = ctx.plugin.state.data;
      await PluginCommands.State.RemoveObject(ctx.plugin, {
        state,
        ref: found.selectionRef,
        removeParentGhosts: true,
      });
      return ApiResult.ok(`Removed the ${found.kind} ${found.label}.`, { nodes: [found.ref] });
    },
  });

  readonly visibility = this.operation({
    id: 'measure.visibility',
    title: 'Show or hide a measurement',
    params: {
      visible: PD.Boolean(true),
      measurement: PD.Text('', { description: 'A measurement by state ref or label; empty: the only one' }),
    },
    async run({ visible, measurement }, ctx) {
      const found = oneMeasurement(ctx, measurement);
      if (ApiResult.isFailure(found)) return found;
      if (found.hidden === !visible) {
        return ApiResult.already(`The ${found.kind} ${found.label} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      await PluginCommands.State.ToggleVisibility(ctx.plugin, {
        state: ctx.plugin.state.data,
        ref: found.selectionRef,
      });
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} the ${found.kind} ${found.label}.`, { nodes: [found.ref] });
    },
  });

  readonly options = this.operation({
    id: 'measure.options',
    title: 'Set the measurement options',
    params: {
      distanceUnitLabel: PD.Text('', { description: 'The unit written after distances (Å); empty: keep' }),
      textColor: PD.Text('', { description: 'The text color as #rrggbb; empty: keep' }),
    },
    async run({ distanceUnitLabel, textColor }, ctx) {
      const manager = ctx.plugin.managers.structure.measurement;
      const current = manager.state.options;
      let color = current.textColor;
      if (textColor) {
        if (!/^#[0-9a-f]{6}$/i.test(textColor.trim())) {
          return ApiResult.fail('invalid', `'${textColor}' is not a color.`, 'Give it as #rrggbb, e.g. #ff8800.');
        }
        color = Color.fromHexStyle(textColor.trim());
      }
      const next = { distanceUnitLabel: distanceUnitLabel || current.distanceUnitLabel, textColor: color };
      if (next.distanceUnitLabel === current.distanceUnitLabel && next.textColor === current.textColor) {
        return ApiResult.already('The measurement options are already so.');
      }
      await manager.setOptions(next);
      return ApiResult.ok(
        `Measurements: unit ${next.distanceUnitLabel}, text color ${Color.toHexStyle(next.textColor)}.`,
        { nodes: readMeasurements(ctx.plugin).map((m) => m.ref) },
      );
    },
  });

  readonly get = this.operation({
    id: 'measure.get',
    title: 'Describe the measurements',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = readMeasurements(ctx.plugin).map(infoOf);
      if (!list.length) return ApiResult.ok('There are no measurements.', {}, { value: list });
      const lines = list.map(
        (m) => `${m.kind} ${m.label}${m.value ? `: ${m.value}` : ''}${m.hidden ? ' (hidden)' : ''} (${m.ref})`,
      );
      return ApiResult.ok(lines.join('\n'), {}, { value: list });
    },
  });
}

/** Said to callers: the measurement without the cell. */
const infoOf = (m: MeasuredCell): MeasurementInfo => ({
  ref: m.ref,
  kind: m.kind,
  label: m.label,
  value: m.value,
  atoms: m.atoms,
  hidden: m.hidden,
});

/** One measurement by state ref (its representation's or its selection's) or label; empty: the only one. */
function oneMeasurement(ctx: ApiContext, name: string | undefined): MeasuredCell | ApiFailure {
  const list = readMeasurements(ctx.plugin);
  return pickNamed(list, name?.trim() || undefined, {
    what: 'measurement',
    ids: (m) => [m.ref, m.selectionRef],
    names: (m) => [m.label, `${MeasureKinds[m.kind].title} ${m.value}`.trim()],
    describe: (m) => `${m.kind} ${m.label} (${m.ref})`,
    none: { message: 'There are no measurements.', hint: 'Measure something first (measure.distance, …).' },
    listed: 'Measurements',
    several: 'measurements exist',
  });
}

// --- points ---------------------------------------------------------------------------------------------------------

interface PointOptions {
  readonly structure: string;
  readonly language: 'pymol' | 'vmd' | 'jmol' | 'molscript';
  readonly text: string;
}

/** The loci a target gives: one structure's, and not empty. */
async function lociOf(
  ctx: ApiContext,
  target: TargetValue,
  o: Pick<PointOptions, 'structure' | 'language'>,
  what: string,
): Promise<StructureElement.Loci | ApiFailure> {
  if (isLociTarget(target)) {
    return StructureElement.Loci.is(target.loci) && !StructureElement.Loci.isEmpty(target.loci)
      ? target.loci
      : ApiResult.fail('invalid', `The ${what} point selects nothing.`);
  }
  const resolved = await ctx.targets.resolve(target, {
    ...(o.structure ? { structure: o.structure } : {}),
    ...(o.language ? { language: o.language } : {}),
  });
  if (ApiResult.isFailure(resolved)) return resolved;
  const parts = resolved.parts.filter((p) => p.selected.atoms > 0);
  if (!parts.length) {
    return ApiResult.fail(
      'invalid',
      `The ${what} point selects nothing: ${resolved.summary}`,
      'Name an atom or a residue that exists, e.g. A:10@CA.',
    );
  }
  if (parts.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `The ${what} point is on ${parts.length} structures.`,
      `Name the structure: ${parts.map((p) => `${p.structure.label} (${p.structure.ref})`).join(', ')}.`,
    );
  }
  return parts[0].loci;
}

/** The latest picks of the selection, newest first (the classic "Distance from the last two picks"). */
function picks(plugin: PluginContext): StructureElement.Loci[] {
  return plugin.managers.structure.selection.additionsHistory.map((e) => e.loci);
}

/** A measurement between points: the targets in order, else the latest picks. */
async function measurePoints(
  ctx: ApiContext,
  kind: keyof typeof PointCount,
  targets: readonly TargetValue[],
  o: PointOptions,
): Promise<ApiResult<MeasurementInfo>> {
  const n = PointCount[kind];
  const given = targets.filter(isGiven).length;
  let locis: StructureElement.Loci[];
  if (given === 0) {
    const latest = picks(ctx.plugin);
    if (latest.length < n) {
      return ApiResult.fail(
        'invalid',
        `A ${kind} needs ${n} ${word(n)}: no targets were given and ${latest.length === 0 ? 'nothing is' : `only ${latest.length} ${word(latest.length)} ${latest.length === 1 ? 'is' : 'are'}`} picked.`,
        `Give ${n} ${n === 1 ? 'target' : 'targets'} (${'abcd'.slice(0, n).split('').join(', ')}), or pick atoms or residues in Select mode (Shift adds): the latest picks become the points.`,
      );
    }
    locis = latest.slice(0, n);
  } else if (given < n) {
    return ApiResult.fail(
      'invalid',
      `A ${kind} needs ${n} ${word(n)}, ${given} given.`,
      `Give all of ${'abcd'.slice(0, n).split('').join(', ')}, or none to use the latest picks.`,
    );
  } else {
    locis = [];
    const names = ['first', 'second', 'third', 'fourth'];
    for (let i = 0; i < n; i++) {
      const l = await lociOf(ctx, targets[i], o, n === 1 ? 'labelled' : names[i]);
      if (ApiResult.isFailure(l)) return l;
      locis.push(l);
    }
  }
  return addMeasurement(ctx, kind, locis, o.text.trim());
}

/** An orientation or a plane: a target's loci, else the current selection. */
async function measureWhole(
  ctx: ApiContext,
  kind: 'orientation' | 'plane',
  o: { target: TargetValue; structure: string; language: PointOptions['language'] },
): Promise<ApiResult<MeasurementInfo>> {
  let locis: StructureElement.Loci[];
  if (isGiven(o.target)) {
    const l = isLociTarget(o.target) ? await lociOf(ctx, o.target, o, 'measured') : await partsOf(ctx, o);
    if (ApiResult.isFailure(l)) return l;
    locis = Array.isArray(l) ? l : [l];
  } else {
    locis = [...ctx.plugin.managers.structure.selection.entries.values()]
      .map((e) => e.selection)
      .filter((l) => !StructureElement.Loci.isEmpty(l));
    if (!locis.length) {
      return ApiResult.fail(
        'invalid',
        `A ${kind} needs something to measure: no target was given and nothing is selected.`,
        'Give a target (a ligand, A:10-40, …) or select something first.',
      );
    }
  }
  return addMeasurement(ctx, kind, locis, '');
}

async function partsOf(
  ctx: ApiContext,
  o: { target: TargetValue; structure: string; language: PointOptions['language'] },
): Promise<StructureElement.Loci[] | ApiFailure> {
  const resolved = await ctx.targets.resolve(o.target, {
    ...(o.structure ? { structure: o.structure } : {}),
    ...(o.language ? { language: o.language } : {}),
  });
  if (ApiResult.isFailure(resolved)) return resolved;
  const parts = resolved.parts.filter((p) => p.selected.atoms > 0);
  if (!parts.length) {
    return ApiResult.fail(
      'invalid',
      `The target selects nothing: ${resolved.summary}`,
      'Name something that exists, e.g. ligand or A:10-40.',
    );
  }
  return parts.map((p) => p.loci);
}

async function addMeasurement(
  ctx: ApiContext,
  kind: MeasureKind,
  locis: StructureElement.Loci[],
  text: string,
): Promise<ApiResult<MeasurementInfo>> {
  const plugin = ctx.plugin;
  const manager = plugin.managers.structure.measurement;
  if (!ctx.copy) {
    const same = findMeasurement(plugin, kind, locis, text);
    if (same) {
      return ApiResult.already(`The ${kind} ${same.label} is already measured${same.value ? `: ${same.value}` : ''}.`, {
        value: infoOf(same),
      });
    }
  }
  const customText = text ? { customText: text } : undefined;
  let made;
  switch (kind) {
    case 'distance':
      made = await manager.addDistance(locis[0], locis[1], customText);
      break;
    case 'angle':
      made = await manager.addAngle(locis[0], locis[1], locis[2], customText);
      break;
    case 'dihedral':
      made = await manager.addDihedral(locis[0], locis[1], locis[2], locis[3], customText);
      break;
    case 'label':
      made = await manager.addLabel(locis[0], customText ? { visualParams: customText } : undefined);
      break;
    case 'orientation':
      made = await manager.addOrientation(locis);
      break;
    case 'plane':
      made = await manager.addPlane(locis);
      break;
  }
  if (!made) {
    return ApiResult.fail(
      'failed',
      `The ${kind} could not be made.`,
      'Its structure is not part of the scene: measure atoms of a structure that is loaded.',
    );
  }
  const cell = plugin.state.data.cells.get(made.representation.ref);
  if (!cell || cell.status === 'error') {
    return ApiResult.fail('failed', `The ${kind} could not be made: ${cell?.errorText ?? 'its cell was not made'}.`);
  }
  const info = readMeasurements(plugin).find((m) => m.ref === made.representation.ref);
  if (!info) return ApiResult.fail('failed', `The ${kind} could not be read back.`);
  return ApiResult.ok(
    `Added the ${kind} ${info.label}${info.value ? `: ${info.value}` : ''}.`,
    { nodes: [info.ref] },
    { value: infoOf(info) },
  );
}
