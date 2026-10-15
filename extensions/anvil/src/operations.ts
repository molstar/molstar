/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { StateSelection } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Structure, StructureSelection } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from '@molstar/prism-viewer/api/area';
import { StructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, oneStructure } from '@molstar/prism-viewer/api/scope';
import { ANVILParams, type ANVILProps } from './algorithm.js';
import { isTransmembrane, MembraneOrientation3D } from './behavior.js';
import { MembraneOrientation, MembraneOrientationProvider } from './prop.js';
import type { MembraneOrientationProps as LookProps } from './representation.js';

// The `membrane` area (operations.md §2, parity.md §6): ANVIL's membrane orientation as operations. `show` places the
// membrane (computed on first use, locally) and adds its planes and rims, `hide` removes them, `fit` places it again
// with ANVIL's settings, `options` sets how it looks (v5's generic representation control), `get` says where it is.
// ANVIL's settings live in the structure's custom-properties node, as v5's custom-property controls kept them, so Undo
// and sessions carry them.

const Tag = MembraneOrientation.Tag.Representation;
const PropName = MembraneOrientationProvider.descriptor.name;

/** Where the membrane is, as `get` says it. */
export interface MembraneInfo {
  readonly structure: string;
  /** Its planes and rims are in the scene. */
  readonly shown: boolean;
  /** The distance between the two planes, in Å. */
  readonly thickness: number;
  /** The planes' radius (without the look's radius factor), in Å. */
  readonly radius: number;
  /** The membrane's normal, a unit vector. */
  readonly normal: readonly [number, number, number];
  readonly center: readonly [number, number, number];
  /** Residues with an atom between the planes. */
  readonly embeddedResidues: number;
  /** ANVIL's settings it was placed with. */
  readonly settings: ANVILProps;
}

const Switch = (what: string) =>
  PD.Select<'keep' | 'on' | 'off'>(
    'keep',
    [
      ['keep', 'Keep'],
      ['on', 'On'],
      ['off', 'Off'],
    ],
    { description: `${what}: on, off, or keep` },
  );

const LookParams = {
  planes: Switch('The two planes'),
  rims: Switch('The rims around them'),
  dashed: Switch('Dashed rims'),
  color: PD.Text('', { description: 'A CSS color for the planes and rims (lightgrey, #cccccc); empty: keep' }),
  opacity: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { description: "The planes' opacity; 0: keep" }),
  radius: PD.Numeric(0, { min: 0, max: 3, step: 0.01 }, { description: "The planes' radius as a factor; 0: keep" }),
  structure: StructureParam(),
};

export class MembraneApi extends ApiArea {
  readonly id = 'membrane';

  readonly get = this.operation({
    id: 'membrane.get',
    title: 'Describe the membrane',
    params: { structure: StructureParam() },
    read: true,
    async run({ structure }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      const placed = await place(ctx, s);
      if (isFailure(placed)) return placed;
      const info = await describe(ctx, s);
      return ApiResult.ok(summary(info), {}, { value: info });
    },
  });

  readonly show = this.operation({
    id: 'membrane.show',
    title: 'Show the membrane',
    params: { structure: StructureParam() },
    async run({ structure }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      if (membraneRef(ctx.plugin, s)) return ApiResult.already(`${labelOf(ctx, s)} already shows its membrane.`);
      const failed = await apply(ctx, s, undefined);
      if (failed) return failed;
      const info = await describe(ctx, s);
      return ApiResult.ok(`Showed the membrane of ${info.structure}: ${placement(info)}.`, changed(s), { value: info });
    },
  });

  readonly hide = this.operation({
    id: 'membrane.hide',
    title: 'Hide the membrane',
    params: { structure: StructureParam() },
    async run({ structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const ref = membraneRef(ctx.plugin, s);
      if (!ref) return ApiResult.already(`${labelOf(ctx, s)} shows no membrane.`);
      await ctx.plugin.state.data.build().delete(ref).commit();
      return ApiResult.ok(`Hid the membrane of ${labelOf(ctx, s)}.`, changed(s));
    },
  });

  readonly fit = this.operation({
    id: 'membrane.fit',
    title: 'Fit the membrane',
    params: { ...ANVILParams, structure: StructureParam() },
    async run({ structure, ...settings }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      const root = s.cell.obj!.data.root;
      const now = MembraneOrientationProvider.get(root).value ? MembraneOrientationProvider.props(root) : undefined;
      if (membraneRef(ctx.plugin, s) && now && PD.areEqual(ANVILParams, now, settings)) {
        return ApiResult.already(`The membrane of ${labelOf(ctx, s)} is already fitted with these settings.`);
      }
      const failed = await apply(ctx, s, settings);
      if (failed) return failed;
      const info = await describe(ctx, s);
      return ApiResult.ok(`Fitted the membrane of ${info.structure}: ${placement(info)}.`, changed(s), {
        value: info,
      });
    },
  });

  readonly options = this.operation({
    id: 'membrane.options',
    title: 'Set how the membrane looks',
    params: LookParams,
    async run({ structure, planes, rims, dashed, color, opacity, radius }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const ref = membraneRef(ctx.plugin, s);
      if (!ref) {
        return ApiResult.fail('blocked', `${labelOf(ctx, s)} shows no membrane.`, 'Show it first (membrane.show).');
      }
      const painted = color ? Color.fromHexStyle(cssHex(color) ?? '') : undefined;
      if (color && !cssHex(color)) {
        return ApiResult.fail('invalid', `'${color}' is not a color.`, 'Give a color name (lightgrey) or #rrggbb.');
      }
      const old = ctx.plugin.state.data.tree.transforms.get(ref).params as LookProps;
      const visuals = new Set(old.visuals);
      if (planes !== 'keep') visuals[planes === 'on' ? 'add' : 'delete']('bilayer-planes');
      if (rims !== 'keep') visuals[rims === 'on' ? 'add' : 'delete']('bilayer-rims');
      if (!visuals.size) {
        return ApiResult.fail(
          'invalid',
          'That would show neither planes nor rims.',
          'Hide it instead (membrane.hide).',
        );
      }
      const next: LookProps = {
        ...old,
        visuals: (['bilayer-planes', 'bilayer-rims'] as const).filter((x) => visuals.has(x)),
        dashedLines: dashed === 'keep' ? old.dashedLines : dashed === 'on',
        color: painted ?? old.color,
        sectorOpacity: opacity || old.sectorOpacity,
        radiusFactor: radius || old.radiusFactor,
      };
      const keys = ['visuals', 'dashedLines', 'color', 'sectorOpacity', 'radiusFactor'] as const;
      if (keys.every((k) => JSON.stringify(next[k]) === JSON.stringify(old[k]))) {
        return ApiResult.already(`The membrane of ${labelOf(ctx, s)} already looks so.`);
      }
      await ctx.plugin.state.data.build().to(ref).update(next).commit();
      return ApiResult.ok(`Changed how the membrane of ${labelOf(ctx, s)} looks.`, { nodes: [ref] });
    },
  });
}

// --- helpers ----------------------------------------------------------------------------------------------------

const changed = (s: StructureRef) => ({ nodes: [s.cell.transform.ref] });

function applicable(ctx: ApiContext, name: string): StructureRef | ApiFailure {
  const s = oneStructure(ctx, name || undefined);
  if (isFailure(s)) return s;
  const data = s.cell.obj?.data;
  if (!data || !MembraneOrientationProvider.isApplicable(data)) {
    return ApiResult.fail(
      'blocked',
      `ANVIL places a membrane around atomic protein chains of 15 residues or more; ${labelOf(ctx, s)} has none.`,
    );
  }
  return s;
}

function membraneRef(plugin: PluginContext, s: StructureRef): string | undefined {
  return StateSelection.findTagInSubtree(plugin.state.data.tree, s.cell.transform.ref, Tag);
}

/** Places the membrane with the settings the structure has (computed once; local, but may take a moment). */
async function place(ctx: ApiContext, s: StructureRef): Promise<MembraneOrientation | ApiFailure> {
  const data = s.cell.obj!.data;
  const plugin = ctx.plugin;
  await plugin.runTask(
    Task.create('Membrane orientation', (runtime) =>
      MembraneOrientationProvider.attach(
        { runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext },
        data,
      ),
    ),
  );
  return MembraneOrientationProvider.get(data.root).value ?? cannotPlace(ctx, s);
}

const cannotPlace = (ctx: ApiContext, s: StructureRef) =>
  ApiResult.fail(
    'failed',
    `ANVIL could not place a membrane for ${labelOf(ctx, s)}.`,
    'It needs protein chains with enough residues; short peptides fail.',
  );

/**
 * Sets ANVIL's settings, in the structure's custom-properties node where it has one (so the state holds them), and
 * adds the membrane if it is not shown; one state update.
 */
async function apply(
  ctx: ApiContext,
  s: StructureRef,
  settings: ANVILProps | undefined,
): Promise<ApiFailure | undefined> {
  const plugin = ctx.plugin;
  const data = s.cell.obj!.data;
  const props = s.properties?.cell;
  const build = plugin.state.data.build();
  if (settings && props) {
    const p = props.transform.params as { autoAttach: string[]; properties: Record<string, unknown> };
    build.to(props.transform.ref).update({
      autoAttach: p.autoAttach.includes(PropName) ? p.autoAttach : [...p.autoAttach, PropName],
      properties: { ...p.properties, [PropName]: settings },
    });
  } else if (settings) {
    await plugin.runTask(
      Task.create('Membrane orientation', (runtime) =>
        MembraneOrientationProvider.attach(
          { runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext },
          data,
          settings,
        ),
      ),
    );
    // Without the node, the planes read the property when they are made: made again for other settings.
    const old = membraneRef(plugin, s);
    if (old) build.delete(old);
  }
  build.to(s.cell.transform.ref).applyOrUpdateTagged(Tag, MembraneOrientation3D);
  try {
    await build.commit({ revertOnError: true });
  } catch {
    return cannotPlace(ctx, s);
  }
  return MembraneOrientationProvider.get(data.root).value ? undefined : cannotPlace(ctx, s);
}

async function describe(ctx: ApiContext, s: StructureRef): Promise<MembraneInfo> {
  const data = s.cell.obj!.data as Structure;
  const m = MembraneOrientationProvider.get(data.root).value!;
  const between = Vec3.sub(Vec3(), m.planePoint1, m.planePoint2);
  const selection = await ctx.plugin.runTask(
    Task.create('Membrane residues', (runtime) => isTransmembrane.getSelection(ctx.plugin, runtime, data)),
  );
  const round = (x: number, d = 1) => Math.round(x * 10 ** d) / 10 ** d;
  const v3 = (v: Vec3, d: number) => [round(v[0], d), round(v[1], d), round(v[2], d)] as const;
  return {
    structure: labelOf(ctx, s),
    shown: !!membraneRef(ctx.plugin, s),
    thickness: round(Math.abs(Vec3.dot(between, m.normalVector))),
    radius: round(m.radius),
    normal: v3(m.normalVector, 3),
    center: v3(m.centroid, 1),
    embeddedResidues: StructureSelection.isEmpty(selection)
      ? 0
      : StructureSelection.unionStructure(selection).polymerResidueCount,
    settings: { ...MembraneOrientationProvider.props(data.root) },
  };
}

const placement = (info: MembraneInfo) =>
  `${info.thickness} Å thick, normal (${info.normal.join(', ')}), ${info.embeddedResidues} residue${info.embeddedResidues === 1 ? '' : 's'} in it`;

function summary(info: MembraneInfo): string {
  return `The membrane of ${info.structure}: ${placement(info)}; ${info.shown ? 'shown' : 'not shown'}.`;
}

/** A CSS color name or #rrggbb, as #rrggbb. */
function cssHex(color: string): string | undefined {
  const c = color.trim().toLowerCase();
  if (/^#[0-9a-f]{6}$/.test(c)) return c;
  const named = (ColorNames as Record<string, Color | undefined>)[c];
  return named === undefined ? undefined : Color.toHexStyle(named);
}
