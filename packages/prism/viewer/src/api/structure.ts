/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from '@molstar/core/state';
import { SpacegroupCell } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { decodeColor } from '@molstar/core/util/color/utils';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import type { Structure } from '@molstar/model/model/structure';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { getTransformFromParams, type TransformParam } from '@molstar/plugin/state/transforms/helpers';
import {
  StructureFromModel,
  StructureFromTrajectory,
  TransformStructureConformation,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';
import { deepEqual } from '@molstar/core/util';
import { ModelUnitcell3D } from '@molstar/plugin/state/transforms/structure/unitcell';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { StructureInfo } from '../sources/viewer-sources.js';
import { type ApiContext, ApiArea, ApiResult } from './area.js';
import { reapplyContents } from './contents.js';
import { EveryStructureParam, OptionalParam, StructureParam } from './params.js';
import { entryRoot, isFailure, labelOf, oneStructure, someStructures } from './scope.js';
import { defaultPreset } from './styles.js';

// `structure` (api.md §2): the scene outline (`get`), representation presets by id, the model of a
// trajectory, the assembly, a structure's visibility and removal. The first scene area: the others follow its shape
// (a `structure` field naming one by ref or label, else the one the person works with; registries as arguments, their
// ids in the hint when one is wrong; "already so" when nothing would change). Parity (parity.md §2): the unit
// cell (`unit-cell`, v5's model-level Unit Cell representation) and the Conformation Transform (`transform`,
// `transform-reset`: v5's `TransformStructureConformation` decorator, the one superposition sets too).

export class StructureApi extends ApiArea {
  readonly id = 'structure';

  /** The loaded structures, as a state to subscribe to. */
  get loaded(): ReadonlyState<readonly StructureInfo[]> {
    return this.ctx.sources.structures;
  }

  readonly list = this.operation({
    id: 'structure.list',
    aliases: ['structure.get'],
    title: 'Describe the scene',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = ctx.sources.structures.value;
      if (!list.length) return ApiResult.ok('Nothing is loaded.', {}, { value: list });
      const lines = list.map(
        (s) =>
          `${s.label} (${s.ref}): ${s.atoms.toLocaleString('en')} atoms, ${s.residues.toLocaleString('en')} residues, ${s.chains} chain${s.chains === 1 ? '' : 's'}${s.components.length ? `; ${s.components.join(', ')}` : ''}`,
      );
      return ApiResult.ok(lines.join('\n'), {}, { value: list });
    },
  });

  /** @deprecated The former name of `list` (api.md §1 "Words"), until the release. */
  readonly get = this.list;

  readonly preset = this.operation({
    id: 'structure.preset',
    title: 'Apply a representation preset',
    params: {
      preset: PD.Text('auto', { description: 'A representation preset by id: auto, polymer-cartoon, illustrative, …' }),
      structure: EveryStructureParam(),
    },
    async run({ preset, structure }, ctx) {
      const provider = ctx.plugin.builders.structure.representation.resolveProvider(preset);
      if (!provider) {
        return ApiResult.fail('invalid', `There is no preset '${preset}'.`, `Presets: ${presetIds(ctx).join(', ')}.`);
      }
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const applicable = targets.filter(
        (s) => !s.cell.obj || !provider.isApplicable || provider.isApplicable(s.cell.obj, ctx.plugin),
      );
      if (!applicable.length) {
        return ApiResult.fail(
          'blocked',
          `${provider.display.name} does not apply to ${targets.map((x) => labelOf(ctx, x)).join(', ')}.`,
        );
      }
      await ctx.plugin.managers.structure.component.applyPreset(applicable, provider);
      // What the person hid (entities, chains) stays hidden in the preset's components (decision 18).
      await reapplyContents(ctx, freshStructures(ctx, applicable));
      return ApiResult.ok(`Applied ${provider.display.name} to ${applicable.map((x) => labelOf(ctx, x)).join(', ')}.`, {
        nodes: applicable.map((s) => s.cell.transform.ref),
      });
    },
  });

  readonly assembly = this.operation({
    id: 'structure.assembly',
    title: 'Show an assembly',
    params: {
      id: PD.Text('', { description: "An assembly id ('1', '2', …); empty: the model's deposited coordinates" }),
      structure: StructureParam(),
    },
    async run({ id, structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const model = s.model?.cell.obj?.data;
      if (!model) return ApiResult.fail('blocked', `${labelOf(ctx, s)} is not built from a model.`);
      const assemblies = ModelSymmetry.Provider.get(model)?.assemblies ?? [];
      const current = s.cell.transform.params as { type: { name: string; params?: { id?: string } } };
      if (id && !assemblies.some((a) => a.id === id)) {
        return ApiResult.fail(
          'invalid',
          `${labelOf(ctx, s)} has no assembly '${id}'.`,
          assemblies.length ? `Assemblies: ${assemblies.map((a) => a.id).join(', ')}.` : 'It has no assemblies.',
        );
      }
      // The structure's common params (bond computation, …) carry over, as in `symmetry`.
      const common = commonTypeParams(current.type.params);
      const want = id ? { name: 'assembly', params: { ...common, id } } : { name: 'model', params: common };
      if (current.type.name === want.name && (want.name === 'model' || current.type.params?.id === id)) {
        return ApiResult.already(
          id ? `${labelOf(ctx, s)} already shows assembly ${id}.` : `${labelOf(ctx, s)} already shows the model.`,
        );
      }
      await rebuildStructure(ctx, s.cell.transform.ref, { ...current, type: want });
      return ApiResult.ok(id ? `${labelOf(ctx, s)} shows assembly ${id}.` : `${labelOf(ctx, s)} shows the model.`, {
        nodes: [s.cell.transform.ref],
      });
    },
  });

  // bond-providers plan §4.1: where a structure's bonds come from, one of the plugin's bond providers.
  readonly bonds = this.operation({
    id: 'structure.bonds',
    title: 'Choose where bonds come from',
    params: {
      provider: PD.Text(BondProvider.DefaultName, {
        description: "A bond provider by name ('default': the built-in bonds; 'none': only the file's)",
      }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The provider's params; missing ones get defaults" },
      ),
      structure: StructureParam(),
    },
    async run({ provider, params, structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const registry = ctx.plugin.bondProviders;
      const transform = s.cell.transform;
      const fromTrajectory = transform.transformer === StructureFromTrajectory;
      const model = fromTrajectory
        ? (s.cell.obj?.data as Structure | undefined)?.representativeModel
        : s.model?.cell.obj?.data;
      if (!model || (!fromTrajectory && transform.transformer !== StructureFromModel)) {
        return ApiResult.fail('blocked', `${labelOf(ctx, s)} is not built from a model.`);
      }
      const options = registry.options(model);
      const p = registry.get(provider);
      if (!p || !options.some((o) => o[0] === provider)) {
        return ApiResult.fail(
          'invalid',
          `No bond provider '${provider}' for ${labelOf(ctx, s)}.`,
          `Providers: ${options.map((o) => o[0]).join(', ')}.`,
        );
      }
      const bonds = { name: provider, params: PD.normalizeParams(p.getParams(model, registry), params ?? {}, 'all') };
      const current = transform.params as any;
      const now = (fromTrajectory ? current?.bonds : current?.type?.params?.bonds) ?? {
        name: BondProvider.DefaultName,
        params: {},
      };
      if (deepEqual(now, bonds))
        return ApiResult.already(`${labelOf(ctx, s)} already takes its bonds from ${p.label}.`);
      // a new structure with the same components (StructureFromModel recreates on a param change)
      const next = fromTrajectory
        ? { ...current, bonds }
        : { ...current, type: { ...current.type, params: { ...current.type.params, bonds } } };
      await ctx.plugin.state.updateTransform(ctx.plugin.state.data, transform.ref, next, 'Bonds');
      return ApiResult.ok(`${labelOf(ctx, s)} takes its bonds from ${p.label}.`, { nodes: [transform.ref] });
    },
  });

  // parity.md §2 (the structure's source): crystal symmetry, v5's Symmetry Mates and Symmetry (indices) types.
  readonly symmetry = this.operation({
    id: 'structure.symmetry',
    // a crystal's symmetry, or a simulation box's periodic images (its P1 cell from the box)
    title: 'Show symmetry or periodic images',
    params: {
      kind: PD.Select<'mates' | 'cells' | 'off'>(
        'mates',
        [
          ['mates', 'Symmetry mates within a radius'],
          ['cells', 'Unit cells around the model'],
          ['off', 'Off: the model alone'],
        ],
        { description: 'mates: the copies within the radius; cells: the copies in the cells around; off: none' },
      ),
      radius: PD.Numeric(5, { min: 0, max: 50, step: 1 }, { description: 'Å around the model (mates)' }),
      cells: PD.Numeric(1, { min: 0, max: 3, step: 1 }, { description: 'Cells on each side (cells): 1 is -1 to 1' }),
      structure: StructureParam(),
    },
    async run({ kind, radius, cells, structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const name = labelOf(ctx, s);
      const model = s.model?.cell.obj?.data;
      if (!model) return ApiResult.fail('blocked', `${name} is not built from a model.`);
      const symmetry = ModelSymmetry.Provider.get(model);
      if (kind !== 'off' && symmetry && SpacegroupCell.isZero(symmetry.spacegroup.cell)) {
        return ApiResult.fail(
          'blocked',
          `${name} has no crystal cell or simulation box.`,
          'Only a crystal or a simulation in a periodic box has copies.',
        );
      }
      const current = s.cell.transform.params as { type: { name: string; params?: Record<string, unknown> } };
      // The structure's common params (bond computation, …) carry over; the type's own are set.
      const common = commonTypeParams(current.type.params);
      const n = Math.round(cells);
      const want =
        kind === 'mates'
          ? { name: 'symmetry-mates', params: { ...common, radius } }
          : kind === 'cells'
            ? { name: 'symmetry', params: { ...common, ijkMin: Vec3.create(-n, -n, -n), ijkMax: Vec3.create(n, n, n) } }
            : { name: 'model', params: common };
      const same =
        current.type.name === want.name &&
        (kind === 'off' ||
          (kind === 'mates'
            ? current.type.params?.radius === radius
            : Vec3.exactEquals(current.type.params?.ijkMin as Vec3, Vec3.create(-n, -n, -n)) &&
              Vec3.exactEquals(current.type.params?.ijkMax as Vec3, Vec3.create(n, n, n))));
      if (same) return ApiResult.already(`${name} already shows that.`);
      await rebuildStructure(ctx, s.cell.transform.ref, { ...current, type: want });
      const did =
        kind === 'mates'
          ? `${name} shows its symmetry mates within ${radius} Å.`
          : kind === 'cells'
            ? `${name} shows the unit cells ${-n} to ${n} around it.`
            : `${name} shows the model alone.`;
      return ApiResult.ok(did, { nodes: [s.cell.transform.ref] });
    },
  });

  // parity.md §2: v5's "Unit Cell" row (`ModelUnitcell3D` on the structure's model; the default preset makes it hidden).
  readonly unitCell = this.operation({
    id: 'structure.unit-cell',
    title: 'Show the unit cell',
    params: {
      visible: OptionalParam(
        PD.Boolean(true, {
          description: "Show or hide the crystal's unit cell; absent: shown when nothing else is given",
        }),
      ),
      color: OptionalParam(PD.Text('', { description: "The cell's color, a CSS color (orange, #ff8800)" })),
      opacity: OptionalParam(
        PD.Numeric(1, { min: 0.05, max: 1, step: 0.05 }, { description: 'Opacity from 0.05 to 1' }),
      ),
      structure: StructureParam(),
    },
    async run(v, ctx) {
      // Absent: kept (api.md §1); called with nothing, it shows the cell (what the operation is for).
      const { color, opacity, structure } = v as { color?: string; opacity?: number; structure: string };
      const given = v as { visible?: boolean };
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const name = labelOf(ctx, s);
      const model = s.model;
      if (!model?.cell.obj) return ApiResult.fail('blocked', `${name} is not built from a model.`);
      const symmetry = ModelSymmetry.Provider.get(model.cell.obj.data);
      if (!symmetry || SpacegroupCell.isZero(symmetry.spacegroup.cell)) {
        return ApiResult.fail('blocked', `${name} has no crystal cell.`, 'Only crystal structures have a unit cell.');
      }
      let cellColor: Color | undefined;
      if (color?.trim()) {
        cellColor = decodeColor(color.trim());
        if (cellColor === undefined) {
          return ApiResult.fail('invalid', `'${color}' is not a color.`, 'Use a CSS color: orange, #ff8800.');
        }
      }
      const existing = model.unitcell?.cell;
      const current = existing?.transform.params as { cellColor: Color; alpha: number } | undefined;
      const style = {
        ...(cellColor !== undefined && cellColor !== current?.cellColor ? { cellColor } : {}),
        ...(opacity !== undefined && opacity !== current?.alpha ? { alpha: opacity } : {}),
      };
      const restyle = Object.keys(style).length > 0;
      const shown = !!existing && !existing.state.isHidden;
      const visible = given.visible ?? (color === undefined && opacity === undefined ? true : shown || !existing);
      if (!restyle && shown === visible) {
        return ApiResult.already(`The unit cell of ${name} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      const plugin = ctx.plugin;
      const state = plugin.state.data;
      if (!existing) {
        await plugin.runTask(
          state.updateTree(
            state
              .build()
              .to(model.cell)
              .apply(ModelUnitcell3D, style, { state: { isHidden: !visible } }),
          ),
        );
      } else {
        if (restyle) {
          await plugin.runTask(
            state.updateTree(
              state
                .build()
                .to(existing)
                .update({ ...current, ...style } as never),
            ),
          );
        }
        if (shown !== visible) {
          plugin.managers.structure.hierarchy.toggleVisibility([model.unitcell!], visible ? 'show' : 'hide');
        }
      }
      const ref = plugin.managers.structure.hierarchy.current.refs.get(model.cell.transform.ref);
      const cellRef = ref?.kind === 'model' ? ref.unitcell?.cell.transform.ref : existing?.transform.ref;
      const did = !visible
        ? `Hid the unit cell of ${name}.`
        : `${name} shows its unit cell (${symmetry.spacegroup.name})${restyle ? `, ${styleText(style)}` : ''}.`;
      return ApiResult.ok(did, { nodes: cellRef ? [cellRef] : [] });
    },
  });

  // parity.md §2: v5's Conformation Transform (a structure moved and rotated by hand). Its values are the transform from
  // where the structure was loaded, as v5's form has them; `compose` applies them on top of the current one instead.
  readonly transform = this.operation({
    id: 'structure.transform',
    title: 'Move or rotate a structure',
    params: {
      translation: OptionalParam(PD.Vec3(Vec3.create(0, 0, 0), {}, { description: 'Translation in Å: [x, y, z]' })),
      axis: OptionalParam(
        PD.Vec3(
          Vec3.create(0, 0, 1),
          {},
          { description: 'Rotation axis: [x, y, z], with an angle ([0, 0, 1] when absent)' },
        ),
      ),
      angle: OptionalParam(
        PD.Numeric(0, { min: -360, max: 360, step: 1 }, { description: 'Rotation angle in degrees' }),
      ),
      center: PD.Select<'centroid' | 'origin'>(
        'centroid',
        [
          ['centroid', 'The structure’s center'],
          ['origin', 'The origin'],
        ],
        { description: 'What the rotation turns about: the structure’s center or the origin' },
      ),
      matrix: OptionalParam(
        PD.Mat4(Mat4.identity(), {
          description: 'A 4×4 matrix (column-major, 16 numbers) instead of the axis, angle and translation',
        }),
      ),
      compose: PD.Boolean(false, { description: 'Apply on top of the current transform instead of replacing it' }),
      structure: StructureParam(),
    },
    async run(v, ctx) {
      const { center, compose, structure } = v;
      // What was given (api.md §1): a matrix, or parts (a translation, an angle about an axis); not both, not none.
      const {
        translation: t,
        axis: a,
        angle: g,
        matrix,
      } = v as {
        translation?: Vec3;
        axis?: Vec3;
        angle?: number;
        matrix?: Mat4;
      };
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const name = labelOf(ctx, s);
      const data = s.cell.obj?.data;
      if (!data) return ApiResult.fail('blocked', `${name} is not built yet.`);
      const byMatrix = matrix !== undefined;
      const byParts = t !== undefined || g !== undefined || a !== undefined;
      if (byMatrix && byParts) {
        return ApiResult.fail(
          'invalid',
          'Give a matrix, or an axis, angle and translation, not both.',
          'Leave the matrix out to move by parts.',
        );
      }
      if (!byMatrix && !byParts) {
        return ApiResult.fail(
          'invalid',
          'Move how?',
          'Give a translation, an angle (about an axis), or a matrix; structure.transform-reset puts it back.',
        );
      }
      if (a !== undefined && g === undefined) {
        return ApiResult.fail('invalid', 'An axis needs an angle.', 'Give the angle in degrees with the axis.');
      }
      const translation = t ?? Vec3.zero();
      const axis = a ?? Vec3.create(0, 0, 1);
      const angle = g ?? 0;
      if (!byMatrix && angle !== 0 && Vec3.magnitude(axis) < 1e-6) {
        return ApiResult.fail('invalid', 'The rotation axis is zero.', 'Give an axis such as [0, 0, 1].');
      }
      const wanted: TransformParam = byMatrix
        ? { name: 'matrix', params: { data: Mat4.clone(matrix!), transpose: false } }
        : {
            name: 'components',
            params: {
              translation: Vec3.clone(translation),
              axis: Vec3.normalize(Vec3(), Vec3.magnitude(axis) < 1e-6 ? Vec3.create(0, 0, 1) : axis),
              angle,
              rotationCenter:
                center === 'centroid'
                  ? { name: 'centroid', params: {} }
                  : { name: 'point', params: { point: Vec3.zero() } },
            },
          };
      const existing = s.transform?.cell;
      const now = currentTransform(s, data);
      const next = compose
        ? {
            name: 'matrix' as const,
            params: { data: Mat4.mul(Mat4(), matrixOf(wanted, data, now), now), transpose: false },
          }
        : wanted;
      if (Mat4.areEqual(matrixOf(next, data), now, 1e-3)) {
        return ApiResult.already(
          existing ? `${name} is already placed so.` : `${name} is already where it was loaded.`,
        );
      }
      const state = ctx.plugin.state.data;
      const build = state.build();
      if (existing) build.to(existing).update({ transform: next });
      else build.to(s.cell).insert(TransformStructureConformation, { transform: next });
      await ctx.plugin.runTask(state.updateTree(build));
      return ApiResult.ok(`Moved ${name}.`, { nodes: [s.cell.transform.ref] });
    },
  });

  readonly transformReset = this.operation({
    id: 'structure.transform-reset',
    title: 'Reset a structure’s transform',
    params: { structure: StructureParam() },
    async run({ structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const name = labelOf(ctx, s);
      const data = s.cell.obj?.data;
      const existing = s.transform?.cell;
      if (!existing || !data || Mat4.isIdentity(currentTransform(s, data), 1e-6)) {
        return ApiResult.already(`${name} is already where it was loaded.`);
      }
      // Removing the decorator would take its subtree (the components) with it: it is set to identity instead.
      const identity: TransformParam = {
        name: 'components',
        params: {
          translation: Vec3.zero(),
          axis: Vec3.create(0, 0, 1),
          angle: 0,
          rotationCenter: { name: 'centroid', params: {} },
        },
      };
      const state = ctx.plugin.state.data;
      await ctx.plugin.runTask(state.updateTree(state.build().to(existing).update({ transform: identity })));
      return ApiResult.ok(`${name} is back where it was loaded.`, { nodes: [s.cell.transform.ref] });
    },
  });

  readonly visibility = this.operation({
    id: 'structure.visibility',
    title: 'Show or hide a structure',
    params: { visible: PD.Boolean(true), structure: EveryStructureParam() },
    async run({ visible, structure }, ctx) {
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const change = targets.filter((s) => !!s.cell.state.isHidden === visible);
      if (!change.length) {
        return ApiResult.already(
          `${targets.map((x) => labelOf(ctx, x)).join(', ')} ${targets.length === 1 ? 'is' : 'are'} already ${visible ? 'shown' : 'hidden'}.`,
        );
      }
      ctx.plugin.managers.structure.hierarchy.toggleVisibility(change, visible ? 'show' : 'hide');
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${change.map((x) => labelOf(ctx, x)).join(', ')}.`, {
        nodes: change.map((s) => s.cell.transform.ref),
      });
    },
  });

  readonly remove = this.operation({
    id: 'structure.remove',
    title: 'Remove a structure',
    params: { structure: StructureParam('A structure by label or state ref; empty: the only one (several: name it)') },
    approval: 'destructive',
    async run({ structure }, ctx) {
      // Destructive: never the focused or selected one by default (targets.md §1).
      const s = oneStructure(ctx, structure || undefined, { destructive: true });
      if (isFailure(s)) return s;
      // The whole entry goes: the model, trajectory and data it came from too, unless something else hangs off them.
      await ctx.plugin.managers.structure.hierarchy.remove([entryRoot(ctx, s.cell.transform.ref)], false);
      return ApiResult.ok(`Removed ${labelOf(ctx, s)}.`, { nodes: [s.cell.transform.ref] });
    },
  });
}

/** The current hierarchy's refs of these structures (their components rebuilt). */
function freshStructures(ctx: ApiContext, list: readonly StructureRef[]): StructureRef[] {
  const refs = new Set(list.map((s) => s.cell.transform.ref));
  return ctx.plugin.managers.structure.hierarchy.current.structures.filter((s) => refs.has(s.cell.transform.ref));
}

/**
 * The matrix of a transform param, rotating about the structure's own center where it says so (as the decorator).
 * `placed` is a transform already on the structure: its center is where the center now is (a rotation composed on top
 * of it turns the structure in place).
 */
function matrixOf(t: TransformParam, s: Structure, placed?: Mat4): Mat4 {
  const center = placed ? Vec3.transformMat4(Vec3(), s.boundary.sphere.center, placed) : s.boundary.sphere.center;
  return getTransformFromParams(t, center);
}

/** A structure's conformation transform now (identity without one); `data` is the structure before it. */
export function currentTransform(s: StructureRef, data: Structure): Mat4 {
  const t = s.transform?.cell.transform.params?.transform;
  return t ? matrixOf(t as TransformParam, data) : Mat4.identity();
}

function styleText(style: { cellColor?: Color; alpha?: number }): string {
  const parts: string[] = [];
  if (style.cellColor !== undefined) parts.push(Color.toHexStyle(style.cellColor));
  if (style.alpha !== undefined) parts.push(`opacity ${style.alpha}`);
  return parts.join(', ');
}

function presetIds(ctx: ApiContext): string[] {
  return ctx.plugin.builders.structure.representation.providers.map((p) => p.id);
}

/** A structure type's params without its own (an assembly's id, symmetry's radius and cells): the common ones. */
function commonTypeParams(params: Record<string, unknown> | undefined): Record<string, unknown> {
  const { id: _, radius: __, ijkMin: ___, ijkMax: ____, ...common } = params ?? {};
  return common;
}

/**
 * A structure's params changed in a way that drops its components (v5's "Structure Type"): the default preset again,
 * with what the person hid left out of it (decision 18, as `structure.preset`).
 */
async function rebuildStructure(ctx: ApiContext, ref: string, params: unknown) {
  const plugin = ctx.plugin;
  const state = plugin.state.data;
  const root = StateTree.getDecoratorRoot(state.tree, ref);
  const children = state.tree.children.get(root).toArray();
  await plugin.managers.structure.hierarchy.remove(children, false);
  await plugin.state.updateTransform(state, ref, params, 'Structure');
  await plugin.builders.structure.representation.applyPreset(ref, defaultPreset(ctx) ?? 'auto');
  await reapplyContents(
    ctx,
    plugin.managers.structure.hierarchy.current.structures.filter((s) => s.cell.transform.ref === ref),
  );
}
