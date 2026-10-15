/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from '@molstar/core/state';
import { SpacegroupCell } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { StructureInfo } from '../sources/viewer-sources.js';
import { type ApiContext, ApiArea, ApiResult } from './area.js';
import { reapplyContents } from './contents.js';
import { EveryStructureParam, StructureParam } from './params.js';
import { entryRoot, isFailure, labelOf, oneStructure, someStructures } from './scope.js';

// `structure` (operations.md §2, plan 2.2): the scene outline (`get`), representation presets by id, the model of a
// trajectory, the assembly, a structure's visibility and removal. The first scene area: the others follow its shape
// (a `structure` field naming one by ref or label, else the one the person works with; registries as arguments, their
// ids in the hint when one is wrong; "already so" when nothing would change).

export class StructureApi extends ApiArea {
  readonly id = 'structure';

  /** The loaded structures. */
  get list(): ReadonlyState<readonly StructureInfo[]> {
    return this.ctx.sources.structures;
  }

  readonly get = this.operation({
    id: 'structure.get',
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

  readonly model = this.operation({
    id: 'structure.model',
    title: 'Show another model',
    params: {
      index: PD.Numeric(0, { min: 0, step: 1 }, { description: 'The model, from 0' }),
      structure: StructureParam(),
    },
    async run({ index, structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const model = s.model;
      const trajectory = model?.trajectory?.cell.obj?.data;
      if (!model || !trajectory) return ApiResult.fail('blocked', `${labelOf(ctx, s)} has no models to choose from.`);
      const count = trajectory.frameCount;
      if (index >= count) {
        return ApiResult.fail(
          'invalid',
          `${labelOf(ctx, s)} has ${count} model${count === 1 ? '' : 's'}.`,
          `Give 0 to ${count - 1}.`,
        );
      }
      const params = model.cell.transform.params as { modelIndex: number };
      if (params.modelIndex === index) return ApiResult.already(`${labelOf(ctx, s)} already shows model ${index + 1}.`);
      await ctx.plugin.state.updateTransform(ctx.plugin.state.data, model.cell.transform.ref, {
        ...params,
        modelIndex: index,
      });
      return ApiResult.ok(`${labelOf(ctx, s)} shows model ${index + 1} of ${count}.`, {
        nodes: [s.cell.transform.ref],
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
      const want = id ? { name: 'assembly', params: { id } } : { name: 'model', params: {} };
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

  // parity.md §2 (the structure's source): crystal symmetry, v5's Symmetry Mates and Symmetry (indices) types.
  readonly symmetry = this.operation({
    id: 'structure.symmetry',
    title: 'Show crystal symmetry',
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
        return ApiResult.fail('blocked', `${name} has no crystal cell.`, 'Only crystal structures have symmetry.');
      }
      const current = s.cell.transform.params as { type: { name: string; params?: Record<string, unknown> } };
      // The structure's common params (bond computation, …) carry over; the type's own are set.
      const { id: _, radius: __, ijkMin: ___, ijkMax: ____, ...common } = current.type.params ?? {};
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
    params: { structure: StructureParam() },
    approval: 'destructive',
    async run({ structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
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

function presetIds(ctx: ApiContext): string[] {
  return ctx.plugin.builders.structure.representation.providers.map((p) => p.id);
}

/** A structure's params changed in a way that drops its components (v5's "Structure Type"): the default preset again. */
async function rebuildStructure(ctx: ApiContext, ref: string, params: unknown) {
  const plugin = ctx.plugin;
  const state = plugin.state.data;
  const root = StateTree.getDecoratorRoot(state.tree, ref);
  const children = state.tree.children.get(root).toArray();
  await plugin.managers.structure.hierarchy.remove(children, false);
  await plugin.state.updateTransform(state, ref, params, 'Structure');
  await plugin.builders.structure.representation.applyPreset(ref, 'auto');
}
