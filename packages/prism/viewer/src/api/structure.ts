/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform, StateTree } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { StructureInfo } from '../sources/viewer-sources.js';
import { type ApiContext, ApiArea, ApiResult } from './area.js';
import { isFailure, labelOf, oneStructure, someStructures } from './scope.js';

// `structure` (operations.md §2, plan 2.2): the scene outline (`get`), representation presets by id, the model of a
// trajectory, the assembly, a structure's visibility and removal. The first scene area: the others follow its shape
// (a `structure` field naming one by ref or label, else the one the person works with; registries as arguments, their
// ids in the hint when one is wrong; "already so" when nothing would change).

const StructureParam = () =>
  PD.Text('', { description: 'A structure by label or state ref; empty: the one in view (or all, where it says)' });

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
      structure: StructureParam(),
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

  readonly visibility = this.operation({
    id: 'structure.visibility',
    title: 'Show or hide a structure',
    params: { visible: PD.Boolean(true), structure: StructureParam() },
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

/** The highest ancestor that holds only this structure (its model, trajectory, data), below the root. */
function entryRoot(ctx: ApiContext, ref: string): string {
  const tree = ctx.plugin.state.data.tree;
  let current = ref;
  for (;;) {
    const parent = tree.transforms.get(current)?.parent;
    if (!parent || parent === StateTransform.RootRef || tree.children.get(parent).size !== 1) return current;
    current = parent;
  }
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
