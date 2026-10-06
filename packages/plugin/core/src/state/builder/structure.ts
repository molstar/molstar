/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import {
  StateObjectRef,
  StateObjectSelector,
  StateTransformer,
  StateTransform,
  StateObjectCell,
} from '@molstar/core/state';
import { PluginStateObject as SO } from '../objects.js';
import { ParseBlob } from '@molstar/plugin/state/formats/cif';
import { TrajectoryFromBlob } from '@molstar/plugin/state/formats/trajectory/mmcif';
import {
  CustomModelProperties,
  CustomStructureProperties,
  ModelFromTrajectory,
  StructureFromModel,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { ModelUnitcell3D } from '@molstar/plugin/state/transforms/structure/unitcell';
import { StructureComponent } from '@molstar/plugin/state/transforms/structure/selection';
import type { RootStructureDefinition } from '../helpers/root-structure.js';
import type { StructureComponentParams, StaticStructureComponentType } from '../helpers/structure-component.js';
import type { BuiltInTrajectoryFormat } from '@molstar/plugin/state/formats/trajectory/catalog';
import type { TrajectoryFormatProvider } from '@molstar/plugin/state/formats/trajectory/provider';
import { StructureRepresentationBuilder } from './structure/representation.js';
import type { StructureSelectionQuery } from '@molstar/plugin/state/queries/structure/query';
import { Task } from '@molstar/core/task';
import { StructureElement } from '@molstar/model/model/structure';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import { SpacegroupCell } from '@molstar/core/math/geometry';
import type { Expression } from '@molstar/model/script/language/expression';
import { TrajectoryHierarchyBuilder } from './structure/hierarchy.js';

export class StructureBuilder {
  private get dataState() {
    return this.plugin.state.data;
  }

  private async parseTrajectoryData(
    data: StateObjectRef<SO.Data.Binary | SO.Data.String>,
    format: BuiltInTrajectoryFormat | TrajectoryFormatProvider,
  ) {
    const provider =
      typeof format === 'string' ? (this.plugin.dataFormats.get(format) as TrajectoryFormatProvider) : format;
    if (!provider) throw new Error(`'${format}' is not a supported data format.`);
    const { trajectory } = await provider.parse(this.plugin, data);
    return trajectory;
  }

  private parseTrajectoryBlob(data: StateObjectRef<SO.Data.Blob>, params: StateTransformer.Params<typeof ParseBlob>) {
    const state = this.dataState;
    const trajectory = state
      .build()
      .to(data)
      .apply(ParseBlob, params, { state: { isGhost: true } })
      .apply(TrajectoryFromBlob, void 0);
    return trajectory.commit({ revertOnError: true });
  }

  readonly hierarchy = new TrajectoryHierarchyBuilder(this.plugin);
  readonly representation = new StructureRepresentationBuilder(this.plugin);

  parseTrajectory(
    data: StateObjectRef<SO.Data.Binary | SO.Data.String>,
    format: BuiltInTrajectoryFormat | TrajectoryFormatProvider,
  ): Promise<StateObjectSelector<SO.Molecule.Trajectory>>;
  parseTrajectory(
    blob: StateObjectRef<SO.Data.Blob>,
    params: StateTransformer.Params<typeof ParseBlob>,
  ): Promise<StateObjectSelector<SO.Molecule.Trajectory>>;
  parseTrajectory(data: StateObjectRef, params: any) {
    const cell = StateObjectRef.resolveAndCheck(this.dataState, data as StateObjectRef);
    if (!cell) throw new Error('Invalid data cell.');

    if (SO.Data.Blob.is(cell.obj)) {
      return this.parseTrajectoryBlob(data, params);
    } else {
      return this.parseTrajectoryData(data, params);
    }
  }

  createModel(
    trajectory: StateObjectRef<SO.Molecule.Trajectory>,
    params?: StateTransformer.Params<typeof ModelFromTrajectory>,
    initialState?: Partial<StateTransform.State>,
  ) {
    const state = this.dataState;
    const model = state
      .build()
      .to(trajectory)
      .apply(ModelFromTrajectory, params || { modelIndex: 0 }, { state: initialState });

    return model.commit({ revertOnError: true });
  }

  insertModelProperties(
    model: StateObjectRef<SO.Molecule.Model>,
    params?: StateTransformer.Params<typeof CustomModelProperties>,
    initialState?: Partial<StateTransform.State>,
  ) {
    const state = this.dataState;
    const props = state.build().to(model).apply(CustomModelProperties, params, { state: initialState });
    return props.commit({ revertOnError: true });
  }

  tryCreateUnitcell(
    model: StateObjectRef<SO.Molecule.Model>,
    params?: StateTransformer.Params<typeof ModelUnitcell3D>,
    initialState?: Partial<StateTransform.State>,
  ) {
    const state = this.dataState;
    const m = StateObjectRef.resolveAndCheck(state, model)?.obj?.data;
    if (!m) return;
    const cell = ModelSymmetry.Provider.get(m)?.spacegroup.cell;
    if (SpacegroupCell.isZero(cell)) return;

    const unitcell = state.build().to(model).apply(ModelUnitcell3D, params, { state: initialState });
    return unitcell.commit({ revertOnError: true });
  }

  createStructure(
    modelRef: StateObjectRef<SO.Molecule.Model>,
    params?: RootStructureDefinition.Params,
    initialState?: Partial<StateTransform.State>,
    tags?: string | string[],
  ) {
    const state = this.dataState;

    if (!params) {
      const model = StateObjectRef.resolveAndCheck(state, modelRef);
      if (model) {
        const symm = ModelSymmetry.Provider.get(model.obj?.data!);
        if (!symm || symm?.assemblies.length === 0) params = { name: 'model', params: {} };
      }
    }

    const structure = state
      .build()
      .to(modelRef)
      .apply(StructureFromModel, { type: params || { name: 'assembly', params: {} } }, { state: initialState, tags });

    return structure.commit({ revertOnError: true });
  }

  insertStructureProperties(
    structure: StateObjectRef<SO.Molecule.Structure>,
    params?: StateTransformer.Params<typeof CustomStructureProperties>,
  ) {
    const state = this.dataState;
    const props = state.build().to(structure).apply(CustomStructureProperties, params);
    return props.commit({ revertOnError: true });
  }

  isComponentTransform(cell: StateObjectCell) {
    return cell.transform.transformer === StructureComponent;
  }

  /** returns undefined if the component is empty/null */
  async tryCreateComponent(
    structure: StateObjectRef<SO.Molecule.Structure>,
    params: StructureComponentParams,
    key: string,
    tags?: string[],
  ): Promise<StateObjectSelector<SO.Molecule.Structure> | undefined> {
    const state = this.dataState;

    const root = state.build().to(structure);

    const keyTag = `structure-component-${key}`;
    const component = root.applyOrUpdateTagged(keyTag, StructureComponent, params, {
      tags: tags ? [...tags, keyTag] : [keyTag],
    });

    await component.commit();

    const selector = component.selector;

    if (!selector.isOk || selector.cell?.obj?.data.elementCount === 0) {
      await state.build().delete(selector.ref).commit();
      return;
    }

    return selector;
  }

  tryCreateComponentFromExpression(
    structure: StateObjectRef<SO.Molecule.Structure>,
    expression: Expression,
    key: string,
    params?: { label?: string; tags?: string[] },
  ) {
    return this.tryCreateComponent(
      structure,
      {
        type: { name: 'expression', params: expression },
        nullIfEmpty: true,
        label: (params?.label || '').trim(),
      },
      key,
      params?.tags,
    );
  }

  tryCreateComponentStatic(
    structure: StateObjectRef<SO.Molecule.Structure>,
    type: StaticStructureComponentType,
    params?: { label?: string; tags?: string[] },
  ) {
    return this.tryCreateComponent(
      structure,
      {
        type: { name: 'static', params: type },
        nullIfEmpty: true,
        label: (params?.label || '').trim(),
      },
      `static-${type}`,
      params?.tags,
    );
  }

  tryCreateComponentFromSelection(
    structure: StateObjectRef<SO.Molecule.Structure>,
    selection: StructureSelectionQuery,
    key: string,
    params?: { label?: string; tags?: string[] },
  ): Promise<StateObjectSelector<SO.Molecule.Structure> | undefined> {
    return this.plugin.runTask(
      Task.create('Query Component', async (taskCtx) => {
        let { label, tags } = params || {};
        label = (label || '').trim();

        const structureData = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;

        if (!structureData) return;

        const transformParams: StructureComponentParams = selection.referencesCurrent
          ? {
              type: {
                name: 'bundle',
                params: StructureElement.Bundle.fromSelection(
                  await selection.getSelection(this.plugin, taskCtx, structureData),
                ),
              },
              nullIfEmpty: true,
              label: label || selection.label,
            }
          : {
              type: { name: 'expression', params: selection.expression },
              nullIfEmpty: true,
              label: label || selection.label,
            };

        if (selection.ensureCustomProperties) {
          await selection.ensureCustomProperties(
            { runtime: taskCtx, assetManager: this.plugin.managers.asset, errorContext: this.plugin.errorContext },
            structureData,
          );
        }

        return this.tryCreateComponent(structure, transformParams, key, tags);
      }),
    );
  }

  constructor(public plugin: PluginContext) {}
}
