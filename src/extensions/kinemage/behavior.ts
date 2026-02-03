/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { KinemageDataProvider, KinemageData } from './prop';
import { StateObjectRef, StateTransformer, StateTransform } from '../../mol-state';
import { Task } from '../../mol-task';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginDragAndDropHandler } from '../../mol-plugin-state/manager/drag-and-drop';
import { KinemageDataRepresentationProvider, KinemageDataParams, KinemageDataRepresentation } from './representation';
import { HydrophobicityColorThemeProvider } from '../../mol-theme/color/hydrophobicity';
import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { DefaultQueryRuntimeTable } from '../../mol-script/runtime/query/compiler';
import { StructureSelectionQuery, StructureSelectionCategory } from '../../mol-plugin-state/helpers/structure-selection-query';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { GenericRepresentationRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { Vec3 } from '../../mol-math/linear-algebra';

const Tag = KinemageData.Tag;

/// @todo Remove once we store into the proper location in the drag and drop handler
/** Global KinemageInfo that is used to display */
let g_kinemageInfo: KinemageData = {
  kinemages: [], activeKinemage: -1
  /// @todo Remove these when we no longer need them.
  , planePoint1: Vec3(), planePoint2: Vec3(), normalVector: Vec3(), radius: 0, centroid: Vec3()
};

export const KinemageExtension = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'kinemage-data-prop',
    category: 'custom-props',
    display: {
        name: 'Kinemage data',
        description: 'Data loaded from Kinemage.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = KinemageDataProvider;

        register(): void {
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.representation.structure.registry.add(KinemageDataRepresentationProvider);
            this.ctx.query.structure.registry.add(isTransmembrane);

            this.ctx.genericRepresentationControls.set(Tag.Representation, selection => {
                const refs: GenericRepresentationRef[] = [];
                selection.structures.forEach(structure => {
                    const memRepr = structure.genericRepresentations?.filter(r => r.cell.transform.transformer.id === KinemageData3D.id)[0];
                    if (memRepr) refs.push(memRepr);
                });
                return [refs, 'Membrane Orientation'];
            });
            this.ctx.builders.structure.representation.registerPreset(KinemageDataPreset);

            this.ctx.managers.dragAndDrop.addHandler(KinemageDragAndDropHandler.name, KinemageDragAndDropHandler.handle);
        }

        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            DefaultQueryRuntimeTable.removeCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

            this.ctx.representation.structure.registry.remove(KinemageDataRepresentationProvider);
            this.ctx.query.structure.registry.remove(isTransmembrane);

            this.ctx.genericRepresentationControls.delete(Tag.Representation);
            this.ctx.builders.structure.representation.unregisterPreset(KinemageDataPreset);

            this.ctx.managers.dragAndDrop.removeHandler(KinemageDragAndDropHandler.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});

//

export const isTransmembrane = StructureSelectionQuery('Residues Embedded in Membrane', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': KinemageData.symbols.isTransmembrane.symbol(),
            })
        ])
    ])
]), {
    description: 'Select residues that are embedded between the membrane layers.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return KinemageDataProvider.attach(ctx, structure);
    }
});

//

export { KinemageData3D };

type KinemageData3D = typeof KinemageData3D
const KinemageData3D = PluginStateTransform.BuiltIn({
    name: 'kinemage-3d',
    display: {
        name: 'Kinemage 3D Data',
        description: '3D Data loaded from Kinemage.'
    },
    from: PluginStateObject.Molecule.Structure,
    to: PluginStateObject.Shape.Representation3D,
    params: (a) => {
        return {
          ...KinemageDataParams,
        };
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Membrane Orientation', async ctx => {
            await KinemageDataProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, a.data);
            const repr = KinemageDataRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => KinemageDataParams);
            await repr.createOrUpdate(params, a.data).runInContext(ctx);
            return new PluginStateObject.Shape.Representation3D({ repr, sourceData: a.data }, { label: 'Membrane Orientation' });
        });
    },
    update({ a, b, newParams }, plugin: PluginContext) {
        return Task.create('Membrane Orientation', async ctx => {
            await KinemageDataProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, a.data);
            const props = { ...b.data.repr.props, ...newParams };
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            b.data.sourceData = a.data;
            return StateTransformer.UpdateResult.Updated;
        });
    },
    isApplicable(a) {
        return KinemageDataProvider.isApplicable(a.data);
    }
});

export const KinemageDataPreset = StructureRepresentationPresetProvider({
    id: 'preset-kinemage',
    display: {
        name: 'Kinemage data', group: 'Annotation',
        description: 'Shows data loaded from Kinemage file.'
    },
    isApplicable(a) {
        return KinemageDataProvider.isApplicable(a.data);
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        if (!KinemageDataProvider.get(structure).value) {
            await plugin.runTask(Task.create('Membrane Orientation', async runtime => {
                await KinemageDataProvider.attach({ runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, structure);
            }));
        }

        const KinemageData = await tryCreateKinemageData(plugin, structureCell);
        const colorTheme = HydrophobicityColorThemeProvider.name as any;
        const preset = await PresetStructureRepresentations.auto.apply(ref, { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme } } }, plugin);

        return { components: preset.components, representations: { ...preset.representations, KinemageData } };
    }
});

export function tryCreateKinemageData(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, params?: StateTransformer.Params<KinemageData3D>, initialState?: Partial<StateTransform.State>) {
    const state = plugin.state.data;
    const KinemageData = state.build().to(structure)
        .applyOrUpdateTagged('kinemage-3d', KinemageData3D, params, { state: initialState });
    return KinemageData.commit({ revertOnError: true });
}

/** Registerable method for handling dragged-and-dropped files */
interface DragAndDropHandler {
  name: string,
  handle: PluginDragAndDropHandler,
}

/** DragAndDropHandler handler for `.kin` files */
const KinemageDragAndDropHandler: DragAndDropHandler = {
  name: 'kin',
  /** Load .kin files. Append to previous plugin state.
   * If multiple files are provided, append them all.
   * Select the last-loaded one from the list.
   * Return `true` if at least one file has been loaded. */
  async handle(files: File[], plugin: PluginContext): Promise<boolean> {
    let applied = false;
    for (const file of files) {
      if (file.name.toLowerCase().endsWith('.kin')) {
        const task = Task.create('Load KIN file', async ctx => {
          console.log('XXX loading KIN file ', file.name);  /// @todo Remove when done debugging
          const kinInfo = await KinemageData.open(file);
          for (const kinData of kinInfo.kinemages) {
            g_kinemageInfo.kinemages.push(kinData);
            g_kinemageInfo.activeKinemage = g_kinemageInfo.kinemages.length - 1;
          }
          console.log('XXX accumulated Kinemages size ', g_kinemageInfo.kinemages.length, ', active is ', g_kinemageInfo.activeKinemage);
          /// @todo See what loadMVS() ... LoadMolstarTree() ... MolstarLoadingActions().primitives() ... applyPrimitiveVisuals() does
        });
        console.log('XXX plugin.runTask');
        await plugin.runTask(task);
        applied = true;
      }
    }
    console.log('XXX KinemageDragAndDropHandler applied=', applied);
    return applied;
  },
};
