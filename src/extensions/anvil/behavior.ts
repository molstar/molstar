/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { MembraneOrientationProvider, MembraneOrientation } from './prop';
import { StateObjectRef, StateTransformer, StateTransform } from '../../mol-state';
import { Task } from '../../mol-task';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { MembraneOrientationRepresentationProvider, MembraneOrientationParams, MembraneOrientationRepresentation } from './representation';
import { HydrophobicityColorThemeProvider } from '../../mol-theme/color/hydrophobicity';
import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { DefaultQueryRuntimeTable } from '../../mol-script/runtime/query/compiler';
import { StructureSelectionQuery, StructureSelectionCategory } from '../../mol-plugin-state/helpers/structure-selection-query';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { GenericRepresentationRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';

const Tag = MembraneOrientation.Tag;

export const ANVILMembraneOrientation = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'anvil-membrane-orientation-prop',
    category: 'custom-props',
    display: {
        name: 'Membrane Orientation',
        description: 'Data calculated with ANVIL algorithm.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = MembraneOrientationProvider;

        register(): void {
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.representation.structure.registry.add(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.add(isTransmembrane);

            this.ctx.genericRepresentationControls.set(Tag.Representation, selection => {
                const refs: GenericRepresentationRef[] = [];
                selection.structures.forEach(structure => {
                    const memRepr = structure.genericRepresentations?.filter(r => r.cell.transform.transformer.id === MembraneOrientation3D.id)[0];
                    if (memRepr) refs.push(memRepr);
                });
                return [refs, 'Membrane Orientation'];
            });
            this.ctx.builders.structure.representation.registerPreset(MembraneOrientationPreset);
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

            this.ctx.representation.structure.registry.remove(MembraneOrientationRepresentationProvider);
            this.ctx.query.structure.registry.remove(isTransmembrane);

            this.ctx.genericRepresentationControls.delete(Tag.Representation);
            this.ctx.builders.structure.representation.unregisterPreset(MembraneOrientationPreset);
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
                'atom-test': MembraneOrientation.symbols.isTransmembrane.symbol(),
            })
        ])
    ])
]), {
    description: 'Select residues that are embedded between the membrane layers.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return MembraneOrientationProvider.attach(ctx, structure);
    }
});

//

export { MembraneOrientation3D };

type MembraneOrientation3D = typeof MembraneOrientation3D
const MembraneOrientation3D = PluginStateTransform.BuiltIn({
    name: 'membrane-orientation-3d',
    display: {
        name: 'Membrane Orientation',
        description: 'Membrane Orientation planes and rims. Data calculated with ANVIL algorithm.'
    },
    from: PluginStateObject.Molecule.Structure,
    to: PluginStateObject.Shape.Representation3D,
    params: (a) => {
        return {
            ...MembraneOrientationParams,
        };
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Membrane Orientation', async ctx => {
            await MembraneOrientationProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, a.data);
            const repr = MembraneOrientationRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => MembraneOrientationParams);
            await repr.createOrUpdate(params, a.data).runInContext(ctx);
            return new PluginStateObject.Shape.Representation3D({ repr, sourceData: a.data }, { label: 'Membrane Orientation' });
        });
    },
    update({ a, b, newParams }, plugin: PluginContext) {
        return Task.create('Membrane Orientation', async ctx => {
            await MembraneOrientationProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, a.data);
            const props = { ...b.data.repr.props, ...newParams };
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            b.data.sourceData = a.data;
            return StateTransformer.UpdateResult.Updated;
        });
    },
    isApplicable(a) {
        return MembraneOrientationProvider.isApplicable(a.data);
    }
});

export const MembraneOrientationPreset = StructureRepresentationPresetProvider({
    id: 'preset-membrane-orientation',
    display: {
        name: 'Membrane Orientation', group: 'Annotation',
        description: 'Shows orientation of membrane layers. Data calculated with ANVIL algorithm.' // TODO add ' or obtained via RCSB PDB'
    },
    isApplicable(a) {
        return MembraneOrientationProvider.isApplicable(a.data);
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        if (!MembraneOrientationProvider.get(structure).value) {
            await plugin.runTask(Task.create('Membrane Orientation', async runtime => {
                await MembraneOrientationProvider.attach({ runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, structure);
            }));
        }

        const membraneOrientation = await tryCreateMembraneOrientation(plugin, structureCell);
        const colorTheme = HydrophobicityColorThemeProvider.name as any;
        const preset = await PresetStructureRepresentations.auto.apply(ref, { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme } } }, plugin);

        return { components: preset.components, representations: { ...preset.representations, membraneOrientation } };
    }
});

export function tryCreateMembraneOrientation(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, params?: StateTransformer.Params<MembraneOrientation3D>, initialState?: Partial<StateTransform.State>) {
    const state = plugin.state.data;
    const membraneOrientation = state.build().to(structure)
        .applyOrUpdateTagged('membrane-orientation-3d', MembraneOrientation3D, params, { state: initialState });
    return membraneOrientation.commit({ revertOnError: true });
}
