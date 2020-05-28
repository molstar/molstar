/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { AssemblySymmetryProvider, AssemblySymmetry, AssemblySymmetryDataProvider } from './prop';
import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { AssemblySymmetryParams, AssemblySymmetryRepresentation } from './representation';
import { AssemblySymmetryClusterColorThemeProvider } from './color';
import { PluginStateTransform, PluginStateObject } from '../../../mol-plugin-state/objects';
import { Task } from '../../../mol-task';
import { PluginContext } from '../../../mol-plugin/context';
import { StateTransformer, StateAction, StateObject, StateTransform, StateObjectRef } from '../../../mol-state';
import { GenericRepresentationRef } from '../../../mol-plugin-state/manager/structure/hierarchy-state';
import { AssemblySymmetryControls } from './ui';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../../mol-plugin-state/builder/structure/representation-preset';

const Tag = AssemblySymmetry.Tag;

export const RCSBAssemblySymmetry = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'rcsb-assembly-symmetry-prop',
    category: 'custom-props',
    display: {
        name: 'Assembly Symmetry',
        description: 'Assembly Symmetry data calculated with BioJava, obtained via RCSB PDB.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = AssemblySymmetryProvider

        register(): void {
            this.ctx.state.data.actions.add(InitAssemblySymmetry3D);
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(AssemblySymmetryClusterColorThemeProvider);

            this.ctx.genericRepresentationControls.set(Tag.Representation, selection => {
                const refs: GenericRepresentationRef[] = [];
                selection.structures.forEach(structure => {
                    const symmRepr = structure.genericRepresentations?.filter(r => r.cell.transform.transformer.id === AssemblySymmetry3D.id)[0];
                    if (symmRepr) refs.push(symmRepr);
                });
                return [refs, 'Symmetries'];
            });
            this.ctx.customStructureControls.set(Tag.Representation, AssemblySymmetryControls as any);
            this.ctx.builders.structure.representation.registerPreset(AssemblySymmetryPreset);
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            this.ctx.state.data.actions.remove(InitAssemblySymmetry3D);
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(AssemblySymmetryClusterColorThemeProvider);

            this.ctx.genericRepresentationControls.delete(Tag.Representation);
            this.ctx.customStructureControls.delete(Tag.Representation);
            this.ctx.builders.structure.representation.unregisterPreset(AssemblySymmetryPreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl)
    })
});

//

export const InitAssemblySymmetry3D = StateAction.build({
    display: {
        name: 'Assembly Symmetry',
        description: 'Initialize Assembly Symmetry axes and cage. Data calculated with BioJava, obtained via RCSB PDB.'
    },
    from: PluginStateObject.Molecule.Structure,
    isApplicable: (a) => AssemblySymmetry.isApplicable(a.data)
})(({ a, ref, state }, plugin: PluginContext) => Task.create('Init Assembly Symmetry', async ctx => {
    try {
        const propCtx = { runtime: ctx, assetManager: plugin.managers.asset };
        await AssemblySymmetryDataProvider.attach(propCtx, a.data);
        const assemblySymmetryData = AssemblySymmetryDataProvider.get(a.data).value;
        const symmetryIndex = assemblySymmetryData ? AssemblySymmetry.firstNonC1(assemblySymmetryData) : -1;
        await AssemblySymmetryProvider.attach(propCtx, a.data, { symmetryIndex });
    } catch(e) {
        plugin.log.error(`Assembly Symmetry: ${e}`);
        return;
    }
    const tree = state.build().to(ref)
        .applyOrUpdateTagged(AssemblySymmetry.Tag.Representation, AssemblySymmetry3D);
    await state.updateTree(tree).runInContext(ctx);
}));

export { AssemblySymmetry3D };

type AssemblySymmetry3D = typeof AssemblySymmetry3D
const AssemblySymmetry3D = PluginStateTransform.BuiltIn({
    name: Tag.Representation,
    display: {
        name: 'Assembly Symmetry',
        description: 'Assembly Symmetry axes and cage. Data calculated with BioJava, obtained via RCSB PDB.'
    },
    from: PluginStateObject.Molecule.Structure,
    to: PluginStateObject.Shape.Representation3D,
    params: (a) => {
        return {
            ...AssemblySymmetryParams,
        };
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Assembly Symmetry', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset }, a.data);
            const assemblySymmetry = AssemblySymmetryProvider.get(a.data).value;
            if (!assemblySymmetry || assemblySymmetry.symbol === 'C1') {
                return StateObject.Null;
            }
            const repr = AssemblySymmetryRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => AssemblySymmetryParams);
            await repr.createOrUpdate(params, a.data).runInContext(ctx);
            const { type, kind, symbol } = assemblySymmetry;
            return new PluginStateObject.Shape.Representation3D({ repr, source: a }, { label: kind, description: `${type} (${symbol})` });
        });
    },
    update({ a, b, newParams }, plugin: PluginContext) {
        return Task.create('Assembly Symmetry', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, assetManager: plugin.managers.asset }, a.data);
            const assemblySymmetry = AssemblySymmetryProvider.get(a.data).value;
            if (!assemblySymmetry || assemblySymmetry.symbol === 'C1') {
                // this should NOT be StateTransformer.UpdateResult.Null
                // because that keeps the old object
                return StateTransformer.UpdateResult.Recreate;
            }
            const props = { ...b.data.repr.props, ...newParams };
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            const { type, kind, symbol } = assemblySymmetry;
            b.label = kind;
            b.description = `${type} (${symbol})`;
            return StateTransformer.UpdateResult.Updated;
        });
    },
    isApplicable(a) {
        return AssemblySymmetry.isApplicable(a.data);
    }
});

//

export const AssemblySymmetryPresetParams = {
    ...StructureRepresentationPresetProvider.CommonParams,
};

export const AssemblySymmetryPreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-rcsb-assembly-symmetry',
    display: {
        name: 'Assembly Symmetry', group: 'Annotation',
        description: 'Shows Assembly Symmetry axes and cage; colors structure according to assembly symmetry cluster membership. Data calculated with BioJava, obtained via RCSB PDB.'
    },
    isApplicable(a) {
        return AssemblySymmetry.isApplicable(a.data);
    },
    params: () => AssemblySymmetryPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        if (!AssemblySymmetryDataProvider.get(structure).value) {
            await plugin.runTask(Task.create('Assembly Symmetry', async runtime => {
                const propCtx = { runtime, assetManager: plugin.managers.asset };
                await AssemblySymmetryDataProvider.attach(propCtx, structure);
                const assemblySymmetryData = AssemblySymmetryDataProvider.get(structure).value;
                const symmetryIndex = assemblySymmetryData ? AssemblySymmetry.firstNonC1(assemblySymmetryData) : -1;
                await AssemblySymmetryProvider.attach(propCtx, structure, { symmetryIndex });
            }));
        }

        const assemblySymmetry = await tryCreateAssemblySymmetry(plugin, structureCell);
        const colorTheme = assemblySymmetry.isOk ? Tag.Cluster as any : undefined;
        const preset = await PresetStructureRepresentations.auto.apply(ref, { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme } } }, plugin);

        return { components: preset.components, representations: { ...preset.representations, assemblySymmetry } };
    }
});

export function tryCreateAssemblySymmetry(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, params?: StateTransformer.Params<AssemblySymmetry3D>, initialState?: Partial<StateTransform.State>) {
    const state = plugin.state.data;
    const assemblySymmetry = state.build().to(structure)
        .applyOrUpdateTagged(AssemblySymmetry.Tag.Representation, AssemblySymmetry3D, params, { state: initialState });
    return assemblySymmetry.commit({ revertOnError: true });
}