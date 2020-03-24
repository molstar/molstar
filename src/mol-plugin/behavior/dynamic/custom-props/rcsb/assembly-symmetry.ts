/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../../mol-util/param-definition'
import { AssemblySymmetryProvider, AssemblySymmetry, getSymmetrySelectParam } from '../../../../../mol-model-props/rcsb/assembly-symmetry';
import { PluginBehavior } from '../../../behavior';
import { AssemblySymmetryParams, AssemblySymmetryRepresentation } from '../../../../../mol-model-props/rcsb/representations/assembly-symmetry';
import { AssemblySymmetryClusterColorThemeProvider } from '../../../../../mol-model-props/rcsb/themes/assembly-symmetry-cluster';
import { PluginStateTransform, PluginStateObject } from '../../../../../mol-plugin-state/objects';
import { Task } from '../../../../../mol-task';
import { PluginContext } from '../../../../context';
import { StateTransformer, StateAction, StateObject } from '../../../../../mol-state';
import { GenericRepresentationRef } from '../../../../../mol-plugin-state/manager/structure/hierarchy-state';

const Tag = AssemblySymmetry.Tag

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
            this.ctx.state.data.actions.add(InitAssemblySymmetry3D)
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(AssemblySymmetryClusterColorThemeProvider)

            this.ctx.customSourceControls.set(Tag.Representation, selection => {
                const refs: GenericRepresentationRef[] = []
                selection.structures.forEach(structure => {
                    const symmRepr = structure.genericRepresentations?.filter(r => r.cell.transform.transformer.id === AssemblySymmetry3D.id)[0]
                    if (symmRepr) refs.push(symmRepr)
                })
                return [refs, 'Symmetries']
            })

            // TODO this should probably be done via a hierarchy preset
            this.subscribeObservable(this.ctx.managers.structure.hierarchy.behaviors.selection, selection => {
                for (const structure of selection.structures) {
                    const symmRepr = structure.genericRepresentations?.filter(r => r.cell.transform.transformer.id === AssemblySymmetry3D.id)[0]

                    if (!symmRepr) {
                        const state = this.ctx.state.data;
                        const symmReprBuilder = state.build().to(structure.cell.transform.ref)
                            .apply(AssemblySymmetry3D, { symmetryIndex: 0 }, { state: {} });
                        this.ctx.updateDataState(symmReprBuilder, { revertOnError: true });
                    }
                }
            })
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            this.ctx.state.data.actions.remove(InitAssemblySymmetry3D)
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(AssemblySymmetryClusterColorThemeProvider)
            this.ctx.customSourceControls.delete(Tag.Representation)
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl)
    })
});

//

const InitAssemblySymmetry3D = StateAction.build({
    display: {
        name: 'Assembly Symmetry',
        description: 'Initialize Assembly Symmetry axes and cage. Data calculated with BioJava, obtained via RCSB PDB.'
    },
    from: PluginStateObject.Molecule.Structure,
    isApplicable: (a) => AssemblySymmetry.isApplicable(a.data)
})(({ a, ref, state }, plugin: PluginContext) => Task.create('Init Assembly Symmetry', async ctx => {
    try {
        await AssemblySymmetryProvider.attach({ runtime: ctx, fetch: plugin.fetch }, a.data)
    } catch(e) {
        plugin.log.error(`Assembly Symmetry: ${e}`)
        return
    }
    const tree = state.build().to(ref).apply(AssemblySymmetry3D);
    await state.updateTree(tree).runInContext(ctx);
}));

export { AssemblySymmetry3D }

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
            symmetryIndex: getSymmetrySelectParam(a?.data),
        }
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Assembly Symmetry', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, fetch: plugin.fetch }, a.data)
            const assemblySymmetry = AssemblySymmetryProvider.get(a.data).value
            if (!assemblySymmetry || assemblySymmetry.length === 0) {
                return StateObject.Null;
            }
            const repr = AssemblySymmetryRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => AssemblySymmetryParams)
            await repr.createOrUpdate(params, a.data).runInContext(ctx);
            const { type, kind, symbol } = assemblySymmetry![params.symmetryIndex]
            return new PluginStateObject.Shape.Representation3D({ repr, source: a }, { label: kind, description: `${type} (${symbol})` });
        });
    },
    update({ a, b, newParams }, plugin: PluginContext) {
        return Task.create('Assembly Symmetry', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, fetch: plugin.fetch }, a.data)
            const assemblySymmetry = AssemblySymmetryProvider.get(a.data).value
            if (!assemblySymmetry || assemblySymmetry.length === 0) {
                return StateTransformer.UpdateResult.Recreate
            }
            const props = { ...b.data.repr.props, ...newParams }
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            const { type, kind, symbol } = assemblySymmetry![newParams.symmetryIndex]
            b.label = kind
            b.description = `${type} (${symbol})`
            return StateTransformer.UpdateResult.Updated;
        });
    },
    isApplicable(a) {
        return AssemblySymmetry.isApplicable(a.data)
    }
});