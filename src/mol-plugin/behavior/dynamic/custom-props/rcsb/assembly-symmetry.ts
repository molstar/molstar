/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../../mol-util/param-definition'
import { AssemblySymmetryProvider, AssemblySymmetry, getSymmetrySelectParam } from '../../../../../mol-model-props/rcsb/assembly-symmetry';
import { PluginBehavior } from '../../../behavior';
import { getAssemblySymmetryAxesRepresentation, AssemblySymmetryAxesParams } from '../../../../../mol-model-props/rcsb/representations/assembly-symmetry-axes';
import { AssemblySymmetryClusterColorThemeProvider } from '../../../../../mol-model-props/rcsb/themes/assembly-symmetry-cluster';
import { PluginStateTransform, PluginStateObject } from '../../../../state/objects';
import { Task } from '../../../../../mol-task';
import { PluginSpec } from '../../../../spec';
import { PluginContext } from '../../../../context';
import { StateTransformer } from '../../../../../mol-state';

export const RCSBAssemblySymmetry = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'rcsb-assembly-symmetry-prop',
    category: 'custom-props',
    display: { name: 'RCSB Assembly Symmetry' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = AssemblySymmetryProvider

        register(): void {
            this.ctx.state.dataState.actions.add(PluginSpec.Action(AssemblySymmetryAxes3D).action)
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.add('rcsb-assembly-symmetry-cluster', AssemblySymmetryClusterColorThemeProvider)
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            // TODO remove `AssemblySymmetryAxes3D` from `this.ctx.state.dataState.actions`
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.remove('rcsb-assembly-symmetry-cluster')
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl)
    })
});

type AssemblySymmetryAxes3D = typeof AssemblySymmetryAxes3D
const AssemblySymmetryAxes3D = PluginStateTransform.BuiltIn({
    name: 'rcsb-assembly-symmetry-axes-3d',
    display: 'RCSB Assembly Symmetry Axes',
    from: PluginStateObject.Molecule.Structure,
    to: PluginStateObject.Shape.Representation3D,
    params: (a, ctx: PluginContext) => {
        return {
            ...AssemblySymmetryAxesParams,
            symmetryIndex: getSymmetrySelectParam(a?.data),
        }
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('RCSB Assembly Symmetry Axes', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, fetch: plugin.fetch }, a.data)
            const repr = await getAssemblySymmetryAxesRepresentation(ctx, a.data, params)
            const { symbol, kind } = AssemblySymmetryProvider.get(a.data).value![params.symmetryIndex]
            return new PluginStateObject.Shape.Representation3D({ repr, source: a }, { label: `Axes`, description: `${symbol} ${kind}` });
        });
    },
    update({ a, b, newParams }, plugin: PluginContext) {
        return Task.create('RCSB Assembly Symmetry Axes', async ctx => {
            await AssemblySymmetryProvider.attach({ runtime: ctx, fetch: plugin.fetch }, a.data)
            await getAssemblySymmetryAxesRepresentation(ctx, a.data, newParams, b.data.repr);
            const { symbol, kind } = AssemblySymmetryProvider.get(a.data).value![newParams.symmetryIndex]
            b.description = `${symbol} ${kind}`
            return StateTransformer.UpdateResult.Updated;
        });
    },
    isApplicable(a) {
        return AssemblySymmetry.isApplicable(a.data)
    }
});