/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from 'mol-plugin/behavior';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { AssemblySymmetry } from 'mol-model-props/rcsb/assembly-symmetry';
import { CustomPropertyRegistry } from 'mol-plugin/util/custom-prop-registry';
import { AssemblySymmetryClusterColorThemeProvider } from 'mol-model-props/rcsb/themes/assembly-symmetry-cluster';
import { AssemblySymmetryAxesRepresentationProvider } from 'mol-model-props/rcsb/representations/assembly-symmetry-axes';

export const RCSBAssemblySymmetry = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'rcsb-assembly-symmetry-prop',
    display: { name: 'RCSB Assembly Symmetry', group: 'Custom Props' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private attach = AssemblySymmetry.createAttachTask(this.ctx.fetch);

        private provider: CustomPropertyRegistry.Provider = {
            option: [AssemblySymmetry.Descriptor.name, 'RCSB Assembly Symmetry'],
            descriptor: AssemblySymmetry.Descriptor,
            defaultSelected: this.params.autoAttach,
            attachableTo: () => true,
            attach: this.attach
        }

        register(): void {
            this.ctx.customModelProperties.register(this.provider);

            // TODO: support filtering of themes and representations based on the input structure
            // in this case, it would check structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.add('rcsb-assembly-symmetry-cluster', AssemblySymmetryClusterColorThemeProvider)
            this.ctx.structureRepresentation.registry.add('rcsb-assembly-symmetry-axes', AssemblySymmetryAxesRepresentationProvider)
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach
            this.params.autoAttach = p.autoAttach;
            this.provider.defaultSelected = p.autoAttach;
            return updated;
        }

        unregister() {
            this.ctx.customModelProperties.unregister(AssemblySymmetry.Descriptor.name);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.remove('rcsb-assembly-symmetry-cluster')
            this.ctx.structureRepresentation.registry.remove('rcsb-assembly-symmetry-axes')
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});