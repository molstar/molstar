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
import { Loci, isDataLoci } from 'mol-model/loci';
import { OrderedSet } from 'mol-data/int';
import { Table } from 'mol-data/db';

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
            this.ctx.lociLabels.addProvider(labelAssemblySymmetryAxes);
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
            this.ctx.lociLabels.removeProvider(labelAssemblySymmetryAxes);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.remove('rcsb-assembly-symmetry-cluster')
            this.ctx.structureRepresentation.registry.remove('rcsb-assembly-symmetry-axes')
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});

function labelAssemblySymmetryAxes(loci: Loci): string | undefined {
    if (isDataLoci(loci) && AssemblySymmetry.is(loci.data) && loci.tag === 'axes') {
        const { rcsb_assembly_symmetry_axis: axis, rcsb_assembly_symmetry: sym } = loci.data.db
        const labels: string[] = []
        OrderedSet.forEach(loci.indices, v => {
            const symmetryId = axis.symmetry_id.value(v)
            const symmetry = Table.pickRow(sym, i => sym.id.value(i) === symmetryId)
            if (symmetry) {
                labels.push(`Axis of order ${axis.order.value(v)} for ${symmetry.kind} ${symmetry.type.toLowerCase()} symmetry`)
            }
        })
        return labels.length ? labels.join(', ') : undefined
    }
    return undefined
}