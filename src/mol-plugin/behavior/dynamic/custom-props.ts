/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { CustomPropertyRegistry } from 'mol-plugin/util/custom-prop-registry';
import { Loci } from 'mol-model/loci';
import { StructureElement } from 'mol-model/structure';
import { OrderedSet } from 'mol-data/int';

// TODO: make auto attach working better for "initial state" by supporting default props in state updates

export const PDBeStructureQualityReport = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'pdbe-structure-quality-report-prop',
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private attach = StructureQualityReport.createAttachTask(
            m => `https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${m.label.toLowerCase()}`,
            this.ctx.fetch
        );

        private provider: CustomPropertyRegistry.Provider = {
            option: [StructureQualityReport.Descriptor.name, 'PDBe Structure Quality Report'],
            descriptor: StructureQualityReport.Descriptor,
            defaultSelected: false,
            attachableTo: () => true,
            attach: this.attach
        }

        register(): void {
            this.ctx.customModelProperties.register(this.provider);
            this.ctx.lociLabels.addProvider(labelPDBeValidation);
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach
            this.params.autoAttach = p.autoAttach;
            this.provider.defaultSelected = p.autoAttach;
            return updated;
        }

        unregister() {
            this.ctx.customModelProperties.unregister(StructureQualityReport.Descriptor.name);
            this.ctx.lociLabels.removeProvider(labelPDBeValidation);
        }
    },
    params: () => ({
        autoAttach: ParamDefinition.Boolean(false)
    }),
    display: { name: 'Focus Loci on Select', group: 'Camera' }
});

function labelPDBeValidation(loci: Loci): string | undefined {
    switch (loci.kind) {
        case 'element-loci':
            const e = loci.elements[0];
            const u = e.unit;
            if (!u.model.customProperties.has(StructureQualityReport.Descriptor)) return void 0;

            const se = StructureElement.create(u, u.elements[OrderedSet.getAt(e.indices, 0)]);
            const issues = StructureQualityReport.getIssues(se);
            if (issues.length === 0) return 'PDBe Validation: No Issues';
            return `PDBe Validation: ${issues.join(', ')}`;

        default: return void 0;
    }
}