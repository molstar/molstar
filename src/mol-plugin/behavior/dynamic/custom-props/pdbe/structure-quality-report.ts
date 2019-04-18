/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { StructureQualityReportColorTheme } from 'mol-model-props/pdbe/themes/structure-quality-report';
import { Loci } from 'mol-model/loci';
import { StructureElement } from 'mol-model/structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../../../behavior';
import { ThemeDataContext } from 'mol-theme/theme';
import { CustomPropertyRegistry } from 'mol-model-props/common/custom-property-registry';

export const PDBeStructureQualityReport = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'pdbe-structure-quality-report-prop',
    category: 'custom-props',
    display: { name: 'PDBe Structure Quality Report' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private attach = StructureQualityReport.createAttachTask(
            m => `https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${m.label.toLowerCase()}`,
            this.ctx.fetch
        );

        private provider: CustomPropertyRegistry.ModelProvider = {
            option: [StructureQualityReport.Descriptor.name, 'PDBe Structure Quality Report'],
            descriptor: StructureQualityReport.Descriptor,
            defaultSelected: this.params.autoAttach,
            attachableTo: () => true,
            attach: this.attach
        }

        register(): void {
            this.ctx.customModelProperties.register(this.provider);
            this.ctx.lociLabels.addProvider(labelPDBeValidation);

            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.add('pdbe-structure-quality-report', {
                label: 'PDBe Structure Quality Report',
                factory: StructureQualityReportColorTheme,
                getParams: () => ({}),
                defaultValues: {},
                isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && !ctx.structure.isEmpty && ctx.structure.models[0].customProperties.has(StructureQualityReport.Descriptor)
            })
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
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.remove('pdbe-structure-quality-report')
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
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