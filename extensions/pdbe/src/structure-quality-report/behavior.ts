/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { StructureQualityReport, StructureQualityReportProvider } from './prop.js';
import { StructureQualityReportColorThemeProvider } from './color.js';
import { Loci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export const PDBeStructureQualityReport = PluginBehavior.create<{ autoAttach: boolean; showTooltip: boolean }>({
  name: 'pdbe-structure-quality-report-prop',
  category: 'custom-props',
  display: {
    name: 'Structure Quality Report',
    description: 'Data from wwPDB Validation Report, obtained via PDBe.',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean; showTooltip: boolean }> {
    private provider = StructureQualityReportProvider;

    private labelPDBeValidation = {
      label: (loci: Loci): string | undefined => {
        if (!this.params.showTooltip) return void 0;

        switch (loci.kind) {
          case 'element-loci':
            if (loci.elements.length === 0) return void 0;
            const e = loci.elements[0];
            const u = e.unit;
            if (!u.model.customProperties.hasReference(StructureQualityReportProvider.descriptor)) return void 0;

            const se = StructureElement.Location.create(loci.structure, u, u.elements[OrderedSet.getAt(e.indices, 0)]);
            const issues = StructureQualityReport.getIssues(se);
            if (issues.length === 0) return 'Validation: No Issues';
            return `Validation: ${issues.join(', ')}`;

          default:
            return void 0;
        }
      },
    };

    private unregisterEntry: (() => void) | undefined;

    register(): void {
      this.ctx.customModelProperties.register(this.provider, this.params.autoAttach);

      const entry: PluginRegistryEntry = {
        structure: { themes: { color: [StructureQualityReportColorThemeProvider] } },
        lociLabels: [this.labelPDBeValidation],
      };
      this.unregisterEntry = this.ctx.register(entry);
    }

    update(p: { autoAttach: boolean; showTooltip: boolean }) {
      const updated = this.params.autoAttach !== p.autoAttach;
      this.params.autoAttach = p.autoAttach;
      this.params.showTooltip = p.showTooltip;
      this.ctx.customModelProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
      return updated;
    }

    unregister() {
      this.ctx.customModelProperties.unregister(StructureQualityReportProvider.descriptor.name);
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
    }
  },
  params: () => ({
    autoAttach: PD.Boolean(false),
    showTooltip: PD.Boolean(true),
  }),
});
