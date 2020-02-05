/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { AccessibleSurfaceAreaProvider } from '../../../../../mol-model-props/computed/accessible-surface-area';
import { Loci } from '../../../../../mol-model/loci';
import { AccessibleSurfaceAreaColorThemeProvider } from '../../../../../mol-model-props/computed/themes/accessible-surface-area';
import { OrderedSet } from '../../../../../mol-data/int';
import { arraySum } from '../../../../../mol-util/array';

export const AccessibleSurfaceArea = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'computed-accessible-surface-area-prop',
    category: 'custom-props',
    display: { name: 'Accessible Surface Area' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        private provider = AccessibleSurfaceAreaProvider

        private label = (loci: Loci): string | undefined => {
            if (!this.params.showTooltip) return

            const { granularity } = this.ctx.interactivity.props
            if (granularity === 'element' || granularity === 'elementInstances') return

            if(loci.kind === 'element-loci') {
                if (loci.elements.length === 0) return;

                const accessibleSurfaceArea = AccessibleSurfaceAreaProvider.get(loci.structure).value
                if (!accessibleSurfaceArea) return;

                const { getSerialIndex } = loci.structure.root.serialMapping
                const { area, serialResidueIndex } = accessibleSurfaceArea
                const seen = new Set<number>()
                let cummulativeArea = 0

                for (const { indices, unit } of loci.elements) {
                    OrderedSet.forEach(indices, idx => {
                        const rSI = serialResidueIndex[getSerialIndex(unit, unit.elements[idx])]
                        if (rSI !== -1 && !seen.has(rSI)) {
                            cummulativeArea += area[rSI]
                            seen.add(rSI)
                        }
                    })
                }
                if (seen.size === 0) return

                return `Accessible Surface Area: ${cummulativeArea.toFixed(2)} \u212B<sup>3</sup>`;

            } else if(loci.kind === 'structure-loci') {
                const accessibleSurfaceArea = AccessibleSurfaceAreaProvider.get(loci.structure).value
                if (!accessibleSurfaceArea) return;

                return `Accessible Surface Area: ${arraySum(accessibleSurfaceArea.area).toFixed(2)} \u212B<sup>3</sup>`;
            }
        }

        update(p: { autoAttach: boolean, showTooltip: boolean }) {
            let updated = (
                this.params.autoAttach !== p.autoAttach ||
                this.params.showTooltip !== p.showTooltip
            )
            this.params.autoAttach = p.autoAttach;
            this.params.showTooltip = p.showTooltip;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        register(): void {
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.add('accessible-surface-area', AccessibleSurfaceAreaColorThemeProvider)
            this.ctx.lociLabels.addProvider(this.label);
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.structureRepresentation.themeCtx.colorThemeRegistry.remove('accessible-surface-area')
            this.ctx.lociLabels.removeProvider(this.label);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        showTooltip: PD.Boolean(true)
    })
});