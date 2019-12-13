/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { AccessibleSurfaceAreaProvider } from '../../../../../mol-model-props/computed/accessible-surface-area';
import { Loci } from '../../../../../mol-model/loci';

export const AccessibleSurfaceArea = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'computed-accessible-surface-area-prop',
    category: 'custom-props',
    display: { name: 'Accessible Surface Area' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        private provider = AccessibleSurfaceAreaProvider

        private label = (loci: Loci): string | undefined => {
            if (!this.params.showTooltip) return void 0;

            switch (loci.kind) {
                case 'element-loci':
                    if (loci.elements.length === 0) return void 0;
                    // const e = loci.elements[0];
                    // const u = e.unit;
                    if (!this.provider.getValue(loci.structure).value) return;

                    return `Accessible Surface Area: ${'TODO'} \u212B<sup>3</sup>`;

                default: return void 0;
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
            this.ctx.lociLabels.addProvider(this.label);
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.lociLabels.removeProvider(this.label);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        showTooltip: PD.Boolean(true)
    })
});