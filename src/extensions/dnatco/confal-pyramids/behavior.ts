/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsProvider } from './property';
import { Loci } from '../../../mol-model/loci';
import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const DnatcoConfalPyramids = PluginBehavior.create<{ autoAttach: boolean, showToolTip: boolean }>({
    name: 'dnatco-confal-pyramids-prop',
    category: 'custom-props',
    display: {
        name: 'Confal Pyramids',
        description: 'Schematic depiction of conformer class and confal value.',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showToolTip: boolean }> {

        private provider = ConfalPyramidsProvider;

        private labelConfalPyramids = {
            label: (loci: Loci): string | undefined => {
                if (!this.params.showToolTip) return void 0;

                /* TODO: Implement this */
                return void 0;
            }
        }

        register(): void {
            this.ctx.customModelProperties.register(this.provider, this.params.autoAttach);
            this.ctx.managers.lociLabels.addProvider(this.labelConfalPyramids);

            /* TODO: Add color and visual providers */
        }

        update(p: { autoAttach: boolean, showToolTip: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.params.showToolTip = p.showToolTip;
            this.ctx.customModelProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            this.ctx.customModelProperties.unregister(ConfalPyramidsProvider.descriptor.name);
            this.ctx.managers.lociLabels.removeProvider(this.labelConfalPyramids);

            /* TODO: Unregister color and visual providers */
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
        showToolTip: PD.Boolean(true)
    })
});
