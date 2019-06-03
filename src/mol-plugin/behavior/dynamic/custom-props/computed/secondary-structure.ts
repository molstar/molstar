/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { PluginBehavior } from '../../../behavior';
import { CustomPropertyRegistry } from '../../../../../mol-model-props/common/custom-property-registry';
import { ComputedSecondaryStructure } from '../../../../../mol-model-props/computed/secondary-structure';

export const MolstarSecondaryStructure = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'molstar-computed-secondary-structure-prop',
    category: 'custom-props',
    display: { name: 'Computed Secondary Structure' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private attach = ComputedSecondaryStructure.createAttachTask();

        private provider: CustomPropertyRegistry.StructureProvider = {
            option: [ComputedSecondaryStructure.Descriptor.name, 'Computed Secondary Structure'],
            descriptor: ComputedSecondaryStructure.Descriptor,
            defaultSelected: this.params.autoAttach,
            attachableTo: () => true,
            attach: this.attach
        }

        register(): void {
            this.ctx.customStructureProperties.register(this.provider);
        }

        update(p: { autoAttach: boolean }) {
            let updated = this.params.autoAttach !== p.autoAttach
            this.params.autoAttach = p.autoAttach;
            this.provider.defaultSelected = p.autoAttach;
            return updated;
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(ComputedSecondaryStructure.Descriptor.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});