/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { SecondaryStructureProvider } from '../../../../../mol-model-props/computed/secondary-structure';

export const SecondaryStructure = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'computed-secondary-structure-prop',
    category: 'custom-props',
    display: { name: 'Secondary Structure' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.customStructureProperties.register(SecondaryStructureProvider);
        }
        unregister() {
            this.ctx.customStructureProperties.unregister(SecondaryStructureProvider.descriptor.name);
        }
    }
});