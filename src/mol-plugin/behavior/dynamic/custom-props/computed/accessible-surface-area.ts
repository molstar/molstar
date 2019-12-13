/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { AccessibleSurfaceAreaProvider } from '../../../../../mol-model-props/computed/accessible-surface-area';

export const AccessibleSurfaceArea = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'computed-accessible-surface-area-prop',
    category: 'custom-props',
    display: { name: 'Accessible Surface Area' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.customStructureProperties.register(AccessibleSurfaceAreaProvider);
        }
        unregister() {
            this.ctx.customStructureProperties.unregister(AccessibleSurfaceAreaProvider.descriptor.name);
        }
    }
});