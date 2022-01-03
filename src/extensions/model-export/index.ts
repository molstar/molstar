/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { ModelExportUI } from './ui';

export const ModelExport = PluginBehavior.create<{}>({
    name: 'extension-model-export',
    category: 'misc',
    display: {
        name: 'Model Export'
    },
    ctor: class extends PluginBehavior.Handler<{}> {
        register(): void {
            this.ctx.customStructureControls.set('model-export', ModelExportUI as any);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.customStructureControls.delete('model-export');
        }
    },
    params: () => ({})
});