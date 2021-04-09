/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { GeometryExporterUI } from './ui';

export const GeometryExport = PluginBehavior.create<{ }>({
    name: 'extension-geo-export',
    category: 'misc',
    display: {
        name: 'Geometry Export'
    },
    ctor: class extends PluginBehavior.Handler<{ }> {
        register(): void {
            this.ctx.customStructureControls.set('geo-export', GeometryExporterUI as any);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.customStructureControls.delete('geo-export');
        }
    },
    params: () => ({ })
});