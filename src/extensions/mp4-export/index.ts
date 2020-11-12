/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { Mp4EncoderUI } from './ui';

export const Mp4Export = PluginBehavior.create<{ }>({
    name: 'extension-mp4-export',
    category: 'misc',
    display: {
        name: 'MP4 Animation Export'
    },
    ctor: class extends PluginBehavior.Handler<{ }> {
        register(): void {
            this.ctx.customStructureControls.set('mp4-export', Mp4EncoderUI as any);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.customStructureControls.delete('mp4-export');
        }
    },
    params: () => ({ })
});