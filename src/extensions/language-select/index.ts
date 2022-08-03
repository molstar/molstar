/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { ScriptImportUI } from './ui';

export const ScriptSetting = PluginBehavior.create<{ }>({
    name: 'extension script language',
    category: 'misc',
    display: {
        name: 'Script Language'
    },
    ctor: class extends PluginBehavior.Handler<{ }> {
        register(): void {
            this.ctx.customImportControls.set('script-language', ScriptImportUI as any);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.customImportControls.delete('script-language');
        }
    },
    params: () => ({ })
});
