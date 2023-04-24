/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginContext } from '../../mol-plugin/context';


export const StructConnExtensionFunctions = {
    foo(plugin: PluginContext) {
        console.log('foo:', plugin);
    }
};
