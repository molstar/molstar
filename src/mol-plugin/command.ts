/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as State from './command/state';
import * as Camera from './command/camera';

export * from './command/command';

export const PluginCommands = {
    State,
    Camera
}