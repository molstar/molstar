/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommand } from './command';
import { Camera } from 'mol-canvas3d/camera';

export const Reset = PluginCommand<{}>({ isImmediate: true });
export const SetSnapshot = PluginCommand<{ snapshot: Camera.Snapshot }>({ isImmediate: true });

export const Snapshots = {
    Add: PluginCommand<{ name?: string, description?: string }>({ isImmediate: true }),
    Remove: PluginCommand<{ id: string }>({ isImmediate: true }),
    Apply: PluginCommand<{ id: string }>({ isImmediate: true }),
    Clear: PluginCommand<{ }>({ isImmediate: true }),
}