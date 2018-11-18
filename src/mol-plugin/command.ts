/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from 'mol-canvas3d/camera';
import { PluginCommand } from './command/base';
import { Transform, State } from 'mol-state';
import { StateAction } from 'mol-state/action';

export * from './command/base';

export const PluginCommands = {
    State: {
        SetCurrentObject: PluginCommand<{ state: State, ref: Transform.Ref }>(),
        ApplyAction: PluginCommand<{ state: State, action: StateAction.Instance, ref?: Transform.Ref }>(),
        Update: PluginCommand<{ state: State, tree: State.Tree | State.Builder }>(),

        RemoveObject: PluginCommand<{ state: State, ref: Transform.Ref }>(),

        ToggleExpanded: PluginCommand<{ state: State, ref: Transform.Ref }>({ isImmediate: true }),
        ToggleVisibility: PluginCommand<{ state: State, ref: Transform.Ref }>({ isImmediate: true }),

        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string }>({ isImmediate: true }),
            Remove: PluginCommand<{ id: string }>({ isImmediate: true }),
            Apply: PluginCommand<{ id: string }>({ isImmediate: true }),
            Clear: PluginCommand<{}>({ isImmediate: true }),

            Upload: PluginCommand<{ name?: string, description?: string, serverUrl: string }>({ isImmediate: true }),
            Fetch: PluginCommand<{ url: string }>()
        }
    },
    Camera: {
        Reset: PluginCommand<{}>({ isImmediate: true }),
        SetSnapshot: PluginCommand<{ snapshot: Camera.Snapshot, durationMs?: number }>({ isImmediate: true }),
        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string }>({ isImmediate: true }),
            Remove: PluginCommand<{ id: string }>({ isImmediate: true }),
            Apply: PluginCommand<{ id: string }>({ isImmediate: true }),
            Clear: PluginCommand<{}>({ isImmediate: true }),
        }
    }
}