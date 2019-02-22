/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from 'mol-canvas3d/camera';
import { PluginCommand } from './command/base';
import { StateTransform, State, StateAction } from 'mol-state';
import { Canvas3DProps } from 'mol-canvas3d/canvas3d';
import { PluginLayoutStateProps } from './layout';

export * from './command/base';

export const PluginCommands = {
    State: {
        SetCurrentObject: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ApplyAction: PluginCommand<{ state: State, action: StateAction.Instance, ref?: StateTransform.Ref }>(),
        Update: PluginCommand<{ state: State, tree: State.Tree | State.Builder, doNotLogTiming?: boolean }>(),

        RemoveObject: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),

        ToggleExpanded: PluginCommand<{ state: State, ref: StateTransform.Ref }>({ isImmediate: true }),
        ToggleVisibility: PluginCommand<{ state: State, ref: StateTransform.Ref }>({ isImmediate: true }),
        Highlight: PluginCommand<{ state: State, ref: StateTransform.Ref }>({ isImmediate: true }),
        ClearHighlight: PluginCommand<{ state: State, ref: StateTransform.Ref }>({ isImmediate: true }),

        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string }>({ isImmediate: true }),
            Remove: PluginCommand<{ id: string }>({ isImmediate: true }),
            Apply: PluginCommand<{ id: string }>({ isImmediate: true }),
            Clear: PluginCommand<{}>({ isImmediate: true }),

            Upload: PluginCommand<{ name?: string, description?: string, serverUrl: string }>({ isImmediate: true }),
            Fetch: PluginCommand<{ url: string }>(),

            DownloadToFile: PluginCommand<{ name?: string }>({ isImmediate: true }),
            OpenFile: PluginCommand<{ file: File }>({ isImmediate: true }),
        }
    },
    Layout: {
        Update: PluginCommand<{ state: Partial<PluginLayoutStateProps> }>({ isImmediate: true })
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
    },
    Canvas3D: {
        SetSettings: PluginCommand<{ settings: Partial<Canvas3DProps> }>({ isImmediate: true })
    }
}