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
import { StructureElement } from 'mol-model/structure';
import { PluginState } from './state';

export * from './command/base';

export const PluginCommands = {
    State: {
        SetCurrentObject: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ApplyAction: PluginCommand<{ state: State, action: StateAction.Instance, ref?: StateTransform.Ref }>(),
        Update: PluginCommand<{ state: State, tree: State.Tree | State.Builder, options?: Partial<State.UpdateOptions> }>(),

        RemoveObject: PluginCommand<{ state: State, ref: StateTransform.Ref, removeParentGhosts?: boolean }>(),

        ToggleExpanded: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ToggleVisibility: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        Highlight: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ClearHighlight: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),

        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string, params?: PluginState.GetSnapshotParams }>(),
            Replace: PluginCommand<{ id: string, params?: PluginState.GetSnapshotParams }>(),
            Move: PluginCommand<{ id: string, dir: -1 | 1 }>(),
            Remove: PluginCommand<{ id: string }>(),
            Apply: PluginCommand<{ id: string }>(),
            Clear: PluginCommand<{}>(),

            Upload: PluginCommand<{ name?: string, description?: string, serverUrl: string }>(),
            Fetch: PluginCommand<{ url: string }>(),

            DownloadToFile: PluginCommand<{ name?: string }>(),
            OpenFile: PluginCommand<{ file: File }>(),
        }
    },
    Interactivity: {
        Structure: {
            Highlight: PluginCommand<{ loci: StructureElement.Loci, isOff?: boolean }>(),
            Select: PluginCommand<{ loci: StructureElement.Loci, isOff?: boolean }>()
        }
    },
    Layout: {
        Update: PluginCommand<{ state: Partial<PluginLayoutStateProps> }>()
    },
    Camera: {
        Reset: PluginCommand<{}>(),
        SetSnapshot: PluginCommand<{ snapshot: Camera.Snapshot, durationMs?: number }>(),
        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string }>(),
            Remove: PluginCommand<{ id: string }>(),
            Apply: PluginCommand<{ id: string }>(),
            Clear: PluginCommand<{}>(),
        }
    },
    Canvas3D: {
        SetSettings: PluginCommand<{ settings: Partial<Canvas3DProps> }>()
    }
}