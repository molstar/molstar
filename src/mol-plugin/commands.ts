/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Camera } from '../mol-canvas3d/camera';
import { PluginCommand } from './command';
import { StateTransform, State, StateAction } from '../mol-state';
import { Canvas3DProps } from '../mol-canvas3d/canvas3d';
import { PluginLayoutStateProps } from './layout';
import { StructureElement } from '../mol-model/structure';
import { PluginState } from './state';
import { PluginToast } from './util/toast';
import { Vec3 } from '../mol-math/linear-algebra';

export const PluginCommands = {
    State: {
        SetCurrentObject: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ApplyAction: PluginCommand<{ state: State, action: StateAction.Instance, ref?: StateTransform.Ref }>(),
        Update: PluginCommand<{ state: State, tree: State.Tree | State.Builder, options?: Partial<State.UpdateOptions> }>(),

        RemoveObject: PluginCommand<{ state: State, ref: StateTransform.Ref, removeParentGhosts?: boolean }>(),

        ToggleExpanded: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),
        ToggleVisibility: PluginCommand<{ state: State, ref: StateTransform.Ref }>(),

        Snapshots: {
            Add: PluginCommand<{ name?: string, description?: string, params?: PluginState.SnapshotParams }>(),
            Replace: PluginCommand<{ id: string, params?: PluginState.SnapshotParams }>(),
            Move: PluginCommand<{ id: string, dir: -1 | 1 }>(),
            Remove: PluginCommand<{ id: string }>(),
            Apply: PluginCommand<{ id: string }>(),
            Clear: PluginCommand<{}>(),

            Upload: PluginCommand<{ name?: string, description?: string, playOnLoad?: boolean, serverUrl: string, params?: PluginState.SnapshotParams }>(),
            Fetch: PluginCommand<{ url: string }>(),

            DownloadToFile: PluginCommand<{ name?: string, type: PluginState.SnapshotType, params?: PluginState.SnapshotParams }>(),
            OpenFile: PluginCommand<{ file: File }>(),
            OpenUrl: PluginCommand<{ url: string, type: PluginState.SnapshotType }>(),
        }
    },
    Interactivity: {
        Object: {
            Highlight: PluginCommand<{ state: State, ref: StateTransform.Ref | StateTransform.Ref[] }>(),
        },
        Structure: {
            Highlight: PluginCommand<{ loci: StructureElement.Loci, isOff?: boolean }>(),
            Select: PluginCommand<{ loci: StructureElement.Loci, isOff?: boolean }>()
        },
        ClearHighlights: PluginCommand<{}>(),
    },
    Layout: {
        Update: PluginCommand<{ state: Partial<PluginLayoutStateProps> }>()
    },
    Toast: {
        Show: PluginCommand<PluginToast>(),
        Hide: PluginCommand<{ key: string }>()
    },
    Camera: {
        Reset: PluginCommand<{ durationMs?: number, snapshot?: Partial<Camera.Snapshot> }>(),
        SetSnapshot: PluginCommand<{ snapshot: Partial<Camera.Snapshot>, durationMs?: number }>(),
        Focus: PluginCommand<{ center: Vec3, radius: number, durationMs?: number }>()
    },
    Canvas3D: {
        SetSettings: PluginCommand<{ settings: Partial<Canvas3DProps> | ((old: Canvas3DProps) => Partial<Canvas3DProps> | void) }>(),
        ResetSettings: PluginCommand<{ }>()
    }
};