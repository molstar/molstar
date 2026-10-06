/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import type { Camera } from '@molstar/graphics/canvas3d/camera';
import { PluginCommand } from '@molstar/plugin/command';
import type { StateTransform, State, StateAction } from '@molstar/core/state';
import type { Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginLayoutStateProps } from '@molstar/plugin/layout';
import type { Structure, StructureElement } from '@molstar/model/model/structure';
import type { PluginState } from '@molstar/plugin/state';
import type { PluginToast } from '@molstar/plugin/util/toast';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import type { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import type { TransitionTrajectory } from '@molstar/graphics/canvas3d/camera/transition-functions';
import type { EasingFunction } from '@molstar/core/math/easing';

export const PluginCommands = {
  State: {
    SetCurrentObject: PluginCommand<{ state: State; ref: StateTransform.Ref }>(),
    ApplyAction: PluginCommand<{ state: State; action: StateAction.Instance; ref?: StateTransform.Ref }>(),
    Update: PluginCommand<{ state: State; tree: State.Tree | State.Builder; options?: Partial<State.UpdateOptions> }>(),

    RemoveObject: PluginCommand<{ state: State; ref: StateTransform.Ref; removeParentGhosts?: boolean }>(),

    ToggleExpanded: PluginCommand<{ state: State; ref: StateTransform.Ref }>(),
    ToggleVisibility: PluginCommand<{ state: State; ref: StateTransform.Ref }>(),

    Snapshots: {
      Add: PluginCommand<{
        key?: string;
        name?: string;
        description?: string;
        descriptionFormat?: PluginStateSnapshotManager.DescriptionFormat;
        params?: PluginState.SnapshotParams;
      }>(),
      Replace: PluginCommand<{ id: string; params?: PluginState.SnapshotParams }>(),
      Move: PluginCommand<{ id: string; dir: -1 | 1 }>(),
      Remove: PluginCommand<{ id: string }>(),
      Apply: PluginCommand<{ id: string }>(),
      Clear: PluginCommand<{}>(),

      Upload: PluginCommand<{
        name?: string;
        description?: string;
        playOnLoad?: boolean;
        serverUrl: string;
        params?: PluginState.SnapshotParams;
      }>(),
      Fetch: PluginCommand<{ url: string }>(),

      DownloadToFile: PluginCommand<{
        name?: string;
        type: PluginState.SnapshotType;
        params?: PluginState.SnapshotParams;
      }>(),
      OpenFile: PluginCommand<{ file: File }>(),
      OpenUrl: PluginCommand<{ url: string; type: PluginState.SnapshotType }>(),
    },
  },
  Interactivity: {
    Object: {
      Highlight: PluginCommand<{ state: State; ref: StateTransform.Ref | StateTransform.Ref[] }>(),
    },
    Structure: {
      Highlight: PluginCommand<{ loci: StructureElement.Loci; isOff?: boolean }>(),
      Select: PluginCommand<{ loci: StructureElement.Loci; isOff?: boolean }>(),
    },
    ClearHighlights: PluginCommand<{}>(),
  },
  Layout: {
    Update: PluginCommand<{ state: Partial<PluginLayoutStateProps> }>(),
  },
  Toast: {
    Show: PluginCommand<PluginToast>(),
    Hide: PluginCommand<{ key: string }>(),
  },
  Camera: {
    Reset: PluginCommand<{
      durationMs?: number;
      easing?: EasingFunction;
      trajectory?: TransitionTrajectory;
      snapshot?: Partial<Camera.Snapshot>;
    }>(),
    SetSnapshot: PluginCommand<{
      snapshot: Partial<Camera.Snapshot>;
      durationMs?: number;
      easing?: EasingFunction;
      trajectory?: TransitionTrajectory;
    }>(),
    Focus: PluginCommand<{
      center: Vec3;
      radius: number;
      durationMs?: number;
      easing?: EasingFunction;
      trajectory?: TransitionTrajectory;
    }>(),
    FocusObject: PluginCommand<
      PluginState.SnapshotFocusInfo & {
        durationMs?: number;
        easing?: EasingFunction;
        trajectory?: TransitionTrajectory;
      }
    >(),
    OrientAxes: PluginCommand<{
      structures?: Structure[];
      durationMs?: number;
      easing?: EasingFunction;
      trajectory?: TransitionTrajectory;
    }>(),
    ResetAxes: PluginCommand<{ durationMs?: number; easing?: EasingFunction; trajectory?: TransitionTrajectory }>(),
  },
  Canvas3D: {
    SetSettings: PluginCommand<{
      settings: Partial<Canvas3DProps> | ((old: Canvas3DProps) => Partial<Canvas3DProps> | void);
    }>(),
    ResetSettings: PluginCommand<{}>(),
  },
};
