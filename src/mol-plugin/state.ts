/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { produce } from '../mol-util/produce';
import { merge } from 'rxjs';
import { Camera } from '../mol-canvas3d/camera';
import { Canvas3DContext, Canvas3DParams, Canvas3DProps } from '../mol-canvas3d/canvas3d';
import { Vec3 } from '../mol-math/linear-algebra';
import { PluginComponent } from '../mol-plugin-state/component';
import { PluginAnimationManager } from '../mol-plugin-state/manager/animation';
import { InteractivityManager } from '../mol-plugin-state/manager/interactivity';
import { StructureComponentManager } from '../mol-plugin-state/manager/structure/component';
import { StructureFocusSnapshot } from '../mol-plugin-state/manager/structure/focus';
import { StructureSelectionSnapshot } from '../mol-plugin-state/manager/structure/selection';
import { PluginStateObject as SO } from '../mol-plugin-state/objects';
import { State, StateTransform, StateTransformer } from '../mol-state';
import { UUID } from '../mol-util';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { PluginBehavior } from './behavior';
import { PluginCommands } from './commands';
import { PluginConfig } from './config';
import { PluginContext } from './context';
import { AnimateStateSnapshotTransition } from '../mol-plugin-state/animation/built-in/state-snapshots';
import { Scheduler } from '../mol-task';

export { PluginState };

class PluginState extends PluginComponent {
    private get animation() { return this.plugin.managers.animation; }

    readonly data = State.create(new SO.Root({}), { runTask: this.plugin.runTask, globalContext: this.plugin, historyCapacity: this.plugin.config.get(PluginConfig.State.HistoryCapacity) });
    readonly behaviors = State.create(new PluginBehavior.Root({}), { runTask: this.plugin.runTask, globalContext: this.plugin, rootState: { isLocked: true } });

    readonly events = {
        cell: {
            stateUpdated: merge(this.data.events.cell.stateUpdated, this.behaviors.events.cell.stateUpdated),
            created: merge(this.data.events.cell.created, this.behaviors.events.cell.created),
            removed: merge(this.data.events.cell.removed, this.behaviors.events.cell.removed),
        },
        object: {
            created: merge(this.data.events.object.created, this.behaviors.events.object.created),
            removed: merge(this.data.events.object.removed, this.behaviors.events.object.removed),
            updated: merge(this.data.events.object.updated, this.behaviors.events.object.updated)
        }
    } as const;

    readonly snapshotParams = this.ev.behavior<PluginState.SnapshotParams>(PluginState.DefaultSnapshotParams);

    setSnapshotParams = (params?: PluginState.SnapshotParams) => {
        this.snapshotParams.next({ ...PluginState.DefaultSnapshotParams, ...params });
    };

    getSnapshot(params?: PluginState.SnapshotParams): PluginState.Snapshot {
        const p = { ...this.snapshotParams.value, ...params };
        return {
            id: UUID.create22(),
            data: p.data ? this.data.getSnapshot() : void 0,
            behaviour: p.behavior ? this.behaviors.getSnapshot() : void 0,
            animation: p.animation ? this.animation.getSnapshot() : void 0,
            startAnimation: p.startAnimation ? !!p.startAnimation : void 0,
            camera: p.camera ? {
                current: this.plugin.canvas3d!.camera.getSnapshot(),
                transitionStyle: p.cameraTransition!.name,
                transitionDurationInMs: p?.cameraTransition?.name === 'animate' ? p.cameraTransition.params.durationInMs : void 0
            } : void 0,
            canvas3dContext: p.canvas3dContext ? { props: this.plugin.canvas3dContext?.props } : void 0,
            canvas3d: p.canvas3d ? { props: this.plugin.canvas3d?.props } : void 0,
            interactivity: p.interactivity ? { props: this.plugin.managers.interactivity.props } : void 0,
            structureFocus: this.plugin.managers.structure.focus.getSnapshot(),
            structureSelection: p.structureSelection ? this.plugin.managers.structure.selection.getSnapshot() : void 0,
            structureComponentManager: p.componentManager ? {
                options: this.plugin.managers.structure.component.state.options
            } : void 0,
            durationInMs: p?.durationInMs
        };
    }

    async setSnapshot(snapshot: PluginState.Snapshot) {
        await this.animation.stop();

        // this needs to go 1st since these changes are already baked into the behavior and data state
        if (snapshot.structureComponentManager?.options) this.plugin.managers.structure.component._setSnapshotState(snapshot.structureComponentManager?.options);
        if (snapshot.behaviour) await this.plugin.runTask(this.behaviors.setSnapshot(snapshot.behaviour));
        if (snapshot.data) await this.plugin.runTask(this.data.setSnapshot(snapshot.data));
        if (snapshot.canvas3d?.props) {
            const settings = PD.normalizeParams(Canvas3DParams, snapshot.canvas3d.props, 'children');
            await PluginCommands.Canvas3D.SetSettings(this.plugin, { settings });
        }
        if (snapshot.canvas3dContext?.props) {
            const props = PD.normalizeParams(Canvas3DContext.Params, snapshot.canvas3dContext.props, 'children');
            this.plugin.canvas3dContext?.setProps(props);
        }
        if (snapshot.interactivity) {
            if (snapshot.interactivity.props) this.plugin.managers.interactivity.setProps(snapshot.interactivity.props);
        }
        if (snapshot.structureFocus) {
            this.plugin.managers.structure.focus.setSnapshot(snapshot.structureFocus);
        }
        if (snapshot.structureSelection) {
            this.plugin.managers.structure.selection.setSnapshot(snapshot.structureSelection);
        }
        if (snapshot.animation) {
            this.animation.setSnapshot(snapshot.animation);
        }
        if (snapshot.camera?.current) {
            PluginCommands.Camera.Reset(this.plugin, {
                snapshot: snapshot.camera.current,
                durationMs: snapshot.camera.transitionStyle === 'animate' ? snapshot.camera.transitionDurationInMs : undefined,
            });
        } else if (snapshot.camera?.focus) {
            PluginCommands.Camera.FocusObject(this.plugin, {
                ...snapshot.camera.focus,
                durationMs: snapshot.camera.transitionStyle === 'animate' ? snapshot.camera.transitionDurationInMs : undefined,
            });
        }

        if (typeof snapshot?.onLoadMarkdownCommands === 'object' && Object.keys(snapshot.onLoadMarkdownCommands).length > 0) {
            this.plugin.managers.markdownExtensions.tryExecute('click', snapshot.onLoadMarkdownCommands);
        }

        if (snapshot.startAnimation) {
            this.animation.start();
            return;
        }

        if (snapshot.transition?.autoplay) {
            await Scheduler.immediatePromise();
            this.animation.play(AnimateStateSnapshotTransition, {});
        }
    }

    async setAnimationSnapshot(snapshot: PluginState.Snapshot, frameIndex: number) {
        await this.animation.stopStateTransitionAnimation();

        const { transition } = snapshot;
        if (!transition) return;
        const finalIndex = Math.min(frameIndex, transition.frames.length - 1);
        const frame = transition.frames[finalIndex] ?? snapshot.data;
        if (frame.data) await this.plugin.runTask(this.data.setSnapshot(frame.data));
        if (frame.canvas3d?.props) {
            const settings = PD.normalizeParams(Canvas3DParams, frame.canvas3d.props, 'children');
            this.plugin.canvas3d?.setProps(settings);
        }
        if (frame.camera?.current) {
            PluginCommands.Camera.Reset(this.plugin, {
                snapshot: frame.camera.current,
                durationMs: frame.camera.transitionStyle === 'animate' ? frame.camera.transitionDurationInMs : undefined,
            });
        }

        if (typeof snapshot?.onLoadMarkdownCommands === 'object' && Object.keys(snapshot.onLoadMarkdownCommands).length > 0) {
            this.plugin.managers.markdownExtensions.tryExecute('click', snapshot.onLoadMarkdownCommands);
        }
    }

    updateTransform(state: State, a: StateTransform.Ref, params: any, canUndo?: string | boolean) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update(this.plugin, { state, tree, options: { canUndo } });
    }

    hasBehavior(behavior: StateTransformer) {
        return this.behaviors.tree.transforms.has(behavior.id);
    }

    updateBehavior<T extends StateTransformer>(behavior: T, params: (old: StateTransformer.Params<T>) => (void | StateTransformer.Params<T>)) {
        const tree = this.behaviors.build();
        if (!this.behaviors.tree.transforms.has(behavior.id)) {
            const defaultParams = behavior.createDefaultParams(void 0 as any, this.plugin);
            tree.to(PluginBehavior.getCategoryId(behavior)).apply(behavior, produce(defaultParams, params) as any, { ref: behavior.id });
        } else {
            tree.to(behavior.id).update(params);
        }
        return this.plugin.runTask(this.behaviors.updateTree(tree));
    }

    dispose() {
        this.behaviors.cells.forEach(cell => {
            if (PluginBehavior.Behavior.is(cell.obj)) {
                cell.obj.data.unregister?.();
                cell.obj.data.dispose?.();
            }
        });

        super.dispose();
        this.data.dispose();
        this.behaviors.dispose();
        this.animation.dispose();
    }

    constructor(private plugin: PluginContext) {
        super();
    }
}

namespace PluginState {
    export type CameraTransitionStyle = 'instant' | 'animate'
    export const SnapshotParams = {
        durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Duration in ms' }),
        data: PD.Boolean(true),
        behavior: PD.Boolean(false),
        structureSelection: PD.Boolean(false),
        componentManager: PD.Boolean(true),
        animation: PD.Boolean(true),
        startAnimation: PD.Boolean(false),
        canvas3d: PD.Boolean(true),
        canvas3dContext: PD.Boolean(true),
        interactivity: PD.Boolean(true),
        camera: PD.Boolean(true),
        cameraTransition: PD.MappedStatic('animate', {
            animate: PD.Group({
                durationInMs: PD.Numeric(250, { min: 100, max: 5000, step: 500 }, { label: 'Duration in ms' }),
            }),
            instant: PD.Group({})
        }, { options: [['animate', 'Animate'], ['instant', 'Instant']] }),
        image: PD.Boolean(false),
    };
    export type SnapshotParams = Partial<PD.Values<typeof SnapshotParams>>
    export const DefaultSnapshotParams = PD.getDefaultValues(SnapshotParams);

    export interface Snapshot {
        id: UUID,
        data?: State.Snapshot,
        behaviour?: State.Snapshot,
        animation?: PluginAnimationManager.Snapshot,
        startAnimation?: boolean,
        camera?: {
            current?: Camera.Snapshot,
            focus?: SnapshotFocusInfo,
            transitionStyle: CameraTransitionStyle,
            transitionDurationInMs?: number
        },
        canvas3d?: {
            props?: Canvas3DProps
        },
        canvas3dContext?: {
            props?: Canvas3DContext.Props
        },
        interactivity?: {
            props?: InteractivityManager.Props
        },
        structureFocus?: StructureFocusSnapshot,
        structureSelection?: StructureSelectionSnapshot,
        structureComponentManager?: {
            options?: StructureComponentManager.Options
        },
        durationInMs?: number,
        transition?: StateTransition,
        onLoadMarkdownCommands?: Record<string, any>
    }

    export interface StateTransition {
        autoplay?: boolean,
        loop?: boolean,
        frames: {
            durationInMs: number,
            data: State.Snapshot,
            camera?: Snapshot['camera'],
            canvas3d?: { props?: Canvas3DProps },
        }[],
    }

    export function getStateTransitionDuration(snapshot: Snapshot): number | undefined {
        const { transition } = snapshot;
        if (!transition) return undefined;
        let totalDuration = 0;
        for (let i = 0; i < transition.frames.length; i++) {
            const frame = transition.frames[i];
            totalDuration += frame.durationInMs;
        }
        return totalDuration;
    }

    export function getStateTransitionFrameIndex(snapshot: Snapshot, timestamp: number): number | undefined {
        const { transition } = snapshot;
        if (!transition) return undefined;

        let t = timestamp;
        if (transition.loop) {
            t %= getStateTransitionDuration(snapshot) ?? 1;
        }

        let currentDuration = 0;
        for (let i = 0; i < transition.frames.length; i++) {
            if (currentDuration >= t) return i;
            const frame = transition.frames[i];
            currentDuration += frame.durationInMs;
        }
        return transition.frames.length - 1;
    }

    export type SnapshotType = 'json' | 'molj' | 'zip' | 'molx'

    export interface SnapshotFocusInfo {
        targets?: SnapshotFocusTargetInfo[],
        direction?: Vec3,
        up?: Vec3,
    }
    /** Final radius to be computed as `radius ?? targetBoundingRadius * radiusFactor + extraRadius` */
    export interface SnapshotFocusTargetInfo {
        /** Reference to plugin state node to focus (undefined means focus whole scene) */
        targetRef?: StateTransform.Ref,
        radius?: number,
        radiusFactor?: number,
        extraRadius?: number,
    }
}
