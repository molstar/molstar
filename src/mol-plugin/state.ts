/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTransform, StateTransformer } from '../mol-state';
import { PluginStateObject as SO } from '../mol-plugin-state/objects';
import { Camera } from '../mol-canvas3d/camera';
import { PluginBehavior } from './behavior';
import { CameraSnapshotManager } from '../mol-plugin-state/camera';
import { PluginStateSnapshotManager } from '../mol-plugin-state/snapshots';
import { RxEventHelper } from '../mol-util/rx-event-helper';
import { Canvas3DProps } from '../mol-canvas3d/canvas3d';
import { PluginCommands } from './commands';
import { PluginAnimationManager } from '../mol-plugin-state/animation/manager';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { UUID } from '../mol-util';
import { InteractivityManager } from '../mol-plugin-state/manager/interactivity';
import { produce } from 'immer';
export { PluginState }

class PluginState {
    private ev = RxEventHelper.create();

    readonly data: State;
    readonly behaviors: State;
    readonly animation: PluginAnimationManager;
    readonly cameraSnapshots = new CameraSnapshotManager();
    readonly snapshots: PluginStateSnapshotManager;

    getSnapshot(params?: PluginState.GetSnapshotParams): PluginState.Snapshot {
        const p = { ...PluginState.DefaultGetSnapshotParams, ...params };
        return {
            id: UUID.create22(),
            data: p.data ? this.data.getSnapshot() : void 0,
            behaviour: p.behavior ? this.behaviors.getSnapshot() : void 0,
            animation: p.animation ? this.animation.getSnapshot() : void 0,
            startAnimation: p.startAnimation ? !!p.startAnimation : void 0,
            camera: p.camera ? {
                current: this.plugin.canvas3d!.camera.getSnapshot(),
                transitionStyle: p.cameraTranstion.name,
                transitionDurationInMs: (params && params.cameraTranstion && params.cameraTranstion.name === 'animate') ? params.cameraTranstion.params.durationInMs : undefined
            } : void 0,
            cameraSnapshots: p.cameraSnapshots ? this.cameraSnapshots.getStateSnapshot() : void 0,
            canvas3d: p.canvas3d ? { props: this.plugin.canvas3d?.props } : void 0,
            interactivity: p.interactivity ? { props: this.plugin.managers.interactivity.props } : void 0,
            durationInMs: params && params.durationInMs
        };
    }

    async setSnapshot(snapshot: PluginState.Snapshot) {
        await this.animation.stop();

        if (snapshot.behaviour) await this.plugin.runTask(this.behaviors.setSnapshot(snapshot.behaviour));
        if (snapshot.data) await this.plugin.runTask(this.data.setSnapshot(snapshot.data));
        if (snapshot.canvas3d) {
            if (snapshot.canvas3d.props) await PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: snapshot.canvas3d.props });
        }
        if (snapshot.interactivity) {
            if (snapshot.interactivity.props) this.plugin.managers.interactivity.setProps(snapshot.interactivity.props);
        }
        if (snapshot.cameraSnapshots) this.cameraSnapshots.setStateSnapshot(snapshot.cameraSnapshots);
        if (snapshot.animation) {
            this.animation.setSnapshot(snapshot.animation);
        }
        if (snapshot.camera) {
            PluginCommands.Camera.Reset(this.plugin, {
                snapshot: snapshot.camera.current,
                durationMs: snapshot.camera.transitionStyle === 'animate'
                    ? snapshot.camera.transitionDurationInMs
                    : void 0
            });
        }
        if (snapshot.startAnimation) {
            this.animation.start();
        }
    }

    applyTransform(state: State, a: StateTransform.Ref, transformer: StateTransformer, params: any) {
        const tree = state.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update(this.plugin, { state, tree });
    }

    updateTransform(state: State, a: StateTransform.Ref, params: any, canUndo?: string | boolean) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update(this.plugin, { state, tree, options: { canUndo } });
    }

    updateBehavior<T extends StateTransformer>(behavior: T, params: (old: StateTransformer.Params<T>) => (void | StateTransformer.Params<T>)) {
        const tree = this.behaviors.build()
        if (!this.behaviors.tree.transforms.has(behavior.id)) {
            const defaultParams = behavior.createDefaultParams(void 0 as any, this.plugin);
            tree.to(PluginBehavior.getCategoryId(behavior)).apply(behavior, produce(defaultParams, params) as any, { ref: behavior.id });
        } else {
            tree.to(behavior.id).update(params);
        }
        return this.plugin.runTask(this.behaviors.updateTree(tree));
    }

    dispose() {
        this.ev.dispose();
        this.data.dispose();
        this.behaviors.dispose();
        this.cameraSnapshots.dispose();
        this.animation.dispose();
    }

    constructor(private plugin: import('./context').PluginContext) {
        this.snapshots = new PluginStateSnapshotManager(plugin);
        this.data = State.create(new SO.Root({ }), { runTask: plugin.runTask, globalContext: plugin });
        this.behaviors = State.create(new PluginBehavior.Root({ }), { runTask: plugin.runTask, globalContext: plugin, rootState: { isLocked: true } });

        this.animation = new PluginAnimationManager(plugin);
    }
}

namespace PluginState {
    export type Kind = 'data' | 'behavior'

    export type CameraTransitionStyle = 'instant' | 'animate'
    export const GetSnapshotParams = {
        durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Duration in ms' }),
        data: PD.Boolean(true),
        behavior: PD.Boolean(false),
        animation: PD.Boolean(true),
        startAnimation: PD.Boolean(false),
        canvas3d: PD.Boolean(true),
        interactivity: PD.Boolean(true),
        camera: PD.Boolean(true),
        // TODO: make camera snapshots same as the StateSnapshots with "child states?"
        cameraSnapshots: PD.Boolean(false),
        cameraTranstion: PD.MappedStatic('animate', {
            animate: PD.Group({
                durationInMs: PD.Numeric(250, { min: 100, max: 5000, step: 500 }, { label: 'Duration in ms' }),
            }),
            instant: PD.Group({ })
        }, { options: [['animate', 'Animate'], ['instant', 'Instant']] })
    };
    export type GetSnapshotParams = Partial<PD.Values<typeof GetSnapshotParams>>
    export const DefaultGetSnapshotParams = PD.getDefaultValues(GetSnapshotParams);

    export interface Snapshot {
        id: UUID,
        data?: State.Snapshot,
        behaviour?: State.Snapshot,
        animation?: PluginAnimationManager.Snapshot,
        startAnimation?: boolean,
        camera?: {
            current: Camera.Snapshot,
            transitionStyle: CameraTransitionStyle,
            transitionDurationInMs?: number
        },
        cameraSnapshots?: CameraSnapshotManager.StateSnapshot,
        canvas3d?: {
            props?: Canvas3DProps
        },
        interactivity?: {
            props?: InteractivityManager.Props
        },
        durationInMs?: number
    }
}
