/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTransform, StateTransformer } from '../mol-state';
import { PluginStateObject as SO } from '../mol-plugin-state/objects';
import { Camera } from '../mol-canvas3d/camera';
import { PluginBehavior } from './behavior';
import { Canvas3DProps } from '../mol-canvas3d/canvas3d';
import { PluginCommands } from './commands';
import { PluginAnimationManager } from '../mol-plugin-state/manager/animation';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { UUID } from '../mol-util';
import { InteractivityManager } from '../mol-plugin-state/manager/interactivity';
import { produce } from 'immer';
import { StructureFocusSnapshot } from '../mol-plugin-state/manager/structure/focus';
import { merge } from 'rxjs';
import { PluginContext } from './context';
import { PluginComponent } from '../mol-plugin-state/component';

export { PluginState };

class PluginState extends PluginComponent {
    private get animation() { return this.plugin.managers.animation; }

    readonly data = State.create(new SO.Root({ }), { runTask: this.plugin.runTask, globalContext: this.plugin });
    readonly behaviors = State.create(new PluginBehavior.Root({ }), { runTask: this.plugin.runTask, globalContext: this.plugin, rootState: { isLocked: true } });

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
    }

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
            canvas3d: p.canvas3d ? { props: this.plugin.canvas3d?.props } : void 0,
            interactivity: p.interactivity ? { props: this.plugin.managers.interactivity.props } : void 0,
            structureFocus: this.plugin.managers.structure.focus.getSnapshot(),
            durationInMs: p?.durationInMs
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
        if (snapshot.structureFocus) {
            this.plugin.managers.structure.focus.setSnapshot(snapshot.structureFocus);
        }
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

    updateTransform(state: State, a: StateTransform.Ref, params: any, canUndo?: string | boolean) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update(this.plugin, { state, tree, options: { canUndo } });
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
                cell.obj.data.unregister();
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
        animation: PD.Boolean(true),
        startAnimation: PD.Boolean(false),
        canvas3d: PD.Boolean(true),
        interactivity: PD.Boolean(true),
        camera: PD.Boolean(true),
        cameraTransition: PD.MappedStatic('animate', {
            animate: PD.Group({
                durationInMs: PD.Numeric(250, { min: 100, max: 5000, step: 500 }, { label: 'Duration in ms' }),
            }),
            instant: PD.Group({ })
        }, { options: [['animate', 'Animate'], ['instant', 'Instant']] })
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
            current: Camera.Snapshot,
            transitionStyle: CameraTransitionStyle,
            transitionDurationInMs?: number
        },
        canvas3d?: {
            props?: Canvas3DProps
        },
        interactivity?: {
            props?: InteractivityManager.Props
        },
        structureFocus?: StructureFocusSnapshot,
        durationInMs?: number
    }

    export type SnapshotType = 'json' | 'molj' | 'zip' | 'molx'
}
