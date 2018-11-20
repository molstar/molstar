/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from 'mol-state';
import { PluginStateObject as SO } from './state/objects';
import { Camera } from 'mol-canvas3d/camera';
import { PluginBehavior } from './behavior';
import { CameraSnapshotManager } from './state/camera';
import { PluginStateSnapshotManager } from './state/snapshots';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { Canvas3DParams } from 'mol-canvas3d/canvas3d';
import { ParamDefinition } from 'mol-util/param-definition';
import { PluginCommands } from './command';
export { PluginState }

class PluginState {
    private ev = RxEventHelper.create();

    readonly dataState: State;
    readonly behaviorState: State;
    readonly cameraSnapshots = new CameraSnapshotManager();

    readonly snapshots = new PluginStateSnapshotManager();

    readonly behavior = {
        kind: this.ev.behavior<PluginState.Kind>('data'),
        currentObject: this.ev.behavior<State.ObjectEvent>({} as any)
    }

    setKind(kind: PluginState.Kind) {
        const current = this.behavior.kind.value;
        if (kind !== current) {
            this.behavior.kind.next(kind);
            this.behavior.currentObject.next(kind === 'data'
                ? this.dataState.behaviors.currentObject.value
                : this.behaviorState.behaviors.currentObject.value)
        }
    }

    getSnapshot(): PluginState.Snapshot {
        return {
            data: this.dataState.getSnapshot(),
            behaviour: this.behaviorState.getSnapshot(),
            cameraSnapshots: this.cameraSnapshots.getStateSnapshot(),
            canvas3d: {
                camera: this.plugin.canvas3d.camera.getSnapshot(),
                viewport: this.plugin.canvas3d.props
            }
        };
    }

    async setSnapshot(snapshot: PluginState.Snapshot) {
        await this.plugin.runTask(this.behaviorState.setSnapshot(snapshot.behaviour));
        await this.plugin.runTask(this.dataState.setSnapshot(snapshot.data));
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: snapshot.canvas3d.viewport || { } });
        this.cameraSnapshots.setStateSnapshot(snapshot.cameraSnapshots);
        this.plugin.canvas3d.camera.setState(snapshot.canvas3d.camera);
        this.plugin.canvas3d.requestDraw(true);
    }

    dispose() {
        this.ev.dispose();
        this.dataState.dispose();
        this.behaviorState.dispose();
        this.cameraSnapshots.dispose();
    }

    constructor(private plugin: import('./context').PluginContext) {
        this.dataState = State.create(new SO.Root({ }), { globalContext: plugin });
        this.behaviorState = State.create(new PluginBehavior.Root({ }), { globalContext: plugin });

        this.dataState.behaviors.currentObject.subscribe(o => {
            if (this.behavior.kind.value === 'data') this.behavior.currentObject.next(o);
        });
        this.behaviorState.behaviors.currentObject.subscribe(o => {
            if (this.behavior.kind.value === 'behavior') this.behavior.currentObject.next(o);
        });

        this.behavior.currentObject.next(this.dataState.behaviors.currentObject.value);
    }
}

namespace PluginState {
    export type Kind = 'data' | 'behavior'

    export interface Snapshot {
        data: State.Snapshot,
        behaviour: State.Snapshot,
        cameraSnapshots: CameraSnapshotManager.StateSnapshot,
        canvas3d: {
            camera: Camera.Snapshot,
            viewport: Canvas3DParams
        }
    }
}
