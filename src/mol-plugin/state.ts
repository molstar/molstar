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

export { PluginState }

class PluginState {
    readonly data: State;
    readonly behavior: State;
    readonly cameraSnapshots: CameraSnapshotManager = new CameraSnapshotManager();

    getSnapshot(): PluginState.Snapshot {
        return {
            data: this.data.getSnapshot(),
            behaviour: this.behavior.getSnapshot(),
            cameraSnapshots: this.cameraSnapshots.getStateSnapshot(),
            canvas3d: {
                camera: this.plugin.canvas3d.camera.getSnapshot()
            }
        };
    }

    async setSnapshot(snapshot: PluginState.Snapshot) {
        await this.plugin.runTask(this.behavior.setSnapshot(snapshot.behaviour));
        await this.plugin.runTask(this.data.setSnapshot(snapshot.data));
        this.cameraSnapshots.setStateSnapshot(snapshot.cameraSnapshots);
        this.plugin.canvas3d.camera.setState(snapshot.canvas3d.camera);
        this.plugin.canvas3d.requestDraw(true);
    }

    dispose() {
        this.data.dispose();
        this.behavior.dispose();
        this.cameraSnapshots.dispose();
    }

    constructor(private plugin: import('./context').PluginContext) {
        this.data = State.create(new SO.Root({ }), { globalContext: plugin });
        this.behavior = State.create(new PluginBehavior.Root({ }), { globalContext: plugin });
    }
}

namespace PluginState {
    export interface Snapshot {
        data: State.Snapshot,
        behaviour: State.Snapshot,
        cameraSnapshots: CameraSnapshotManager.StateSnapshot,
        canvas3d: {
            camera: Camera.Snapshot
        }
    }
}
