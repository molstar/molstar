/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTree } from 'mol-state';
import { PluginStateObjects as SO } from './state/objects';
import { CombinedCamera } from 'mol-canvas3d/camera/combined';

export { PluginState }

class PluginState {
    readonly data: State;
    readonly behavior: State;

    getSnapshot(): PluginState.Snapshot {
        return {
            data: this.data.getSnapshot(),
            behaviour: this.behavior.getSnapshot(),
            canvas3d: {
                camera: { ...this.plugin.canvas3d.camera }
            }
        };
    }

    setSnapshot(snapshot: PluginState.Snapshot) {
        this.behavior.setSnapshot(snapshot.behaviour);
        this.data.setSnapshot(snapshot.data);

        // TODO: handle camera
        // console.log({ old: { ...this.plugin.canvas3d.camera  }, new: snapshot.canvas3d.camera });
        // CombinedCamera.copy(snapshot.canvas3d.camera, this.plugin.canvas3d.camera);
        // CombinedCamera.update(this.plugin.canvas3d.camera);
        // this.plugin.canvas3d.center
        // console.log({ copied: { ...this.plugin.canvas3d.camera  } });
        // this.plugin.canvas3d.requestDraw(true);
        // console.log('updated camera');
    }

    updateData(tree: StateTree) {
        return this.plugin.runTask(this.data.update(tree));
    }

    updateBehaviour(tree: StateTree) {
        return this.plugin.runTask(this.behavior.update(tree));
    }

    dispose() {
        this.data.dispose();
    }

    constructor(private plugin: import('./context').PluginContext) {
        this.data = State.create(new SO.DataRoot({ label: 'Root' }, { }), { globalContext: plugin });
        this.behavior = State.create(new SO.BehaviorRoot({ label: 'Root' }, { }), { globalContext: plugin });
    }
}

namespace PluginState {
    export interface Snapshot {
        data: State.Snapshot,
        behaviour: State.Snapshot,
        canvas3d: {
            camera: CombinedCamera
        }
    }
}
