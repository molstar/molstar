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
    readonly behaviour: State;

    readonly canvas = {
        camera: CombinedCamera.create()
    };

    getSnapshot(): PluginState.Snapshot {
        return {
            data: this.data.getSnapshot(),
            behaviour: this.behaviour.getSnapshot(),
            canvas: {
                camera: { ...this.canvas.camera }
            }
        };
    }

    setSnapshot(snapshot: PluginState.Snapshot) {
        // TODO events
        this.behaviour.setSnapshot(snapshot.behaviour);
        this.data.setSnapshot(snapshot.data);
        this.canvas.camera = { ...snapshot.canvas.camera };
    }

    async updateData(tree: StateTree) {
        // TODO: "task observer"
        await this.data.update(tree).run(p => console.log(p), 250);
    }

    dispose() {
        this.data.dispose();
    }

    constructor(globalContext: unknown) {
        this.data = State.create(new SO.Root({ label: 'Root' }, { }), { globalContext });
        this.behaviour = State.create(new SO.Root({ label: 'Root' }, { }), { globalContext });
    }
}

namespace PluginState {
    export interface Snapshot {
        data: State.Snapshot,
        behaviour: State.Snapshot,
        canvas: PluginState['canvas']
    }
}
