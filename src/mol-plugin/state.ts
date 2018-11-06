/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from 'mol-state';
import { PluginStateObjects as SO } from './state/objects';

export { PluginState }

class PluginState {
    readonly data: State;

    getSnapshot(): PluginState.Snapshot {
        throw 'nyi';
    }

    setSnapshot(snapshot: PluginState.Snapshot) {
        throw 'nyi';
    }

    setDataSnapshot(snapshot: State.Snapshot) {
        throw 'nyi';
    }

    dispose() {
        this.data.dispose();
    }

    constructor(globalContext: unknown) {
        this.data = State.create(new SO.Root({ label: 'Root' }, { }), { globalContext });
    }
}

namespace PluginState {
    export interface Snapshot { }
}
