/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedMap } from 'immutable';
import { UUID } from 'mol-util';
import { PluginState } from '../state';
import { PluginComponent } from 'mol-plugin/component';

export { PluginStateSnapshotManager }

class PluginStateSnapshotManager extends PluginComponent<{ entries: OrderedMap<string, PluginStateSnapshotManager.Entry> }> {
    readonly events = {
        changed: this.ev()
    };

    getEntry(id: string) {
        return this.state.entries.get(id);
    }

    remove(id: string) {
        if (!this.state.entries.has(id)) return;
        this.updateState({ entries: this.state.entries.delete(id) });
        this.events.changed.next();
    }

    add(e: PluginStateSnapshotManager.Entry) {
        this.updateState({ entries: this.state.entries.set(e.id, e) });
        this.events.changed.next();
    }

    clear() {
        if (this.state.entries.size === 0) return;
        this.updateState({ entries: OrderedMap<string, PluginStateSnapshotManager.Entry>() });
        this.events.changed.next();
    }

    constructor() {
        super({ entries: OrderedMap<string, PluginStateSnapshotManager.Entry>() });
    }
}

namespace PluginStateSnapshotManager {
    export interface Entry {
        id: UUID,
        timestamp: string,
        name?: string,
        description?: string,
        snapshot: PluginState.Snapshot
    }

    export function Entry(snapshot: PluginState.Snapshot, name?: string, description?: string): Entry {
        return { id: UUID.create22(), timestamp: new Date().toLocaleString(), name, snapshot, description };
    }

    export interface StateSnapshot {
        entries: Entry[]
    }
}