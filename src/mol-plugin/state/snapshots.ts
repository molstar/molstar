/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedMap } from 'immutable';
import { UUID } from 'mol-util';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { PluginState } from '../state';

export { PluginStateSnapshotManager }

class PluginStateSnapshotManager {
    private ev = RxEventHelper.create();
    private _entries = OrderedMap<string, PluginStateSnapshotManager.Entry>().asMutable();

    readonly events = {
        changed: this.ev()
    };

    get entries() { return this._entries; }

    getEntry(id: string) {
        return this._entries.get(id);
    }

    remove(id: string) {
        if (!this._entries.has(id)) return;
        this._entries.delete(id);
        this.events.changed.next();
    }

    add(e: PluginStateSnapshotManager.Entry) {
        this._entries.set(e.id, e);
        this.events.changed.next();
    }

    clear() {
        if (this._entries.size === 0) return;
        this._entries = OrderedMap<string, PluginStateSnapshotManager.Entry>().asMutable();
        this.events.changed.next();
    }

    dispose() {
        this.ev.dispose();
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