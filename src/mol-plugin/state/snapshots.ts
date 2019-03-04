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

class PluginStateSnapshotManager extends PluginComponent<{ current?: UUID | undefined, entries: OrderedMap<string, PluginStateSnapshotManager.Entry> }> {
    readonly events = {
        changed: this.ev()
    };

    getEntry(id: string) {
        return this.state.entries.get(id);
    }

    remove(id: string) {
        if (!this.state.entries.has(id)) return;
        this.updateState({
            current: this.state.current === id ? void 0 : this.state.current,
            entries: this.state.entries.delete(id)
        });
        this.events.changed.next();
    }

    add(e: PluginStateSnapshotManager.Entry) {
        this.updateState({ current: e.snapshot.id, entries: this.state.entries.set(e.snapshot.id, e) });
        this.events.changed.next();
    }

    clear() {
        if (this.state.entries.size === 0) return;
        this.updateState({ current: void 0, entries: OrderedMap<string, PluginStateSnapshotManager.Entry>() });
        this.events.changed.next();
    }

    setCurrent(id: string) {
        const e = this.getEntry(id);
        if (e) {
            this.updateState({ current: id as UUID });
            this.events.changed.next();
        }
        return e && e.snapshot;
    }

    getNextId(id: string | undefined, dir: -1 | 1) {
        const xs = this.state.entries;
        const keys = xs.keys();
        let k = keys.next();
        let prev = k.value;
        const fst = prev;
        while (!k.done) {
            k = keys.next();
            if (k.value === id && dir === -1) return prev;
            if (!k.done && prev === id && dir === 1) return k.value;
            if (!k.done) prev = k.value;
            else break;
        }
        if (dir === -1) return prev;
        return fst;
    }

    setRemoteSnapshot(snapshot: PluginStateSnapshotManager.RemoteSnapshot): PluginState.Snapshot | undefined {
        this.clear();
        const entries = this.state.entries.withMutations(m => {
            for (const e of snapshot.entries) {
                m.set(e.snapshot.id, e);
            }
        });
        const current = snapshot.current
            ? snapshot.current
            : snapshot.entries.length > 0
            ? snapshot.entries[0].snapshot.id
            : void 0;
        this.updateState({ current, entries });
        this.events.changed.next();
        if (!current) return;
        const ret = this.getEntry(current);
        return ret && ret.snapshot;
    }

    getRemoteSnapshot(options?: { name?: string, description?: string }): PluginStateSnapshotManager.RemoteSnapshot {
        // TODO: diffing and all that fancy stuff
        return {
            timestamp: +new Date(),
            name: options && options.name,
            description: options && options.description,
            current: this.state.current,
            entries: this.state.entries.valueSeq().toArray()
        };
    }

    constructor() {
        super({ current: void 0, entries: OrderedMap<string, PluginStateSnapshotManager.Entry>() });
    }
}

namespace PluginStateSnapshotManager {
    export interface Entry {
        timestamp: number,
        name?: string,
        description?: string,
        snapshot: PluginState.Snapshot
    }

    export function Entry(snapshot: PluginState.Snapshot, name?: string, description?: string): Entry {
        return { timestamp: +new Date(), name, snapshot, description };
    }

    export interface RemoteSnapshot {
        timestamp: number,
        name?: string,
        description?: string,
        current: UUID | undefined,
        entries: Entry[]
    }
}