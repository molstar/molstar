/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from 'mol-canvas3d/camera';
import { OrderedMap } from 'immutable';
import { UUID } from 'mol-util';
import { RxEventHelper } from 'mol-util/rx-event-helper';

export { CameraSnapshotManager }

class CameraSnapshotManager {
    private ev = RxEventHelper.create();
    private _entries = OrderedMap<string, CameraSnapshotManager.Entry>().asMutable();

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

    add(e: CameraSnapshotManager.Entry) {
        this._entries.set(e.id, e);
        this.events.changed.next();
    }

    clear() {
        if (this._entries.size === 0) return;
        this._entries = OrderedMap<string, CameraSnapshotManager.Entry>().asMutable();
        this.events.changed.next();
    }

    getStateSnapshot(): CameraSnapshotManager.StateSnapshot {
        const entries: CameraSnapshotManager.Entry[] = [];
        this._entries.forEach(e => entries.push(e!));
        return { entries };
    }

    setStateSnapshot(state: CameraSnapshotManager.StateSnapshot ) {
        this._entries = OrderedMap<string, CameraSnapshotManager.Entry>().asMutable();
        for (const e of state.entries) {
            this._entries.set(e.id, e);
        }
        this.events.changed.next();
    }

    dispose() {
        this.ev.dispose();
    }
}

namespace CameraSnapshotManager {
    export interface Entry {
        id: UUID,
        name: string,
        description?: string,
        snapshot: Camera.Snapshot
    }

    export function Entry(name: string, snapshot: Camera.Snapshot, description?: string): Entry {
        return { id: UUID.create22(), name, snapshot, description };
    }

    export interface StateSnapshot {
        entries: Entry[]
    }
}