/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from 'mol-canvas3d/camera';
import { OrderedMap } from 'immutable';
import { UUID } from 'mol-util';
import { PluginComponent } from 'mol-plugin/component';
import { PluginContext } from 'mol-plugin/context';

export { CameraSnapshotManager }

class CameraSnapshotManager extends PluginComponent<{ entries: OrderedMap<string, CameraSnapshotManager.Entry> }> {
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

    add(e: CameraSnapshotManager.Entry) {
        this.updateState({ entries: this.state.entries.set(e.id, e) });
        this.events.changed.next();
    }

    clear() {
        if (this.state.entries.size === 0) return;
        this.updateState({ entries: OrderedMap<string, CameraSnapshotManager.Entry>() });
        this.events.changed.next();
    }

    getStateSnapshot(): CameraSnapshotManager.StateSnapshot {
        const entries: CameraSnapshotManager.Entry[] = [];
        this.state.entries.forEach(e => entries.push(e!));
        return { entries };
    }

    setStateSnapshot(state: CameraSnapshotManager.StateSnapshot ) {
        const entries = OrderedMap<string, CameraSnapshotManager.Entry>().asMutable();
        for (const e of state.entries) {
            entries.set(e.id, e);
        }
        this.updateState({ entries: entries.asImmutable() });
        this.events.changed.next();
    }

    constructor(ctx: PluginContext) {
        super(ctx, { entries: OrderedMap<string, CameraSnapshotManager.Entry>() });
    }
}

namespace CameraSnapshotManager {
    export interface Entry {
        id: UUID,
        timestamp: string,
        name?: string,
        description?: string,
        snapshot: Camera.Snapshot
    }

    export function Entry(snapshot: Camera.Snapshot, name?: string, description?: string): Entry {
        return { id: UUID.create22(), timestamp: new Date().toLocaleString(), name, snapshot, description };
    }

    export interface StateSnapshot {
        entries: Entry[]
    }
}