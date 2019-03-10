/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { List } from 'immutable';
import { UUID } from 'mol-util';
import { PluginState } from '../state';
import { PluginComponent } from 'mol-plugin/component';
import { PluginContext } from 'mol-plugin/context';

export { PluginStateSnapshotManager }

class PluginStateSnapshotManager extends PluginComponent<{
    current?: UUID | undefined,
    entries: List<PluginStateSnapshotManager.Entry>,
    isPlaying: boolean,
    nextSnapshotDelayInMs: number
}> {
    static DefaultNextSnapshotDelayInMs = 1500;

    private entryMap = new Map<string, PluginStateSnapshotManager.Entry>();

    readonly events = {
        changed: this.ev()
    };

    currentGetSnapshotParams: PluginState.GetSnapshotParams = PluginState.DefaultGetSnapshotParams as any;

    getIndex(e: PluginStateSnapshotManager.Entry) {
        return this.state.entries.indexOf(e);
    }

    getEntry(id: string | undefined) {
        if (!id) return;
        return this.entryMap.get(id);
    }

    remove(id: string) {
        const e = this.entryMap.get(id);
        if (!e) return;

        this.entryMap.delete(id);
        this.updateState({
            current: this.state.current === id ? void 0 : this.state.current,
            entries: this.state.entries.delete(this.getIndex(e))
        });
        this.events.changed.next();
    }

    add(e: PluginStateSnapshotManager.Entry) {
        this.entryMap.set(e.snapshot.id, e);
        this.updateState({ current: e.snapshot.id, entries: this.state.entries.push(e) });
        this.events.changed.next();
    }

    replace(id: string, snapshot: PluginState.Snapshot) {
        const old = this.getEntry(id);
        if (!old) return;

        const idx = this.getIndex(old);
        // The id changes here!
        const e = PluginStateSnapshotManager.Entry(snapshot, {
            name: old.name,
            description: old.description
        });
        this.entryMap.set(snapshot.id, e);
        this.updateState({ current: e.snapshot.id, entries: this.state.entries.set(idx, e) });
        this.events.changed.next();
    }

    move(id: string, dir: -1 | 1) {
        const len = this.state.entries.size;
        if (len < 2) return;

        const e = this.getEntry(id);
        if (!e) return;
        const from = this.getIndex(e);
        let to = (from + dir) % len;
        if (to < 0) to += len;
        const f = this.state.entries.get(to);

        const entries = this.state.entries.asMutable();
        entries.set(to, e);
        entries.set(from, f);

        this.updateState({ current: e.snapshot.id, entries: entries.asImmutable() });
        this.events.changed.next();
    }

    clear() {
        if (this.state.entries.size === 0) return;
        this.entryMap.clear();
        this.updateState({ current: void 0, entries: List<PluginStateSnapshotManager.Entry>() });
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
        const len = this.state.entries.size;
        if (!id) {
            if (len === 0) return void 0;
            const idx = dir === -1 ? len - 1 : 0;
            return this.state.entries.get(idx).snapshot.id;
        }

        const e = this.getEntry(id);
        if (!e) return;
        let idx = this.getIndex(e);
        if (idx < 0) return;

        idx = (idx + dir) % len;
        if (idx < 0) idx += len;

        return this.state.entries.get(idx).snapshot.id;
    }

    async setRemoteSnapshot(snapshot: PluginStateSnapshotManager.RemoteSnapshot): Promise<PluginState.Snapshot | undefined> {
        this.clear();
        const entries = List<PluginStateSnapshotManager.Entry>().asMutable()
        for (const e of snapshot.entries) {
            this.entryMap.set(e.snapshot.id, e);
            entries.push(e);
        }
        const current = snapshot.current
            ? snapshot.current
            : snapshot.entries.length > 0
            ? snapshot.entries[0].snapshot.id
            : void 0;
        this.updateState({
            current,
            entries: entries.asImmutable(),
            isPlaying: false,
            nextSnapshotDelayInMs: snapshot.playback ? snapshot.playback.nextSnapshotDelayInMs : PluginStateSnapshotManager.DefaultNextSnapshotDelayInMs
        });
        this.events.changed.next();
        if (!current) return;
        const entry = this.getEntry(current);
        const next = entry && entry.snapshot;
        if (!next) return;
        await this.plugin.state.setSnapshot(next);
        if (snapshot.playback && snapshot.playback.isPlaying) this.play(true);
        return next;
    }

    getRemoteSnapshot(options?: { name?: string, description?: string, playOnLoad?: boolean }): PluginStateSnapshotManager.RemoteSnapshot {
        // TODO: diffing and all that fancy stuff
        return {
            timestamp: +new Date(),
            name: options && options.name,
            description: options && options.description,
            current: this.state.current,
            playback: {
                isPlaying: !!(options && options.playOnLoad),
                nextSnapshotDelayInMs: this.state.nextSnapshotDelayInMs
            },
            entries: this.state.entries.valueSeq().toArray()
        };
    }

    private timeoutHandle: any = void 0;
    private next = async () => {
        this.timeoutHandle = void 0;
        const next = this.getNextId(this.state.current, 1);
        if (!next || next === this.state.current) {
            this.stop();
            return;
        }
        const snapshot = this.setCurrent(next)!;
        await this.plugin.state.setSnapshot(snapshot);
        const delay = typeof snapshot.durationInMs !== 'undefined' ? snapshot.durationInMs : this.state.nextSnapshotDelayInMs;
        if (this.state.isPlaying) this.timeoutHandle = setTimeout(this.next, delay);
    };

    play(delayFirst: boolean = false) {
        this.updateState({ isPlaying: true });

        if (delayFirst) {
            const e = this.getEntry(this.state.current);
            if (!e) {
                this.next();
                return;
            }
            this.events.changed.next();
            const snapshot = e.snapshot;
            const delay = typeof snapshot.durationInMs !== 'undefined' ? snapshot.durationInMs : this.state.nextSnapshotDelayInMs;
            this.timeoutHandle = setTimeout(this.next, delay);
        } else {
            this.next();
        }
    }

    stop() {
        this.updateState({ isPlaying: false });
        if (typeof this.timeoutHandle !== 'undefined') clearTimeout(this.timeoutHandle);
        this.timeoutHandle = void 0;
        this.events.changed.next();
    }

    togglePlay() {
        if (this.state.isPlaying) {
            this.stop();
            this.plugin.state.animation.stop();
        }
        else this.play();
    }

    constructor(private plugin: PluginContext) {
        super({
            current: void 0,
            entries: List(),
            isPlaying: false,
            nextSnapshotDelayInMs: PluginStateSnapshotManager.DefaultNextSnapshotDelayInMs
        });
        // TODO make nextSnapshotDelayInMs editable
    }
}

namespace PluginStateSnapshotManager {
    export interface Entry {
        timestamp: number,
        name?: string,
        description?: string,
        snapshot: PluginState.Snapshot
    }

    export function Entry(snapshot: PluginState.Snapshot, params: {name?: string, description?: string }): Entry {
        return { timestamp: +new Date(), snapshot, ...params };
    }

    export interface RemoteSnapshot {
        timestamp: number,
        name?: string,
        description?: string,
        current: UUID | undefined,
        playback: {
            isPlaying: boolean,
            nextSnapshotDelayInMs: number,
        },
        entries: Entry[]
    }
}