/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { List } from 'immutable';
import { UUID } from '../../mol-util';
import { PluginState } from '../../mol-plugin/state';
import { StatefulPluginComponent } from '../component';
import { PluginContext } from '../../mol-plugin/context';
import { utf8ByteCount, utf8Write } from '../../mol-io/common/utf8';
import { Asset } from '../../mol-util/assets';
import { zip } from '../../mol-util/zip/zip';
import { readFromFile } from '../../mol-util/data-source';
import { objectForEach } from '../../mol-util/object';
import { PLUGIN_VERSION } from '../../mol-plugin/version';

export { PluginStateSnapshotManager };

class PluginStateSnapshotManager extends StatefulPluginComponent<{
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

    async setStateSnapshot(snapshot: PluginStateSnapshotManager.StateSnapshot): Promise<PluginState.Snapshot | undefined> {
        if (snapshot.version !== PLUGIN_VERSION) {
            // TODO
            // console.warn('state snapshot version mismatch');
        }

        this.clear();
        const entries = List<PluginStateSnapshotManager.Entry>().asMutable();
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

    private syncCurrent(options?: { name?: string, description?: string, params?: PluginState.SnapshotParams }) {
        const snapshot = this.plugin.state.getSnapshot(options?.params);
        if (this.state.entries.size === 0 || !this.state.current) {
            this.add(PluginStateSnapshotManager.Entry(snapshot, { name: options?.name, description: options?.description }));
        } else {
            this.replace(this.state.current, snapshot);
        }
    }

    getStateSnapshot(options?: { name?: string, description?: string, playOnLoad?: boolean, params?: PluginState.SnapshotParams }): PluginStateSnapshotManager.StateSnapshot {
        // TODO: diffing and all that fancy stuff
        this.syncCurrent(options);

        return {
            timestamp: +new Date(),
            version: PLUGIN_VERSION,
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

    async serialize(options?: { type: 'json' | 'molj' | 'zip' | 'molx', params?: PluginState.SnapshotParams }) {
        const json = JSON.stringify(this.getStateSnapshot({ params: options?.params }), null, 2);

        if (!options?.type || options.type === 'json' || options.type === 'molj') {
            return new Blob([json], {type : 'application/json;charset=utf-8'});
        } else {
            const state = new Uint8Array(utf8ByteCount(json));
            utf8Write(state, 0, json);

            const zipDataObj: { [k: string]: Uint8Array } = {
                'state.json': state
            };

            const assets: [UUID, Asset][] = [];

            // TODO: there can be duplicate entries: check for this?
            for (const { asset, file } of this.plugin.managers.asset.assets) {
                assets.push([asset.id, asset]);
                zipDataObj[`assets/${asset.id}`] = new Uint8Array(await file.arrayBuffer());
            }

            if (assets.length > 0) {
                const index = JSON.stringify(assets, null, 2);
                const data = new Uint8Array(utf8ByteCount(index));
                utf8Write(data, 0, index);
                zipDataObj['assets.json'] = data;
            }

            const zipFile = zip(zipDataObj);
            return new Blob([zipFile], {type : 'application/zip'});
        }
    }

    async open(file: File) {
        try {
            const fn = file.name.toLowerCase();
            if (fn.endsWith('json') || fn.endsWith('molj')) {
                const data = await this.plugin.runTask(readFromFile(file, 'string'));
                const snapshot = JSON.parse(data);

                if (PluginStateSnapshotManager.isStateSnapshot(snapshot)) {
                    return this.setStateSnapshot(snapshot);
                } else if (PluginStateSnapshotManager.isStateSnapshot(snapshot.data)) {
                    return this.setStateSnapshot(snapshot.data);
                } else {
                    this.plugin.state.setSnapshot(snapshot);
                }
            } else {
                const data = await this.plugin.runTask(readFromFile(file, 'zip'));
                const assets = Object.create(null);

                objectForEach(data, (v, k) => {
                    if (k === 'state.json' || k === 'assets.json') return;
                    const name = k.substring(k.indexOf('/') + 1);
                    assets[name] = new File([v], name);
                });
                const stateFile = new File([data['state.json']], 'state.json');
                const stateData = await this.plugin.runTask(readFromFile(stateFile, 'string'));

                if (data['assets.json']) {
                    const file = new File([data['assets.json']], 'assets.json');
                    const json = JSON.parse(await this.plugin.runTask(readFromFile(file, 'string')));

                    for (const [id, asset] of json) {
                        this.plugin.managers.asset.set(asset, assets[id]);
                    }
                }

                const snapshot = JSON.parse(stateData);
                return this.setStateSnapshot(snapshot);
            }
        } catch (e) {
            this.plugin.log.error(`Reading state: ${e}`);
        }
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
            this.plugin.managers.animation.stop();
        } else {
            this.play();
        }
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

    export function Entry(snapshot: PluginState.Snapshot, params: { name?: string, description?: string }): Entry {
        return { timestamp: +new Date(), snapshot, ...params };
    }

    export function isStateSnapshot(x?: any): x is StateSnapshot {
        const s = x as StateSnapshot;
        return !!s && !!s.timestamp && !!s.entries;
    }

    export interface StateSnapshot {
        timestamp: number,
        version: string,
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