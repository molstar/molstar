/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { Zip } from '../../mol-util/zip/zip';
import { readFromFile } from '../../mol-util/data-source';
import { objectForEach } from '../../mol-util/object';
import { PLUGIN_VERSION } from '../../mol-plugin/version';
import { canvasToBlob } from '../../mol-canvas3d/util';
import { Task } from '../../mol-task';
import { StringLike } from '../../mol-io/common/string-like';
import { SingleTaskQueue } from '../../mol-util/single-task-queue';

export { PluginStateSnapshotManager };

interface StateManagerState {
    current?: UUID,
    currentAnimationFrame?: number,
    entries: List<PluginStateSnapshotManager.Entry>,
    isPlaying: boolean,
    nextSnapshotDelayInMs: number
}

class PluginStateSnapshotManager extends StatefulPluginComponent<StateManagerState> {
    static DefaultNextSnapshotDelayInMs = 1500;

    private entryMap = new Map<string, PluginStateSnapshotManager.Entry>();
    private defaultSnapshotId: UUID | undefined = undefined;

    protected updateState(state: Partial<StateManagerState>) {
        if ('current' in state && !('curentAnimationFrame' in state)) {
            return super.updateState({ ...state, currentAnimationFrame: 0 });
        } else {
            return super.updateState(state);
        }
    }

    readonly events = {
        changed: this.ev(),
        opened: this.ev(),
    };

    get current() {
        const id = this.state.current;
        return this.state.entries.find(e => e.snapshot.id === id);
    }

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

        if (e?.image) this.plugin.managers.asset.delete(e.image);
        this.entryMap.delete(id);
        this.updateState({
            current: this.state.current === id ? void 0 : this.state.current,
            entries: this.state.entries.delete(this.getIndex(e))
        });
        this.events.changed.next(void 0);
    }

    add(e: PluginStateSnapshotManager.Entry) {
        this.entryMap.set(e.snapshot.id, e);
        this.updateState({ current: e.snapshot.id, entries: this.state.entries.push(e) });
        this.events.changed.next(void 0);
    }

    replace(id: string, snapshot: PluginState.Snapshot, params?: PluginStateSnapshotManager.EntryParams) {
        const old = this.getEntry(id);
        if (!old) return;

        this.defaultSnapshotId = undefined;

        if (old?.image) this.plugin.managers.asset.delete(old.image);
        const idx = this.getIndex(old);
        // The id changes here!
        const e = PluginStateSnapshotManager.Entry(snapshot, {
            key: params?.key ?? old.key,
            name: params?.name ?? old.name,
            description: params?.description ?? old.description,
            descriptionFormat: params?.descriptionFormat ?? old.descriptionFormat,
            image: params?.image,
        });
        this.entryMap.set(snapshot.id, e);
        this.updateState({ current: e.snapshot.id, entries: this.state.entries.set(idx, e) });
        this.events.changed.next(void 0);
    }

    move(id: string, dir: -1 | 1) {
        const len = this.state.entries.size;
        if (len < 2) return;

        const e = this.getEntry(id);
        if (!e) return;
        const from = this.getIndex(e);
        let to = (from + dir) % len;
        if (to < 0) to += len;
        const f = this.state.entries.get(to)!;

        const entries = this.state.entries.asMutable();
        entries.set(to, e);
        entries.set(from, f);

        this.updateState({ current: e.snapshot.id, entries: entries.asImmutable() });
        this.events.changed.next(void 0);
    }

    update(e: PluginStateSnapshotManager.Entry, options: { key?: string, name?: string, description?: string, descriptionFormat?: PluginStateSnapshotManager.DescriptionFormat }) {
        const idx = this.getIndex(e);
        if (idx < 0) return;
        const entries = this.state.entries.set(idx, {
            ...e,
            key: options.key?.trim() || undefined,
            name: options.name?.trim() || undefined,
            description: options.description?.trim() || undefined,
            descriptionFormat: options.descriptionFormat,
        });
        this.updateState({ entries });
        this.entryMap.set(e.snapshot.id, this.state.entries.get(idx)!);
        this.events.changed.next(void 0);
    }

    clear() {
        if (this.state.entries.size === 0) return;

        this.entryMap.forEach(e => {
            if (e?.image) this.plugin.managers.asset.delete(e.image);
        });
        this.entryMap.clear();
        this.updateState({ current: void 0, entries: List<PluginStateSnapshotManager.Entry>() });
        this.events.changed.next(void 0);
    }

    applyKey(key: string) {
        const e = this.state.entries.find(e => e.key === key);
        if (!e) return;

        this.updateState({ current: e.snapshot.id as UUID });
        this.events.changed.next(void 0);
        this.plugin.state.setSnapshot(e.snapshot);
    }

    setCurrent(id: string) {
        const e = this.getEntry(id);
        if (e) {
            this.updateState({ current: id as UUID });
            this.events.changed.next(void 0);
        }
        return e && e.snapshot;
    }

    private animationFrameQueue = new SingleTaskQueue();
    setSnapshotAnimationFrame(frame: number, load = false) {
        if (this.updateState({ currentAnimationFrame: frame })) {
            this.events.changed.next(void 0);
        }

        if (load) {
            this.animationFrameQueue.run(() => {
                const entry = this.getEntry(this.state.current);
                if (!entry) return Promise.resolve();
                return this.plugin.state.setAnimationSnapshot(entry.snapshot, frame);
            });
        }
    }

    getNextId(id: string | undefined, dir: -1 | 1) {
        const len = this.state.entries.size;
        if (!id) {
            if (len === 0) return void 0;
            const idx = dir === -1 ? len - 1 : 0;
            return this.state.entries.get(idx)!.snapshot.id;
        }

        const e = this.getEntry(id);
        if (!e) return;
        let idx = this.getIndex(e);
        if (idx < 0) return;

        idx = (idx + dir) % len;
        if (idx < 0) idx += len;

        return this.state.entries.get(idx)!.snapshot.id;
    }

    applyNext(dir: -1 | 1) {
        const next = this.getNextId(this.state.current, dir);
        if (next) {
            const snapshot = this.setCurrent(next);
            if (snapshot) return this.plugin.state.setSnapshot(snapshot);
        }
    }

    async setStateSnapshot(snapshot: PluginStateSnapshotManager.StateSnapshot): Promise<PluginState.Snapshot | undefined> {
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
        this.events.changed.next(void 0);
        if (!current) return;
        const entry = this.getEntry(current);
        const next = entry && entry.snapshot;
        if (!next) return;
        await this.plugin.state.setSnapshot(next);
        if (snapshot.playback && snapshot.playback.isPlaying) this.play(true);
        return next;
    }

    private async syncCurrent(options?: { name?: string, description?: string, descriptionFormat?: PluginStateSnapshotManager.DescriptionFormat, params?: PluginState.SnapshotParams }) {
        const isEmpty = this.state.entries.size === 0;
        const canReplace = this.state.entries.size === 1 && this.state.current && (!this.defaultSnapshotId || this.state.current === this.defaultSnapshotId);
        if (!isEmpty && !canReplace) return;

        const snapshot = this.plugin.state.getSnapshot(options?.params);
        const image = (options?.params?.image ?? this.plugin.state.snapshotParams.value.image) ? await PluginStateSnapshotManager.getCanvasImageAsset(this.plugin, `${snapshot.id}-image.png`) : undefined;

        if (isEmpty) {
            this.add(PluginStateSnapshotManager.Entry(snapshot, { name: options?.name, description: options?.description, descriptionFormat: options?.descriptionFormat, image }));
        } else if (canReplace) {
            // Replace the current state only if there is a single snapshot that has been created automatically
            const current = this.getEntry(this.state.current);
            if (current?.image) this.plugin.managers.asset.delete(current.image);
            this.replace(this.state.current!, snapshot, { image });
        }

        this.defaultSnapshotId = snapshot.id;
    }

    async getStateSnapshot(options?: { name?: string, description?: string, playOnLoad?: boolean, params?: PluginState.SnapshotParams }): Promise<PluginStateSnapshotManager.StateSnapshot> {
        await this.syncCurrent(options);

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
        const json = JSON.stringify(await this.getStateSnapshot({ params: options?.params }), null, 2);

        if (!options?.type || options.type === 'json' || options.type === 'molj') {
            return new Blob([json], { type: 'application/json;charset=utf-8' });
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

            const zipFile = await this.plugin.runTask(Zip(zipDataObj));
            return new Blob([zipFile], { type: 'application/zip' });
        }
    }

    async open(file: File) {
        try {
            const fn = file.name.toLowerCase();
            if (fn.endsWith('json') || fn.endsWith('molj')) {
                const data = await this.plugin.runTask(readFromFile(file, 'string'));
                const snapshot = JSON.parse(StringLike.toString(data));

                if (PluginStateSnapshotManager.isStateSnapshot(snapshot)) {
                    await this.setStateSnapshot(snapshot);
                } else if (PluginStateSnapshotManager.isStateSnapshot(snapshot.data)) {
                    await this.setStateSnapshot(snapshot.data);
                } else {
                    await this.plugin.state.setSnapshot(snapshot);
                }
            } else {
                const data = await this.plugin.runTask(readFromFile(file, 'zip'));
                const assetData = Object.create(null);

                objectForEach(data, (v, k) => {
                    if (k === 'state.json' || k === 'assets.json') return;
                    const name = k.substring(k.indexOf('/') + 1);
                    assetData[name] = v;
                });
                const stateFile = new File([data['state.json']], 'state.json');
                const snapshot = await this.plugin.runTask(readFromFile(stateFile, 'json'));

                if (data['assets.json']) {
                    const file = new File([data['assets.json']], 'assets.json');
                    const json = await this.plugin.runTask(readFromFile(file, 'json'));

                    for (const [id, asset] of json) {
                        this.plugin.managers.asset.set(asset, new File([assetData[id]], asset.name));
                    }
                }

                await this.setStateSnapshot(snapshot);
            }
            this.events.opened.next(void 0);
        } catch (e) {
            console.error(e);
            this.plugin.log.error('Error reading state');
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

    private async startPlayback() {
        const { current } = this;
        if (!current) return;

        // If there is a transition associated with the current snapshot, replay it
        if (current.snapshot.transition) {
            const snapshot = this.setCurrent(this.state.current!)!;
            await this.plugin.state.setSnapshot(snapshot);
            const delay = typeof snapshot.durationInMs !== 'undefined' ? snapshot.durationInMs : this.state.nextSnapshotDelayInMs;
            if (this.state.isPlaying) this.timeoutHandle = setTimeout(this.next, delay);
        } else {
            return this.next();
        }
    }

    play(delayFirst: boolean = false) {
        this.updateState({ isPlaying: true });

        if (delayFirst) {
            const e = this.getEntry(this.state.current);
            if (!e) {
                this.next();
                return;
            }
            this.events.changed.next(void 0);
            const snapshot = e.snapshot;
            const delay = typeof snapshot.durationInMs !== 'undefined' ? snapshot.durationInMs : this.state.nextSnapshotDelayInMs;
            this.timeoutHandle = setTimeout(this.next, delay);
        } else {
            this.startPlayback();
        }
    }

    stop() {
        this.plugin.managers.animation.stop();
        this.updateState({ isPlaying: false });
        if (typeof this.timeoutHandle !== 'undefined') clearTimeout(this.timeoutHandle);
        this.timeoutHandle = void 0;
        this.events.changed.next(void 0);
    }

    togglePlay() {
        if (this.state.isPlaying) {
            this.stop();
            this.plugin.managers.animation.stop();
        } else {
            this.play();
        }
    }

    dispose() {
        super.dispose();
        this.entryMap.clear();
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
    export type DescriptionFormat = 'markdown' | 'plaintext';

    export interface EntryParams {
        key?: string,
        name?: string,
        /** Information about the snapshot, to be shown in the UI when the snapshot is loaded. */
        description?: string,
        /** Toggle between markdown and plaintext interpretation of `description`. Default is markdown. */
        descriptionFormat?: DescriptionFormat,
        image?: Asset,
    }

    export interface Entry extends EntryParams {
        timestamp: number,
        snapshot: PluginState.Snapshot
    }

    export function Entry(snapshot: PluginState.Snapshot, params: EntryParams): Entry {
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

    export async function getCanvasImageAsset(ctx: PluginContext, name: string): Promise<Asset | undefined> {
        const task = Task.create('Render Screenshot', async runtime => {
            if (!ctx.helpers.viewportScreenshot) return;
            const p = await ctx.helpers.viewportScreenshot.getPreview(runtime, 512);
            if (!p) return;

            const blob = await canvasToBlob(p.canvas, 'png');
            const file = new File([blob], name);
            const image: Asset = { kind: 'file', id: UUID.create22(), name };
            ctx.managers.asset.set(image, file);
            return image;
        });
        return ctx.runTask(task);
    }
}