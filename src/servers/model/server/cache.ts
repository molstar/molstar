/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ConsoleLogger } from '../../../mol-util/console-logger';
import { LinkedList } from '../../../mol-data/generic';
import { ModelServerConfig as ServerConfig } from '../config';

interface CacheEntry<T> {
    key: string,
    approximateSize: number,
    timeoutId: NodeJS.Timer | undefined,
    item: T
}

type CacheNode<T> = LinkedList.Node<CacheEntry<T>>

export interface CacheParams {
    useCache: boolean,
    maxApproximateSizeInBytes: number, // = 2 * 1014 * 1024 * 1024; // 2 GB
    entryTimeoutInMs: number // = 10 * 60 * 1000; // 10 minutes
}

export class Cache<T> {
    private entries = LinkedList<CacheEntry<T>>();
    private entryMap = new Map<string, CacheNode<T>>();
    private approximateSize = 0;

    private clearTimeout(e: CacheNode<T>) {
        if (typeof e.value.timeoutId !== 'undefined') {
            clearTimeout(e.value.timeoutId);
            e.value.timeoutId = void 0;
        }
    }

    private dispose(e: CacheNode<T>) {
        this.clearTimeout(e);

        if (e.inList) {
            this.entries.remove(e);
            this.approximateSize -= e.value.approximateSize;
        }
        this.entryMap.delete(e.value.key);
    }

    private refresh(e: CacheNode<T>) {
        this.clearTimeout(e);

        e.value.timeoutId = setTimeout(() => this.expireNode(e), ServerConfig.cacheEntryTimeoutMs);
        this.entries.remove(e);
        this.entries.addFirst(e.value);
    }

    private expireNode(e: CacheNode<T>, notify = true) {
        if (notify) ConsoleLogger.log('Cache', `${e.value.key} expired.`);
        this.dispose(e);
    }

    expireAll() {
        for (let e = this.entries.first; e; e = e.next) this.expireNode(e, false);
    }

    expire(key: string) {
        const entry = this.entryMap.get(key);
        if (!entry) return;
        this.expireNode(entry);
    }

    add(item: T) {
        const key = this.keyGetter(item);
        const approximateSize = this.sizeGetter(item);

        if (this.entryMap.has(key)) this.dispose(this.entryMap.get(key)!);

        if (ServerConfig.cacheMaxSizeInBytes < this.approximateSize + approximateSize) {
            if (this.entries.last) this.dispose(this.entries.last);
        }

        this.approximateSize += approximateSize;
        const entry: CacheEntry<T> = { key, approximateSize, timeoutId: void 0, item };
        const e = this.entries.addFirst(entry);
        this.entryMap.set(key, e);
        this.refresh(e);
        ConsoleLogger.log('Cache', `${key} added.`);

        return item;
    }

    has(key: string) {
        return this.entryMap.has(key);
    }

    get(key: string) {
        if (!this.entryMap.has(key)) return void 0;
        let e = this.entryMap.get(key)!;
        this.refresh(e);
        ConsoleLogger.log('Cache', `${key} accessed.`);
        return e.value.item;
    }

    constructor(private keyGetter: (i: T) => string, private sizeGetter: (i: T) => number) {
    }
}