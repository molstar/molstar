/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/dsehnal/LiteMol
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Scheduler from './scheduler'

class Computation<A> {
    run(ctx?: Computation.Context) {
        return this.runObservable(ctx).result;
    }

    runObservable(ctx?: Computation.Context): Computation.Running<A> {
        const context = ctx ? ctx as ObservableContext : new ObservableContext();

        return {
            subscribe: (context as ObservableContext).subscribe || NoOpSubscribe,
            result: new Promise<A>(async (resolve, reject) => {
                try {
                    if (context.started) context.started();
                    const result = await this.computation(context);
                    resolve(result);
                } catch (e) {
                    if (Computation.PRINT_CONSOLE_ERROR) console.error(e);
                    reject(e);
                } finally {
                    if (context.finished) context.finished();
                }
            })
        };
    }

    constructor(private computation: (ctx: Computation.Context) => Promise<A>) {

    }
}


namespace Computation {
    export let PRINT_CONSOLE_ERROR = false;

    export function create<A>(computation: (ctx: Context) => Promise<A>) {
        return new Computation(computation);
    }

    export function resolve<A>(a: A) {
        return new Computation<A>(_ => Promise.resolve(a));
    }

    export function reject<A>(reason: any) {
        return new Computation<A>(_ => Promise.reject(reason));
    }

    export interface Params {
        isSynchronous: boolean,
        updateRateMs: number
    }

    export const Aborted = 'Aborted';

    export interface Progress {
        message: string,
        isIndeterminate: boolean,
        current: number,
        max: number,
        elapsedMs: number,
        requestAbort?: () => void
    }

    export interface Context {
        readonly requiresUpdate: boolean,
        requestAbort(): void,
        /**
         * Checks if the computation was aborted. If so, throws.
         * Otherwise, updates the progress.
         *
         * Returns the number of ms since the last update.
         */
        updateProgress(msg: string, abort?: boolean | (() => void), current?: number, max?: number): Promise<void> | void
    }

    export type ProgressObserver = (progress: Readonly<Progress>) => void;

    export interface Running<A> {
        subscribe(onProgress: ProgressObserver): void;
        result: Promise<A>
    }

    export const contextWithoutUpdates: Context = {
        requiresUpdate: false,
        requestAbort() { },
        updateProgress(msg, abort, current, max) { }
    }

    export function observableContext(params?: Partial<Params>) {
        return new ObservableContext(params);
    }



    declare var process: any;
    declare var window: any;

    export const now: () => number = (function () {
        if (typeof window !== 'undefined' && window.performance) {
            const perf = window.performance;
            return function () { return perf.now(); }
        } else if (typeof process !== 'undefined' && process.hrtime !== 'undefined') {
            return function () {
                let t = process.hrtime();
                return t[0] * 1000 + t[1] / 1000000;
            };
        } else {
            return function () { return +new Date(); }
        }
    })();

    export interface Chunked {
        /**
         * Get automatically computed chunk size
         * Or set it a default value.
         */
        chunkSize: number,
        readonly requiresUpdate: boolean,
        updateProgress: Context['updateProgress'],
        context: Context
    }

    export function chunked(ctx: Context, defaultChunkSize: number): Chunked {
        return new ChunkedImpl(ctx, defaultChunkSize);
    }
}

const DefaulUpdateRateMs = 150;
const NoOpSubscribe = () => { }

class ObservableContext implements Computation.Context {
    readonly updateRate: number;
    private isSynchronous: boolean;
    private level = 0;
    private startedTime = 0;
    private abortRequested = false;
    private lastUpdated = 0;
    private observers: Computation.ProgressObserver[] | undefined = void 0;
    private progress: Computation.Progress = { message: 'Working...', current: 0, max: 0, elapsedMs: 0, isIndeterminate: true, requestAbort: void 0 };

    lastDelta = 0;

    private checkAborted() {
        if (this.abortRequested) throw Computation.Aborted;
    }

    private abortRequester = () => { this.abortRequested = true };

    subscribe = (obs: Computation.ProgressObserver) => {
        if (!this.observers) this.observers = [];
        this.observers.push(obs);
    }

    requestAbort() {
        try {
            if (this.abortRequester) {
                this.abortRequester.call(null);
            }
        } catch (e) { }
    }

    updateProgress(msg: string, abort?: boolean | (() => void), current?: number, max?: number): Promise<void> | void {
        this.checkAborted();

        const time = Computation.now();

        if (typeof abort === 'boolean') {
            this.progress.requestAbort = abort ? this.abortRequester : void 0;
        } else {
            if (abort) this.abortRequester = abort;
            this.progress.requestAbort = abort ? this.abortRequester : void 0;
        }

        this.progress.message = msg;
        this.progress.elapsedMs = time - this.startedTime;
        if (isNaN(current!) || isNaN(max!)) {
            this.progress.isIndeterminate = true;
        } else {
            this.progress.isIndeterminate = false;
            this.progress.current = current!;
            this.progress.max = max!;
        }

        if (this.observers) {
            const p = { ...this.progress };
            for (const o of this.observers) Scheduler.immediate(o, p);
        }

        this.lastDelta = time - this.lastUpdated;
        this.lastUpdated = time;

        return Scheduler.immediatePromise();
    }

    get requiresUpdate() {
        this.checkAborted();
        if (this.isSynchronous) return false;
        return Computation.now() - this.lastUpdated > this.updateRate / 2;
    }

    started() {
        if (!this.level) this.startedTime = Computation.now();
        this.level++;
    }

    finished() {
        this.level--;
        if (this.level < 0) {
            throw new Error('Bug in code somewhere, Computation.resolve/reject called too many times.');
        }
        if (!this.level) this.observers = void 0;
    }

    constructor(params?: Partial<Computation.Params>) {
        this.updateRate = (params && params.updateRateMs) || DefaulUpdateRateMs;
        this.isSynchronous = !!(params && params.isSynchronous);
    }
}


class ChunkedImpl implements Computation.Chunked {
    private currentChunkSize: number;

    private computeChunkSize() {
        const lastDelta = (this.context as ObservableContext).lastDelta || 0;
        if (!lastDelta) return this.defaultChunkSize;
        const rate = (this.context as ObservableContext).updateRate || 0;
        return Math.round(this.currentChunkSize * rate / lastDelta + 1);
    }

    get chunkSize() {
        return this.defaultChunkSize;
    }

    set chunkSize(value: number) {
        this.defaultChunkSize = value;
        this.currentChunkSize = value;
    }

    get requiresUpdate() {
        const ret = this.context.requiresUpdate;
        if (!ret) this.currentChunkSize += this.chunkSize;
        return ret;
    }

    async updateProgress(msg: string, abort?: boolean | (() => void), current?: number, max?: number) {
        await this.context.updateProgress(msg, abort, current, max);
        this.defaultChunkSize = this.computeChunkSize();
    }

    constructor(public context: Computation.Context, private defaultChunkSize: number) {
        this.currentChunkSize = defaultChunkSize;
    }
}

export default Computation;