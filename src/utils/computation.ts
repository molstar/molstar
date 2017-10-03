/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/dsehnal/LiteMol
 * @author David Sehnal <david.sehnal@gmail.com>
 */

class Computation<A> {
    run(ctx?: Computation.Context) {
        return this.runWithContext(ctx).result;
    }

    runWithContext(ctx?: Computation.Context): Computation.Running<A> {
        const context = ctx ? ctx as Computation.ObservableContext : new Computation.ObservableContext();

        return {
            subscribe: (context as Computation.ObservableContext).subscribe || Helpers.NoOpSubscribe,
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
    const DefaulUpdateRateMs = 100;

    export let PRINT_CONSOLE_ERROR = false;

    export function resolve<A>(a: A) {
        return new Computation<A>(() => Promise.resolve(a));
    }

    export function reject<A>(reason: any) {
        return new Computation<A>(() => Promise.reject(reason));
    }

    export interface Params {
        isSynchronous: boolean,
        updateRateMs: number
    }

    export const Aborted = 'Aborted';

    export interface Progress {
        message: string;
        isIndeterminate: boolean;
        current: number;
        max: number;
        requestAbort?: () => void;
    }

    export interface Context {
        readonly requiresUpdate: boolean,
        requestAbort(): void,
        /**
         * Checks if the computation was aborted. If so, throws.
         * Otherwise, updates the progress.
         */
        updateProgress(msg: string, abort?: boolean | (() => void), current?: number, max?: number): Promise<void> | void
    }

    type ProgressObserver = (progress: Readonly<Progress>) => void;

    export interface Running<A> {
        subscribe(onProgress: ProgressObserver): void;
        result: Promise<A>
    }

    export const ContextWithoutUpdates: Context = {
        requiresUpdate: false,
        requestAbort() { },
        updateProgress(msg, abort, current, max) { }
    }

    export class ObservableContext implements Context {
        private updateRate: number;
        private isSynchronous: boolean;
        private level = 0;
        private abortRequested = false;
        private lastUpdated = 0;
        private observers: ProgressObserver[] | undefined = void 0;
        private progress: Progress = { message: 'Working...', current: 0, max: 0, isIndeterminate: true, requestAbort: void 0 };

        private checkAborted() {
            if (this.abortRequested) throw Aborted;
        }

        private abortRequester = () => { this.abortRequested = true };

        subscribe = (obs: ProgressObserver) => {
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

            if (typeof abort === 'boolean') {
                this.progress.requestAbort = abort ? this.abortRequester : void 0;
            } else {
                if (abort) this.abortRequester = abort;
                this.progress.requestAbort = abort ? this.abortRequester : void 0;
            }

            this.progress.message = msg;
            if (isNaN(current!)) {
                this.progress.isIndeterminate = true;
            } else {
                this.progress.isIndeterminate = false;
                this.progress.current = current!;
                this.progress.max = max!;
            }

            if (this.observers) {
                const p = { ...this.progress };
                for (const o of this.observers) setTimeout(o, 0, p);
            }

            this.lastUpdated = Helpers.getTime();

            return new Promise<void>(res => setTimeout(res, 0));
        }

        get requiresUpdate() {
            this.checkAborted();
            if (this.isSynchronous) return false;
            return Helpers.getTime() - this.lastUpdated > this.updateRate;
        }

        started() {
            this.level++;
        }

        finished() {
            this.level--;
            if (this.level < 0) {
                throw new Error('Bug in code somewhere, Computation.resolve/reject called too many times.');
            }
            if (!this.level) this.observers = void 0;
        }

        constructor(params?: Params) {
            this.updateRate = (params && params.updateRateMs) || DefaulUpdateRateMs;
            this.isSynchronous = !!(params && params.isSynchronous);
        }
    }
}

namespace Helpers {
    declare var process: any;
    declare var window: any;

    export const getTime: () => number = (function () {
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

    export const NoOpSubscribe = () => {};
}

export default Computation;