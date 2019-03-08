/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { LinkedList } from 'mol-data/generic';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { UUID } from 'mol-util';

export { PluginCommand }

interface PluginCommand<T = unknown> {
    readonly id: UUID,
    dispatch(ctx: PluginContext, params: T, isChild?: boolean): Promise<void>,
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription,
    params: { isImmediate: boolean }
}

/** namespace.id must a globally unique identifier */
function PluginCommand<T>(params?: Partial<PluginCommand<T>['params']>): PluginCommand<T> {
    return new Impl({ isImmediate: false, ...params });
}

class Impl<T> implements PluginCommand<T> {
    dispatch(ctx: PluginContext, params: T, isChild?: boolean): Promise<void> {
        return ctx.commands.dispatch(this, params, isChild);
    }
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription {
        return ctx.commands.subscribe(this, action);
    }
    id = UUID.create22();
    constructor(public params: PluginCommand<T>['params']) {
    }
}

namespace PluginCommand {
    export type Id = string & { '@type': 'plugin-command-id' }

    export interface Subscription {
        unsubscribe(): void
    }

    export type Action<T> = (params: T) => unknown | Promise<unknown>
    type Instance = { cmd: PluginCommand<any>, params: any, isChild: boolean, resolve: () => void, reject: (e: any) => void }

    export class Manager {
        private subs = new Map<string, Action<any>[]>();
        private queue = LinkedList<Instance>();
        private disposing = false;

        private ev = RxEventHelper.create();

        readonly behaviour = {
            locked: this.ev.behavior<boolean>(false)
        };

        lock(locked: boolean = true) {
            this.behaviour.locked.next(locked);
        }

        subscribe<T>(cmd: PluginCommand<T>, action: Action<T>): Subscription {
            let actions = this.subs.get(cmd.id);
            if (!actions) {
                actions = [];
                this.subs.set(cmd.id, actions);
            }
            actions.push(action);

            return {
                unsubscribe: () => {
                    const actions = this.subs.get(cmd.id);
                    if (!actions) return;
                    const idx = actions.indexOf(action);
                    if (idx < 0) return;
                    for (let i = idx + 1; i < actions.length; i++) {
                        actions[i - 1] = actions[i];
                    }
                    actions.pop();
                }
            }
        }


        /** Resolves after all actions have completed */
        dispatch<T>(cmd: PluginCommand<T>, params: T, isChild = false) {
            return new Promise<void>((resolve, reject) => {
                if (this.disposing) {
                    reject('disposed');
                    return;
                }

                const actions = this.subs.get(cmd.id);
                if (!actions) {
                    resolve();
                    return;
                }

                const instance: Instance = { cmd, params, resolve, reject, isChild };

                if (cmd.params.isImmediate || isChild) {
                    this.resolve(instance);
                } else {
                    this.queue.addLast(instance);
                    this.next();
                }
            });
        }

        dispose() {
            this.subs.clear();
            while (this.queue.count > 0) {
                this.queue.removeFirst();
            }
        }

        private async resolve(instance: Instance) {
            const actions = this.subs.get(instance.cmd.id);
            if (!actions) {
                try {
                    instance.resolve();
                } finally {
                    if (!instance.cmd.params.isImmediate && !this.disposing) this.next();
                }
                return;
            }

            try {
                if (!instance.cmd.params.isImmediate && !instance.isChild) this.executing = true;
                // TODO: should actions be called "asynchronously" ("setImmediate") instead?
                for (const a of actions) {
                    await a(instance.params);
                }
                instance.resolve();
            } catch (e) {
                instance.reject(e);
            } finally {
                if (!instance.cmd.params.isImmediate && !instance.isChild) {
                    this.executing = false;
                    if (!this.disposing) this.next();
                }
            }
        }

        private executing = false;
        private async next() {
            if (this.queue.count === 0 || this.executing) return;
            const instance = this.queue.removeFirst()!;
            this.resolve(instance);
        }
    }
}