/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { LinkedList } from 'mol-data/generic';
import { RxEventHelper } from 'mol-util/rx-event-helper';

export { PluginCommand }

interface PluginCommand<T = unknown> {
    readonly id: PluginCommand.Id,
    dispatch(ctx: PluginContext, params: T): Promise<void>,
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription,
    params?: { toJSON(params: T): any, fromJSON(json: any): T }
}

/** namespace.id must a globally unique identifier */
function PluginCommand<T>(namespace: string, id: string, params?: PluginCommand<T>['params']): PluginCommand<T> {
    return new Impl(`${namespace}.${id}` as PluginCommand.Id, params);
}

const cmdRepo = new Map<string, PluginCommand<any>>();
class Impl<T> implements PluginCommand<T> {
    dispatch(ctx: PluginContext, params: T): Promise<void> {
        return ctx.commands.dispatch(this, params)
    }
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription {
        return ctx.commands.subscribe(this, action);
    }
    constructor(public id: PluginCommand.Id, public params: PluginCommand<T>['params']) {
        if (cmdRepo.has(id)) throw new Error(`Command id '${id}' already in use.`);
        cmdRepo.set(id, this);
    }
}

namespace PluginCommand {
    export type Id = string & { '@type': 'plugin-command-id' }

    export interface Subscription {
        unsubscribe(): void
    }

    export type Action<T> = (params: T) => void | Promise<void>
    type Instance = { id: string, params: any, resolve: () => void, reject: (e: any) => void }

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
        dispatch<T>(cmd: PluginCommand<T> | Id, params: T) {
            return new Promise<void>((resolve, reject) => {
                if (this.disposing) {
                    reject('disposed');
                    return;
                }

                const id = typeof cmd === 'string' ? cmd : (cmd as PluginCommand<T>).id;
                const actions = this.subs.get(id);
                if (!actions) {
                    resolve();
                    return;
                }

                this.queue.addLast({ id, params, resolve, reject });
                this.next();
            });
        }

        dispose() {
            this.subs.clear();
            while (this.queue.count > 0) {
                this.queue.removeFirst();
            }
        }

        private executing = false;
        private async next() {
            if (this.queue.count === 0 || this.executing) return;
            const cmd = this.queue.removeFirst()!;

            const actions = this.subs.get(cmd.id);
            if (!actions) {
                return;
            }

            try {
                this.executing = true;
                // TODO: should actions be called "asynchronously" ("setImmediate") instead?
                for (const a of actions) {
                    await a(cmd.params);
                }
                cmd.resolve();
            } catch (e) {
                cmd.reject(e);
            } finally {
                this.executing = false;
                if (!this.disposing) this.next();
            }
        }
    }
}