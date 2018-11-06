/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { LinkedList } from 'mol-data/generic';

export { PluginCommand }

/** namespace.id must a globally unique identifier */
function PluginCommand<T>(namespace: string, id: string, params?: PluginCommand.Descriptor<T>['params']): PluginCommand.Descriptor<T> {
    return new Impl(`${namespace}.${id}` as PluginCommand.Id, params);
}

const cmdRepo = new Map<string, PluginCommand.Descriptor<any>>();
class Impl<T> implements PluginCommand.Descriptor<T> {
    dispatch(ctx: PluginContext, params: T): Promise<void> {
        return ctx.commands.dispatch(this, params)
    }
    constructor(public id: PluginCommand.Id, public params: PluginCommand.Descriptor<T>['params']) {
        if (cmdRepo.has(id)) throw new Error(`Command id '${id}' already in use.`);
        cmdRepo.set(id, this);
    }
}

namespace PluginCommand {
    export type Id = string & { '@type': 'plugin-command-id' }

    export interface Descriptor<T = unknown> {
        readonly id: PluginCommand.Id,
        dispatch(ctx: PluginContext, params: T): Promise<void>,
        params?: { toJSON(params: T): any, fromJSON(json: any): T }
    }

    type Action<T> = (params: T) => void | Promise<void>
    type Instance = { id: string, params: any, resolve: () => void, reject: (e: any) => void }

    export class Manager {
        private subs = new Map<string, Action<any>[]>();
        private queue = LinkedList<Instance>();
        private disposing = false;

        subscribe<T>(cmd: Descriptor<T>, action: Action<T>) {
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
        dispatch<T>(cmd: Descriptor<T> | Id, params: T) {
            return new Promise<void>((resolve, reject) => {
                if (!this.disposing) {
                    reject('disposed');
                    return;
                }

                const id = typeof cmd === 'string' ? cmd : (cmd as Descriptor<T>).id;
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

        private async next() {
            if (this.queue.count === 0) return;
            const cmd = this.queue.removeFirst()!;

            const actions = this.subs.get(cmd.id);
            if (!actions) return;

            try {
                // TODO: should actions be called "asynchronously" ("setImmediate") instead?
                for (const a of actions) {
                    await a(cmd.params);
                }
                cmd.resolve();
            } catch (e) {
                cmd.reject(e);
            } finally {
                if (!this.disposing) this.next();
            }
        }
    }
}


// TODO: command interface and queue.
// How to handle command resolving? Track how many subscriptions a command has?