/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../context';
import { UUID } from '../../mol-util';

export { PluginCommand }

interface PluginCommand<T = unknown> {
    readonly id: UUID,
    dispatch(ctx: PluginContext, params: T): Promise<void>,
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription
}

/** namespace.id must a globally unique identifier */
function PluginCommand<T>(): PluginCommand<T> {
    return new Impl();
}

class Impl<T> implements PluginCommand<T> {
    dispatch(ctx: PluginContext, params: T): Promise<void> {
        return ctx.commands.dispatch(this, params);
    }
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription {
        return ctx.commands.subscribe(this, action);
    }
    id = UUID.create22();
    constructor() {
    }
}

namespace PluginCommand {
    export type Id = string & { '@type': 'plugin-command-id' }

    export interface Subscription {
        unsubscribe(): void
    }

    export type Action<T> = (params: T) => unknown | Promise<unknown>
    type Instance = { cmd: PluginCommand<any>, params: any, resolve: () => void, reject: (e: any) => void }

    export class Manager {
        private subs = new Map<string, Action<any>[]>();
        private disposing = false;

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
        dispatch<T>(cmd: PluginCommand<T>, params: T) {
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

                this.resolve({ cmd, params, resolve, reject });
            });
        }

        dispose() {
            this.subs.clear();
        }

        private async resolve(instance: Instance) {
            const actions = this.subs.get(instance.cmd.id);
            if (!actions) {
                instance.resolve();
                return;
            }

            try {
                // TODO: should actions be called "asynchronously" ("setImmediate") instead?
                for (const a of actions) {
                    await a(instance.params);
                }
                instance.resolve();
            } catch (e) {
                instance.reject(e);
            }
        }
    }
}