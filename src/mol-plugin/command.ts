/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { UUID } from '../mol-util';

export { PluginCommand, PluginCommandManager };

interface PluginCommand<T = unknown> {
    (ctx: PluginContext, params?: T): Promise<void>,
    readonly id: UUID,
    subscribe(ctx: PluginContext, action: PluginCommand.Action<T>): PluginCommand.Subscription
}

function PluginCommand<T>(): PluginCommand<T> {
    const ret: PluginCommand<T> = ((ctx, params) => ctx.commands.dispatch(ret, params || {} as any)) as PluginCommand<T>;
    ret.subscribe = (ctx, action) => ctx.commands.subscribe(ret, action);
    (ret.id as UUID) = UUID.create22();

    return ret;
}

namespace PluginCommand {
    export type Id = string & { '@type': 'plugin-command-id' }

    export interface Subscription {
        unsubscribe(): void
    }

    export type Action<T> = (params: T) => unknown | Promise<unknown>
}

type Instance = { cmd: PluginCommand<any>, params: any, resolve: () => void, reject: (e: any) => void }
class PluginCommandManager {
    private subs = new Map<string, PluginCommand.Action<any>[]>();
    private disposing = false;

    subscribe<T>(cmd: PluginCommand<T>, action: PluginCommand.Action<T>): PluginCommand.Subscription {
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
        };
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