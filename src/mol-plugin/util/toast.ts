/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol (c) David Sehnal
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StatefulPluginComponent } from '../../mol-plugin-state/component';
import { OrderedMap } from 'immutable';
import { PluginContext } from '../context';
import { PluginCommands } from '../commands';

export interface PluginToast {
    title: string,
    /**
     * The message can be either a string, html string, or an arbitrary React component.
     */
    message: string | React.ComponentClass,
    /**
     * Only one message with a given key can be shown.
     */
    key?: string,
    /**
     * Specify a timeout for the message in milliseconds.
     */
    timeoutMs?: number
}

export class PluginToastManager extends StatefulPluginComponent<{
    entries: OrderedMap<number, PluginToastManager.Entry>
}> {
    readonly events = {
        changed: this.ev()
    };

    private serialNumber = 0;
    private serialId = 0;

    private findByKey(key: string): PluginToastManager.Entry | undefined {
        return this.state.entries.find(e => !!e && e.key === key);
    }

    private show(toast: PluginToast) {
        let entries = this.state.entries;
        let e: PluginToastManager.Entry | undefined = void 0;
        const id = ++this.serialId;
        let serialNumber: number;
        if (toast.key && (e = this.findByKey(toast.key))) {
            if (e.timeout !== void 0) clearTimeout(e.timeout);
            serialNumber = e.serialNumber;
            entries = entries.remove(e.id);
        } else {
            serialNumber = ++this.serialNumber;
        }

        e = {
            id,
            serialNumber,
            key: toast.key,
            title: toast.title,
            message: toast.message,
            timeout: this.timeout(id, toast.timeoutMs),
            hide: () => this.hideId(id)
        };

        if (this.updateState({ entries: entries.set(id, e) })) this.events.changed.next();
    }

    private timeout(id: number, delay?: number) {
        if (delay === void 0) return void 0;

        if (delay < 0) delay = 500;
        return <number><any>setTimeout(() => {
            const e = this.state.entries.get(id);
            e.timeout = void 0;
            this.hide(e);
        }, delay);
    }

    private hideId(id: number) {
        this.hide(this.state.entries.get(id));
    }

    private hide(e: PluginToastManager.Entry | undefined) {
        if (!e) return;
        if (e.timeout !== void 0) clearTimeout(e.timeout);
        e.hide = <any>void 0;
        if (this.updateState({ entries: this.state.entries.delete(e.id) })) this.events.changed.next();
    }

    constructor(plugin: PluginContext) {
        super({ entries: OrderedMap<number, PluginToastManager.Entry>() });

        PluginCommands.Toast.Show.subscribe(plugin, e => this.show(e));
        PluginCommands.Toast.Hide.subscribe(plugin, e => this.hide(this.findByKey(e.key)));
    }
}

export namespace PluginToastManager {
    export interface Entry {
        id: number,
        serialNumber: number,
        key?: string,
        title: string,
        message: string | React.ComponentClass,
        hide: () => void,
        timeout?: number
    }
}