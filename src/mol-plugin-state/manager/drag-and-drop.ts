/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { OpenFiles } from '../actions/file';
import { Asset } from '../../mol-util/assets';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';

export type PluginDragAndDropHandler = (files: File[], plugin: PluginContext) => Promise<boolean> | boolean

export class DragAndDropManager {
    private handlers: [name: string, handler: PluginDragAndDropHandler][] = [];

    addHandler(name: string, handler: PluginDragAndDropHandler) {
        const index = this.handlers.findIndex(h => h[0] === name);
        if (index < 0) this.handlers.push([name, handler]);
        else this.handlers[index][1] = handler;
    }

    removeHandler(name: string) {
        const index = this.handlers.findIndex(h => h[0] === name);
        if (index >= 0) this.handlers.splice(index, 1);
    }

    async handle(files: File[]) {
        for (let i = this.handlers.length - 1; i >= 0; i--) {
            const handler = this.handlers[i][1];
            const handled = await handler(files, this.plugin);
            if (handled) return;
        }

        defaultDragAndDropHandler(this.plugin, files);
    }

    dispose() {
        this.handlers.length = 0;
    }

    constructor(public plugin: PluginContext) {
    }
}

function defaultDragAndDropHandler(plugin: PluginContext, files: File[]) {
    const sessions = files.filter(f => {
        const fn = f.name.toLowerCase();
        return fn.endsWith('.molx') || fn.endsWith('.molj');
    });

    if (sessions.length > 0) {
        PluginCommands.State.Snapshots.OpenFile(plugin, { file: sessions[0] });
    } else {
        plugin.runTask(plugin.state.data.applyAction(OpenFiles, {
            files: files.map(f => Asset.File(f)),
            format: { name: 'auto', params: {} },
            visuals: true
        }));
    }
}