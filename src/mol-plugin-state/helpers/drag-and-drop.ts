/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { OpenFiles } from '../actions/file';
import { Asset } from '../../mol-util/assets';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';

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

export async function handleDragAndDrop(plugin: PluginContext, files: File[]) {
    const handlers = Array.from(plugin.customDragAndDropHandlers.values());

    for (const handler of handlers) {
        const handled = await handler(files, plugin);
        if (handled) return;
    }

    defaultDragAndDropHandler(plugin, files);
}