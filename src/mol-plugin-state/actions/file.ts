/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { Task } from '../../mol-task';
import { getFileInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';

export const OpenFiles = StateAction.build({
    display: { name: 'Open Files', description: 'Load one or more files and optionally create default visuals' },
    from: PluginStateObject.Root,
    params: (a, ctx: PluginContext) => {
        const { extensions, options } = ctx.dataFormats;
        return {
            files: PD.FileList({ accept: Array.from(extensions.values()).map(e => `.${e}`).join(',') + ',.gz,.zip', multiple: true }),
            format: PD.MappedStatic('auto', {
                auto: PD.EmptyGroup(),
                specific: PD.Select(options[0][0], options)
            }),
            visuals: PD.Boolean(true, { description: 'Add default visuals' }),
        };
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Open Files', async taskCtx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    await state.transaction(async () => {
        if (params.files === null) {
            plugin.log.error('No file(s) selected');
            return;
        }
        for (const file of params.files) {
            try {
                const info = getFileInfo(file.file!);
                const isBinary = plugin.dataFormats.binaryExtensions.has(info.ext);
                const { data } = await plugin.builders.data.readFile({ file, isBinary });
                const provider = params.format.name === 'auto'
                    ? plugin.dataFormats.auto(info, data.cell?.obj!)
                    : plugin.dataFormats.get(params.format.params);

                if (!provider) {
                    plugin.log.warn(`OpenFiles: could not find data provider for '${info.name}.${info.ext}'`);
                    continue;
                }

                // need to await so that the enclosing Task finishes after the update is done.
                const parsed = await provider.parse(plugin, data);
                if (params.visuals) {
                    await provider.visuals?.(plugin, parsed);
                }
            } catch (e) {
                plugin.log.error(e);
            }
        }
    }).runInContext(taskCtx);
}));

export const DownloadFile = StateAction.build({
    display: { name: 'Download File', description: 'Load one or more file from an URL' },
    from: PluginStateObject.Root,
    params: (a, ctx: PluginContext) => {
        const { options } = ctx.dataFormats;
        return {
            url: PD.Url(''),
            format: PD.Select(options[0][0], options),
            isBinary: PD.Boolean(false),
            visuals: PD.Boolean(true, { description: 'Add default visuals' }),
        };
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Open Files', async taskCtx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    await state.transaction(async () => {
        try {
            const provider = plugin.dataFormats.get(params.format);
            if (!provider) {
                plugin.log.warn(`DownloadFile: could not find data provider for '${params.format}'`);
                return;
            }

            const data = await plugin.builders.data.download({ url: params.url, isBinary: params.isBinary });
            const parsed = await provider.parse(plugin, data);
            if (params.visuals) {
                await provider.visuals?.(plugin, parsed);
            }
        } catch (e) {
            plugin.log.error(e);
        }
    }).runInContext(taskCtx);
}));