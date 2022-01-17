/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8ByteCount, utf8Write } from '../../mol-io/common/utf8';
import { to_mmCIF } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { getFormattedTime } from '../../mol-util/date';
import { download } from '../../mol-util/download';
import { zip } from '../../mol-util/zip/zip';

export async function exportHierarchy(plugin: PluginContext, options?: { format?: 'cif' | 'bcif' }) {
    try {
        await _exportHierarchy(plugin, options);
    } catch (e) {
        console.error(e);
        plugin.log.error(`Model export failed. See console for details.`);
    }
}

async function _exportHierarchy(plugin: PluginContext, options?: { format?: 'cif' | 'bcif' }) {
    const format = options?.format ?? 'cif';
    const { structures } = plugin.managers.structure.hierarchy.current;

    const files: [name: string, data: string | Uint8Array][] = [];
    const entryMap = new Map<string, number>();

    for (const _s of structures) {
        const s = _s.transform?.cell.obj?.data ?? _s.cell.obj?.data;
        if (!s) continue;
        if (s.models.length > 1) {
            plugin.log.warn(`[Export] Skipping ${_s.cell.obj?.label}: Multimodel exports not supported.`);
            continue;
        }

        const name = entryMap.has(s.model.entryId)
            ? `${s.model.entryId}_${entryMap.get(s.model.entryId)! + 1}.${format}`
            : `${s.model.entryId}.${format}`;
        entryMap.set(s.model.entryId, (entryMap.get(s.model.entryId) ?? 0) + 1);
        files.push([name, to_mmCIF(s.model.entryId, s, format === 'bcif', { copyAllCategories: true })]);
    }

    if (files.length === 1) {
        download(new Blob([files[0][1]]), files[0][0]);
    } else if (files.length > 1) {
        const zipData: Record<string, Uint8Array> = {};
        for (const [fn, data] of files) {
            if (data instanceof Uint8Array) {
                zipData[fn] = data;
            } else {
                const bytes = new Uint8Array(utf8ByteCount(data));
                utf8Write(bytes, 0, data);
                zipData[fn] = bytes;
            }
        }
        const task = Task.create('Export Models', async ctx => {
            return zip(ctx, zipData);
        });
        const buffer = await plugin.runTask(task);
        download(new Blob([new Uint8Array(buffer, 0, buffer.byteLength)]), `structures_${getFormattedTime()}.zip`);
    }
}