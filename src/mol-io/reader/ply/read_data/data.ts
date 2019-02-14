/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../../../../mol-plugin/state/objects';
import { PluginStateObject as SO } from '../../../../mol-plugin/state/objects';
import { Task } from 'mol-task';
import PLY from 'mol-io/reader/ply/parse_data/ply_parser'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Transformer } from 'mol-state';
import { readFromFile } from './data-source';

export { ReadFile_ascii }
type ReadFile_ascii = typeof ReadFile_ascii
const ReadFile_ascii = PluginStateTransform.BuiltIn({
    name: 'ReadFile_ascii',
    display: { name: 'ReadFile_ascii', description: 'Read string data from the specified file' },
    from: SO.Root,
    to: [SO.Data.String],
    params: {
        file: PD.File(),
        label: PD.makeOptional(PD.Text('')),
        isBinary: PD.makeOptional(PD.Boolean(false, { description: 'If true, open file as as binary (string otherwise)' }))
    }
})({
    apply({ params: p }) {
        return Task.create('Open File', async ctx => {
            const data = await readFromFile(p.file).runInContext(ctx);
            return  new SO.Data.String(data as string, { label: p.label ? p.label : p.file.name });
        });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.label !== newParams.label) {
            (b.label as string) = newParams.label || oldParams.file.name;
            return Transformer.UpdateResult.Updated;
        }
        return Transformer.UpdateResult.Unchanged;
    },
    isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user loaded files.' })
});


export { ParsePLY }
type ParsePLY = typeof ParsePLY
const ParsePLY = PluginStateTransform.BuiltIn({
    name: 'parse-ply',
    display: { name: 'Parse PLY', description: 'Parse PLY from String' },
    from: [SO.Data.String],
    to: SO.Format.Ply
})({
    apply({ a }) {
        return Task.create('Parse PLY', async ctx => {
            const parsed = await (PLY(a.data).runInContext(ctx));
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Ply(parsed.result);
        });
    }
});