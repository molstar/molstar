/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { Task } from 'mol-task';
import CIF from 'mol-io/reader/cif'
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Transformer } from 'mol-state';
import { readFromFile } from 'mol-util/data-source';
import * as PLY from 'mol-io/reader/ply/parse_data/ply_parser'

export { Download }
type Download = typeof Download
const Download = PluginStateTransform.BuiltIn({
    name: 'download',
    display: { name: 'Download', description: 'Download string or binary data from the specified URL' },
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    params: {
        url: PD.Text('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
        label: PD.makeOptional(PD.Text('')),
        isBinary: PD.makeOptional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' }))
    }
})({
    apply({ params: p }, globalCtx: PluginContext) {
        return Task.create('Download', async ctx => {
            const data = await globalCtx.fetch(p.url, p.isBinary ? 'binary' : 'string').runInContext(ctx);
            return p.isBinary
                ? new SO.Data.Binary(data as Uint8Array, { label: p.label ? p.label : p.url })
                : new SO.Data.String(data as string, { label: p.label ? p.label : p.url });
        });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary) return Transformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            (b.label as string) = newParams.label || newParams.url;
            return Transformer.UpdateResult.Updated;
        }
        return Transformer.UpdateResult.Unchanged;
    }
});

export { ReadFile }
type ReadFile = typeof ReadFile
const ReadFile = PluginStateTransform.BuiltIn({
    name: 'read-file',
    display: { name: 'Read File', description: 'Read string or binary data from the specified file' },
    from: SO.Root,
    to: [SO.Data.String, SO.Data.Binary],
    params: {
        file: PD.File(),
        label: PD.makeOptional(PD.Text('')),
        isBinary: PD.makeOptional(PD.Boolean(false, { description: 'If true, open file as as binary (string otherwise)' }))
    }
})({
    apply({ params: p }) {
        return Task.create('Open File', async ctx => {
            const data = await readFromFile(p.file, p.isBinary ? 'binary' : 'string').runInContext(ctx);
            return p.isBinary
                ? new SO.Data.Binary(data as Uint8Array, { label: p.label ? p.label : p.file.name })
                : new SO.Data.String(data as string, { label: p.label ? p.label : p.file.name });
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

export { ParseCif }
type ParseCif = typeof ParseCif
const ParseCif = PluginStateTransform.BuiltIn({
    name: 'parse-cif',
    display: { name: 'Parse CIF', description: 'Parse CIF from String or Binary data' },
    from: [SO.Data.String, SO.Data.Binary],
    to: SO.Format.Cif
})({
    apply({ a }) {
        return Task.create('Parse CIF', async ctx => {
            const parsed = await (SO.Data.String.is(a) ? CIF.parse(a.data) : CIF.parseBinary(a.data)).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Cif(parsed.result);
        });
    }
});

export { ParsePly }
type ParsePly = typeof ParsePly
const ParsePly = PluginStateTransform.BuiltIn({
    name: 'parse-ply',
    display: { name: 'Parse PLY', description: 'Parse PLY from Binary data' },
    from: [SO.Data.String],
    to: SO.Format.Ply
})({
    apply({ a }) {
        return Task.create('Parse PLY', async ctx => {
            const parsed = await PLY.parse(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Ply(parsed.result, { label: parsed.result.name || 'PLY Data' });
        });
    }
});