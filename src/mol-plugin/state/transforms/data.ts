/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { Task } from 'mol-task';
import CIF from 'mol-io/reader/cif'
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StateTransformer } from 'mol-state';
import { readFromFile, ajaxGetMany } from 'mol-util/data-source';
import * as CCP4 from 'mol-io/reader/ccp4/parser'
import * as DSN6 from 'mol-io/reader/dsn6/parser'
import * as PLY from 'mol-io/reader/ply/parser'

export { Download }
type Download = typeof Download
const Download = PluginStateTransform.BuiltIn({
    name: 'download',
    display: { name: 'Download', description: 'Download string or binary data from the specified URL' },
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    params: {
        url: PD.Text('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
        label: PD.Optional(PD.Text('')),
        isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' }))
    }
})({
    apply({ params: p }, globalCtx: PluginContext) {
        return Task.create('Download', async ctx => {
            const data = await globalCtx.fetch({ url: p.url, type: p.isBinary ? 'binary' : 'string' }).runInContext(ctx);
            return p.isBinary
                ? new SO.Data.Binary(data as Uint8Array, { label: p.label ? p.label : p.url })
                : new SO.Data.String(data as string, { label: p.label ? p.label : p.url });
        });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            (b.label as string) = newParams.label || newParams.url;
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    }
});

export { DownloadBlob }
type DownloadBlob = typeof DownloadBlob
const DownloadBlob = PluginStateTransform.BuiltIn({
    name: 'download-blob',
    display: { name: 'Download Blob', description: 'Download multiple string or binary data from the specified URLs.' },
    from: SO.Root,
    to: SO.Data.Blob,
    params: {
        sources: PD.ObjectList({
            id: PD.Text('', { label: 'Unique ID' }),
            url: PD.Text('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
            isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' })),
            canFail: PD.Optional(PD.Boolean(false, { description: 'Indicate whether the download can fail and not be included in the blob as a result.' }))
        }, e => `${e.id}: ${e.url}`),
        maxConcurrency: PD.Optional(PD.Numeric(4, { min: 1, max: 12, step: 1 }, { description: 'The maximum number of concurrent downloads.' }))
    }
})({
    apply({ params }, plugin: PluginContext) {
        return Task.create('Download Blob', async ctx => {
            const entries: SO.Data.BlobEntry[] = [];
            const data = await ajaxGetMany(ctx, params.sources, params.maxConcurrency || 4);

            for (let i = 0; i < data.length; i++) {
                const r = data[i], src = params.sources[i];
                if (r.kind === 'error') plugin.log.warn(`Download ${r.id} (${src.url}) failed: ${r.error}`);
                else {
                    entries.push(src.isBinary
                        ? { id: r.id, kind: 'binary', data: r.result as Uint8Array }
                        : { id: r.id, kind: 'string', data: r.result as string });
                }
            }
            return new SO.Data.Blob(entries, { label: 'Data Blob', description: `${entries.length} ${entries.length === 1 ? 'entry' : 'entries'}` });
        });
    },
    // TODO: ??
    // update({ oldParams, newParams, b }) {
    //     return 0 as any;
    //     // if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary) return StateTransformer.UpdateResult.Recreate;
    //     // if (oldParams.label !== newParams.label) {
    //     //     (b.label as string) = newParams.label || newParams.url;
    //     //     return StateTransformer.UpdateResult.Updated;
    //     // }
    //     // return StateTransformer.UpdateResult.Unchanged;
    // }
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
        label: PD.Optional(PD.Text('')),
        isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, open file as as binary (string otherwise)' }))
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
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    },
    isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user loaded files.' })
});

export { ParseBlob }
type ParseBlob = typeof ParseBlob
const ParseBlob = PluginStateTransform.BuiltIn({
    name: 'parse-blob',
    display: { name: 'Parse Blob', description: 'Parse multiple data enties' },
    from: SO.Data.Blob,
    to: SO.Format.Blob,
    params: {
        formats: PD.ObjectList({
            id: PD.Text('', { label: 'Unique ID' }),
            format: PD.Select<'cif'>('cif', [['cif', 'cif']])
        }, e => `${e.id}: ${e.format}`)
    }
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Parse Blob', async ctx => {
            const map = new Map<string, string>();
            for (const f of params.formats) map.set(f.id, f.format);

            const entries: SO.Format.BlobEntry[] = [];

            for (const e of a.data) {
                if (!map.has(e.id)) continue;

                const parsed = await (e.kind === 'string' ? CIF.parse(e.data) : CIF.parseBinary(e.data)).runInContext(ctx);
                if (parsed.isError) throw new Error(`${e.id}: ${parsed.message}`);
                entries.push({ id: e.id, kind: 'cif', data: parsed.result });
            }

            return new SO.Format.Blob(entries, { label: 'Format Blob', description: `${entries.length} ${entries.length === 1 ? 'entry' : 'entries'}` });
        });
    },
    // TODO: ??
    // update({ oldParams, newParams, b }) {
    //     return 0 as any;
    //     // if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary) return StateTransformer.UpdateResult.Recreate;
    //     // if (oldParams.label !== newParams.label) {
    //     //     (b.label as string) = newParams.label || newParams.url;
    //     //     return StateTransformer.UpdateResult.Updated;
    //     // }
    //     // return StateTransformer.UpdateResult.Unchanged;
    // }
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

export { ParseCcp4 }
type ParseCcp4 = typeof ParseCcp4
const ParseCcp4 = PluginStateTransform.BuiltIn({
    name: 'parse-ccp4',
    display: { name: 'Parse CCP4/MRC/MAP', description: 'Parse CCP4/MRC/MAP from Binary data' },
    from: [SO.Data.Binary],
    to: SO.Format.Ccp4
})({
    apply({ a }) {
        return Task.create('Parse CCP4/MRC/MAP', async ctx => {
            const parsed = await CCP4.parse(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Ccp4(parsed.result);
        });
    }
});

export { ParseDsn6 }
type ParseDsn6 = typeof ParseDsn6
const ParseDsn6 = PluginStateTransform.BuiltIn({
    name: 'parse-dsn6',
    display: { name: 'Parse DSN6/BRIX', description: 'Parse CCP4/BRIX from Binary data' },
    from: [SO.Data.Binary],
    to: SO.Format.Dsn6
})({
    apply({ a }) {
        return Task.create('Parse DSN6/BRIX', async ctx => {
            const parsed = await DSN6.parse(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Dsn6(parsed.result);
        });
    }
});