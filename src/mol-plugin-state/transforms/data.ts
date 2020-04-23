/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isTypedArray } from '../../mol-data/db/column-helpers';
import * as CCP4 from '../../mol-io/reader/ccp4/parser';
import { CIF } from '../../mol-io/reader/cif';
import * as DSN6 from '../../mol-io/reader/dsn6/parser';
import * as PLY from '../../mol-io/reader/ply/parser';
import { parsePsf } from '../../mol-io/reader/psf/parser';
import { PluginContext } from '../../mol-plugin/context';
import { StateObject, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ajaxGetMany } from '../../mol-util/data-source';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { Asset } from '../../mol-util/assets';
import { parseCube } from '../../mol-io/reader/cube/parser';
import { parseDx } from '../../mol-io/reader/dx/parser';

export { Download };
export { DownloadBlob };
export { RawData };
export { ReadFile };
export { ParseBlob };
export { ParseCif };
export { ParseCube };
export { ParsePsf };
export { ParsePly };
export { ParseCcp4 };
export { ParseDsn6 };
export { ParseDx };
export { ImportString };
export { ImportJson };
export { ParseJson };

type Download = typeof Download
const Download = PluginStateTransform.BuiltIn({
    name: 'download',
    display: { name: 'Download', description: 'Download string or binary data from the specified URL' },
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    params: {
        url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
        label: PD.Optional(PD.Text('')),
        isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' }))
    }
})({
    apply({ params: p, cache }, plugin: PluginContext) {
        return Task.create('Download', async ctx => {
            let url = Asset.getUrlAsset(plugin.managers.asset, p.url);
            const asset = await plugin.managers.asset.resolve(url, p.isBinary ? 'binary' : 'string').runInContext(ctx);
            (cache as any).asset = asset;
            return p.isBinary
                ? new SO.Data.Binary(asset.data as Uint8Array, { label: p.label ? p.label : url.url })
                : new SO.Data.String(asset.data as string, { label: p.label ? p.label : url.url });
        });
    },
    dispose({ cache }) {
        ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            b.label = newParams.label || ((typeof newParams.url === 'string') ? newParams.url : newParams.url.url);
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    }
});

type DownloadBlob = typeof DownloadBlob
const DownloadBlob = PluginStateTransform.BuiltIn({
    name: 'download-blob',
    display: { name: 'Download Blob', description: 'Download multiple string or binary data from the specified URLs.' },
    from: SO.Root,
    to: SO.Data.Blob,
    params: {
        sources: PD.ObjectList({
            id: PD.Text('', { label: 'Unique ID' }),
            url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
            isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' })),
            canFail: PD.Optional(PD.Boolean(false, { description: 'Indicate whether the download can fail and not be included in the blob as a result.' }))
        }, e => `${e.id}: ${e.url}`),
        maxConcurrency: PD.Optional(PD.Numeric(4, { min: 1, max: 12, step: 1 }, { description: 'The maximum number of concurrent downloads.' }))
    }
})({
    apply({ params, cache }, plugin: PluginContext) {
        return Task.create('Download Blob', async ctx => {
            const entries: SO.Data.BlobEntry[] = [];
            const data = await ajaxGetMany(ctx, plugin.managers.asset, params.sources, params.maxConcurrency || 4);

            const assets: Asset.Wrapper[] = [];

            for (let i = 0; i < data.length; i++) {
                const r = data[i], src = params.sources[i];
                if (r.kind === 'error') plugin.log.warn(`Download ${r.id} (${src.url}) failed: ${r.error}`);
                else {
                    assets.push(r.result);
                    entries.push(src.isBinary
                        ? { id: r.id, kind: 'binary', data: r.result.data as Uint8Array }
                        : { id: r.id, kind: 'string', data: r.result.data as string });
                }
            }
            (cache as any).assets = assets;
            return new SO.Data.Blob(entries, { label: 'Data Blob', description: `${entries.length} ${entries.length === 1 ? 'entry' : 'entries'}` });
        });
    },
    dispose({ cache }, plugin: PluginContext) {
        const assets: Asset.Wrapper[] | undefined = (cache as any)?.assets;
        if (!assets) return;
        for (const a of assets) {
            a.dispose();
        }
    }
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

type RawData = typeof RawData
const RawData = PluginStateTransform.BuiltIn({
    name: 'raw-data',
    display: { name: 'Raw Data', description: 'Raw data supplied by value.' },
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    params: {
        data: PD.Value<string | number[]>('', { isHidden: true }),
        label: PD.Optional(PD.Text(''))
    }
})({
    apply({ params: p }) {
        return Task.create('Raw Data', async () => {
            if (typeof p.data !== 'string' && isTypedArray(p.data)) {
                throw new Error('Supplied binary data must be a plain array.');
            }
            return typeof p.data === 'string'
                ? new SO.Data.String(p.data as string, { label: p.label ? p.label : 'String' })
                : new SO.Data.Binary(new Uint8Array(p.data), { label: p.label ? p.label : 'Binary' });
        });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.data !== newParams.data) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            b.label = newParams.label || b.label;
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    }
});

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
    apply({ params: p, cache }, plugin: PluginContext) {
        return Task.create('Open File', async ctx => {
            if (p.file === null) {
                plugin.log.error('No file(s) selected');
                return StateObject.Null;
            }

            const asset = await plugin.managers.asset.resolve(p.file, p.isBinary ? 'binary' : 'string').runInContext(ctx);
            (cache as any).asset = asset;
            const o = p.isBinary
                ? new SO.Data.Binary(asset.data as Uint8Array, { label: p.label ? p.label : p.file.name })
                : new SO.Data.String(asset.data as string, { label: p.label ? p.label : p.file.name });

            return o;
        });
    },
    dispose({ cache }) {
        ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.label !== newParams.label) {
            (b.label as string) = newParams.label || oldParams.file?.name || '';
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    },
    isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user loaded files.' })
});

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

type ParseCube = typeof ParseCube
const ParseCube = PluginStateTransform.BuiltIn({
    name: 'parse-cube',
    display: { name: 'Parse Cube', description: 'Parse Cube from String data' },
    from: SO.Data.String,
    to: SO.Format.Cube
})({
    apply({ a }) {
        return Task.create('Parse Cube', async ctx => {
            const parsed = await parseCube(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Cube(parsed.result);
        });
    }
});

type ParsePsf = typeof ParsePsf
const ParsePsf = PluginStateTransform.BuiltIn({
    name: 'parse-psf',
    display: { name: 'Parse PSF', description: 'Parse PSF from String data' },
    from: [SO.Data.String],
    to: SO.Format.Psf
})({
    apply({ a }) {
        return Task.create('Parse PSF', async ctx => {
            const parsed = await parsePsf(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Psf(parsed.result);
        });
    }
});

type ParsePly = typeof ParsePly
const ParsePly = PluginStateTransform.BuiltIn({
    name: 'parse-ply',
    display: { name: 'Parse PLY', description: 'Parse PLY from String data' },
    from: [SO.Data.String],
    to: SO.Format.Ply
})({
    apply({ a }) {
        return Task.create('Parse PLY', async ctx => {
            const parsed = await PLY.parse(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Ply(parsed.result, { label: parsed.result.comments[0] || 'PLY Data' });
        });
    }
});

type ParseCcp4 = typeof ParseCcp4
const ParseCcp4 = PluginStateTransform.BuiltIn({
    name: 'parse-ccp4',
    display: { name: 'Parse CCP4/MRC/MAP', description: 'Parse CCP4/MRC/MAP from Binary data' },
    from: [SO.Data.Binary],
    to: SO.Format.Ccp4
})({
    apply({ a }) {
        return Task.create('Parse CCP4/MRC/MAP', async ctx => {
            const parsed = await CCP4.parse(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Ccp4(parsed.result);
        });
    }
});

type ParseDsn6 = typeof ParseDsn6
const ParseDsn6 = PluginStateTransform.BuiltIn({
    name: 'parse-dsn6',
    display: { name: 'Parse DSN6/BRIX', description: 'Parse CCP4/BRIX from Binary data' },
    from: [SO.Data.Binary],
    to: SO.Format.Dsn6
})({
    apply({ a }) {
        return Task.create('Parse DSN6/BRIX', async ctx => {
            const parsed = await DSN6.parse(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Dsn6(parsed.result);
        });
    }
});

type ParseDx = typeof ParseDx
const ParseDx = PluginStateTransform.BuiltIn({
    name: 'parse-dx',
    display: { name: 'Parse DX', description: 'Parse DX from Binary/String data' },
    from: [SO.Data.Binary, SO.Data.String],
    to: SO.Format.Dx
})({
    apply({ a }) {
        return Task.create('Parse DX', async ctx => {
            const parsed = await parseDx(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Format.Dx(parsed.result);
        });
    }
});

type ImportString = typeof ImportString
const ImportString = PluginStateTransform.BuiltIn({
    name: 'import-string',
    display: { name: 'Import String', description: 'Import given data as a string' },
    from: SO.Root,
    to: SO.Data.String,
    params: {
        data: PD.Value(''),
        label: PD.Optional(PD.Text('')),
    }
})({
    apply({ params: { data, label } }) {
        return new SO.Data.String(data, { label: label || '' });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.data !== newParams.data) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            b.label = newParams.label || '';
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    },
    isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user imported strings.' })
});

type ImportJson = typeof ImportJson
const ImportJson = PluginStateTransform.BuiltIn({
    name: 'import-json',
    display: { name: 'Import JSON', description: 'Import given data as a JSON' },
    from: SO.Root,
    to: SO.Format.Json,
    params: {
        data: PD.Value<any>({}),
        label: PD.Optional(PD.Text('')),
    }
})({
    apply({ params: { data, label } }) {
        return new SO.Format.Json(data, { label: label || '' });
    },
    update({ oldParams, newParams, b }) {
        if (oldParams.data !== newParams.data) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.label !== newParams.label) {
            b.label = newParams.label || '';
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Unchanged;
    },
    isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user imported JSON.' })
});

type ParseJson = typeof ParseJson
const ParseJson = PluginStateTransform.BuiltIn({
    name: 'parse-json',
    display: { name: 'Parse JSON', description: 'Parse JSON from String data' },
    from: [SO.Data.String],
    to: SO.Format.Json
})({
    apply({ a }) {
        return Task.create('Parse JSON', async ctx => {
            const json = await (new Response(a.data)).json(); // async JSON parsing via fetch API
            return new SO.Format.Json(json);
        });
    }
});