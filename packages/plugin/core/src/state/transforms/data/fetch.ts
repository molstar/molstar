/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { Asset, AssetManager } from '@molstar/core/util/assets';
import { StateTransformer, StateObject } from '@molstar/core/state';
import { fetchDataMany } from '@molstar/core/util/data-source';
import { isGzip, ungzip } from '@molstar/core/util/zip/zip';
import { utf8ReadLong } from '@molstar/core/util/utf8';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { ColorNames } from '@molstar/core/util/color/names';

const AsResponseDescription =
  'If true (the default), hold the bytes and let the parser read them as text or binary (Response Data); isBinary is ignored. False: String or Binary Data, as isBinary says.';

/** `asResponse` unset is `true`: every parser reads a response (`.v6/data-response`) */
const asResponseOf = (p: { asResponse?: boolean }) => p.asResponse ?? true;

const DecompressParam = PD.Optional(
  PD.Select<'auto' | 'none'>(
    'auto',
    [
      ['auto', 'Auto'],
      ['none', 'None'],
    ],
    {
      description:
        'Auto: gzip data (told by its bytes, not the name) is inflated, as a server that sends .gz files without Content-Encoding needs. None: the bytes as received.',
    },
  ),
);

/** `decompress` unset is `auto` (`util/_spec/assets.md` §5) */
const decompressOf = (p: { decompress?: 'auto' | 'none' }) => p.decompress !== 'none';

/**
 * `asResponse`: reads the asset eagerly, as `fetch()` starts its request, and holds one shared view of its bytes
 * (`util/_spec/assets.md`) for the node's lifetime; the reader chooses the form (`SO.Data.ResponseBody`).
 */
async function readResponse(
  ctx: RuntimeContext,
  plugin: PluginContext,
  asset: Asset,
  cache: any,
  label: string,
  decompress: boolean,
) {
  const held = await AssetManager.resolveData(plugin.managers.asset, asset, 'binary', {
    shared: true,
    decompress,
  }).runInContext(ctx);
  cache.asset = held;
  // turned to text, the body holds one copy: the asset lets go of its bytes (`util/_spec/assets.md` §6)
  const onText = () => plugin.managers.asset.letGoOfBytes?.(asset);
  return new SO.Data.Response(new SO.Data.ResponseBody(held.data, { onText }), { label });
}

export { Download };
type Download = typeof Download;
const Download = PluginStateTransform.BuiltIn({
  name: 'download',
  display: { name: 'Download', description: 'Download string or binary data from the specified URL' },
  from: [SO.Root],
  to: [SO.Data.String, SO.Data.Binary, SO.Data.Response],
  params: {
    url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', {
      description: 'Resource URL. Must be the same domain or support CORS.',
    }),
    label: PD.Optional(PD.Text('')),
    isBinary: PD.Optional(
      PD.Boolean(false, {
        description: 'If true, download data as binary (string otherwise). Ignored with asResponse.',
        hideIf: (p) => asResponseOf(p),
      }),
    ),
    asResponse: PD.Optional(PD.Boolean(true, { description: AsResponseDescription })),
    decompress: DecompressParam,
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Download', async (ctx) => {
      const url = Asset.getUrlAsset(plugin.managers.asset, p.url);
      const decompress = decompressOf(p);
      if (asResponseOf(p)) return readResponse(ctx, plugin, url, cache, p.label || url.url, decompress);
      const asset = await AssetManager.resolveData(plugin.managers.asset, url, p.isBinary ? 'binary' : 'string', {
        decompress,
      }).runInContext(ctx);
      (cache as any).asset = asset;
      return p.isBinary
        ? new SO.Data.Binary(asset.data as Uint8Array<ArrayBuffer>, { label: p.label ? p.label : url.url })
        : new SO.Data.String(asset.data as string, { label: p.label ? p.label : url.url });
    });
  },
  dispose({ cache }) {
    ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
  },
  update({ oldParams, newParams, b }) {
    if (
      oldParams.url !== newParams.url ||
      oldParams.isBinary !== newParams.isBinary ||
      asResponseOf(oldParams) !== asResponseOf(newParams) ||
      decompressOf(oldParams) !== decompressOf(newParams)
    )
      return StateTransformer.UpdateResult.Recreate;
    if (oldParams.label !== newParams.label) {
      b.label = newParams.label || (typeof newParams.url === 'string' ? newParams.url : newParams.url.url);
      return StateTransformer.UpdateResult.Updated;
    }
    return StateTransformer.UpdateResult.Unchanged;
  },
});

export { DownloadBlob };
type DownloadBlob = typeof DownloadBlob;
const DownloadBlob = PluginStateTransform.BuiltIn({
  name: 'download-blob',
  display: { name: 'Download Blob', description: 'Download multiple string or binary data from the specified URLs.' },
  from: SO.Root,
  to: [SO.Data.Blob, SO.Data.ResponseBlob],
  params: {
    sources: PD.ObjectList(
      {
        id: PD.Text('', { label: 'Unique ID' }),
        url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', {
          description: 'Resource URL. Must be the same domain or support CORS.',
        }),
        isBinary: PD.Optional(
          PD.Boolean(false, {
            description: 'If true, download data as binary (string otherwise). Ignored with asResponse.',
          }),
        ),
        canFail: PD.Optional(
          PD.Boolean(false, {
            description: 'Indicate whether the download can fail and not be included in the blob as a result.',
          }),
        ),
      },
      (e) => `${e.id}: ${e.url}`,
    ),
    maxConcurrency: PD.Optional(
      PD.Numeric(4, { min: 1, max: 12, step: 1 }, { description: 'The maximum number of concurrent downloads.' }),
    ),
    asResponse: PD.Optional(
      PD.Boolean(true, {
        description:
          'If true (the default), hold every entry as bytes and let the parser read each as text or binary (Response Blob); isBinary is ignored. False: a Data Blob of strings and binaries, as isBinary says.',
      }),
    ),
    decompress: DecompressParam,
  },
})({
  apply({ params, cache }, plugin: PluginContext) {
    return Task.create('Download Blob', async (ctx) => {
      const decompress = decompressOf(params);
      const asResponse = asResponseOf(params);
      const sources = asResponse
        ? params.sources.map((src) => ({ ...src, isBinary: true, shared: true, decompress }))
        : params.sources.map((src) => ({ ...src, decompress }));
      const data = await fetchDataMany(ctx, plugin.managers.asset, sources, params.maxConcurrency || 4);
      const entries: SO.Data.BlobEntry[] = [];
      const responses: SO.Data.ResponseBlobEntry[] = [];

      const assets: Asset.Wrapper[] = [];

      for (let i = 0; i < data.length; i++) {
        const r = data[i],
          src = params.sources[i];
        if (r.kind === 'error') plugin.log.warn(`Download ${r.id} (${src.url}) failed: ${r.error}`);
        else {
          assets.push(r.result);
          if (asResponse) {
            // one shared hold of each entry's bytes for the node's lifetime, as `Download`'s `asResponse`
            const asset = Asset.getUrlAsset(plugin.managers.asset, src.url);
            const onText = () => plugin.managers.asset.letGoOfBytes?.(asset);
            const body = new SO.Data.ResponseBody(r.result.data as Uint8Array<ArrayBuffer>, { onText });
            responses.push({ id: r.id, data: body });
          } else {
            entries.push(
              src.isBinary
                ? { id: r.id, kind: 'binary', data: r.result.data as Uint8Array<ArrayBuffer> }
                : { id: r.id, kind: 'string', data: r.result.data as string },
            );
          }
        }
      }
      (cache as any).assets = assets;
      const count = asResponse ? responses.length : entries.length;
      const props = { description: `${count} ${count === 1 ? 'entry' : 'entries'}` };
      return asResponse
        ? new SO.Data.ResponseBlob(responses, { label: 'Response Blob', ...props })
        : new SO.Data.Blob(entries, { label: 'Data Blob', ...props });
    });
  },
  dispose({ cache }, plugin: PluginContext) {
    const assets: Asset.Wrapper[] | undefined = (cache as any)?.assets;
    if (!assets) return;
    for (const a of assets) {
      a.dispose();
    }
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

export { DeflateData };
type DeflateData = typeof DeflateData;
const DeflateData = PluginStateTransform.BuiltIn({
  name: 'defalate-data',
  display: { name: 'Deflate', description: 'Deflate compressed data' },
  params: {
    method: PD.Select('gzip', [['gzip', 'gzip']]), // later on we might have to add say brotli
    isString: PD.Boolean(false),
    stringEncoding: PD.Optional(PD.Select('utf-8', [['utf-8', 'UTF8']])),
    label: PD.Optional(PD.Text('')),
  },
  from: [SO.Data.Binary, SO.Data.Response],
  to: [SO.Data.Binary, SO.Data.String],
})({
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Gzip', async (ctx) => {
      // bytes that are not gzip pass through, copied (its readers may change them, the input's may not): a `Download`
      // inflates gzip itself (`decompress`, auto by default)
      const bytes = SO.Data.readBytes(a);
      const decompressedData = isGzip(bytes) ? await ungzip(ctx, bytes) : bytes.slice();
      const label = params.label ? params.label : a.label;
      // handle decoding based on stringEncoding param
      if (params.isString) {
        const textData = utf8ReadLong(decompressedData);
        return new SO.Data.String(textData, { label });
      }
      return new SO.Data.Binary(decompressedData as Uint8Array<ArrayBuffer>, { label });
    });
  },
});

export { RawData };
type RawData = typeof RawData;
const RawData = PluginStateTransform.BuiltIn({
  name: 'raw-data',
  display: { name: 'Raw Data', description: 'Raw data supplied by value.' },
  from: [SO.Root],
  to: [SO.Data.String, SO.Data.Binary],
  params: {
    data: PD.Value<string | number[] | ArrayBuffer | Uint8Array<ArrayBuffer>>('', { isHidden: true }),
    label: PD.Optional(PD.Text('')),
  },
})({
  apply({ params: p }) {
    return Task.create('Raw Data', async () => {
      if (typeof p.data === 'string') {
        return new SO.Data.String(p.data as string, { label: p.label ? p.label : 'String' });
      } else if (Array.isArray(p.data)) {
        return new SO.Data.Binary(new Uint8Array(p.data), { label: p.label ? p.label : 'Binary' });
      } else if (p.data instanceof ArrayBuffer) {
        return new SO.Data.Binary(new Uint8Array(p.data), { label: p.label ? p.label : 'Binary' });
      } else if (p.data instanceof Uint8Array) {
        return new SO.Data.Binary(p.data, { label: p.label ? p.label : 'Binary' });
      } else {
        assertUnreachable(p.data);
      }
    });
  },
  update({ oldParams, newParams, b }) {
    if (oldParams.data !== newParams.data) return StateTransformer.UpdateResult.Recreate;
    if (oldParams.label !== newParams.label) {
      b.label = newParams.label || b.label;
      return StateTransformer.UpdateResult.Updated;
    }
    return StateTransformer.UpdateResult.Unchanged;
  },
  customSerialization: {
    toJSON(p) {
      if (typeof p.data === 'string' || Array.isArray(p.data)) {
        return p;
      } else if (p.data instanceof ArrayBuffer) {
        const v = new Uint8Array(p.data);
        const data = new Array(v.length);
        for (let i = 0, _i = v.length; i < _i; i++) data[i] = v[i];
        return { data, label: p.label };
      } else if (p.data instanceof Uint8Array) {
        const data = new Array(p.data.length);
        for (let i = 0, _i = p.data.length; i < _i; i++) data[i] = p.data[i];
        return { data, label: p.label };
      }
    },
    fromJSON(data: any) {
      return data;
    },
  },
});

export { ReadFile };
type ReadFile = typeof ReadFile;
const ReadFile = PluginStateTransform.BuiltIn({
  name: 'read-file',
  display: { name: 'Read File', description: 'Read string or binary data from the specified file' },
  from: SO.Root,
  to: [SO.Data.String, SO.Data.Binary, SO.Data.Response],
  params: {
    file: PD.File(),
    label: PD.Optional(PD.Text('')),
    isBinary: PD.Optional(
      PD.Boolean(false, {
        description: 'If true, open file as as binary (string otherwise). Ignored with asResponse.',
        hideIf: (p) => asResponseOf(p),
      }),
    ),
    asResponse: PD.Optional(PD.Boolean(true, { description: AsResponseDescription })),
    decompress: DecompressParam,
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Open File', async (ctx) => {
      if (p.file === null) {
        plugin.log.error('No file(s) selected');
        return StateObject.Null;
      }

      const decompress = decompressOf(p);
      if (asResponseOf(p)) return readResponse(ctx, plugin, p.file, cache, p.label || p.file.name, decompress);
      const asset = await AssetManager.resolveData(plugin.managers.asset, p.file, p.isBinary ? 'binary' : 'string', {
        decompress,
      }).runInContext(ctx);
      (cache as any).asset = asset;
      const o = p.isBinary
        ? new SO.Data.Binary(asset.data as Uint8Array<ArrayBuffer>, { label: p.label ? p.label : p.file.name })
        : new SO.Data.String(asset.data as string, { label: p.label ? p.label : p.file.name });

      return o;
    });
  },
  dispose({ cache }) {
    ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
  },
  update({ oldParams, newParams, b }) {
    if (
      oldParams.file?.id !== newParams.file?.id ||
      oldParams.isBinary !== newParams.isBinary ||
      asResponseOf(oldParams) !== asResponseOf(newParams) ||
      decompressOf(oldParams) !== decompressOf(newParams)
    )
      return StateTransformer.UpdateResult.Recreate;
    if (oldParams.label !== newParams.label) {
      (b.label as string) = newParams.label || oldParams.file?.name || '';
      return StateTransformer.UpdateResult.Updated;
    }
    return StateTransformer.UpdateResult.Unchanged;
  },
  isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user loaded files.' }),
});

export { ReadResponse };
type ReadResponse = typeof ReadResponse;
/** Response Data as String or Binary Data, for transforms that do not read a response. */
const ReadResponse = PluginStateTransform.BuiltIn({
  name: 'read-response',
  display: { name: 'Read Response', description: 'Read Response Data as string or binary data' },
  from: SO.Data.Response,
  to: [SO.Data.String, SO.Data.Binary],
  params: {
    isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, read as binary (string otherwise)' })),
  },
})({
  apply({ a, params: p }) {
    return Task.create('Read Response', async (ctx) => {
      const props = { label: a.label, description: a.description };
      // a copy: Binary Data may be changed by its readers, the held bytes may not
      return p.isBinary
        ? new SO.Data.Binary(a.data.bytes.slice(), props)
        : new SO.Data.String(await a.data.text(ctx), props);
    });
  },
  update({ oldParams, newParams }) {
    return oldParams.isBinary !== newParams.isBinary
      ? StateTransformer.UpdateResult.Recreate
      : StateTransformer.UpdateResult.Unchanged;
  },
});

export { ImportString };
type ImportString = typeof ImportString;
const ImportString = PluginStateTransform.BuiltIn({
  name: 'import-string',
  display: { name: 'Import String', description: 'Import given data as a string' },
  from: SO.Root,
  to: SO.Data.String,
  params: {
    data: PD.Value(''),
    label: PD.Optional(PD.Text('')),
  },
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
  isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user imported strings.' }),
});

export { LazyVolume };
type LazyVolume = typeof LazyVolume;
const LazyVolume = PluginStateTransform.BuiltIn({
  name: 'lazy-volume',
  display: { name: 'Lazy Volume', description: 'A placeholder for lazy loaded volume representation' },
  from: SO.Root,
  to: SO.Volume.Lazy,
  params: {
    url: PD.Url(''),
    isBinary: PD.Boolean(false),
    format: PD.Text('ccp4'), // TODO: use Select based on available formats
    entryId: PD.Value<string | string[]>('', { isHidden: true }),
    isovalues: PD.ObjectList(
      {
        type: PD.Text<'absolute' | 'relative'>('relative'), // TODO: Select
        value: PD.Numeric(0),
        color: PD.Color(ColorNames.black),
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
        volumeIndex: PD.Numeric(0),
      },
      (e) => `${e.type} ${e.value}`,
    ),
  },
})({
  apply({ a, params }) {
    return Task.create('Lazy Volume', async (ctx) => {
      const entryId = Array.isArray(params.entryId) ? params.entryId.join(', ') : params.entryId;
      return new SO.Volume.Lazy(params, { label: `${entryId || params.url}`, description: 'Lazy Volume' });
    });
  },
});
