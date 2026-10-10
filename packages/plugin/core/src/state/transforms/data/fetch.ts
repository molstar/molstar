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
import { Task } from '@molstar/core/task';
import { Asset } from '@molstar/core/util/assets';
import { StateTransformer, StateObject } from '@molstar/core/state';
import { ajaxGetMany } from '@molstar/core/util/data-source';
import { ungzip } from '@molstar/core/util/zip/zip';
import { utf8ReadLong } from '@molstar/core/util/utf8';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { ColorNames } from '@molstar/core/util/color/names';

export { Download };
type Download = typeof Download;
const Download = PluginStateTransform.BuiltIn({
  name: 'download',
  display: { name: 'Download', description: 'Download string or binary data from the specified URL' },
  from: [SO.Root],
  to: [SO.Data.String, SO.Data.Binary],
  params: {
    url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', {
      description: 'Resource URL. Must be the same domain or support CORS.',
    }),
    label: PD.Optional(PD.Text('')),
    isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' })),
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Download', async (ctx) => {
      const url = Asset.getUrlAsset(plugin.managers.asset, p.url);
      const asset = await plugin.managers.asset.resolve(url, p.isBinary ? 'binary' : 'string').runInContext(ctx);
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
    if (oldParams.url !== newParams.url || oldParams.isBinary !== newParams.isBinary)
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
  to: SO.Data.Blob,
  params: {
    sources: PD.ObjectList(
      {
        id: PD.Text('', { label: 'Unique ID' }),
        url: PD.Url('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', {
          description: 'Resource URL. Must be the same domain or support CORS.',
        }),
        isBinary: PD.Optional(
          PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)' }),
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
  },
})({
  apply({ params, cache }, plugin: PluginContext) {
    return Task.create('Download Blob', async (ctx) => {
      const entries: SO.Data.BlobEntry[] = [];
      const data = await ajaxGetMany(ctx, plugin.managers.asset, params.sources, params.maxConcurrency || 4);

      const assets: Asset.Wrapper[] = [];

      for (let i = 0; i < data.length; i++) {
        const r = data[i],
          src = params.sources[i];
        if (r.kind === 'error') plugin.log.warn(`Download ${r.id} (${src.url}) failed: ${r.error}`);
        else {
          assets.push(r.result);
          entries.push(
            src.isBinary
              ? { id: r.id, kind: 'binary', data: r.result.data as Uint8Array<ArrayBuffer> }
              : { id: r.id, kind: 'string', data: r.result.data as string },
          );
        }
      }
      (cache as any).assets = assets;
      return new SO.Data.Blob(entries, {
        label: 'Data Blob',
        description: `${entries.length} ${entries.length === 1 ? 'entry' : 'entries'}`,
      });
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
  from: [SO.Data.Binary],
  to: [SO.Data.Binary, SO.Data.String],
})({
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Gzip', async (ctx) => {
      const decompressedData = await ungzip(ctx, a.data);
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
  to: [SO.Data.String, SO.Data.Binary],
  params: {
    file: PD.File(),
    label: PD.Optional(PD.Text('')),
    isBinary: PD.Optional(PD.Boolean(false, { description: 'If true, open file as as binary (string otherwise)' })),
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Open File', async (ctx) => {
      if (p.file === null) {
        plugin.log.error('No file(s) selected');
        return StateObject.Null;
      }

      const asset = await plugin.managers.asset.resolve(p.file, p.isBinary ? 'binary' : 'string').runInContext(ctx);
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
    if (oldParams.label !== newParams.label) {
      (b.label as string) = newParams.label || oldParams.file?.name || '';
      return StateTransformer.UpdateResult.Updated;
    }
    return StateTransformer.UpdateResult.Unchanged;
  },
  isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user loaded files.' }),
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
