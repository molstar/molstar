/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { Asset } from '@molstar/core/util/assets';
import { unzip } from '@molstar/core/util/zip/zip';
import { CIF } from '@molstar/io/reader/cif';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';

// AlphaFold DB's PAE as a state cell of the model it belongs to: a download and a JSON parse in one transform. The
// built-in `download` applies to the root only, so its cell could not go with the structure; this one is a child of the
// model cell, so removing or reloading the model removes the PAE (and its parsed arrays) with it, and a saved session
// restores it with the model (its address is the parameter, fetched again on load).

export { FetchPaeJson };
type FetchPaeJson = typeof FetchPaeJson;
const FetchPaeJson = PluginStateTransform.BuiltIn({
  name: 'model-archive-fetch-pae-json',
  display: { name: 'AlphaFold DB PAE', description: "Download and parse a model's PAE JSON" },
  from: SO.Molecule.Model,
  to: SO.Format.Json,
  params: {
    url: PD.Url('', { description: 'The PAE JSON on AlphaFold DB.' }),
    label: PD.Optional(PD.Text('')),
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Fetch PAE', async (ctx) => {
      const url = Asset.getUrlAsset(plugin.managers.asset, p.url);
      const asset = await plugin.managers.asset.resolve(url, 'string').runInContext(ctx);
      (cache as { asset?: Asset.Wrapper }).asset = asset;
      const json = await new Response(String(asset.data)).json();
      return new SO.Format.Json(json, { label: p.label || url.url });
    });
  },
  dispose({ cache }) {
    (cache as { asset?: Asset.Wrapper } | undefined)?.asset?.dispose();
  },
  update({ oldParams, newParams }) {
    return oldParams.url !== newParams.url
      ? StateTransformer.UpdateResult.Recreate
      : StateTransformer.UpdateResult.Unchanged;
  },
});

// Model Archive's PAE as a state cell of the model, as AlphaFold DB's above: the associated file the model's own file
// names (`associated.ts`), a CIF or an entry of a zip, downloaded, unpacked and parsed in one transform.

/** What an associated archive may hold before anything is inflated (Model Archive's: a PAE and an MSA, tens of MB). */
const ArchiveLimits = { maxEntries: 256, maxBytes: 1 << 30 };

export { FetchAssociatedPae };
type FetchAssociatedPae = typeof FetchAssociatedPae;
const FetchAssociatedPae = PluginStateTransform.BuiltIn({
  name: 'model-archive-fetch-associated-pae',
  display: { name: 'Model Archive PAE', description: "Download and parse a model's associated PAE file" },
  from: SO.Molecule.Model,
  to: SO.Format.Cif,
  params: {
    url: PD.Url('', { description: 'The associated file: a CIF, or a zip with the CIF in it.' }),
    path: PD.Optional(PD.Text('', { description: "The CIF's path in the zip." })),
    label: PD.Optional(PD.Text('')),
  },
})({
  apply({ params: p, cache }, plugin: PluginContext) {
    return Task.create('Fetch PAE', async (ctx) => {
      const url = Asset.getUrlAsset(plugin.managers.asset, p.url);
      const asset = await plugin.managers.asset.resolve(url, 'binary').runInContext(ctx);
      let bytes: Uint8Array = asset.data;
      if (p.path) {
        // The archive is not needed once its entry is out (a PAE's zip is MBs, kept beside the parse otherwise).
        try {
          const buf = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
          const entries = await unzip(ctx, buf, false, ArchiveLimits, (name) => name === p.path);
          const entry = entries[p.path];
          if (!(entry instanceof Uint8Array)) throw new Error(`the archive has no ${p.path}`);
          bytes = entry;
        } finally {
          asset.dispose();
        }
      } else {
        (cache as { asset?: Asset.Wrapper }).asset = asset;
      }
      const parsed = await CIF.parse(bytes).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      // An empty file is an error, not a Null cell: the fetch reports it and leaves nothing under the model.
      if (parsed.result.blocks.length === 0) throw new Error('the file has no data block');
      return new SO.Format.Cif(parsed.result, { label: p.label || p.path || url.url });
    });
  },
  dispose({ cache }) {
    (cache as { asset?: Asset.Wrapper } | undefined)?.asset?.dispose();
  },
  update({ oldParams, newParams }) {
    return oldParams.url !== newParams.url || oldParams.path !== newParams.path
      ? StateTransformer.UpdateResult.Recreate
      : StateTransformer.UpdateResult.Unchanged;
  },
});
