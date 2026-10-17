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
