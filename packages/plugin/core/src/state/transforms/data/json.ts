/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateTransformer } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { StringLike } from '@molstar/core/util/string-like';
import { decodeData } from '@molstar/core/util/data-source';

export { ImportJson };
type ImportJson = typeof ImportJson;
const ImportJson = PluginStateTransform.BuiltIn({
  name: 'import-json',
  display: { name: 'Import JSON', description: 'Import given data as a JSON' },
  from: SO.Root,
  to: SO.Format.Json,
  params: {
    data: PD.Value<any>({}),
    label: PD.Optional(PD.Text('')),
  },
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
  isSerializable: () => ({ isSerializable: false, reason: 'Cannot serialize user imported JSON.' }),
});

export { ParseJson };
type ParseJson = typeof ParseJson;
const ParseJson = PluginStateTransform.BuiltIn({
  name: 'parse-json',
  display: { name: 'Parse JSON', description: 'Parse JSON from String, Binary or Response data' },
  from: [SO.Data.String, SO.Data.Binary, SO.Data.Response],
  to: SO.Format.Json,
})({
  apply({ a }) {
    return Task.create('Parse JSON', async (ctx) => {
      const json =
        a instanceof SO.Data.String
          ? await new Response(StringLike.toString(a.data)).json() // async JSON parsing via fetch API
          : await decodeData(ctx, SO.Data.readBytes(a), 'json');
      return new SO.Format.Json(json);
    });
  },
});
