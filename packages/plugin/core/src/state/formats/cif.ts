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
import { CIF } from '@molstar/io/reader/cif';
import { StateObject } from '@molstar/core/state';

/**
 * A response's bytes go through `CIF.parse`, which tells text CIF from BinaryCIF (as `ParseCif`); a binary entry stays
 * strict BinaryCIF (`io/_spec/cif.md` §2).
 */
function parseBlobEntry(e: SO.Data.BlobEntry | SO.Data.ResponseBlobEntry) {
  if (!('kind' in e)) return CIF.parse(e.data.bytes);
  return e.kind === 'string' ? CIF.parse(e.data) : CIF.parseBinary(e.data);
}

export { ParseBlob };
type ParseBlob = typeof ParseBlob;
const ParseBlob = PluginStateTransform.BuiltIn({
  name: 'parse-blob',
  display: { name: 'Parse Blob', description: 'Parse multiple data enties' },
  from: [SO.Data.Blob, SO.Data.ResponseBlob],
  to: SO.Format.Blob,
  params: {
    formats: PD.ObjectList(
      {
        id: PD.Text('', { label: 'Unique ID' }),
        format: PD.Select<'cif'>('cif', [['cif', 'cif']]),
      },
      (e) => `${e.id}: ${e.format}`,
    ),
  },
})({
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Parse Blob', async (ctx) => {
      const map = new Map<string, string>();
      for (const f of params.formats) map.set(f.id, f.format);

      const entries: SO.Format.BlobEntry[] = [];

      for (const e of a.data) {
        if (!map.has(e.id)) continue;

        const parsed = await parseBlobEntry(e).runInContext(ctx);
        if (parsed.isError) throw new Error(`${e.id}: ${parsed.message}`);
        entries.push({ id: e.id, kind: 'cif', data: parsed.result });
      }

      return new SO.Format.Blob(entries, {
        label: 'Format Blob',
        description: `${entries.length} ${entries.length === 1 ? 'entry' : 'entries'}`,
      });
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

export { ParseCif };
type ParseCif = typeof ParseCif;
const ParseCif = PluginStateTransform.BuiltIn({
  name: 'parse-cif',
  display: { name: 'Parse CIF', description: 'Parse CIF from String, Binary or Response data' },
  from: [SO.Data.String, SO.Data.Binary, SO.Data.Response],
  to: SO.Format.Cif,
})({
  apply({ a }) {
    return Task.create('Parse CIF', async (ctx) => {
      // a response's bytes as they are held: text CIF is tokenized from them, never decoded to a string, and the
      // parsed file shares them (`io/_spec/cif.md` §7)
      const data = a instanceof SO.Data.Response ? a.data.bytes : a.data;
      const parsed = await CIF.parse(data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      if (parsed.result.blocks.length === 0) return StateObject.Null;
      return new SO.Format.Cif(parsed.result);
    });
  },
});
