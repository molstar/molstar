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

export { Download }
namespace Download { export interface Params { url: string, isBinary?: boolean, label?: string } }
const Download = PluginStateTransform.Create<SO.Root, SO.Data.String | SO.Data.Binary, Download.Params>({
    name: 'download',
    display: {
        name: 'Download',
        description: 'Download string or binary data from the specified URL'
    },
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    params: () => ({
        url: PD.Text('https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif', { description: 'Resource URL. Must be the same domain or support CORS.' }),
        label: PD.Text('', { isOptional: true }),
        isBinary: PD.Boolean(false, { description: 'If true, download data as binary (string otherwise)', isOptional: true })
    }),
    apply({ params: p }, globalCtx: PluginContext) {
        return Task.create('Download', async ctx => {
            // TODO: track progress
            const data = await globalCtx.fetch(p.url, p.isBinary ? 'binary' : 'string');
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

export { ParseCif }
namespace ParseCif { export interface Params { } }
const ParseCif = PluginStateTransform.Create<SO.Data.String | SO.Data.Binary, SO.Data.Cif, ParseCif.Params>({
    name: 'parse-cif',
    display: {
        name: 'Parse CIF',
        description: 'Parse CIF from String or Binary data'
    },
    from: [SO.Data.String, SO.Data.Binary],
    to: [SO.Data.Cif],
    apply({ a }) {
        return Task.create('Parse CIF', async ctx => {
            const parsed = await (SO.Data.String.is(a) ? CIF.parse(a.data) : CIF.parseBinary(a.data)).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Data.Cif(parsed.result);
        });
    }
});