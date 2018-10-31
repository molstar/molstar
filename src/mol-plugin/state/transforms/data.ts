/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
import { Task } from 'mol-task';
import CIF from 'mol-io/reader/cif'

export const Download = PluginStateTransform.Create<SO.Root, SO.Data.String | SO.Data.Binary, { url: string, isBinary?: boolean, label?: string }>({
    name: 'download',
    from: [SO.Root],
    to: [SO.Data.String, SO.Data.Binary],
    apply({ params: p }) {
        return Task.create('Download', async ctx => {
            // TODO: track progress
            const req = await fetch(p.url);
            return p.isBinary
                ? new SO.Data.Binary({ label: p.label ? p.label : p.url }, new Uint8Array(await req.arrayBuffer()))
                : new SO.Data.String({ label: p.label ? p.label : p.url }, await req.text());
        });
    }
});

export const ParseCif = PluginStateTransform.Create<SO.Data.String | SO.Data.Binary, SO.Data.Cif, { }>({
    name: 'parse-cif',
    from: [SO.Data.String, SO.Data.Binary],
    to: [SO.Data.Cif],
    apply({ a }) {
        return Task.create('Parse CIF', async ctx => {
            const parsed = await (SO.Data.String.is(a) ? CIF.parse(a.data) : CIF.parseBinary(a.data)).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            return new SO.Data.Cif({ label: 'CIF File' }, parsed.result);
        });
    }
});