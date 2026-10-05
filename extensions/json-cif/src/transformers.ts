/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '@molstar/plugin/state/objects';
import { StateTransformer } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { ParamDefinition } from '@molstar/core/util/param-definition';
import type { JSONCifFile } from '@molstar/json-cif-extension/model';
import { parseJSONCif } from '@molstar/json-cif-extension/parser';

const Transform = StateTransformer.builderFactory('json-cif');

export const ParseJSONCifFileData = Transform({
    name: 'parse-json-cif-data',
    from: PluginStateObject.Root,
    to: PluginStateObject.Format.Cif,
    params: {
        data: ParamDefinition.Value<JSONCifFile>(undefined as any, { isHidden: true }),
    }
})({
    apply({ params }) {
        return Task.create('Parse JSON Cif', async ctx => {
            const parsed = parseJSONCif(params.data);
            return new PluginStateObject.Format.Cif(parsed, { label: 'CIF Data' });
        });
    }
});