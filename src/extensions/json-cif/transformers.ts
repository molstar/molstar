/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';
import { JSONCifFile } from './model';
import { parseJSONCif } from './parser';

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