/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from '../../mol-state';
import { shallowEqualObjects } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';

export { CreateGroup };
type CreateGroup = typeof CreateGroup
const CreateGroup = PluginStateTransform.BuiltIn({
    name: 'create-group',
    display: { name: 'Group' },
    from: [],
    to: SO.Group,
    params: {
        label: PD.Text('Group'),
        description: PD.Optional(PD.Text(''))
    }
})({
    apply({ params }) {
        return new SO.Group({}, params);
    },
    update({ oldParams, newParams, b }) {
        if (shallowEqualObjects(oldParams, newParams)) return StateTransformer.UpdateResult.Unchanged;
        b.label = newParams.label;
        b.description = newParams.description;
        return StateTransformer.UpdateResult.Updated;
    }
});