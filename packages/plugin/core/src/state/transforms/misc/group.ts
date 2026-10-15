/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shallowEqualObjects } from '@molstar/core/util';
import { StateTransformer } from '@molstar/core/state';

export { CreateGroup };
type CreateGroup = typeof CreateGroup;
const CreateGroup = PluginStateTransform.BuiltIn({
  name: 'create-group',
  display: { name: 'Group' },
  from: [],
  to: SO.Group,
  params: {
    label: PD.Text('Group'),
    description: PD.Optional(PD.Text('')),
  },
})({
  apply({ params }) {
    return new SO.Group({}, params);
  },
  update({ oldParams, newParams, b }) {
    if (shallowEqualObjects(oldParams, newParams)) return StateTransformer.UpdateResult.Unchanged;
    b.label = newParams.label;
    b.description = newParams.description;
    return StateTransformer.UpdateResult.Updated;
  },
});
