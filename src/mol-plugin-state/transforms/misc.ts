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

// export { ValueRefTest };
// type ValueRefTest = typeof ValueRefTest
// const ValueRefTest = PluginStateTransform.BuiltIn({
//     name: 'value-ref-test',
//     display: { name: 'ValueRef Test' },
//     from: SO.Root,
//     to: SO.Data.String,
//     params: (_, ctx: PluginContext) => {
//         const getOptions = () => ctx.state.data.selectQ(q => q.rootsOfType(SO.Molecule.Model)).map(m => [m.transform.ref, m.obj?.label || m.transform.ref] as [string, string]);
//         return {
//             ref: PD.ValueRef<SO.Molecule.Model>(getOptions, ctx.state.data.tryGetCellData, { defaultRef: getOptions()[0]?.[0] })
//         };
//     }
// })({
//     apply({ params }) {
//         const model = params.ref.getValue();
//         console.log(model);
//         return new SO.Data.String(`Model: ${model.label}`, { label: model.label });
//     }
// });