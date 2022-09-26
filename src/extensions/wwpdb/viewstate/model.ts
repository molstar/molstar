/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { wwPDBViewState } from './state';

export class wwPDBViewStateObject extends PluginStateObject.Create<{ state: wwPDBViewState }>({ name: 'wwPDB View State', typeClass: 'Object' }) { }

export const LoadWwPDBViewState = StateAction.build({
    display: { name: 'Load wwPDB View State' },
    params: {
        file: ParamDefinition.File({ accept: '.json' })
    },
    from: PluginStateObject.Root
})(({ state, params }, plugin: PluginContext) => Task.create('Load wwPDB View State', async taskCtx => {
    try {
        const json = JSON.parse(await params.file?.file?.text() ?? '');
        if (!json) return;
        return await state.build().toRoot().apply(CreateWwPDBViewState, { state: json }).commit();
    } catch (err) {
        plugin.log.error(`wwPDB View State: ${err}`);
    }
}));

const WwPDBTransformerFactory = StateTransformer.builderFactory('wwpdb');

export const CreateWwPDBViewState = WwPDBTransformerFactory({
    name: 'create-viewstate',
    display: 'Create wwPDB View State',
    from: PluginStateObject.Root,
    to: wwPDBViewStateObject,
    params: {
        state: ParamDefinition.Value<wwPDBViewState>({} as any, { isHidden: true }),
        label: ParamDefinition.Text(),
    },
})({
    apply({ params }) {
        return new wwPDBViewStateObject({ state: params.state }, { label: params.label ?? 'wwPDB View State' });
    }
});
