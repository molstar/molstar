/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { VolsegEntry } from './entry-root';
import { isDefined } from './helpers';


export const VolsegGlobalStateParams = {
    tryUseGpu: PD.Boolean(true, { description: 'Attempt using GPU for faster rendering. \nCaution: with some hardware setups, this might render some objects incorrectly or not at all.' }),
    selectionMode: PD.Boolean(true, { description: 'Allow selecting/deselecting a segment by clicking on it.' }),
};
export type VolsegGlobalStateParamValues = PD.Values<typeof VolsegGlobalStateParams>;


export class VolsegGlobalState extends PluginStateObject.CreateBehavior<VolsegGlobalStateData>({ name: 'Vol & Seg Global State' }) { }

export class VolsegGlobalStateData extends PluginBehavior.WithSubscribers<VolsegGlobalStateParamValues> {
    private ref: string;
    currentState = new BehaviorSubject(PD.getDefaultValues(VolsegGlobalStateParams));

    constructor(plugin: PluginContext, params: VolsegGlobalStateParamValues) {
        super(plugin, params);
        this.currentState.next(params);
    }

    register(ref: string) {
        this.ref = ref;
    }
    unregister() {
        this.ref = '';
    }
    isRegistered() {
        return this.ref !== '';
    }
    async updateState(plugin: PluginContext, state: Partial<VolsegGlobalStateParamValues>) {
        const oldState = this.currentState.value;

        const promises = [];
        const allEntries = plugin.state.data.selectQ(q => q.ofType(VolsegEntry)).map(cell => cell.obj?.data).filter(isDefined);
        if (state.tryUseGpu !== undefined && state.tryUseGpu !== oldState.tryUseGpu) {
            for (const entry of allEntries) {
                promises.push(entry.setTryUseGpu(state.tryUseGpu));
            }
        }
        if (state.selectionMode !== undefined && state.selectionMode !== oldState.selectionMode) {
            for (const entry of allEntries) {
                promises.push(entry.setSelectionMode(state.selectionMode));
            }
        }
        await Promise.all(promises);
        await plugin.build().to(this.ref).update(state).commit();
    }

    static getGlobalState(plugin: PluginContext): VolsegGlobalStateParamValues | undefined {
        return plugin.state.data.selectQ(q => q.ofType(VolsegGlobalState))[0]?.obj?.data.currentState.value;
    }
}
