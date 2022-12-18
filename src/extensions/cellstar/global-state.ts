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
import { CellstarEntry } from './entry-root';
import { isDefined } from './helpers';


export const CellstarGlobalStateParams = {
    tryUseGpu: PD.Boolean(true, { description: 'Attempt using GPU for faster rendering. \nCaution: with some hardware setups, this might render some objects incorrectly or not at all.' }),
    selectionMode: PD.Boolean(true, { description: 'Allow selecting/deselecting a segment by clicking on it.' }),
};
export type CellstarGlobalStateParamValues = PD.Values<typeof CellstarGlobalStateParams>;


export class CellstarGlobalState extends PluginStateObject.CreateBehavior<CellstarGlobalStateData>({ name: 'Vol & Seg Global State' }) { }

export class CellstarGlobalStateData extends PluginBehavior.WithSubscribers<CellstarGlobalStateParamValues> {
    private ref: string;
    currentState = new BehaviorSubject(PD.getDefaultValues(CellstarGlobalStateParams));

    constructor(plugin: PluginContext, params: CellstarGlobalStateParamValues) {
        super(plugin, params);
        this.currentState.next(params);
    }

    register(ref: string) {
        this.ref = ref;
    }
    async updateState(plugin: PluginContext, state: Partial<CellstarGlobalStateParamValues>) {
        const oldState = this.currentState.value;

        const promises = [];
        const allEntries = plugin.state.data.selectQ(q => q.ofType(CellstarEntry)).map(cell => cell.obj?.data).filter(isDefined);
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

    static getGlobalState(plugin: PluginContext): CellstarGlobalStateParamValues | undefined {
        return plugin.state.data.selectQ(q => q.ofType(CellstarGlobalState))[0]?.obj?.data.currentState.value;
    }
}
