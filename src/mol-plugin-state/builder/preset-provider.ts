/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectRef } from '../../mol-state';
import { PluginContext } from '../../mol-plugin/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

export interface PresetProvider<O extends StateObject = StateObject, P = any, S = {}> {
    id: string,
    display: { name: string, group?: string, description?: string },
    isApplicable?(a: O, plugin: PluginContext): boolean,
    params?(a: O | undefined, plugin: PluginContext): PD.For<P>,
    apply(a: StateObjectRef<O>, params: P, plugin: PluginContext): Promise<S> | S
}