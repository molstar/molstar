/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO: primites

import { StateObject, State, StateObjectCell, StateBuilder, StateTransformer, StateTransform } from '../../../mol-state';
import { RuntimeContext } from '../../../mol-task';
import { PluginContext } from '../../context';

export { StateAction, BuilderAction }

type StateAction<P = any, O extends StateObject = StateObject, R = {}> =
    (cell: StateObjectCell<O>, params: P, ctx: { ctx: RuntimeContext, state: State, plugin: PluginContext }) => Promise<R> | R;
function StateAction<P = any, O extends StateObject = StateObject, R = {}>(action: StateAction<P, O, R>) { return action; }

type BuilderAction<P = any, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer, R = {}> =
    (builder: StateBuilder.To<O, T>, params: P, ctx: { options?: Partial<StateTransform.Options>, plugin: PluginContext }) => R;
function BuilderAction<P = any, O extends StateObject = StateObject, T extends StateTransformer = StateTransformer, R = {}>(action: BuilderAction<P, O, T, R>) { return action; }