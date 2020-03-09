/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { State, StateObjectCell } from '../../../mol-state';
import { RuntimeContext } from '../../../mol-task';
import { Structure } from '../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';

export interface StructureRepresentationProvider<P = any, S = {}> {
    id: string,
    display: { name: string, group: string, description?: string },
    isApplicable?(structure: Structure, plugin: PluginContext): boolean,
    params?(structure: Structure | undefined, plugin: PluginContext): PD.Def<P>,
    apply(ctx: RuntimeContext, state: State, structure: StateObjectCell<PluginStateObject.Molecule.Structure>, params: P, plugin: PluginContext): Promise<S> | S,
    // TODO: Custom remove function for more complicated things
    // remove?(state: State, ref: string, plugin: PluginContext): void
}

export namespace StructureRepresentationProvider {
    export type Params<P extends StructureRepresentationProvider> = P extends StructureRepresentationProvider<infer T> ? T : never;
    export type State<P extends StructureRepresentationProvider> = P extends StructureRepresentationProvider<infer _, infer S> ? S : never;
}

export const enum RepresentationProviderTags {
    Representation = 'preset-structure-representation',
    Component = 'preset-structure-component'
}

export function StructureRepresentationProvider<P, S>(repr: StructureRepresentationProvider<P, S>) { return repr; }