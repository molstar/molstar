/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../context';
import { State, StateObjectCell } from '../../../../mol-state';
import { Task } from '../../../../mol-task';
import { Structure } from '../../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';

export interface StructureRepresentationProvider<P = any> {
    id: string,
    display: { name: string, group: string, description?: string },
    isApplicable?(structure: Structure, plugin: PluginContext): boolean,
    params?(structure: Structure | undefined, plugin: PluginContext): PD.Def<P>,
    // TODO: have create return a "representation structure object" that allows modifications
    apply(state: State, structure: StateObjectCell<PluginStateObject.Molecule.Structure>, params: P, plugin: PluginContext): Task<any> | Promise<void> | void
}

export function StructureRepresentationProvider<P>(repr: StructureRepresentationProvider<P>) { return repr; }