/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommand } from './command';
import { Transform, StateTree, Transformer } from 'mol-state';

export const SetCurrentObject = PluginCommand<{ ref: Transform.Ref }>('ms-data', 'set-current-object');
export const Update = PluginCommand<{ tree: StateTree }>('ms-data', 'update');
export const UpdateObject = PluginCommand<{ ref: Transform.Ref, params: any }>('ms-data', 'update-object');
export const CreateObject = PluginCommand<{ parentRef?: Transform.Ref, transformer: Transformer, params: any }>('ms-data', 'create-object');