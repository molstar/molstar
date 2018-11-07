/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, Transformer } from 'mol-state';

export type TypeClass = 'root' | 'data' | 'prop'

export namespace PluginStateObject {
    export type TypeClass = 'Root' | 'Group' | 'Data' | 'Object' | 'Representation' | 'Behavior'
    export interface TypeInfo { name: string, shortName: string, description: string, typeClass: TypeClass }
    export interface Props { label: string, description?: string }

    export const Create = StateObject.factory<TypeInfo, Props>();
}

export namespace PluginStateTransform {
    export const Create = Transformer.factory('ms-plugin');
}