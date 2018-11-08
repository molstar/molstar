/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, Transformer } from 'mol-state';
import { Representation } from 'mol-repr';
import { PluginBehavior } from 'mol-plugin/behavior/behavior';

export type TypeClass = 'root' | 'data' | 'prop'

export namespace PluginStateObject {
    export type Any = StateObject<Props, any, TypeInfo>

    export type TypeClass = 'Root' | 'Group' | 'Data' | 'Object' | 'Representation3D' | 'Behavior'
    export interface TypeInfo { name: string, shortName: string, description: string, typeClass: TypeClass }
    export interface Props { label: string, description?: string }

    export const Create = StateObject.factory<TypeInfo, Props>();

    export function isRepresentation3D(o?: Any): o is StateObject<Props, Representation.Any, TypeInfo> {
        return !!o && o.type.typeClass === 'Representation3D';
    }

    export function isBehavior(o?: Any): o is StateObject<Props, PluginBehavior, TypeInfo> {
        return !!o && o.type.typeClass === 'Behavior';
    }

    export function CreateRepresentation3D<T extends Representation.Any>(type: { name: string, shortName: string, description: string }) {
        return Create<T>({ ...type, typeClass: 'Representation3D' })
    }

    export function CreateBehavior<T extends PluginBehavior>(type: { name: string, shortName: string, description: string }) {
        return Create<T>({ ...type, typeClass: 'Behavior' })
    }
}

export namespace PluginStateTransform {
    export const Create = Transformer.factory('ms-plugin');
}