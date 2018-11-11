/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { Transform } from './transform';

export { StateObject, StateObjectCell }

interface StateObject<D = any, T extends StateObject.Type = { name: string, typeClass: any }> {
    readonly id: UUID,
    readonly type: T,
    readonly data: D,
    readonly label: string,
    readonly description?: string,
}

namespace StateObject {
    export function factory<T extends Type>() {
        return <D = { }>(type: T) => create<D, T>(type);
    }

    export type Type<Cls extends string = string> = { name: string, typeClass: Cls }
    export type Ctor = { new(...args: any[]): StateObject, type: any }

    export function create<Data, T extends Type>(type: T) {
        return class implements StateObject<Data, T> {
            static type = type;
            static is(obj?: StateObject): obj is StateObject<Data, T> { return !!obj && type === obj.type; }
            id = UUID.create();
            type = type;
            label: string;
            description?: string;
            constructor(public data: Data, props?: { label: string, description?: string }) {
                this.label = props && props.label || type.name;
                this.description = props && props.description;
            }
        }
    }
}

interface StateObjectCell {
    transform: Transform,

    // Which object was used as a parent to create data in this cell
    sourceRef: Transform.Ref | undefined,

    version: string
    status: StateObjectCell.Status,

    errorText?: string,
    obj?: StateObject
}

namespace StateObjectCell {
    export type Status = 'ok' | 'error' | 'pending' | 'processing'

    export interface State {
        isObjectHidden: boolean,
        isTransformHidden: boolean,
        isBinding: boolean,
        isCollapsed: boolean
    }

    export const DefaultState: State = {
        isObjectHidden: false,
        isTransformHidden: false,
        isBinding: false,
        isCollapsed: false
    };
}