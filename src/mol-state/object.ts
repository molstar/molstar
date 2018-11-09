/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { Transform } from './transform';

export { StateObject, StateObjectCell }

interface StateObject<P = any, D = any, T = any> {
    readonly id: UUID,
    readonly type: StateObject.Type<T>,
    readonly props: P,
    readonly data: D
}

namespace StateObject {
    export function factory<Type, CommonProps>() {
        return <D = { }, P = {}>(type: Type) => create<P & CommonProps, D, Type>(type);
    }

    export type Type<I = unknown> = I
    export type Ctor = { new(...args: any[]): StateObject, type: any }

    export function create<Props, Data, Type>(type: Type) {
        return class implements StateObject<Props, Data, Type> {
            static type = type;
            static is(obj?: StateObject): obj is StateObject<Props, Data, Type> { return !!obj && type === obj.type; }
            id = UUID.create();
            type = type;
            constructor(public props: Props, public data: Data) { }
        }
    }
}

interface StateObjectCell {
    ref: Transform.Ref,
    version: string
    status: StateObjectCell.Status,

    state: unknown,

    errorText?: string,
    obj?: StateObject
}

namespace StateObjectCell {
    export type Status = 'ok' | 'error' | 'pending' | 'processing'

    export interface State {
        isObjectHidden: boolean,
        isHidden: boolean,
        isBinding: boolean,
        isCollapsed: boolean
    }

    export const DefaultState: State = {
        isObjectHidden: false,
        isHidden: false,
        isBinding: false,
        isCollapsed: false
    };
}