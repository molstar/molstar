/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { Transform } from './transform';

export { StateObject, StateObjectBox }

interface StateObject<P = any, D = any, Type = { }> {
    readonly id: UUID,
    readonly type: StateObject.Type,
    readonly props: P,
    readonly data: D
}

namespace StateObject {
    export interface Type<Info = any> {
        info: Info
    }

    export function factory<TypeInfo, CommonProps>() {
        return <D = { }, P = {}>(typeInfo: TypeInfo) => create<P & CommonProps, D, TypeInfo>(typeInfo);
    }

    export type Ctor = { new(...args: any[]): StateObject, type: Type }

    export function create<Props, Data, TypeInfo>(typeInfo: TypeInfo) {
        const dataType: Type<TypeInfo> = { info: typeInfo };
        return class implements StateObject<Props, Data, Type<TypeInfo>> {
            static type = dataType;
            static is(obj?: StateObject): obj is StateObject<Props, Data> { return !!obj && dataType === obj.type; }
            id = UUID.create();
            type = dataType;
            constructor(public props: Props, public data: Data) { }
        }
    }
}

interface StateObjectBox {
    ref: Transform.Ref,
    props: unknown,

    status: StateObjectBox.Status,
    errorText?: string,
    obj?: StateObject,
    version: string
}

namespace StateObjectBox {
    export type Status = 'ok' | 'error' | 'pending' | 'processing'

    export interface Props {
        isVisible: boolean,
        isHidden: boolean,
        isBound: boolean,
        isExpanded: boolean
    }
}