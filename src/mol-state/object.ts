
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';
import { UUID } from 'mol-util';

/** A mutable state object */
export interface StateObject<P = unknown, D = unknown> {
    readonly id: UUID,
    readonly type: StateObject.Type,
    readonly props: P,
    readonly data: D
}

export namespace StateObject {
    export enum StateType {
        // The object has been successfully created
        Ok,
        // An error occured during the creation of the object
        Error,
        // The object is queued to be created
        Pending,
        // The object is currently being created
        Processing
    }

    export interface Type<Info = any> {
        info: Info
    }

    export function factory<TypeInfo, CommonProps>() {
        return <D = { }, P = {}>(typeInfo: TypeInfo) => create<P & CommonProps, D, TypeInfo>(typeInfo);
    }

    export type Ctor = { new(...args: any[]): StateObject, type: Type }

    export function create<Props, Data, TypeInfo>(typeInfo: TypeInfo) {
        const dataType: Type<TypeInfo> = { info: typeInfo };
        return class implements StateObject<Props, Data> {
            static type = dataType;
            static is(obj?: StateObject): obj is StateObject<Props, Data> { return !!obj && dataType === obj.type; }
            id = UUID.create();
            type = dataType;
            ref = 'not set' as Transform.Ref;
            constructor(public props: Props, public data: Data) { }
        }
    }

    export interface Node {
        ref: Transform.Ref,
        state: StateType,
        props: unknown,
        errorText?: string,
        obj?: StateObject,
        version: string
    }
}