
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';

/** A mutable state object */
export interface StateObject<P = unknown, D = unknown> {
    ref: Transform.Ref,
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
        kind: string,
        info: Info
    }

    export function factory<TypeInfo, CommonProps>() {
        return <D = { }, P = {}>(kind: string, info: TypeInfo) => create<P & CommonProps, D, TypeInfo>(kind, info);
    }

    export function create<Props, Data, TypeInfo>(kind: string, typeInfo: TypeInfo) {
        const dataType: Type<TypeInfo> = { kind, info: typeInfo };
        return class implements StateObject<Props, Data> {
            static type = dataType;
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