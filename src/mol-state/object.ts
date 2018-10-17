/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A mutable state object */
export interface StateObject<T extends StateObject.Type = any> {
    '@type': T,
    label: string,
    version: number
}

export namespace StateObject {
    export type TypeOf<T>
        = T extends StateObject<infer X> ? [X]
        : T extends [StateObject<infer X>] ? [X]
        : T extends [StateObject<infer X>, StateObject<infer Y>] ? [X, Y]
        : unknown[];

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

    export type Type = string & { '@type': 'state-object-type' }
}