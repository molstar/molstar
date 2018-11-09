/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { Transformer } from './transformer';
import { UUID } from 'mol-util';

export interface Transform<A extends StateObject = StateObject, B extends StateObject = StateObject, P = unknown> {
    readonly transformer: Transformer<A, B, P>,
    readonly params: P,
    readonly ref: Transform.Ref,
    readonly version: string,
    readonly defaultProps?: unknown
}

export namespace Transform {
    export type Ref = string /* & { '@type': 'transform-ref' } */

    export const RootRef = '-=root=-' as Ref;

    export interface Options { ref?: Ref, defaultProps?: unknown }

    export function create<A extends StateObject, B extends StateObject, P>(transformer: Transformer<A, B, P>, params?: P, options?: Options): Transform<A, B, P> {
        const ref = options && options.ref ? options.ref : UUID.create() as string as Ref;
        return {
            transformer,
            params: params || {} as any,
            ref,
            version: UUID.create(),
            defaultProps: options && options.defaultProps
        }
    }

    export function updateParams<T>(t: Transform, params: any): Transform {
        return { ...t, params, version: UUID.create() };
    }

    export function createRoot(ref: Ref): Transform {
        return create(Transformer.ROOT, {}, { ref });
    }

    export interface Serialized {
        transformer: string,
        params: any,
        ref: string,
        version: string,
        defaultProps?: unknown
    }

    function _id(x: any) { return x; }
    export function toJSON(t: Transform): Serialized {
        const pToJson = t.transformer.definition.customSerialization
            ? t.transformer.definition.customSerialization.toJSON
            : _id;
        return {
            transformer: t.transformer.id,
            params: pToJson(t.params),
            ref: t.ref,
            version: t.version,
            defaultProps: t.defaultProps
        };
    }

    export function fromJSON(t: Serialized): Transform {
        const transformer = Transformer.get(t.transformer);
        const pFromJson = transformer.definition.customSerialization
            ? transformer.definition.customSerialization.toJSON
            : _id;
        return {
            transformer,
            params: pFromJson(t.params),
            ref: t.ref as Ref,
            version: t.version,
            defaultProps: t.defaultProps
        };
    }
}