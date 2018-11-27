/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { Transformer } from './transformer';
import { UUID } from 'mol-util';

export interface Transform<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
    readonly parent: Transform.Ref,
    readonly transformer: Transformer<A, B, P>,
    readonly props: Transform.Props,
    readonly ref: Transform.Ref,
    readonly params: P,
    readonly version: string
}

export namespace Transform {
    export type Ref = string

    export const RootRef = '-=root=-' as Ref;

    export interface Props {
        tag?: string
        isGhost?: boolean,
        isBinding?: boolean
    }

    export interface Options {
        ref?: string,
        props?: Props
    }

    export function create<A extends StateObject, B extends StateObject, P extends {} = {}>(parent: Ref, transformer: Transformer<A, B, P>, params?: P, options?: Options): Transform<A, B, P> {
        const ref = options && options.ref ? options.ref : UUID.create22() as string as Ref;
        return {
            parent,
            transformer,
            props: (options && options.props) || { },
            ref,
            params: params as any,
            version: UUID.create22()
        }
    }

    export function withParams<T>(t: Transform, params: any): Transform {
        return { ...t, params, version: UUID.create22() };
    }

    export function createRoot(): Transform {
        return create(RootRef, Transformer.ROOT, {}, { ref: RootRef });
    }

    export interface Serialized {
        parent: string,
        transformer: string,
        params: any,
        props: Props,
        ref: string,
        version: string
    }

    function _id(x: any) { return x; }
    export function toJSON(t: Transform): Serialized {
        const pToJson = t.transformer.definition.customSerialization
            ? t.transformer.definition.customSerialization.toJSON
            : _id;
        return {
            parent: t.parent,
            transformer: t.transformer.id,
            params: pToJson(t.params),
            props: t.props,
            ref: t.ref,
            version: t.version
        };
    }

    export function fromJSON(t: Serialized): Transform {
        const transformer = Transformer.get(t.transformer);
        const pFromJson = transformer.definition.customSerialization
            ? transformer.definition.customSerialization.toJSON
            : _id;
        return {
            parent: t.parent as Ref,
            transformer,
            params: pFromJson(t.params),
            props: t.props,
            ref: t.ref as Ref,
            version: t.version
        };
    }
}