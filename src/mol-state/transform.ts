/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from './transformer';
import { UUID } from 'mol-util';

export { Transform as StateTransform }

interface Transform<T extends StateTransformer = StateTransformer> {
    readonly parent: Transform.Ref,
    readonly transformer: T,
    readonly props: Transform.Props,
    readonly ref: Transform.Ref,
    readonly params?: StateTransformer.Params<T>,
    readonly version: string
}

namespace Transform {
    export type Ref = string
    export type Transformer<T extends Transform> = T extends Transform<infer S> ? S : never

    export const RootRef = '-=root=-' as Ref;

    export interface Props {
        tag?: string
        isGhost?: boolean,
        // determine if the corresponding cell can be deleted by the user.
        isLocked?: boolean
    }

    export interface Options {
        ref?: string,
        props?: Props
    }

    export function create<T extends StateTransformer>(parent: Ref, transformer: T, params?: StateTransformer.Params<T>, options?: Options): Transform<T> {
        const ref = options && options.ref ? options.ref : UUID.create22() as string as Ref;
        return {
            parent,
            transformer,
            props: (options && options.props) || { },
            ref,
            params,
            version: UUID.create22()
        }
    }

    export function withParams(t: Transform, params: any): Transform {
        return { ...t, params, version: UUID.create22() };
    }

    export function withParent(t: Transform, parent: Ref): Transform {
        return { ...t, parent, version: UUID.create22() };
    }

    export function withNewVersion(t: Transform): Transform {
        return { ...t, version: UUID.create22() };
    }

    export function createRoot(props?: Props): Transform {
        return create(RootRef, StateTransformer.ROOT, {}, { ref: RootRef, props });
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
            params: t.params ? pToJson(t.params) : void 0,
            props: t.props,
            ref: t.ref,
            version: t.version
        };
    }

    export function fromJSON(t: Serialized): Transform {
        const transformer = StateTransformer.get(t.transformer);
        const pFromJson = transformer.definition.customSerialization
            ? transformer.definition.customSerialization.toJSON
            : _id;
        return {
            parent: t.parent as Ref,
            transformer,
            params: t.params ? pFromJson(t.params) : void 0,
            props: t.props,
            ref: t.ref as Ref,
            version: t.version
        };
    }
}