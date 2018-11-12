/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell } from './object';
import { Transformer } from './transformer';
import { UUID } from 'mol-util';

export interface Transform<A extends StateObject = StateObject, B extends StateObject = StateObject, P = unknown> {
    readonly parent: Transform.Ref,
    readonly transformer: Transformer<A, B, P>,
    readonly params: P,
    readonly props: Transform.Props,
    readonly ref: Transform.Ref,
    readonly version: string,
    readonly cellState: StateObjectCell.State
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
        props?: Props,
        cellState?: Partial<StateObjectCell.State>
    }

    export function create<A extends StateObject, B extends StateObject, P>(parent: Ref, transformer: Transformer<A, B, P>, params?: P, options?: Options): Transform<A, B, P> {
        const ref = options && options.ref ? options.ref : UUID.create22() as string as Ref;
        return {
            parent,
            transformer,
            params: params || {} as any,
            props: (options && options.props) || { },
            ref,
            version: UUID.create22(),
            cellState: { ...StateObjectCell.DefaultState, ...(options && options.cellState) }
        }
    }

    export function withParams<T>(t: Transform, params: any): Transform {
        return { ...t, params, version: UUID.create22() };
    }

    export function withCellState<T>(t: Transform, state: Partial<StateObjectCell.State>): Transform {
        return { ...t, cellState: { ...t.cellState, ...state } };
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
        version: string,
        cellState: StateObjectCell.State,
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
            version: t.version,
            cellState: t.cellState,
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
            version: t.version,
            cellState: t.cellState,
        };
    }
}