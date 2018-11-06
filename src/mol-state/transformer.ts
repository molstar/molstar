/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task';
import { StateObject } from './object';
import { Transform } from './transform';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export interface Transformer<A extends StateObject = StateObject, B extends StateObject = StateObject, P = unknown> {
    apply(params?: P, props?: Partial<Transform.Options>): Transform<A, B, P>,
    readonly namespace: string,
    readonly id: Transformer.Id,
    readonly definition: Transformer.Definition<A, B, P>
}

export namespace Transformer {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends Transformer<any, any, any>> = T extends Transformer<any, any, infer P> ? P : unknown;
    export type To<T extends Transformer<any, any, any>> = T extends Transformer<any, infer B, any> ? B : unknown;
    export type ControlsFor<A extends StateObject, Props> = { [P in keyof Props]?: PD.Any }

    export interface ApplyParams<A extends StateObject = StateObject, P = unknown> {
        a: A,
        params: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown
    }

    export interface UpdateParams<A extends StateObject = StateObject, B extends StateObject = StateObject, P = unknown> {
        a: A,
        b: B,
        oldParams: P,
        newParams: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown
    }

    export enum UpdateResult { Unchanged, Updated, Recreate }

    export interface Definition<A extends StateObject = StateObject, B extends StateObject = StateObject, P = unknown> {
        readonly name: string,
        readonly from: { type: StateObject.Type }[],
        readonly to: { type: StateObject.Type }[],
        readonly display?: { readonly name: string, readonly description?: string },

        /**
         * Apply the actual transformation. It must be pure (i.e. with no side effects).
         * Returns a task that produces the result of the result directly.
         */
        apply(params: ApplyParams<A, P>, globalCtx: unknown): Task<B> | B,

        /**
         * Attempts to update the entity in a non-destructive way.
         * For example changing a color scheme of a visual does not require computing new geometry.
         * Return/resolve to undefined if the update is not possible.
         */
        update?(params: UpdateParams<A, B, P>, globalCtx: unknown): Task<UpdateResult> | UpdateResult,

        params?: {
            /** Check the parameters and return a list of errors if the are not valid. */
            default?(a: A, globalCtx: unknown): P,
            /** Specify default control descriptors for the parameters */
            controls?(a: A, globalCtx: unknown): ControlsFor<A, P>,
            /** Check the parameters and return a list of errors if the are not valid. */
            validate?(a: A, params: P, globalCtx: unknown): string[] | undefined,
            /** Optional custom parameter equality. Use deep structural equal by default. */
            areEqual?(oldParams: P, newParams: P): boolean
        }

        /** Test if the transform can be applied to a given node */
        isApplicable?(a: A, globalCtx: unknown): boolean,

        /** By default, returns true */
        isSerializable?(params: P): { isSerializable: true } | { isSerializable: false; reason: string },

        /** Custom conversion to and from JSON */
        customSerialization?: { toJSON(params: P, obj?: B): any, fromJSON(data: any): P }
    }

    const registry = new Map<Id, Transformer<any, any>>();
    const fromTypeIndex: Map<StateObject.Type, Transformer[]> = new Map();

    function _index(tr: Transformer) {
        for (const t of tr.definition.from) {
            if (fromTypeIndex.has(t.type)) {
                fromTypeIndex.get(t.type)!.push(tr);
            } else {
                fromTypeIndex.set(t.type, [tr]);
            }
        }
    }

    export function get(id: string): Transformer {
        const t = registry.get(id as Id);
        if (!t) {
            throw new Error(`A transformer with signature '${id}' is not registered.`);
        }
        return t;
    }

    export function fromType(type: StateObject.Type): ReadonlyArray<Transformer> {
        return fromTypeIndex.get(type) || [];
    }

    export function create<A extends StateObject, B extends StateObject, P>(namespace: string, definition: Definition<A, B, P>) {
        const { name } = definition;
        const id = `${namespace}.${name}` as Id;

        if (registry.has(id)) {
            throw new Error(`A transform with id '${name}' is already registered. Please pick a unique identifier for your transforms and/or register them only once. This is to ensure that transforms can be serialized and replayed.`);
        }

        const t: Transformer<A, B, P> = {
            apply(params, props) { return Transform.create<A, B, P>(t as any, params, props); },
            namespace,
            id,
            definition
        };
        registry.set(id, t);
        _index(t);

        return t;
    }

    export function factory(namespace: string) {
        return <A extends StateObject, B extends StateObject, P>(definition: Definition<A, B, P>) => create(namespace, definition);
    }

    export const ROOT = create<any, any, any>('build-in', {
        name: 'root',
        from: [],
        to: [],
        apply() { throw new Error('should never be applied'); }
    })
}