/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task';
import { StateObject } from './object';
import { Transform } from './transform';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StateAction } from './action';

export interface Transformer<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
    apply(parent: Transform.Ref, params?: P, props?: Partial<Transform.Options>): Transform<A, B, P>,
    toAction(): StateAction<A, void, P>,
    readonly namespace: string,
    readonly id: Transformer.Id,
    readonly definition: Transformer.Definition<A, B, P>
}

export namespace Transformer {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends Transformer<any, any, any>> = T extends Transformer<any, any, infer P> ? P : unknown;
    export type From<T extends Transformer<any, any, any>> = T extends Transformer<infer A, any, any> ? A : unknown;
    export type To<T extends Transformer<any, any, any>> = T extends Transformer<any, infer B, any> ? B : unknown;

    export function is(obj: any): obj is Transformer {
        return !!obj && typeof (obj as Transformer).toAction === 'function' && typeof (obj as Transformer).apply === 'function';
    }

    export interface ApplyParams<A extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        params: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown
    }

    export interface UpdateParams<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        b: B,
        oldParams: P,
        newParams: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown
    }

    export interface AutoUpdateParams<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        b: B,
        oldParams: P,
        newParams: P
    }

    export enum UpdateResult { Unchanged, Updated, Recreate }

    /** Specify default control descriptors for the parameters */
    // export type ParamsDefinition<A extends StateObject = StateObject, P = any> = (a: A, globalCtx: unknown) => { [K in keyof P]: PD.Any }

    export interface DefinitionBase<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
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

        /** Determine if the transformer can be applied automatically on UI change. Default is false. */
        canAutoUpdate?(params: AutoUpdateParams<A, B, P>, globalCtx: unknown): boolean,

        /** Test if the transform can be applied to a given node */
        isApplicable?(a: A, globalCtx: unknown): boolean,

        /** By default, returns true */
        isSerializable?(params: P): { isSerializable: true } | { isSerializable: false; reason: string },

        /** Custom conversion to and from JSON */
        readonly customSerialization?: { toJSON(params: P, obj?: B): any, fromJSON(data: any): P }
    }

    export interface Definition<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> extends DefinitionBase<A, B, P> {
        readonly name: string,
        readonly from: StateObject.Ctor[],
        readonly to: StateObject.Ctor[],
        params?(a: A, globalCtx: unknown): { [K in keyof P]: PD.Any },
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

    export function create<A extends StateObject, B extends StateObject, P extends {} = {}>(namespace: string, definition: Definition<A, B, P>) {
        const { name } = definition;
        const id = `${namespace}.${name}` as Id;

        if (registry.has(id)) {
            throw new Error(`A transform with id '${name}' is already registered. Please pick a unique identifier for your transforms and/or register them only once. This is to ensure that transforms can be serialized and replayed.`);
        }

        const t: Transformer<A, B, P> = {
            apply(parent, params, props) { return Transform.create<A, B, P>(parent, t, params, props); },
            toAction() { return StateAction.fromTransformer(t); },
            namespace,
            id,
            definition
        };
        registry.set(id, t);
        _index(t);

        return t;
    }

    export function factory(namespace: string) {
        return <A extends StateObject, B extends StateObject, P extends {} = {}>(definition: Definition<A, B, P>) => create(namespace, definition);
    }

    export function factory1(namespace: string) {
        return Builder.build(namespace);
    }

    export namespace Builder {
        type ParamDefinition<P> = { [K in keyof P]-?: PD.Base<P[K]> }

        export interface Type<A extends StateObject.Ctor, B extends StateObject.Ctor> {
            name: string,
            from: A | A[],
            to: B | B[]
        }

        export interface TypeAndParams<A extends StateObject.Ctor, B extends StateObject.Ctor, P> extends Type<A, B> {
            params: ParamDefinition<P>
        }

        export interface TypeAndParamProvider<A extends StateObject.Ctor, B extends StateObject.Ctor, P> extends Type<A, B> {
            paramProvider(a: A, globalCtx: unknown): ParamDefinition<P>
        }

        export interface Root {
            <A extends StateObject.Ctor, B extends StateObject.Ctor>(info: Type<A, B>): Define<StateObject.From<A>, StateObject.From<B>, {}>,
            <A extends StateObject.Ctor, B extends StateObject.Ctor, P>(info: TypeAndParams<A, B, P>): Define<StateObject.From<A>, StateObject.From<B>, Params<P>>,
            <A extends StateObject.Ctor, B extends StateObject.Ctor, P>(info: TypeAndParamProvider<A, B, P>): Define<StateObject.From<A>, StateObject.From<B>, Params<P>>        
        }

        type Optionals<P> = { [K in keyof P]-?: undefined extends P[K] ? K : never }[keyof P]
        type NonOptionals<P> = { [K in keyof P]-?: undefined extends P[K] ? never: K }[keyof P]
        type Params<P> = Pick<P, NonOptionals<P>> & Partial<Pick<P, Optionals<P>>>

        export interface Define<A extends StateObject, B extends StateObject, P> {
            (def: DefinitionBase<A, B, P>): Transformer<A, B, P>
        }

        function root(namespace: string, info: Type<any, any> & TypeAndParams<any, any, any> & TypeAndParamProvider<any, any, any>): Define<any, any, any> {
            return def => create(namespace, {
                name: info.name,
                from: info.from instanceof Array ? info.from : [info.from],
                to: info.to instanceof Array ? info.to : [info.to],
                params: info.paramProvider
                    ? info.paramProvider as any
                    : info.params
                    ? () => info.params
                    : void 0,
                ...def
            });
        }

        export function build(namespace: string): Root {
            return (info: any) => root(namespace, info);
        }
    }

    export function build(namespace: string): Builder.Root {
        return Builder.build(namespace);
    }

    export const ROOT = create<any, any, {}>('build-in', {
        name: 'root',
        from: [],
        to: [],
        apply() { throw new Error('should never be applied'); },
        update() { return UpdateResult.Unchanged; }
    })
}