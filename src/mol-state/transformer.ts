/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../mol-task';
import { StateObject, StateObjectCell } from './object';
import { StateTransform } from './transform';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { StateAction } from './action';
import { capitalize } from '../mol-util/string';
import { StateTreeSpine } from './tree/spine';

export { Transformer as StateTransformer };

interface Transformer<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = any> {
    apply(parent: StateTransform.Ref, params?: P, props?: Partial<StateTransform.Options>): StateTransform<this>,
    toAction(): StateAction<A, void, P>,
    readonly namespace: string,
    readonly id: Transformer.Id,
    readonly definition: Transformer.Definition<A, B, P>,
    /** create a fresh copy of the params which can be edited in place */
    createDefaultParams(a: A, globalCtx: unknown): P
}

namespace Transformer {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends Transformer<any, any, any>> = T extends Transformer<any, any, infer P> ? P : unknown;
    export type From<T extends Transformer<any, any, any>> = T extends Transformer<infer A, any, any> ? A : unknown;
    export type To<T extends Transformer<any, any, any>> = T extends Transformer<any, infer B, any> ? B : unknown;
    export type Cell<T extends Transformer<any, any, any>> = T extends Transformer<any, infer B, any> ? StateObjectCell<B> : unknown;

    export function getParamDefinition<T extends Transformer>(t: T, a: From<T> | undefined, globalCtx: unknown): PD.For<Params<T>> {
        return t.definition.params ? t.definition.params(a, globalCtx) as any : { } as any;
    }

    export function is(obj: any): obj is Transformer {
        return !!obj && typeof (obj as Transformer).toAction === 'function' && typeof (obj as Transformer).apply === 'function';
    }

    export interface ApplyParams<A extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        params: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown,
        spine: StateTreeSpine,
        dependencies?: { [k: string]: StateObject<unknown> }
    }

    export interface UpdateParams<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        b: B,
        oldParams: P,
        newParams: P,
        /** A cache object that is purged each time the corresponding StateObject is removed or recreated. */
        cache: unknown,
        spine: StateTreeSpine,
        dependencies?: { [k: string]: StateObject<unknown> }
    }

    export interface AutoUpdateParams<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
        a: A,
        b: B,
        oldParams: P,
        newParams: P
    }

    export interface DisposeParams<B extends StateObject = StateObject, P extends {} = {}> {
        b: B | undefined,
        params: P | undefined,
        cache: unknown
    }

    export enum UpdateResult { Unchanged, Updated, Recreate, Null }

    /** Specify default control descriptors for the parameters */
    // export type ParamsDefinition<A extends StateObject = StateObject, P = any> = (a: A, globalCtx: unknown) => { [K in keyof P]: PD.Any }

    export interface DefinitionBase<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> {
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

        /** Parameter interpolation */
        interpolate?(src: P, target: P, t: number, globalCtx: unknown): P

        /**
         * Cleanup resources
         *
         * Automatically called on deleting an object and on recreating it
         * (i.e. when update returns UpdateResult.Recreate or UpdateResult.Null)
         *
         * Not called on UpdateResult.Updated because the resources might not
         * have been invalidated. In this case, the possible cleanup has to be handled
         * manually.
         */
        dispose?(params: DisposeParams<B, P>, globalCtx: unknown): void

        /** Custom conversion to and from JSON */
        readonly customSerialization?: { toJSON(params: P, obj?: B): any, fromJSON(data: any): P }
    }

    export interface Definition<A extends StateObject = StateObject, B extends StateObject = StateObject, P extends {} = {}> extends DefinitionBase<A, B, P> {
        readonly name: string,
        readonly from: StateObject.Ctor[],
        readonly to: StateObject.Ctor[],
        readonly display: { readonly name: string, readonly description?: string },
        params?(a: A | undefined, globalCtx: unknown): { [K in keyof P]: PD.Any },

        /**
         * Decorators are special Transformers mapping the object to the same type.
         *
         * Special rules apply:
         * - applying decorator always "inserts" it instead
         * - applying to a decorated Transform is applied to the decorator instead (transitive)
         */
        readonly isDecorator?: boolean
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

    export function getAll() {
        return Array.from(registry.values());
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
            apply(parent, params, props) { return StateTransform.create<Transformer<A, B, P>>(parent, t, params, props); },
            toAction() { return StateAction.fromTransformer(t); },
            namespace,
            id,
            definition,
            createDefaultParams(a, globalCtx) { return definition.params ? PD.getDefaultValues( definition.params(a, globalCtx)) : {} as any; }
        };
        registry.set(id, t);
        _index(t);

        return t;
    }

    export function factory(namespace: string) {
        return <A extends StateObject, B extends StateObject, P extends {} = {}>(definition: Definition<A, B, P>) => create(namespace, definition);
    }

    export function builderFactory(namespace: string) {
        return Builder.build(namespace);
    }

    export namespace Builder {
        export interface Type<A extends StateObject.Ctor, B extends StateObject.Ctor, P extends { }> {
            name: string,
            from: A | A[],
            to: B | B[],
            /** The source StateObject can be undefined: used for generating docs. */
            params?: PD.For<P> | ((a: StateObject.From<A> | undefined, globalCtx: any) => PD.For<P>),
            display?: string | { name: string, description?: string },
            isDecorator?: boolean
        }

        export interface Root {
            <A extends StateObject.Ctor, B extends StateObject.Ctor, P extends { }>(info: Type<A, B, P>): Define<StateObject.From<A>, StateObject.From<B>, PD.Normalize<P>>
        }

        export interface Define<A extends StateObject, B extends StateObject, P> {
            (def: DefinitionBase<A, B, P>): Transformer<A, B, P>
        }

        function root(namespace: string, info: Type<any, any, any>): Define<any, any, any> {
            return def => create(namespace, {
                name: info.name,
                from: info.from instanceof Array ? info.from : [info.from],
                to: info.to instanceof Array ? info.to : [info.to],
                display: typeof info.display === 'string'
                    ? { name: info.display }
                    : !!info.display
                        ? info.display
                        : { name: capitalize(info.name.replace(/[-]/g, ' ')) },
                params: typeof info.params === 'object'
                    ? () => info.params as any
                    : !!info.params
                        ? info.params as any
                        : void 0,
                isDecorator: info.isDecorator,
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
        display: { name: 'Root', description: 'For internal use.' },
        apply() { throw new Error('should never be applied'); },
        update() { return UpdateResult.Unchanged; }
    });
}