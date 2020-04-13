/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../mol-task';
import { UUID } from '../mol-util';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { StateObject, StateObjectCell } from './object';
import { State } from './state';
import { StateTransformer } from './transformer';
import { StateTransform } from './transform';

export { StateAction };

interface StateAction<A extends StateObject = StateObject, T = any, P extends {} = {}> {
    create(params: P): StateAction.Instance,
    readonly id: UUID,
    readonly definition: StateAction.Definition<A, T, P>,
    /** create a fresh copy of the params which can be edited in place */
    createDefaultParams(a: A, globalCtx: unknown): P
}

namespace StateAction {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends StateAction<any, any, any>> = T extends StateAction<any, any, infer P> ? P : unknown;
    export type ReType<T extends StateAction<any, any, any>> = T extends StateAction<any, infer T, any> ? T : unknown;
    export type ControlsFor<Props> = { [P in keyof Props]?: PD.Any }

    export interface Instance {
        action: StateAction,
        params: any
    }

    export interface ApplyParams<A extends StateObject = StateObject, P extends {} = {}> {
        ref: string,
        cell: StateObjectCell,
        a: A,
        state: State,
        params: P
    }

    export interface DefinitionBase<A extends StateObject = StateObject, T = any, P extends {} = {}> {
        /**
         * Apply an action that modifies the State specified in Params.
         */
        run(params: ApplyParams<A, P>, globalCtx: unknown): T | Task<T>,

        /** Test if the transform can be applied to a given node */
        isApplicable?(a: A, aTransform: StateTransform<StateTransformer<any, A, any>>, globalCtx: unknown): boolean
    }

    export interface Definition<A extends StateObject = StateObject, T = any, P extends {} = {}> extends DefinitionBase<A, T, P> {
        readonly from: StateObject.Ctor[],
        readonly display: { readonly name: string, readonly description?: string },
        params?(a: A, globalCtx: unknown): { [K in keyof P]: PD.Any }
    }

    export function create<A extends StateObject, T, P extends {} = {}>(definition: Definition<A, T, P>): StateAction<A, T, P> {
        const action: StateAction<A, T, P> = {
            create(params) { return { action, params }; },
            id: UUID.create22(),
            definition,
            createDefaultParams(a, globalCtx) { return definition.params ? PD.getDefaultValues( definition.params(a, globalCtx)) : {} as any; }
        };
        return action;
    }

    export function fromTransformer<T extends StateTransformer>(transformer: T) {
        const def = transformer.definition;
        return create<StateTransformer.From<T>, void, StateTransformer.Params<T>>({
            from: def.from,
            display: def.display,
            params: def.params as StateTransformer.Definition<StateTransformer.From<T>, any, StateTransformer.Params<T>>['params'],
            isApplicable: transformer.definition.isApplicable
                ? (a, t, ctx) => transformer.definition.isApplicable!(a, ctx)
                : void 0,
            run({ cell, state, params }) {
                const tree = state.build().to(cell.transform.ref).apply(transformer, params);
                return state.updateTree(tree) as unknown as Task<void>;
            }
        });
    }

    export namespace Builder {
        export interface Type<A extends StateObject.Ctor, P extends { }> {
            from?: A | A[],
            params?: PD.For<P> | ((a: StateObject.From<A>, globalCtx: any) => PD.For<P>),
            display?: string | { name: string, description?: string },
            isApplicable?: DefinitionBase<StateObject.From<A>, any, P>['isApplicable']
        }

        export interface Root {
            <A extends StateObject.Ctor, P extends { }>(info: Type<A, P>): Define<StateObject.From<A>, PD.Normalize<P>>
        }

        export interface Define<A extends StateObject, P> {
            <T>(def: DefinitionBase<A, T, P> | DefinitionBase<A, T, P>['run']): StateAction<A, T, P>,
        }

        function root(info: Type<any, any>): Define<any, any> {
            return def => create({
                from: info.from instanceof Array
                    ? info.from
                    : !!info.from ? [info.from] : [],
                display: typeof info.display === 'string'
                    ? { name: info.display }
                    : !!info.display
                        ? info.display
                        : { name: 'Unnamed State Action' },
                params: typeof info.params === 'object'
                    ? () => info.params as any
                    : !!info.params
                        ? info.params as any
                        : void 0,
                isApplicable: info.isApplicable,
                ...(typeof def === 'function'
                    ? { run: def }
                    : def)
            });
        }

        export const build: Root = (info: any) => root(info);
    }

    export const build = Builder.build;
}