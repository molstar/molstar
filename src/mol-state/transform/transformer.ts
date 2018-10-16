/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task';
import { EventDispatcher } from '../context/event';
import { ModelNode } from '../model/node';
import { ModelTree } from '../model/tree';

export interface Transformer<A extends ModelNode, B extends ModelNode, P = any> {
    readonly id: Transformer.Id,
    readonly definition: Transformer.Definition<A, B, P>
}

export namespace Transformer {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends Transformer<any, any, any>> = T extends Transformer<any, any, infer P> ? P : unknown;

    export interface Definition<A extends ModelNode, B extends ModelNode, P> {
        readonly name: string,
        readonly namespace: string,
        readonly description?: string,
        readonly from: ModelNode.TypeOf<A>[],
        readonly to: ModelNode.TypeOf<B>[]

        /**
         * Apply the actual transformation. It must be pure (i.e. with no side effects).
         * Returns a task that produces the result of the result directly.
         */
        apply(a: A, params: P, context: TransformContext): Task<B> | B,

        /**
         * Attempts to update the entity in a non-destructive way.
         * For example changing a color scheme of a visual does not require computing new geometry.
         * Return/resolve to undefined if the update is not possible.
         *
         * The ability to resolve the task to undefined is present for "async updates" (i.e. containing an ajax call).
         */
        update?(a: A, b: B, newParams: P, context: TransformContext): Task<B | undefined> | B | undefined,

        /** Check the parameters and return a list of errors if the are not valid. */
        defaultParams?(a: A, context: TransformContext): P,

        /**  */
        defaultControls?(a: A, context: TransformContext): ControlsFor<P>,

        /** Check the parameters and return a list of errors if the are not valid. */
        validateParams?(a: A, params: P, context: TransformContext): string[] | undefined,

        /** Test if the transform can be applied to a given node */
        isApplicable?(a: A, context: TransformContext): boolean,

        /** By default, returns true */
        isSerializable?(params: P): { isSerializable: true } | { isSerializable: false; reason: string },
    }

    export type ControlsFor<Props> = { [P in keyof Props]: any }

    /** A tree context constructed dynamically duing application of transforms. */
    export interface TransformContext {
        /** An event dispatcher for executing child tasks. */
        dispatcher: EventDispatcher,

        globalContext: any,
        tree: ModelTree
    }
}