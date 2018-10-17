/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task';
import { StateObject } from '../model/object';
import { TransformContext } from './context';

export interface Transformer<A extends StateObject, B extends StateObject, P = any> {
    readonly id: Transformer.Id,
    readonly name: string,
    readonly namespace: string,
    readonly description?: string,
    readonly from: StateObject.TypeOf<A>[],
    readonly to: StateObject.TypeOf<B>[],

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

    /** Specify default control descriptors for the parameters */
    defaultControls?(a: A, context: TransformContext): Transformer.ControlsFor<P>,

    /** Check the parameters and return a list of errors if the are not valid. */
    validateParams?(a: A, params: P, context: TransformContext): string[] | undefined,

    /** Test if the transform can be applied to a given node */
    isApplicable?(a: A, context: TransformContext): boolean,

    /** By default, returns true */
    isSerializable?(params: P): { isSerializable: true } | { isSerializable: false; reason: string },

    /** Custom conversion to and from JSON */
    customSerialization?: { toJSON(params: P, obj?: B): any, fromJSON(data: any): P }
}

export namespace Transformer {
    export type Id = string & { '@type': 'transformer-id' }
    export type Params<T extends Transformer<any, any, any>> = T extends Transformer<any, any, infer P> ? P : unknown;
    export type ControlsFor<Props> = { [P in keyof Props]?: any }
}