/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { PathReporter } from 'io-ts/PathReporter';
import { onelinerJsonString } from '../../util/json';


/** All types that can be used in tree node params.
 * Can be extended, this is just to list them all in one place and possibly catch some typing errors */
type AllowedValueTypes = string | number | boolean | null | [number, number, number] | string[] | number[] | {};

/** Type definition for a string  */
export const str = iots.string;
/** Type definition for an integer  */
export const int = iots.Integer;
/** Type definition for a float or integer number  */
export const float = iots.number;
/** Type definition for a boolean  */
export const bool = iots.boolean;
/** Type definition for a tuple, e.g. `tuple([str, int, int])`  */
export const tuple = iots.tuple;
/** Type definition for a list/array, e.g. `list(str)`  */
export const list = iots.array;
/** Type definition for union types, e.g. `union([str, int])` means string or integer  */
export const union = iots.union;
/** Type definition used to create objects */
export const obj = iots.type;
/** Type definition used to create partial objects */
export const partial = iots.partial;

/** Type definition for nullable types, e.g. `nullable(str)` means string or `null`  */
export function nullable<T extends iots.Type<any>>(type: T) {
    return union([type, iots.null]);
}
/** Type definition for literal types, e.g. `literal('red', 'green', 'blue')` means 'red' or 'green' or 'blue'  */
export function literal<V extends string | number | boolean>(...values: V[]) {
    if (values.length === 0) {
        throw new Error(`literal type must have at least one value`);
    }
    const typeName = `(${values.map(v => onelinerJsonString(v)).join(' | ')})`;
    return new iots.Type<V>(
        typeName,
        ((value: any) => values.includes(value)) as any,
        (value, ctx) => values.includes(value as any) ? { _tag: 'Right', right: value as any } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid value for literal type ${typeName}` }] },
        value => value
    );
}
/** Type definition for mapping between two types, e.g. `mapping(str, float)` means type `{ [key in string]: number }` */
export function mapping<A extends iots.Type<any>, B extends iots.Type<any>>(from: A, to: B) {
    return iots.record(from, to);
}


interface FieldBase<V extends AllowedValueTypes = any, R extends boolean = boolean> {
    /** Definition of allowed types for the field */
    type: iots.Type<V>,
    /** If `required===true`, the value must always be defined in molviewspec format (can be `null` if `type` allows it).
     * If `required===false`, the value can be ommitted (meaning that a default should be used).
     * If `type` allows `null`, the default must be `null`. */
    required: R,
    /** Description of what the field value means */
    description: string,
}

/** Schema for param field which must always be provided (has no default value) */
export interface RequiredField<V extends AllowedValueTypes = any> extends FieldBase<V> {
    required: true,
}
export function RequiredField<V extends AllowedValueTypes>(type: iots.Type<V>, description: string): RequiredField<V> {
    return { type, required: true, description };
}

/** Schema for param field which can be dropped (meaning that a default value will be used) */
export interface OptionalField<V extends AllowedValueTypes = any> extends FieldBase<V> {
    required: false,
    /** Default value for optional field.
     * If field type allows `null`, default must be `null` (this is to avoid issues in languages that do not distinguish `null` and `undefined`). */
    default: DefaultValue<V>,
}
export function OptionalField<V extends AllowedValueTypes>(type: iots.Type<V>, defaultValue: DefaultValue<V>, description: string): OptionalField<V> {
    return { type, required: false, description, default: defaultValue };
}

/** Schema for one field in params (i.e. a value in a top-level key-value pair) */
export type Field<V extends AllowedValueTypes = any> = RequiredField<V> | OptionalField<V>;

/** Type of valid default value for value type `V` (if the type allows `null`, the default must be `null`) */
type DefaultValue<V extends AllowedValueTypes> = null extends V ? null : V;

/** Type of valid value for field of type `F` (never includes `undefined`, even if field is optional) */
export type ValueFor<F extends Field | iots.Any> = F extends Field<infer V> ? V : F extends iots.Any ? iots.TypeOf<F> : never;

/** Return `undefined` if `value` has correct type for `field`, regardsless of if required or optional.
 * Return description of validation issues, if `value` has wrong type. */
export function fieldValidationIssues<F extends Field>(field: F, value: any): string[] | undefined {
    const validation = field.type.decode(value);
    if (validation._tag === 'Right') {
        return undefined;
    } else {
        return PathReporter.report(validation);
    }
}
