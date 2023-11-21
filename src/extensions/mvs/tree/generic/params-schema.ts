/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as iots from 'io-ts';
import { PathReporter } from 'io-ts/PathReporter';
import { isPlainObject, mapObjectMap } from '../../../../mol-util/object';


/** All types that can be used in tree node params.
 * Can be extended, this is just to list them all in one place and possibly catch some typing errors */
type AllowedValueTypes = string | number | boolean | null | [number, number, number] | string[] | number[] | {}

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
/** Type definition for nullable types, e.g. `nullable(str)` means string or `null`  */
export function nullable<T extends iots.Type<any>>(type: T) {
    return union([type, iots.null]);
}
/** Type definition for literal types, e.g. `literal('red', 'green', 'blue')` means 'red' or 'green' or 'blue'  */
export function literal<V extends string | number | boolean>(v1: V, v2?: V, ...others: V[]) {
    if (v2 === undefined) {
        return iots.literal(v1);
    } else {
        return union([iots.literal(v1), iots.literal(v2), ...others.map(v => iots.literal(v))]);
    }
}


/** Schema for one field in params (i.e. a value in a top-level key-value pair) */
interface Field<V extends AllowedValueTypes = any, R extends boolean = boolean> {
    /** Definition of allowed types for the field */
    type: iots.Type<V>,
    /** If `required===true`, the value must always be defined in molviewspec format (can be `null` if `type` allows it).
     * If `required===false`, the value can be ommitted (meaning that a default should be used).
     * If `type` allows `null`, the default must be `null`. */
    required: R,
    /** Description of what the field value means */
    description?: string,
}
/** Schema for param field which must always be provided (has no default value) */
export interface RequiredField<V extends AllowedValueTypes = any> extends Field<V> {
    required: true,
}
export function RequiredField<V extends AllowedValueTypes>(type: iots.Type<V>, description?: string): RequiredField<V> {
    return { type, required: true, description };
}

/** Schema for param field which can be dropped (meaning that a default value will be used) */
export interface OptionalField<V extends AllowedValueTypes = any> extends Field<V> {
    required: false,
}
export function OptionalField<V extends AllowedValueTypes>(type: iots.Type<V>, description?: string): OptionalField<V> {
    return { type, required: false, description };
}

/** Type of valid value for field of type `F` (never includes `undefined`, even if field is optional) */
export type ValueFor<F extends Field | iots.Any> = F extends Field<infer V> ? V : F extends iots.Any ? iots.TypeOf<F> : never

/** Type of valid default value for field of type `F` (if the field's type allows `null`, the default must be `null`) */
export type DefaultFor<F extends Field> = F extends Field<infer V> ? (null extends V ? null : V) : never

/** Return `undefined` if `value` has correct type for `field`, regardsless of if required or optional.
 * Return description of validation issues, if `value` has wrong type. */
export function fieldValidationIssues<F extends Field, V>(field: F, value: V): V extends ValueFor<F> ? undefined : string[] {
    const validation = field.type.decode(value);
    if (validation._tag === 'Right') {
        return undefined as any;
    } else {
        return PathReporter.report(validation) as any;
    }
}


/** Schema for "params", i.e. a flat collection of key-value pairs */
export type ParamsSchema<TKey extends string = string> = { [key in TKey]: Field }

/** Variation of a params schema where all fields are required */
export type AllRequired<TParamsSchema extends ParamsSchema> = { [key in keyof TParamsSchema]: TParamsSchema[key] extends Field<infer V> ? RequiredField<V> : never }
export function AllRequired<TParamsSchema extends ParamsSchema>(paramsSchema: TParamsSchema): AllRequired<TParamsSchema> {
    return mapObjectMap(paramsSchema, field => RequiredField(field.type, field.description)) as AllRequired<TParamsSchema>;
}

/** Type of values for a params schema (optional fields can be missing) */
export type ValuesFor<P extends ParamsSchema> =
    { [key in keyof P as (P[key] extends RequiredField<any> ? key : never)]: ValueFor<P[key]> }
    & { [key in keyof P as (P[key] extends OptionalField<any> ? key : never)]?: ValueFor<P[key]> }

/** Type of full values for a params schema, i.e. including all optional fields */
export type FullValuesFor<P extends ParamsSchema> = { [key in keyof P]: ValueFor<P[key]> }

/** Type of default values for a params schema, i.e. including only optional fields */
export type DefaultsFor<P extends ParamsSchema> = { [key in keyof P as (P[key] extends Field<any, false> ? key : never)]: ValueFor<P[key]> }


/** Return `undefined` if `values` contains correct value types for `schema`,
 * return description of validation issues, if `values` have wrong type.
 * If `options.requireAll`, all parameters (including optional) must have a value provided.
 * If `options.noExtra` is true, presence of any extra parameters is treated as an issue.
 */
export function paramsValidationIssues<P extends ParamsSchema, V extends { [k: string]: any }>(schema: P, values: V, options: { requireAll?: boolean, noExtra?: boolean } = {}): string[] | undefined {
    if (!isPlainObject(values)) return [`Parameters must be an object, not ${values}`];
    for (const key in schema) {
        const paramDef = schema[key];
        if (Object.hasOwn(values, key)) {
            const value = values[key];
            const issues = fieldValidationIssues(paramDef, value);
            if (issues) return [`Invalid type for parameter "${key}":`, ...issues.map(s => '  ' + s)];
        } else {
            if (paramDef.required) return [`Missing required parameter "${key}".`];
            if (options.requireAll) return [`Missing optional parameter "${key}".`];
        }
    }
    if (options.noExtra) {
        for (const key in values) {
            if (!Object.hasOwn(schema, key)) return [`Unknown parameter "${key}".`];
        }
    }
    return undefined;
}
