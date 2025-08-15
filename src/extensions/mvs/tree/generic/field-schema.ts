/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { onelinerJsonString } from '../../../../mol-util/json';

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
/** Type definition for a dictionary/mapping/record, e.g. `dict(str, float)` means type `{ [K in string]: number }` */
export const dict = iots.record;

/** Type definition used to create objects, e.g. `object({ name: str, age: float }, { address: str })` means type `{ name: string, age: number, address?: string }` */
export function object<P extends iots.Props, Q extends iots.Props>(props: P, optionalProps: undefined, name?: string): iots.TypeC<P>;
export function object<P extends iots.Props, Q extends iots.Props>(props: P, optionalProps: Q, name?: string): iots.IntersectionC<[iots.TypeC<P>, iots.PartialC<Q>]>;
export function object<P extends iots.Props, Q extends iots.Props>(props: P, optionalProps?: Q, name?: string) {
    if (!optionalProps) {
        return iots.type(props, name);
    }

    if (name === undefined) {
        const nameChunks = [];
        for (const key in props) {
            nameChunks.push(`${key}: ${props[key].name}`);
        }
        for (const key in optionalProps) {
            nameChunks.push(`${key}?: ${optionalProps[key].name}`);
        }
        name = `{ ${nameChunks.join(', ')} }`;
    }
    return iots.intersection([iots.type(props), iots.partial(optionalProps)], name);
}

/** Type definition used to create partial objects, e.g. `partial({ name: str, age: float })` means type `{ name?: string, age?: number }` */
export function partial<P extends iots.Props>(props: P, name?: string) {
    if (name === undefined) {
        const nameChunks = [];
        for (const key in props) {
            nameChunks.push(`${key}?: ${props[key].name}`);
        }
        name = `{ ${nameChunks.join(', ')} }`;
    }
    return iots.partial(props, name);
}

/** Type definition for union types, e.g. `union(str, int)` means string or integer */
export function union<T1 extends iots.Mixed, T2 extends iots.Mixed, TOthers extends iots.Mixed[]>(first: T1, second: T2, ...others: TOthers): iots.UnionC<[T1, T2, ...TOthers]> {
    const baseTypes: iots.Mixed[] = [];
    for (const type of [first, second, ...others]) {
        if (type instanceof iots.UnionType) {
            baseTypes.push(...type.types);
        } else {
            baseTypes.push(type);
        }
    }
    return iots.union(baseTypes as any);
}

/** Type definition for nullable types, e.g. `nullable(str)` means string or `null`  */
export function nullable<V>(type: iots.Type<V>): iots.Type<V | null> {
    return union(type, iots.null);
}

/** Type definition for literal types, e.g. `literal('red', 'green', 'blue')` means 'red' or 'green' or 'blue'  */
export function literal<V extends string | number | boolean>(...values: V[]) {
    if (values.length === 0) {
        throw new Error(`literal type must have at least one value`);
    }
    const typeName = values.length === 1 ? onelinerJsonString(values[0]) : `(${values.map(v => onelinerJsonString(v)).join(' | ')})`;
    const valueSet = new Set(values);
    return new iots.Type<V>(
        typeName,
        ((value: any) => valueSet.has(value)) as any,
        (value, ctx) => valueSet.has(value as any) ? { _tag: 'Right', right: value as any } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid value for literal type ${typeName}` }] },
        value => value
    );
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
        return reportErrors(validation.left);
    }
}

// Inlining `reportErrors` instead of `import { PathReporter } from 'io-ts/PathReporter'`;
// because it breaks Deno usage.

function reportErrors(errors: iots.Errors): string[] | undefined {
    if (errors.length === 0) return undefined;
    return errors.map(getMessage);
}

function getMessage(e: iots.ValidationError) {
    return e.message !== undefined
        ? e.message
        : `Invalid value ${stringifyError(e.value)} supplied to ${getContextPath(e.context)}`;
}

function getContextPath(context: iots.ValidationError['context']) {
    return context.map(a => `${a.key}: ${a.type.name}`).join('/');
}

function getFunctionName(f: Function & { displayName?: string }) {
    return f.displayName || f.name || `<function ${f.length}>`;
}

function stringifyError(v: any) {
    if (typeof v === 'function') {
        return getFunctionName(v);
    }
    if (typeof v === 'number' && !isFinite(v)) {
        if (isNaN(v)) {
            return 'NaN';
        }
        return v > 0 ? 'Infinity' : '-Infinity';
    }
    return JSON.stringify(v);
}