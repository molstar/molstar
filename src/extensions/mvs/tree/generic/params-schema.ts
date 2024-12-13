/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { PathReporter } from 'io-ts/PathReporter';
import { onelinerJsonString } from '../../../../mol-util/json';
import { isPlainObject, mapObjectMap, omitObjectKeys } from '../../../../mol-util/object';


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
/** Mapping between two types */
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
    default: V, // TODO enforce default null for nullable types
}
export function OptionalField<V extends AllowedValueTypes>(type: iots.Type<V>, defaultValue: V, description: string): OptionalField<V> {
    return { type, required: false, description, default: defaultValue };
}

/** Schema for one field in params (i.e. a value in a top-level key-value pair) */
type Field<V extends AllowedValueTypes = any> = RequiredField<V> | OptionalField<V>;


/** Type of valid value for field of type `F` (never includes `undefined`, even if field is optional) */
export type ValueFor<F extends Field | iots.Any> = F extends Field<infer V> ? V : F extends iots.Any ? iots.TypeOf<F> : never

/** Type of valid default value for field of type `F` (if the field's type allows `null`, the default must be `null`) */
export type DefaultFor<F extends Field> = F extends Field<infer V> ? (null extends V ? null : V) : never

/** Return `undefined` if `value` has correct type for `field`, regardsless of if required or optional.
 * Return description of validation issues, if `value` has wrong type. */
export function fieldValidationIssues<F extends Field, V>(field: F, value: V): string[] | undefined {
    const validation = field.type.decode(value);
    if (validation._tag === 'Right') {
        return undefined;
    } else {
        return PathReporter.report(validation);
    }
}


/** Flat collection of key-value pairs with field schemas */
type Fields = { [key in string]: Field }

export interface SimpleParamsSchema<TFields extends Fields = Fields> {
    _type_: 'simple',
    fields: TFields,
}
export function SimpleParamsSchema<TFields extends Fields>(fields: TFields): SimpleParamsSchema<TFields> {
    return { _type_: 'simple', fields };
}

/** Type of values for an object of fields */
type ValuesForFields<F extends Fields> =
    { [key in keyof F as (F[key] extends RequiredField<any> ? key : never)]: ValueFor<F[key]> }
    & { [key in keyof F as (F[key] extends OptionalField<any> ? key : never)]?: ValueFor<F[key]> };

type ValuesForSimpleParamsSchema<TSchema extends SimpleParamsSchema> = ValuesForFields<TSchema['fields']>;
type AllRequiredSimple<TSchema extends SimpleParamsSchema> = SimpleParamsSchema<AllRequiredFields<TSchema['fields']>>;




type Cases = { [c in string]: SimpleParamsSchema }

export interface UnionParamsSchema<TDiscriminator extends string = string, TCases extends Cases = Cases> {
    _type_: 'union',
    discriminator: TDiscriminator,
    cases: TCases,
}
export function UnionParamsSchema<TDiscriminator extends string, TCases extends Cases>(discriminator: TDiscriminator, cases: TCases): UnionParamsSchema<TDiscriminator, TCases> {
    return { _type_: 'union', discriminator, cases };
}
type ValuesForUnionParamsSchema<TSchema extends UnionParamsSchema, TCase extends keyof TSchema['cases'] = keyof TSchema['cases']>
    = TCase extends keyof TSchema['cases'] ? { [disc in TSchema['discriminator']]: TCase } & ValuesFor_new<TSchema['cases'][TCase]> : never;
// `extends` clause seems superfluous here, but is needed to properly create discriminated union type

type AllRequiredUnion<TSchema extends UnionParamsSchema>
    = UnionParamsSchema<TSchema['discriminator'], { [c in keyof TSchema['cases']]: AllRequired_new<TSchema['cases'][c]> }>;


/** Schema for "params", i.e. a flat collection of key-value pairs */
export type ParamsSchema_new = SimpleParamsSchema | UnionParamsSchema;

/** Type of values for a params schema (optional fields can be missing) */
export type ValuesFor_new<T extends ParamsSchema_new>
    = T extends SimpleParamsSchema ? ValuesForSimpleParamsSchema<T> : T extends UnionParamsSchema ? ValuesForUnionParamsSchema<T> : never;
export type AllRequired_new<T extends ParamsSchema_new>
    = T extends SimpleParamsSchema ? AllRequiredSimple<T> : T extends UnionParamsSchema ? AllRequiredUnion<T> : never;

function AllRequired_simple<TSchema extends SimpleParamsSchema>(schema: TSchema): AllRequired_new<TSchema> {
    const newFields = mapObjectMap(schema.fields, field => RequiredField(field.type, field.description));
    return SimpleParamsSchema(newFields) as AllRequired_new<TSchema>;
}

function AllRequired_union<TSchema extends UnionParamsSchema>(schema: TSchema): AllRequired_new<TSchema> {
    const newCases = mapObjectMap(schema.cases, c => AllRequired_new(c));
    return UnionParamsSchema(schema.discriminator, newCases) as AllRequired_new<TSchema>;
}

/** Variation of a params schema where all fields are required */
export function AllRequired_new<TSchema extends ParamsSchema_new>(schema: TSchema): AllRequired_new<TSchema> {
    if (schema._type_ === 'simple') {
        return AllRequired_simple(schema) as AllRequired_new<TSchema>;
    } else {
        return AllRequired_union(schema) as AllRequired_new<TSchema>;
    }
}


/** Make all fields are required */
type AllRequiredFields<TParamsSchema extends Fields>
    = { [key in keyof TParamsSchema]: TParamsSchema[key] extends Field<infer V> ? RequiredField<V> : never }
function AllRequiredFields<TParamsSchema extends Fields>(paramsSchema: TParamsSchema): AllRequiredFields<TParamsSchema> {
    return mapObjectMap(paramsSchema, field => RequiredField(field.type, field.description)) as AllRequiredFields<TParamsSchema>;
}

/** Type of full values for a params schema, i.e. including all optional fields */
export type FullValuesFor_new<P extends ParamsSchema_new> = ValuesFor_new<AllRequired_new<P>>;


interface ValidationOptions {
    /** Check that all parameters (including optional) have a value provided. */
    requireAll?: boolean,
    /** Check there are extra parameters other that those defined in the schema. */
    noExtra?: boolean,
}


function paramsValidationIssues_simple(schema: SimpleParamsSchema, values: { [k: string]: any }, options: ValidationOptions): string[] | undefined {
    for (const key in schema.fields) {
        const fieldSchema = schema.fields[key];

        if (Object.hasOwn(values, key)) {
            const value = values[key];
            const issues = fieldValidationIssues(fieldSchema, value);
            if (issues) return [`Invalid value for parameter "${key}":`, ...issues.map(s => '  ' + s)];
        } else {
            if (fieldSchema.required) return [`Missing required parameter "${key}".`];
            if (options.requireAll) return [`Missing optional parameter "${key}".`];
        }
    }
    if (options.noExtra) {
        for (const key in values) {
            if (!Object.hasOwn(schema.fields, key)) return [`Unknown parameter "${key}".`];
        }
    }
    return undefined;
}

function paramsValidationIssues_union(schema: UnionParamsSchema, values: { [k: string]: any }, options: ValidationOptions): string[] | undefined {
    if (!Object.hasOwn(values, schema.discriminator)) {
        return [`Missing required parameter "${schema.discriminator}".`];
    }
    const case_ = values[schema.discriminator];
    const subschema = schema.cases[case_];
    if (subschema === undefined) {
        const allowedCases = Object.keys(schema.cases).map(x => `"${x}"`).join(' | ');
        return [
            `Invalid value for parameter "${schema.discriminator}":`,
            `"${case_}" is not a valid value for literal type (${allowedCases})`,
        ];
    }
    const issues = paramsValidationIssues_new(subschema, omitObjectKeys(values, [schema.discriminator]), options);
    if (issues) {
        issues.unshift(`(case "${schema.discriminator}": "${case_}")`);
        return issues.map(s => '  ' + s);
    }
    return undefined;
}

/** Return `undefined` if `values` contains correct value types for `schema`,
 * return description of validation issues, if `values` have wrong type.
 * If `options.requireAll`, all parameters (including optional) must have a value provided.
 * If `options.noExtra` is true, presence of any extra parameters is treated as an issue. */
export function paramsValidationIssues_new<P extends ParamsSchema_new>(schema: P, values: { [k: string]: any }, options: ValidationOptions = {}): string[] | undefined {
    if (!isPlainObject(values)) return [`Parameters must be an object, not ${values}`];

    if (schema._type_ === 'simple') {
        return paramsValidationIssues_simple(schema, values, options);
    } else {
        return paramsValidationIssues_union(schema, values, options);
    }
}


function addParamDefaults_simple<P extends SimpleParamsSchema>(schema: P, values: any): any {
    const out = { ...values };
    for (const key in schema.fields) {
        const field = schema.fields[key];
        if (!field.required && out[key] === undefined) {
            out[key] = field.default;
        }
    }
    return out;
}

function addParamDefaults_union<P extends UnionParamsSchema>(schema: P, values: any): any {
    const case_ = values[schema.discriminator];
    const subschema = schema.cases[case_];
    return addParamDefaults_new(subschema, values);
}

export function addParamDefaults_new<P extends ParamsSchema_new>(schema: P, values: ValuesFor_new<P>): FullValuesFor_new<P> {
    if (schema._type_ === 'simple') {
        return addParamDefaults_simple(schema, values);
    } else {
        return addParamDefaults_union(schema, values);
    }
}
