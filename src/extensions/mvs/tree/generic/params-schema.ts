/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isPlainObject, mapObjectMap, omitObjectKeys } from '../../../../mol-util/object';
import { Field, fieldValidationIssues, OptionalField, RequiredField, ValueFor } from './field-schema';


type Fields = { [key in string]: Field };

/** Type of `ParamsSchema` where all fields are completely independent */
export interface SimpleParamsSchema<TFields extends Fields = Fields> {
    type: 'simple',
    /** Parameter fields */
    fields: TFields,
}
export function SimpleParamsSchema<TFields extends Fields>(fields: TFields): SimpleParamsSchema<TFields> {
    return { type: 'simple', fields };
}

type ValuesForFields<F extends Fields> =
    { [key in keyof F as (F[key] extends RequiredField<any> ? key : never)]: ValueFor<F[key]> }
    & { [key in keyof F as (F[key] extends OptionalField<any> ? key : never)]?: ValueFor<F[key]> };

type ValuesForSimpleParamsSchema<TSchema extends SimpleParamsSchema> = ValuesForFields<TSchema['fields']>;

type AllRequiredFields<F extends Fields>
    = { [key in keyof F]: F[key] extends Field<infer V> ? RequiredField<V> : never };

type AllRequiredSimple<TSchema extends SimpleParamsSchema> = SimpleParamsSchema<AllRequiredFields<TSchema['fields']>>;


type Cases = { [case_ in string]: SimpleParamsSchema };
// Tried to have this recursive ({ [case_ in string]: ParamsSchema }) but ran into "ts(2589) Type instantiation is excessively deep..."

/** Type of `ParamsSchema` where one field (discriminator) determines what other fields are allowed (i.e. discriminated union type) */
export interface UnionParamsSchema<TDiscriminator extends string = string, TCases extends Cases = Cases> {
    type: 'union',
    /** Name of parameter field that determines the rest (allowed values are defined by keys of `cases`) */
    discriminator: TDiscriminator,
    /** Description for the discriminator parameter field */
    discriminatorDescription: string,
    /** `ParamsSchema` for the rest, for each case of discriminator value */
    cases: TCases,
}
export function UnionParamsSchema<TDiscriminator extends string, TCases extends Cases>(discriminator: TDiscriminator, discriminatorDescription: string, cases: TCases): UnionParamsSchema<TDiscriminator, TCases> {
    return { type: 'union', discriminator, discriminatorDescription, cases };
}

type ValuesForUnionParamsSchema<TSchema extends UnionParamsSchema, TCase extends keyof TSchema['cases'] = keyof TSchema['cases']>
    = TCase extends keyof TSchema['cases'] ? { [discriminator in TSchema['discriminator']]: TCase } & ValuesFor<TSchema['cases'][TCase]> : never;
// `extends` clause seems superfluous here, but is needed to properly create discriminated union type

type AllRequiredUnion<TSchema extends UnionParamsSchema>
    = UnionParamsSchema<TSchema['discriminator'], { [case_ in keyof TSchema['cases']]: AllRequired<TSchema['cases'][case_]> }>;


/** Schema for "params", i.e. a flat collection of key-value pairs */
export type ParamsSchema = SimpleParamsSchema | UnionParamsSchema;

/** Type of values for a params schema (optional fields can be missing) */
export type ValuesFor<P extends ParamsSchema>
    = P extends SimpleParamsSchema ? ValuesForSimpleParamsSchema<P> : P extends UnionParamsSchema ? ValuesForUnionParamsSchema<P> : never;

/** Variation of a params schema where all fields are required */
export type AllRequired<P extends ParamsSchema>
    = P extends SimpleParamsSchema ? AllRequiredSimple<P> : P extends UnionParamsSchema ? AllRequiredUnion<P> : never;

function AllRequiredSimple<TSchema extends SimpleParamsSchema>(schema: TSchema): AllRequired<TSchema> {
    const newFields = mapObjectMap(schema.fields, field => RequiredField(field.type, field.description));
    return SimpleParamsSchema(newFields) as AllRequired<TSchema>;
}
function AllRequiredUnion<TSchema extends UnionParamsSchema>(schema: TSchema): AllRequired<TSchema> {
    const newCases = mapObjectMap(schema.cases, AllRequired);
    return UnionParamsSchema(schema.discriminator, schema.discriminatorDescription, newCases) as AllRequired<TSchema>;
}
export function AllRequired<TSchema extends ParamsSchema>(schema: TSchema): AllRequired<TSchema> {
    if (schema.type === 'simple') {
        return AllRequiredSimple(schema) as AllRequired<TSchema>;
    } else {
        return AllRequiredUnion(schema) as AllRequired<TSchema>;
    }
}

/** Type of full values for a params schema, i.e. including all optional fields */
export type FullValuesFor<P extends ParamsSchema> = ValuesFor<AllRequired<P>>;


interface ValidationOptions {
    /** Check that all parameters (including optional) have a value provided. */
    requireAll?: boolean,
    /** Check there are extra parameters other that those defined in the schema. */
    noExtra?: boolean,
}

/** Return `undefined` if `values` contains correct value types for `schema`,
 * return description of validation issues, if `values` have wrong type.
 * If `options.requireAll`, all parameters (including optional) must have a value provided.
 * If `options.noExtra` is true, presence of any extra parameters is treated as an issue. */
export function paramsValidationIssues<P extends ParamsSchema>(schema: P, values: { [k: string]: any }, options: ValidationOptions = {}): string[] | undefined {
    if (!isPlainObject(values)) return [`Parameters must be an object, not ${values}`];

    if (schema.type === 'simple') {
        return simpleParamsValidationIssue(schema, values, options);
    } else {
        return unionParamsValidationIssues(schema, values, options);
    }
}

function simpleParamsValidationIssue(schema: SimpleParamsSchema, values: { [k: string]: any }, options: ValidationOptions): string[] | undefined {
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

function unionParamsValidationIssues(schema: UnionParamsSchema, values: { [k: string]: any }, options: ValidationOptions): string[] | undefined {
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
    const issues = paramsValidationIssues(subschema, omitObjectKeys(values, [schema.discriminator]), options);
    if (issues) {
        issues.unshift(`(case "${schema.discriminator}": "${case_}")`);
        return issues.map(s => '  ' + s);
    }
    return undefined;
}

/** Add default parameter values to `values` based on a parameter schema (only for optional parameters) */
export function addParamDefaults<P extends ParamsSchema>(schema: P, values: ValuesFor<P>): FullValuesFor<P> {
    if (schema.type === 'simple') {
        return addSimpleParamsDefaults(schema, values);
    } else {
        return addUnionParamsDefaults(schema, values);
    }
}

function addSimpleParamsDefaults(schema: SimpleParamsSchema, values: any): any {
    const out = { ...values };
    for (const key in schema.fields) {
        const field = schema.fields[key];
        if (!field.required && out[key] === undefined) {
            out[key] = field.default;
        }
    }
    return out;
}

function addUnionParamsDefaults(schema: UnionParamsSchema, values: any): any {
    const case_ = values[schema.discriminator];
    const subschema = schema.cases[case_];
    return addParamDefaults(subschema, values);
}
