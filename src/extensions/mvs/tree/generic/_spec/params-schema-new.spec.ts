/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as iots from 'io-ts';

import { OptionalField, RequiredField, SimpleParamsSchema, UnionParamsSchema, bool, fieldValidationIssues, float, int, literal, nullable, paramsValidationIssues_new, str, union } from '../params-schema';


describe('fieldValidationIssues', () => {
    it('fieldValidationIssues string', async () => {
        const stringField = RequiredField(str, 'Testing required field stringField');
        expect(fieldValidationIssues(stringField, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringField, '')).toBeUndefined();
        expect(fieldValidationIssues(stringField, 5)).toBeTruthy();
        expect(fieldValidationIssues(stringField, null)).toBeTruthy();
        expect(fieldValidationIssues(stringField, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues string choice', async () => {
        const colorParam = RequiredField(literal('red', 'green', 'blue', 'yellow'), 'Testing required field colorParam');
        expect(fieldValidationIssues(colorParam, 'red')).toBeUndefined();
        expect(fieldValidationIssues(colorParam, 'green')).toBeUndefined();
        expect(fieldValidationIssues(colorParam, 'blue')).toBeUndefined();
        expect(fieldValidationIssues(colorParam, 'yellow')).toBeUndefined();
        expect(fieldValidationIssues(colorParam, 'banana')).toBeTruthy();
        expect(fieldValidationIssues(colorParam, 5)).toBeTruthy();
        expect(fieldValidationIssues(colorParam, null)).toBeTruthy();
        expect(fieldValidationIssues(colorParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues number choice', async () => {
        const numberParam = RequiredField(literal(1, 2, 3, 4), 'Testing required field numberParam');
        expect(fieldValidationIssues(numberParam, 1)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 2)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 3)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 4)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 5)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, '1')).toBeTruthy();
        expect(fieldValidationIssues(numberParam, null)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues int', async () => {
        const numberParam = RequiredField(int, 'Testing required field numberParam');
        expect(fieldValidationIssues(numberParam, 1)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 0)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 0.5)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, '1')).toBeTruthy();
        expect(fieldValidationIssues(numberParam, null)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues union', async () => {
        const stringOrNumberParam = RequiredField(union([str, float]), 'Testing required field stringOrNumberParam');
        expect(fieldValidationIssues(stringOrNumberParam, 1)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, 2)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, '')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, true)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNumberParam, null)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNumberParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues nullable', async () => {
        const stringOrNullParam = RequiredField(nullable(str), 'Testing required field stringOrNullParam');
        expect(fieldValidationIssues(stringOrNullParam, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, '')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, null)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, 1)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNullParam, true)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNullParam, undefined)).toBeTruthy();
    });
});

const schema = SimpleParamsSchema({
    name: OptionalField(str, 'Anonymous', 'Testing optional field name'),
    surname: RequiredField(str, 'Testing optional field surname'),
    lunch: RequiredField(bool, 'Testing optional field lunch'),
    age: OptionalField(int, 0, 'Testing optional field age'),
});

describe('validateParams', () => {
    it('validateParams', async () => {
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(schema, {}, { noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', age: 29 }, { noExtra: true })).toBeTruthy(); // missing `lunch`
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true, married: false }, { noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true, married: false })).toBeUndefined(); // extra param `married`
    });
});


describe('validateFullParams', () => {
    it('validateFullParams', async () => {
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(schema, { surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(schema, {}, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { requireAll: true, noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues_new(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: false })).toBeUndefined(); // extra param `married`
    });
});


const unionSchema = UnionParamsSchema('kind', {
    person: SimpleParamsSchema({
        name: OptionalField(str, 'Anonymous', 'Testing optional field name'),
        surname: RequiredField(str, 'Testing optional field surname'),
        lunch: RequiredField(bool, 'Testing optional field lunch'),
        age: OptionalField(int, 0, 'Testing optional field age'),
    }),
    object: SimpleParamsSchema({
        weight: RequiredField(float, 'Testing optional field weight'),
        color: OptionalField(str, 'colorless', 'Testing optional field color'),
    }),
});

describe('validateUnionParams', () => {
    it('validateUnionParams', async () => {
        expect(paramsValidationIssues_new(unionSchema, { surname: 'Doe', lunch: true }, { noExtra: true })).toBeTruthy(); // missing discriminator param `kind`

        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person' }, { noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', age: 29 }, { noExtra: true })).toBeTruthy(); // missing `lunch`
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, married: false }, { noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues_new(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, married: false })).toBeUndefined(); // extra param `married`

        expect(paramsValidationIssues_new(unionSchema, { kind: 'object', weight: 42, color: 'black' }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'object', weight: 42 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues_new(unionSchema, { kind: 'object', color: 'black' }, { noExtra: true })).toBeTruthy(); // missing param `weight`
        expect(paramsValidationIssues_new(unionSchema, { kind: 'object', weight: 42, name: 'John' }, { noExtra: true })).toBeTruthy(); // extra param `name`

        expect(paramsValidationIssues_new(unionSchema, { kind: 'spanish_inquisition' }, { noExtra: true })).toBeTruthy(); // unexpected value for discriminator param `kind`
    });
});
