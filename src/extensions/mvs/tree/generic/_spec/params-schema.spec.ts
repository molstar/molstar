/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as iots from 'io-ts';

import { fieldValidationIssues, RequiredField, literal, nullable, paramsValidationIssues, OptionalField } from '../params-schema';


describe('fieldValidationIssues', () => {
    it('fieldValidationIssues string', async () => {
        const stringField = RequiredField(iots.string);
        expect(fieldValidationIssues(stringField, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringField, '')).toBeUndefined();
        expect(fieldValidationIssues(stringField, 5)).toBeTruthy();
        expect(fieldValidationIssues(stringField, null)).toBeTruthy();
        expect(fieldValidationIssues(stringField, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues string choice', async () => {
        const colorParam = RequiredField(literal('red', 'green', 'blue', 'yellow'));
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
        const numberParam = RequiredField(literal(1, 2, 3, 4));
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
        const numberParam = RequiredField(iots.Integer);
        expect(fieldValidationIssues(numberParam, 1)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 0)).toBeUndefined();
        expect(fieldValidationIssues(numberParam, 0.5)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, '1')).toBeTruthy();
        expect(fieldValidationIssues(numberParam, null)).toBeTruthy();
        expect(fieldValidationIssues(numberParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues union', async () => {
        const stringOrNumberParam = RequiredField(iots.union([iots.string, iots.number]));
        expect(fieldValidationIssues(stringOrNumberParam, 1)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, 2)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, '')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNumberParam, true)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNumberParam, null)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNumberParam, undefined)).toBeTruthy();
    });
    it('fieldValidationIssues nullable', async () => {
        const stringOrNullParam = RequiredField(nullable(iots.string));
        expect(fieldValidationIssues(stringOrNullParam, 'hello')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, '')).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, null)).toBeUndefined();
        expect(fieldValidationIssues(stringOrNullParam, 1)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNullParam, true)).toBeTruthy();
        expect(fieldValidationIssues(stringOrNullParam, undefined)).toBeTruthy();
    });
});

const schema = {
    name: OptionalField(iots.string),
    surname: RequiredField(iots.string),
    lunch: RequiredField(iots.boolean),
    age: OptionalField(iots.number),
};

describe('validateParams', () => {
    it('validateParams', async () => {
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(schema, {}, { noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', age: 29 }, { noExtra: true })).toBeTruthy(); // missing `lunch`
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true, married: false }, { noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true, married: false })).toBeUndefined(); // extra param `married`
    });
});


describe('validateFullParams', () => {
    it('validateFullParams', async () => {
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(schema, { surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(schema, {}, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { requireAll: true, noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues(schema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: false })).toBeUndefined(); // extra param `married`
    });
});

