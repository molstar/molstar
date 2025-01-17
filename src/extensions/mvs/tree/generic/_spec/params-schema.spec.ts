/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { OptionalField, RequiredField, bool, float, int, str } from '../field-schema';
import { SimpleParamsSchema, UnionParamsSchema, paramsValidationIssues } from '../params-schema';


const simpleSchema = SimpleParamsSchema({
    name: OptionalField(str, 'Anonymous', 'Testing optional field name'),
    surname: RequiredField(str, 'Testing optional field surname'),
    lunch: RequiredField(bool, 'Testing optional field lunch'),
    age: OptionalField(int, 0, 'Testing optional field age'),
});

describe('validateParams', () => {
    it('validateParams', async () => {
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(simpleSchema, {}, { noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', age: 29 }, { noExtra: true })).toBeTruthy(); // missing `lunch`
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true, married: false }, { noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true, married: false })).toBeUndefined(); // extra param `married`
    });
});


describe('validateFullParams', () => {
    it('validateFullParams', async () => {
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(simpleSchema, { surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 29 }, { requireAll: true, noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(simpleSchema, {}, { requireAll: true, noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { requireAll: true, noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues(simpleSchema, { name: 'John', surname: 'Doe', lunch: true, age: 29, married: true }, { requireAll: true, noExtra: false })).toBeUndefined(); // extra param `married`
    });
});


const unionSchema = UnionParamsSchema(
    'kind',
    'Description for "kind"',
    {
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
    },
);

describe('validateUnionParams', () => {
    it('validateUnionParams', async () => {
        expect(paramsValidationIssues(unionSchema, { surname: 'Doe', lunch: true }, { noExtra: true })).toBeTruthy(); // missing discriminator param `kind`

        expect(paramsValidationIssues(unionSchema, { kind: 'person', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true, age: 29 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'person' }, { noExtra: true })).toBeTruthy();
        expect(paramsValidationIssues(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', age: 29 }, { noExtra: true })).toBeTruthy(); // missing `lunch`
        expect(paramsValidationIssues(unionSchema, { kind: 'person', name: 'John', surname: 'Doe', lunch: true, age: 'old' }, { noExtra: true })).toBeTruthy(); // wrong type of `age`
        expect(paramsValidationIssues(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, married: false }, { noExtra: true })).toBeTruthy(); // extra param `married`
        expect(paramsValidationIssues(unionSchema, { kind: 'person', surname: 'Doe', lunch: true, married: false })).toBeUndefined(); // extra param `married`

        expect(paramsValidationIssues(unionSchema, { kind: 'object', weight: 42, color: 'black' }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'object', weight: 42 }, { noExtra: true })).toBeUndefined();
        expect(paramsValidationIssues(unionSchema, { kind: 'object', color: 'black' }, { noExtra: true })).toBeTruthy(); // missing param `weight`
        expect(paramsValidationIssues(unionSchema, { kind: 'object', weight: 42, name: 'John' }, { noExtra: true })).toBeTruthy(); // extra param `name`

        expect(paramsValidationIssues(unionSchema, { kind: 'spanish_inquisition' }, { noExtra: true })).toBeTruthy(); // unexpected value for discriminator param `kind`
    });
});
