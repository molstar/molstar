/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { RequiredField, fieldValidationIssues, float, int, literal, nullable, str, union } from '../field-schema';


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
        const stringOrNumberParam = RequiredField(union(str, float), 'Testing required field stringOrNumberParam');
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
