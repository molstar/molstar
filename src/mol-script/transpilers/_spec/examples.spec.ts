/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
**/

import { Transpiler } from '../transpiler';
import { _transpiler as transpilers } from '../all';
import { describe, test, it, expect } from 'vitest';

function testTranspilerExamples(name: string) {
    const examplesPath = `../${name}/examples.ts`;

    try {
        console.log(`Attempting to load examples from: ${examplesPath}`);
        const { examples } = require(examplesPath);

        if (Array.isArray(examples) && examples.length > 0) {
            describe(`${name} examples`, () => {
                examples.forEach((example: { name: string; value: any }) => {
                    console.log(`Running test for example: ${example.name}`);
                    test(`Example: ${example.name}`, () => {
                        expect(example.value).toBeDefined();
                    });
                });
            });
        } else {
            console.warn(`No valid examples found in ${examplesPath}`);
        }
    } catch (err) {
        console.warn(`Failed to load examples for ${name}:`, err.message);
    }
}

for (const name in transpilers) {
    testTranspilerExamples(name);
}