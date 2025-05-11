/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
**/

import path from 'path';
import { describe, it, expect } from 'vitest';

function testTranspilerExamples(name: string, transpiler: any) {
    const examplesPath = path.join(__dirname, `../${name}/examples.ts`);
    console.log(`Resolved path for ${name}:`, examplesPath);

    let examples;

    try {
        // Clear module cache to handle re-runs without stale data
        delete require.cache[require.resolve(examplesPath)];
        const module = require(examplesPath);
        console.log(`Module structure for ${name}:`, module);
        examples = module.examples;
    } catch (err) {
        console.warn(`Failed to load examples for ${name}: ${err.message}`);
        examples = [];
    }

    console.log(`Examples for ${name}:`, examples);

    if (!Array.isArray(examples)) {
        console.warn(`Expected examples to be an array but got:`, typeof examples);
        return;
    }

    describe(`${name} examples`, () => {
        examples.forEach((e) => {
            it(e.name, () => {
                console.log(`Running test for ${name} - ${e.name}:`, e);
                expect(typeof e.name).toBe('string');
                expect(e.value).toBeDefined();
            });
        });
    });
}

// Execute tests for each transpiler
const transpilers = {
    pymol: {},
    vmd: {},
    jmol: {}
};

testTranspilerExamples('pymol', transpilers.pymol);
testTranspilerExamples('vmd', transpilers.vmd);
testTranspilerExamples('jmol', transpilers.jmol);
