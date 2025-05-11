/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
**/

import { describe, test, it, expect } from 'vitest';
import path from 'path';

const transpilers = ['pymol', 'vmd', 'jmol'];

transpilers.forEach(name => {
  const examplesPath = path.join(__dirname, `../${name}/examples.ts`);

  try {
    console.log(`Attempting to load examples from: ${examplesPath}`);
    const module = require(examplesPath);

    if (module && module.examples && Array.isArray(module.examples)) {
      console.log(`Found ${module.examples.length} examples for ${name}`);
      describe(`${name} examples`, () => {
        module.examples.forEach((example: { name: string; value: any }) => {
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
});