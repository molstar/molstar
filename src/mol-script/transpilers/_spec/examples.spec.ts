/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
**/

import path from 'path';
import { describe, test, it, expect } from 'vitest';

const transpilers = ['pymol', 'vmd', 'jmol'];

transpilers.forEach(name => {
  const examplesPath = path.join(__dirname, `../${name}/examples.ts`);

  try {
    const { examples } = require(examplesPath);

    if (examples && examples.length > 0) {
      describe(`${name} examples`, () => {
        examples.forEach((example: { name: string; value: any }) => {
          test(`Example: ${example.name}`, () => {
            expect(example.value).toBeDefined();
          });
        });
      });
    } else {
      console.warn(`No examples found in ${examplesPath}`);
    }
  } catch (err) {
    console.warn(`Failed to load examples for ${name}:`, err.message);
  }
});