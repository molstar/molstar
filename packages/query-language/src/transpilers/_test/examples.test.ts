/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
 **/

import { describe, it } from 'vitest';
import type { Transpiler } from '../transpiler.js';
import { transpiler as jmol } from '../jmol/parser.js';
import { transpiler as pymol } from '../pymol/parser.js';
import { transpiler as vmd } from '../vmd/parser.js';
import { examples as jmolExamples } from '../jmol/examples.js';
import { examples as pymolExamples } from '../pymol/examples.js';
import { examples as vmdExamples } from '../vmd/examples.js';

function testTranspilerExamples(name: string, examples: { name: string; value: string }[], transpiler: Transpiler) {
  describe(`${name} examples`, () => {
    for (const e of examples) {
      it(e.name, () => {
        // check if it transpiles and compiles/typechecks.
        transpiler(e.value);
      });
    }
  });
}

testTranspilerExamples('pymol', pymolExamples, pymol);
testTranspilerExamples('vmd', vmdExamples, vmd);
testTranspilerExamples('jmol', jmolExamples, jmol);
