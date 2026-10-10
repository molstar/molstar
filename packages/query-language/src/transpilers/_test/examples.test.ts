/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * Adapted from MolQL project
 **/

import type { Transpiler } from '../transpiler.js';
import { transpiler as jmol } from '../jmol/parser.js';
import { transpiler as pymol } from '../pymol/parser.js';
import { transpiler as vmd } from '../vmd/parser.js';

function testTranspilerExamples(name: string, transpiler: Transpiler) {
  describe(`${name} examples`, () => {
    const examples = require(`../${name}/examples`).examples;
    //        console.log(examples);
    for (const e of examples) {
      it(e.name, () => {
        // check if it transpiles and compiles/typechecks.
        transpiler(e.value);
      });
    }
  });
}

testTranspilerExamples('pymol', pymol);
testTranspilerExamples('vmd', vmd);
testTranspilerExamples('jmol', jmol);
