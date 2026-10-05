/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { transpileMolScript } from './script/mol-script/symbols.js';
import { parseMolScript } from './language/parser.js';
import { parse } from './transpile.js';
import type { Expression } from './language/expression.js';
import {
  type StructureElement,
  QueryContext,
  StructureSelection,
  type Structure,
  type QueryFn,
  type QueryContextOptions,
} from '@molstar/model/model/structure';
import { compile } from './runtime/query/compiler.js';
import { MolScriptBuilder } from './language/builder.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import type { Script } from '@molstar/core/util/script';

export { ScriptImpl as Script };

type ScriptImpl = Script;

function ScriptImpl(expression: string, language: Script['language']): Script {
  return { expression, language };
}

namespace ScriptImpl {
  export const Info: { [k in Script['language']]: string } = {
    'mol-script': 'Mol-Script',
    pymol: 'PyMOL',
    vmd: 'VMD',
    jmol: 'Jmol',
  };
  export type Language = Script['language'];

  export function is(x: any): x is Script {
    return !!x && typeof (x as Script).expression === 'string' && !!(x as Script).language;
  }

  export function areEqual(a: Script, b: Script) {
    return a.language === b.language && a.expression === b.expression;
  }

  export function toExpression(script: Script): Expression {
    switch (script.language) {
      case 'mol-script':
        const parsed = parseMolScript(script.expression);
        if (parsed.length === 0) throw new Error('No query');
        return transpileMolScript(parsed[0]);
      case 'pymol':
      case 'jmol':
      case 'vmd':
        return parse(script.language, script.expression);
      default:
        assertUnreachable(script.language);
    }
  }

  export function toQuery(script: Script): QueryFn<StructureSelection> {
    const expression = toExpression(script);
    return compile<StructureSelection>(expression);
  }

  export function toLoci(script: Script, structure: Structure): StructureElement.Loci {
    const query = toQuery(script);
    const result = query(new QueryContext(structure));
    return StructureSelection.toLociWithSourceUnits(result);
  }

  export function getStructureSelection(
    expr: Expression | ((builder: typeof MolScriptBuilder) => Expression),
    structure: Structure,
    options?: QueryContextOptions,
  ) {
    const e = typeof expr === 'function' ? expr(MolScriptBuilder) : expr;
    const query = compile<StructureSelection>(e);
    return query(new QueryContext(structure, options));
  }
}
