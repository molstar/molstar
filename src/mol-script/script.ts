/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { transpileMolScript } from './script/mol-script/symbols';
import { parseMolScript } from './language/parser';
import { parse } from './transpile';
import { Expression } from './language/expression';
import { StructureElement, QueryContext, StructureSelection, Structure, QueryFn, QueryContextOptions } from '../mol-model/structure';
import { compile } from './runtime/query/compiler';
import { MolScriptBuilder } from './language/builder';
import { assertUnreachable } from '../mol-util/type-helpers';

export { Script };

interface Script { expression: string, language: Script.Language }

function Script(expression: string, language: Script.Language): Script {
    return { expression, language };
}

namespace Script {
    export const Info = {
        'mol-script': 'Mol-Script',
        'pymol': 'PyMOL',
        'vmd': 'VMD',
        'jmol': 'Jmol',
    };
    export type Language = keyof typeof Info;

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

    export function getStructureSelection(expr: Expression | ((builder: typeof MolScriptBuilder) => Expression), structure: Structure, options?: QueryContextOptions) {
        const e = typeof expr === 'function' ? expr(MolScriptBuilder) : expr;
        const query = compile<StructureSelection>(e);
        return query(new QueryContext(structure, options));
    }
}
