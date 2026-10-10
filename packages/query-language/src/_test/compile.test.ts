/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { compileScript } from '../compile.js';
import { expressionValidationIssues } from '../language/validation.js';
import { examples as pymol } from '../transpilers/pymol/examples.js';
import { examples as vmd } from '../transpilers/vmd/examples.js';
import { examples as jmol } from '../transpilers/jmol/examples.js';
import { Script } from '../script.js';

describe('Text-to-MolQL compilation without runtime registration', () => {
  it('translates MolScript aliases', () => {
    const expression = compileScript('mol-script', '(sel.atom.all)');
    expect(expressionValidationIssues(expression)).toBeUndefined();
    expect(JSON.stringify(expression)).toContain('structure-query.generator.all');
  });

  for (const [language, examples] of [
    ['pymol', pymol],
    ['vmd', vmd],
    ['jmol', jmol],
  ] as const) {
    for (const example of examples) {
      it(`${language}: ${example.name}`, () => {
        const expression = compileScript(language, example.value);
        expect(expressionValidationIssues(JSON.parse(JSON.stringify(expression)))).toBeUndefined();
      });
    }
    it(`${language}: reports parser errors`, () => expect(() => compileScript(language, '(')).toThrow());
  }
  it('reports unsupported languages and empty MolScript', () => {
    expect(() => compileScript('unknown' as any, 'all')).toThrow('Unsupported script language');
    expect(() => compileScript('mol-script', '')).toThrow('No query');
  });

  it('keeps plugin language registration independent of the authoring facade', () => {
    expect(Script.getAvailableLanguages()).toEqual(['mol-script']);
    compileScript('pymol', 'all');
    expect(Script.getAvailableLanguages()).toEqual(['mol-script']);
    expect(() => Script.toExpression({ language: 'pymol', expression: 'all' })).toThrow('not available');
  });
});
