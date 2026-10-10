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

  it.each([
    '(sel.atom.res :typo true)',
    '(sel.atom.res :typo (not-a-function))',
    '(sel.atom.atoms true false)',
    '(sel.atom.atom-groups :atom-test (sel.atom.res :typo true))',
    '(sel.atom.res (not-a-function))',
    '(bond.is :typo true)',
  ])('rejects invalid arguments and calls before expanding macros: %s', (source) => {
    expect(() => compileScript('mol-script', source)).toThrow(/Unknown argument|Unknown callable symbol/);
  });

  it.each([
    '(sel.atom.res true)',
    '(sel.atom.chains (= atom.chain `A`))',
    '(sel.atom.atoms :0 true)',
    '(structure-query.generator.all)',
  ])('accepts valid macros and canonical calls: %s', (source) => {
    expect(expressionValidationIssues(compileScript('mol-script', source))).toBeUndefined();
  });

  it.each([
    '(sel.atom.all) (not-a-function)',
    '(sel.atom.all) (sel.atom.atom-groups :atom-tset true)',
    '(sel.atom.all) (sel.atom.all)',
  ])('rejects multiple top-level expressions: %s', (source) => {
    expect(() => compileScript('mol-script', source)).toThrow('Expected exactly one MolScript expression');
  });

  it('accepts comments around a single expression', () => {
    expect(compileScript('mol-script', '; before\n(sel.atom.all)\n; after\n')).toEqual(
      compileScript('mol-script', '(sel.atom.all)'),
    );
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
