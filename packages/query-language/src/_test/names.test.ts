/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { compileScript } from '../compile.js';
import { MolScriptBuilder as B } from '../language/builder.js';
import { expressionValidationIssues } from '../language/validation.js';
import { transpiler } from '../transpilers/pymol/parser.js';

// PyMOL's objects and named selections through `resolveName` (TranspileOptions): what a command layer knows.

const lig = B.struct.generator.atomGroups({
  'residue-test': B.core.rel.eq([B.ammp('label_comp_id'), 'HEM']),
});
const all = B.struct.generator.all();
const names: Record<string, ReturnType<typeof B.struct.generator.all>> = { lig, '1abc': all, sele: lig };
const options = { resolveName: (name: string) => names[name] };

describe('PyMOL names', () => {
  it('fails on a name without a resolver, as before', () => {
    expect(() => transpiler('foobar')).toThrow(/object' not supported/);
    expect(() => transpiler('/pept//a/10')).toThrow(/object' not supported/);
  });

  it('fails on a name the resolver does not know', () => {
    expect(() => transpiler('foobar', options)).toThrow(/no object or selection named 'foobar'/);
  });

  it('reads a named selection as its expression', () => {
    expect(transpiler('lig', options)).toEqual(lig);
    expect(transpiler('sele', options)).toEqual(lig);
    const e = transpiler('lig and chain A', options);
    expect(JSON.stringify(e)).toContain('HEM');
    expect(expressionValidationIssues(e)).toBeUndefined();
  });

  it('reads an object in a macro, with and without the rest of the macro', () => {
    expect(transpiler('/1abc', options)).toEqual(all);
    const e = transpiler('/1abc//A/10/CA', options);
    expect(JSON.stringify(e)).toContain('structure-query.modifier.intersect-by');
    expect(expressionValidationIssues(e)).toBeUndefined();
  });

  it('passes the options through compileScript, and leaves the next call without them', () => {
    expect(expressionValidationIssues(compileScript('pymol', 'not lig', options))).toBeUndefined();
    expect(() => compileScript('pymol', 'lig')).toThrow(/object' not supported/);
  });
});
