/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { MolScriptBuilder as B } from '../../language/builder.js';
import { transpiler as pymol } from '../pymol/parser.js';
import { transpiler as vmd } from '../vmd/parser.js';

// Text the transpilers cannot read fails as a parse error that says where, never as a stack overflow or a TypeError:
// what a person types is unfinished most of the time (a command line, a target field).

const Unfinished = ['', '(', '((', 'chain', 'chain ', 'resn', 'b', 'not', 'not ', 'byres', 'byres ', 'name CA and'];

function failsReadably(transpile: (s: string) => unknown, text: string) {
  try {
    transpile(text);
  } catch (e) {
    expect(e, `'${text}'`).not.toBeInstanceOf(RangeError);
    expect(e, `'${text}'`).not.toBeInstanceOf(TypeError);
    return;
  }
  expect.fail(`'${text}' read`);
}

describe('unreadable text', () => {
  it.each(Unfinished)("PyMOL fails readably on '%s'", (text) => failsReadably(pymol, text));
  it.each([...Unfinished, 'resn ', 'b2', 'chain A and resn'])("VMD fails readably on '%s'", (text) =>
    failsReadably(vmd, text),
  );
});

// Unclosed groups, deep: each level is read once, not once per alternative (three a level was exponential: 30 levels
// did not finish in minutes). No clock: an exponential grammar times the test out.
describe('deep groups', () => {
  const deep = (n: number, closed: boolean) => `${'('.repeat(n)}chain A${closed ? ')'.repeat(n) : ''}`;

  it('PyMOL reads 40 closed levels and fails readably on 40 unclosed', () => {
    expect(pymol(deep(40, true))).toBeTruthy();
    failsReadably(pymol, deep(40, false));
  });

  it('VMD reads 20 closed levels and fails readably on 20 unclosed', () => {
    expect(vmd(deep(20, true))).toBeTruthy();
    failsReadably(vmd, deep(20, false));
  });
});

describe('PyMOL grammar', () => {
  const all = B.struct.generator.all();

  it('reads a name that starts with a reserved word (b2, both)', () => {
    for (const name of ['b2', 'both', 'notes', 'within_a']) {
      expect(pymol(`${name} and name CA`, { resolveName: (n) => (n === name ? all : undefined) })).toBeTruthy();
    }
  });

  it('reads a prefix operator whose operand is another prefix operator', () => {
    expect(JSON.stringify(pymol('byres bound_to resn CYS'))).toContain('query-in-selection');
    expect(pymol('not not chain A')).toBeTruthy();
  });

  it('reads a group holding a group and an operator', () => {
    expect(pymol('byres ((chain A and polymer) within 4 of (chain A and resn HEM))')).toBeTruthy();
  });
});
