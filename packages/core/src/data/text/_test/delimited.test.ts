/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type DelimitedOptions, readDelimited } from '../delimited.js';

function records(text: string, options?: DelimitedOptions) {
  const r = readDelimited(text, options);
  if (!r.ok) throw new Error(r.message);
  return r.records;
}

function error(text: string, options?: DelimitedOptions) {
  const r = readDelimited(text, options);
  if (r.ok) throw new Error('expected an error');
  return r;
}

describe('readDelimited', () => {
  it('reads records of cells, spaces kept', () => {
    expect(records('a,b\n  x \t, -1 \n')).toEqual([
      ['a', 'b'],
      ['  x \t', ' -1 '],
    ]);
  });

  it('quoted cells: delimiters, line breaks, doubled quotes', () => {
    expect(records('"Int Col",b\r\n"1,2","x\r\ny"\r\n"say ""hi""",""\r\n""""," "\r\n')).toEqual([
      ['Int Col', 'b'],
      ['1,2', 'x\r\ny'],
      ['say "hi"', ''],
      ['"', ' '],
    ]);
  });

  it('a quote inside an unquoted cell is text', () => {
    expect(records('a\nx"y\n')).toEqual([['a'], ['x"y']]);
  });

  it('empty cells, also the last of a record', () => {
    expect(records('a,b,c\n,,\n1,,\n,2,')).toEqual([
      ['a', 'b', 'c'],
      ['', '', ''],
      ['1', '', ''],
      ['', '2', ''],
    ]);
  });

  it('line breaks: CRLF, LF and CR; the last optional; empty lines skipped; lines counted', () => {
    for (const text of ['a,b\r\n1,2\r\n3,4', 'a,b\n1,2\n3,4\n', 'a,b\r1,2\r3,4\r']) {
      const r = readDelimited(text);
      expect(r.ok && r.records).toEqual([
        ['a', 'b'],
        ['1', '2'],
        ['3', '4'],
      ]);
      expect(r.ok && r.lines).toEqual([1, 2, 3]);
    }
    const r = readDelimited('a,b\n\n1,2\r\n\r\n"3\n",4\n\n5,6');
    expect(r.ok && r.lines).toEqual([1, 3, 5, 8]);
  });

  it('ragged records are as written', () => {
    expect(records('a,b\n1\n2,3,4\n')).toEqual([['a', 'b'], ['1'], ['2', '3', '4']]);
  });

  it('comments only when set, and then only at the start of a line', () => {
    expect(records('a\n#x\n')).toEqual([['a'], ['#x']]);
    expect(records('a,b\n# comment\n1,#2\n', { comment: '#' })).toEqual([
      ['a', 'b'],
      ['1', '#2'],
    ]);
  });

  it('skips a byte order mark', () => {
    expect(records('﻿a,b\n1,2\n')).toEqual([
      ['a', 'b'],
      ['1', '2'],
    ]);
  });

  it('tabs, and quoted cells with a tab or space delimiter', () => {
    expect(records('a\tb\n"x"\ty\n', { delimiter: '\t' })).toEqual([
      ['a', 'b'],
      ['x', 'y'],
    ]);
    expect(records('a b\n"x" y\n', { delimiter: ' ' })).toEqual([
      ['a', 'b'],
      ['x', 'y'],
    ]);
  });

  it('empty input', () => {
    expect(records('')).toEqual([]);
    expect(records('\n\r\n')).toEqual([]);
  });

  it('errors: an unterminated quote, text after a closing quote', () => {
    expect(error('a\n"x\ny\n').line).toBe(2);
    expect(error('a,b\n1,2\n"x"y,3\n').line).toBe(3);
    expect(records('a,b\n"x"  ,1\n')).toEqual([
      ['a', 'b'],
      ['x', '1'],
    ]);
  });

  it('refuses bad options', () => {
    expect(() => readDelimited('a', { comment: '"' })).toThrow();
    expect(() => readDelimited('a', { comment: ',' })).toThrow();
    expect(() => readDelimited('a', { delimiter: '\n' })).toThrow();
    expect(() => readDelimited('a', { delimiter: ',,' })).toThrow();
  });

  it('round-trips random tables written by the RFC (seeded)', () => {
    let seed = 7;
    const rand = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const alphabet = ['a', 'b', ' ', ',', '"', '\n', '\r\n', 'é', '\t', '1'];
    const cell = () => {
      let s = '';
      const n = Math.floor(rand() * 5);
      for (let i = 0; i < n; i++) s += alphabet[Math.floor(rand() * alphabet.length)];
      return s;
    };
    const write = (c: string) => (/[",\r\n]/.test(c) || c === '' ? `"${c.replace(/"/g, '""')}"` : c);
    for (let t = 0; t < 50; t++) {
      const columns = 1 + Math.floor(rand() * 4);
      const table = Array.from({ length: 1 + Math.floor(rand() * 6) }, () => Array.from({ length: columns }, cell));
      const text = table.map((r) => r.map(write).join(',')).join(rand() < 0.5 ? '\n' : '\r\n');
      expect(records(text)).toEqual(table);
    }
  });
});
