/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { utf8Encode } from '@molstar/core/util/utf8';
import { Tokenizer } from '../common/text/tokenizer.js';
import { TextByteTokenizer, trimBytesStr } from '../common/text/text-byte-tokenizer.js';
import { TextByteTokenColumnProvider } from '../common/text/column/text-byte-token.js';
import { TextByteFixedColumnProvider } from '../common/text/column/text-byte-fixed.js';
import { CifField } from '../cif/data-model.js';
import { parseCsv } from '../csv/parser.js';

const text = 'line one\r\n\tvalue  two \rthree\n\n  xé 中  y\nlast';

function lines(t: Tokenizer | TextByteTokenizer, read: (t: any) => string) {
  const out: string[] = [];
  while (t.position < t.length) out.push(read(t));
  return { out, lineNumber: t.lineNumber };
}

describe('text byte tokenizer', () => {
  it('reads lines as the string tokenizer does', () => {
    const bytes = utf8Encode(text);
    expect(lines(TextByteTokenizer(bytes), TextByteTokenizer.readLine)).toEqual(
      lines(Tokenizer(text), Tokenizer.readLine),
    );
    expect(lines(TextByteTokenizer(bytes), TextByteTokenizer.readLineTrim)).toEqual(
      lines(Tokenizer(text), Tokenizer.readLineTrim),
    );
  });

  it('reads values as the string tokenizer does', () => {
    function values(t: any, T: typeof Tokenizer | typeof TextByteTokenizer) {
      const out: [string, number][] = [];
      while (true) {
        T.skipWhitespace(t);
        if (t.position >= t.length) break;
        T.markStart(t);
        T.eatValue(t);
        out.push([T.getTokenString(t), t.lineNumber]);
      }
      return out;
    }
    expect(values(TextByteTokenizer(utf8Encode(text)), TextByteTokenizer)).toEqual(values(Tokenizer(text), Tokenizer));
  });

  it('readLines and markLines', () => {
    const t = TextByteTokenizer(utf8Encode('a\nb\nc'));
    expect(TextByteTokenizer.readLines(t, 2)).toEqual(['a', 'b']);
    const tokens = TextByteTokenizer.markLines(TextByteTokenizer(utf8Encode('a\nbc\nd')), 2);
    expect(tokens.count).toBe(2);
    expect(Array.from(tokens.indices).slice(0, 4)).toEqual([0, 1, 2, 4]);
  });

  it('readAllLinesAsync and readLinesAsync', async () => {
    const data = utf8Encode(Array.from({ length: 1000 }, (_, i) => `${i}`).join('\n'));
    const all = await TextByteTokenizer.readAllLinesAsync(data, SyncRuntimeContext, 7);
    expect(all.count).toBe(1000);
    const col = TextByteTokenColumnProvider(all)(Column.Schema.int);
    expect(col.value(999)).toBe(999);

    const t = TextByteTokenizer(data);
    const some = await TextByteTokenizer.readLinesAsync(t, 500, SyncRuntimeContext, 7);
    expect(some.count).toBe(500);
    expect(TextByteTokenizer.readLine(t)).toBe('500');
  });

  it('bytesOf encodes strings and passes bytes', () => {
    const bytes = new Uint8Array([97]);
    expect(TextByteTokenizer.bytesOf(bytes)).toBe(bytes);
    expect(Array.from(TextByteTokenizer.bytesOf('é'))).toEqual([0xc3, 0xa9]);
    // without a leading byte order mark, as a decoded string was
    expect(Array.from(TextByteTokenizer.bytesOf(new Uint8Array([0xef, 0xbb, 0xbf, 97])))).toEqual([97]);
    expect(Array.from(TextByteTokenizer.bytesOf('\ufeffa'))).toEqual([97]);
  });

  it('trimBytesStr', () => {
    const data = utf8Encode(' \té x \t');
    expect(trimBytesStr(data, 0, data.length)).toBe('é x');
    expect(trimBytesStr(data, 0, 2)).toBe('');
  });
});

describe('text byte columns', () => {
  it('token column', () => {
    const data = utf8Encode('321 é');
    const col = TextByteTokenColumnProvider({ data, indices: [0, 1, 1, 2, 2, 3, 4, 6], count: 4 });
    expect(col(Column.Schema.int).toArray()).toEqual([3, 2, 1, 0]);
    expect(col(Column.Schema.str).value(3)).toBe('é');
    expect(col(Column.Schema.str).areValuesEqual(0, 1)).toBe(false);
  });

  it('fixed column', () => {
    const lines = ['1.123 abc', '1.00  a', '1.1    bcd   ', '', ' 5'];
    const data = utf8Encode(lines.join('\n'));
    const indices: number[] = [];
    let last = 0;
    for (const l of lines) {
      indices.push(last, Math.min(last + l.length, data.length));
      last += l.length + 1;
    }
    const col = TextByteFixedColumnProvider({ data, indices, count: lines.length });
    expect(col(0, 5, Column.Schema.float).toArray()).toEqual([1.123, 1.0, 1.1, 0, 5]);
    expect(col(5, 4, Column.Schema.str).toArray()).toEqual(['abc', 'a', 'bc', '', '']);
    // a width past the end of the data reads the number to the line's end
    expect(col(0, 8, Column.Schema.int).value(4)).toBe(5);
  });

  it('CIF field', () => {
    const data = utf8Encode('. ? 12 é');
    const field = CifField.ofTextByteTokens({ data, indices: [0, 1, 2, 3, 4, 6, 7, 9], count: 4 });
    expect(field.toStringArray()).toEqual(['', '', '12', 'é']);
    expect(field.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
    expect(field.valueKind(1)).toBe(Column.ValueKinds.Unknown);
    expect(field.int(2)).toBe(12);
  });
});

describe('string cache', () => {
  // 6000 values of 0 to 10 bytes (past the cache's 4096-lookup trial), 1 in 5 with a non-ASCII char, from a small alphabet: many repeats and slot collisions
  const values: string[] = [];
  let seed = 7;
  const rand = (n: number) => (seed = (seed * 48271) % 0x7fffffff) % n;
  for (let i = 0; i < 6000; i++) {
    let v = '';
    for (let j = 0, l = rand(11); j < l; j++) v += 'abcAB12'[rand(7)];
    if (rand(5) === 0) v += '\u00e9';
    values.push(v);
  }

  it('token column and CIF field read every value', () => {
    const data = utf8Encode(values.join(' '));
    const indices: number[] = [];
    let p = 0;
    for (const v of values) {
      const l = utf8Encode(v).length;
      indices.push(p, p + l);
      p += l + 1;
    }
    const tokens = { data, indices, count: values.length };
    const col = TextByteTokenColumnProvider(tokens)(Column.Schema.str);
    const field = CifField.ofTextByteTokens(tokens);
    for (let pass = 0; pass < 2; pass++) {
      for (let i = 0; i < values.length; i++) {
        expect(col.value(i)).toBe(values[i]);
        expect(field.str(i)).toBe(values[i] === '.' || values[i] === '?' ? '' : values[i]);
      }
    }
  });

  it('CSV column, with quoted values', async () => {
    const text = 'a\n' + values.map((v, i) => (i % 3 === 0 ? `"${v}""q"` : v)).join('\n') + '\n';
    const parsed = await parseCsv(text).run();
    if (parsed.isError) throw new Error(parsed.message);
    const a = parsed.result.table.getColumn('a')!;
    // empty lines are skipped: the empty values not quoted are not records
    const expected = values.map((v, i) => (i % 3 === 0 ? `${v}"q` : v)).filter((v, i) => i % 3 === 0 || v !== '');
    expect(a.rowCount).toBe(expected.length);
    for (let i = 0; i < expected.length; i++) expect(a.str(i)).toBe(expected[i]);
  });
});
