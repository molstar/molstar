/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parseCsv } from '../csv/parser.js';

async function parse(data: string | Uint8Array, opts?: Parameters<typeof parseCsv>[1]) {
  const parsed = await parseCsv(data, opts).run();
  if (parsed.isError) throw new Error(parsed.toString());
  return parsed.result.table;
}

async function parseError(data: string) {
  const parsed = await parseCsv(data).run();
  if (!parsed.isError) throw new Error('expected an error');
  return parsed;
}

function strings(table: Awaited<ReturnType<typeof parse>>) {
  return table.columnNames.map((n) => table.getColumn(n)!.toStringArray());
}

describe('csv reader', () => {
  it('basic', async () => {
    const table = await parse('StrCol,IntCol,FloatCol\nstring1,-1,-0.34e3\nstring2,42,2.44');
    expect(table.rowCount).toBe(2);
    expect(table.getColumn('StrCol')!.toStringArray()).toEqual(['string1', 'string2']);
    expect(table.getColumn('IntCol')!.toIntArray()).toEqual([-1, 42]);
    expect(table.getColumn('FloatCol')!.toFloatArray()).toEqual([-340.0, 2.44]);
  });

  it('reads bytes as strings do', async () => {
    const text = 'a,b\nxé,1\n"y,中",2\n';
    const fromBytes = await parse(utf8Encode(text));
    expect(strings(fromBytes)).toEqual([
      ['xé', 'y,中'],
      ['1', '2'],
    ]);
    expect(strings(fromBytes)).toEqual(strings(await parse(text)));
  });

  it('keeps spaces in fields; numbers parse past them', async () => {
    const table = await parse('a,b\n  x \t, -1 \n');
    expect(table.getColumn('a')!.str(0)).toBe('  x \t');
    expect(table.getColumn('b')!.str(0)).toBe(' -1 ');
    expect(table.getColumn('b')!.int(0)).toBe(-1);
    expect(table.getColumn('b')!.float(0)).toBe(-1);
  });

  it('quoted fields: delimiters, line breaks, doubled quotes', async () => {
    const table = await parse('"Int Col",b\r\n"1,2","x\r\ny"\r\n"say ""hi""",""\r\n""""," "\r\n');
    expect(table.columnNames).toEqual(['Int Col', 'b']);
    expect(strings(table)).toEqual([
      ['1,2', 'say "hi"', '"'],
      ['x\r\ny', '', ' '],
    ]);
    const b = table.getColumn('b')!;
    expect(b.valueKind(1)).toBe(Column.ValueKinds.NotPresent);
    expect(b.valueKind(2)).toBe(Column.ValueKinds.Present);
  });

  it('a quote inside an unquoted field is text', async () => {
    const table = await parse('a\nx"y\n');
    expect(table.getColumn('a')!.str(0)).toBe('x"y');
  });

  it('empty fields, also the last of a record', async () => {
    const table = await parse('a,b,c\n,,\n1,,\n,2,\n');
    expect(table.rowCount).toBe(3);
    expect(strings(table)).toEqual([
      ['', '1', ''],
      ['', '', '2'],
      ['', '', ''],
    ]);
    expect(table.getColumn('a')!.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
  });

  it('line breaks: CRLF, LF and CR; the last is optional; empty lines skipped', async () => {
    for (const text of ['a,b\r\n1,2\r\n3,4', 'a,b\n1,2\n3,4\n', 'a,b\r1,2\r3,4\r', 'a,b\n\n1,2\r\n\r\n3,4\n\n']) {
      const table = await parse(text);
      expect(table.rowCount).toBe(2);
      expect(strings(table)).toEqual([
        ['1', '3'],
        ['2', '4'],
      ]);
    }
  });

  it('short records read as empty, long ones are cut to the header', async () => {
    const table = await parse('a,b\n1\n2,3,4\n');
    expect(strings(table)).toEqual([
      ['1', '2'],
      ['', '3'],
    ]);
    expect(table.getColumn('b')!.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
  });

  it('. and ? are values', async () => {
    const table = await parse('a\n.\n?\n');
    expect(table.getColumn('a')!.toStringArray()).toEqual(['.', '?']);
    expect(table.getColumn('a')!.valueKind(0)).toBe(Column.ValueKinds.Present);
  });

  it('# is text unless comment is set, and then only at the start of a line', async () => {
    expect(strings(await parse('a,b\n#x,1\n'))).toEqual([['#x'], ['1']]);
    const table = await parse('a,b\n# comment\n1,#2\n', { comment: '#' });
    expect(strings(table)).toEqual([['1'], ['#2']]);
  });

  it('byte order mark', async () => {
    const table = await parse('﻿a,b\n1,2\n');
    expect(table.columnNames).toEqual(['a', 'b']);
  });

  it('no column names', async () => {
    const table = await parse('\n1,2\n3,4\n', { noColumnNames: true });
    expect(table.columnNames).toEqual(['0', '1']);
    expect(strings(table)).toEqual([
      ['1', '3'],
      ['2', '4'],
    ]);
  });

  it('tabs', async () => {
    const table = await parse('StrCol\tIntCol\tFloatCol\nstring1\t-1\t-0.34e3\nstring2\t42\t2.44', { delimiter: '\t' });
    expect(table.getColumn('StrCol')!.toStringArray()).toEqual(['string1', 'string2']);
    expect(table.getColumn('IntCol')!.toIntArray()).toEqual([-1, 42]);
    expect(table.getColumn('FloatCol')!.toFloatArray()).toEqual([-340.0, 2.44]);
  });

  it('quoted fields with a tab or space delimiter', async () => {
    expect(strings(await parse('a\tb\n"x"\ty\n', { delimiter: '\t' }))).toEqual([['x'], ['y']]);
    expect(strings(await parse('a b\n"x" y\n', { delimiter: ' ' }))).toEqual([['x'], ['y']]);
  });

  it('comment must differ from the quote and the delimiter', async () => {
    await expect(parseCsv('a\n', { comment: '"' }).run()).rejects.toThrow();
    await expect(parseCsv('a\n', { comment: ',' }).run()).rejects.toThrow();
  });

  it('empty input', async () => {
    const table = await parse('');
    expect(table.rowCount).toBe(0);
    expect(table.columnNames).toEqual([]);
  });

  it('header only', async () => {
    const table = await parse('a,b\n');
    expect(table.rowCount).toBe(0);
    expect(table.getColumn('a')!.rowCount).toBe(0);
  });

  it('errors: unterminated quote, text after a closing quote', async () => {
    const unterminated = await parseError('a\n"x\ny\n');
    expect(unterminated.line).toBe(2);
    const after = await parseError('a,b\n1,2\n"x"y,3\n');
    expect(after.line).toBe(3);
    // spaces after a closing quote are fine
    expect(strings(await parse('a,b\n"x"  ,1\n'))).toEqual([['x'], ['1']]);
  });

  it('equal values', async () => {
    const table = await parse('a\n"x""y"\nx"y\n"x""y"\n');
    const a = table.getColumn('a')!;
    expect(a.areValuesEqual(0, 2)).toBe(true);
    // the same value, quoted and not
    expect(a.str(1)).toBe('x"y');
    expect(a.areValuesEqual(0, 1)).toBe(true);
  });
});
