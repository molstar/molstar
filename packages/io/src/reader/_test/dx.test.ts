/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parseDx } from '../dx/parser.js';

const header = `# Data from APBS, Ångström
object 1 class gridpositions counts 2 2 2
origin -1.0 -2.0 0.5
delta 0.5 0 0
delta 0 0.25 0
delta 0 0 2.0
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
`;

const dx = `${header}1.0 2.0 -3.5
4.0e-1 5 6
7 8
attribute "dep" string "positions"
`;

describe('dx reader', () => {
  it('parses text DX', async () => {
    const parsed = await parseDx(dx, 'test').run();
    if (parsed.isError) throw new Error(parsed.message);
    const { header, values } = parsed.result;
    expect(Array.from(header.dim)).toEqual([2, 2, 2]);
    expect(Array.from(header.min)).toEqual([-1, -2, 0.5]);
    expect(Array.from(header.h)).toEqual([0.5, 0.25, 2]);
    expect(Array.from(values)).toEqual([1, 2, -3.5, 0.4, 5, 6, 7, 8]);
  });

  it('parses text DX the same from a string and from its UTF-8 bytes', async () => {
    for (const text of [dx, dx.replace(/\n/g, '\r\n')]) {
      const fromString = await parseDx(text, 'test').run();
      const fromBytes = await parseDx(utf8Encode(text), 'test', { binary: false }).run();
      if (fromString.isError || fromBytes.isError) throw new Error('parse failed');
      expect(fromBytes.result).toEqual(fromString.result);
    }
  });

  it('skips a leading UTF-8 byte order mark', async () => {
    const text = dx.slice(dx.indexOf('object 1'));
    const withBom = new Uint8Array([0xef, 0xbb, 0xbf, ...utf8Encode(text)]);
    const fromBom = await parseDx(withBom, 'test', { binary: false }).run();
    const fromString = await parseDx(text, 'test').run();
    if (fromString.isError || fromBom.isError) throw new Error('parse failed');
    expect(fromBom.result).toEqual(fromString.result);
  });

  it('parses binary DX after a non-ASCII header', async () => {
    const head = utf8Encode(header);
    const bytes = new Uint8Array(head.length + 8 * 8);
    bytes.set(head);
    const dv = new DataView(bytes.buffer, head.length);
    for (let i = 0; i < 8; i++) dv.setFloat64(8 * i, i - 2.5, true);
    const parsed = await parseDx(bytes, 'test').run();
    if (parsed.isError) throw new Error(parsed.message);
    expect(Array.from(parsed.result.values)).toEqual([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]);
  });

  it('fails without an object 3 record', async () => {
    await expect(parseDx('object 1 class gridpositions counts 2 2 2\n', 'test').run()).rejects.toThrow(
      "no 'object 3' record found",
    );
  });
});
