/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { decodeLongPath, decodePax, encodePax, encodeTarHeader } from '../tar';

describe('model server tar helpers', () => {
  it('round-trips PAX path, link path, and custom fields', () => {
    const name = `nested/${'long-segment/'.repeat(10)}file.cif`;
    const encoded = encodePax({ name, linkname: 'target.cif', pax: { comment: 'fixture' } });

    expect(decodePax(encoded)).toEqual({ path: name, linkpath: 'target.cif', comment: 'fixture' });
  });

  it('decodes a GNU long path and emits a USTAR header', () => {
    const name = 'nested/structure.cif';
    const header = encodeTarHeader({ name });

    expect(decodeLongPath(Buffer.from(name))).toBe(name);
    expect(header).not.toBeNull();
    expect(header!.toString('utf8', 0, name.length)).toBe(name);
    expect(header!.toString('ascii', 257, 262)).toBe('ustar');
  });
});
