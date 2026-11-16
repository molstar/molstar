/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { namedCatalog } from '../named-catalog.js';

describe('namedCatalog', () => {
  it('returns the catalog as given, its keys checked against the names at compile time', () => {
    const catalog = { a: { name: 'a' as const, v: 1 }, b: { name: 'b' as const, v: 2 } };
    expect(namedCatalog(catalog)).toBe(catalog);
    // a key that is not the name of its provider does not compile
    // @ts-expect-error
    namedCatalog({ a: { name: 'b' as const } });
  });
});
