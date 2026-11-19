/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { dataFilter } from '../behaviors/theming.js';

// The viewer's filter is a layer's trailing `(not …)` only (files.md, One track a feature type): what leads it, the
// annotation and a leading item form, is the track's base, kept when the filter is set, toggled or cleared.

describe('dataFilter', () => {
  it('a trailing `(not …)` is the user’s filter: the annotation before it is the base', () => {
    expect(dataFilter('/data/f/(not :of (item :category linker))', 'category')).toEqual({
      base: '/data/f',
      hidden: ['linker'],
    });
    expect(dataFilter('/data/f/(not :of (item :category [a b]))', 'category')).toEqual({
      base: '/data/f',
      hidden: ['a', 'b'],
    });
  });

  it('a leading item form is the layer’s base: kept with the annotation, never listed as hidden', () => {
    expect(dataFilter('/data/f/(item :category Domain)/(not :of (item :category linker))', 'category')).toEqual({
      base: '/data/f/(item :category Domain)',
      hidden: ['linker'],
    });
    // a track of one feature type with no user filter: the base is the whole data, nothing hidden
    expect(dataFilter('/data/f/(item :category Domain)', 'category')).toEqual({
      base: '/data/f/(item :category Domain)',
      hidden: [],
    });
    // a quoted value with a space
    expect(dataFilter('/data/f/(item :category "DNA binding")/(not :of (item :id x))', 'id')).toEqual({
      base: '/data/f/(item :category "DNA binding")',
      hidden: ['x'],
    });
  });

  it('data with no filter is its own base', () => {
    expect(dataFilter('/data/f', 'category')).toEqual({ base: '/data/f', hidden: [] });
  });
});
