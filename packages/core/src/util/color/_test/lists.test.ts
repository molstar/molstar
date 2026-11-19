/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import * as Lists from '../lists.js';

/** A list's constant by its key: `red-blue` is `ColorListRedBlue`. */
const constantOf = (key: string) =>
  `ColorList${key
    .split('-')
    .map((p) => p[0].toUpperCase() + p.slice(1))
    .join('')}`;

describe('the color lists', () => {
  it('each list is a named constant, and `ColorLists` holds exactly those constants, by name', () => {
    const keys = Object.keys(Lists.ColorLists);
    expect(keys.length).toBe(59);
    for (const key of keys)
      expect((Lists as Record<string, unknown>)[constantOf(key)], key).toBe(
        Lists.ColorLists[key as Lists.ColorListName],
      );
    // no named list outside the table (the names and options are not lists)
    const named = Object.keys(Lists).filter(
      (k) => /^ColorList[A-Z0-9]/.test(k) && !/^ColorList(Names|Options)/.test(k),
    );
    expect(named.length).toBe(keys.length);
    expect(Lists.ColorListNames).toEqual(keys);
  });
});
