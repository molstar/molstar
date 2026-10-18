/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type NameRule, namedPicker, pickNamed } from '../scope.js';

// `namedPicker` (api/scope.ts) answers as `pickNamed` does, for every kind of name: a ref, a key, a name of each tier,
// one that matches several, none, an empty one, and over an empty list.

interface Item {
  ref: string;
  key?: string;
  names: (string | undefined)[];
}

const rule: NameRule<Item> = {
  what: 'thing',
  ids: (x) => (x.key ? [x.ref, x.key] : [x.ref]),
  names: (x) => x.names,
  describe: (x) => `${x.names[0] ?? '?'} (${x.ref})`,
  none: { message: 'Nothing.', hint: 'Make one.' },
  listed: 'Things',
  more: ' or name it',
};

const items: Item[] = [
  { ref: 'r1', key: 'k1', names: ['Alpha', 'a', 'shared'] },
  { ref: 'r2', key: 'r2', names: ['Beta', 'b'] },
  { ref: 'r3', names: ['Gamma', undefined, 'shared'] },
  { ref: 'r4', key: 'k1', names: ['alpha two', 'b', 'tier3'] },
  { ref: 'r5', names: ['', 'r1'] },
];

describe('namedPicker', () => {
  it('answers as pickNamed does', () => {
    const names = [
      undefined,
      '',
      'r1',
      'r2',
      'k1',
      'Alpha',
      'ALPHA',
      ' alpha ',
      'a',
      'b',
      'shared',
      'tier3',
      'gamma',
      'nothing',
      '   ',
      'r5',
    ];
    for (const list of [items, items.slice(0, 1), []]) {
      const pick = namedPicker(list, rule);
      for (const name of names) expect(pick(name), `${name} of ${list.length}`).toEqual(pickNamed(list, name, rule));
    }
  });
});
