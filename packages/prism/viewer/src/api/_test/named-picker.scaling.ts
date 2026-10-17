/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { it } from 'vitest';
import { checkLinear, scaling } from '../../../../ui/src/_test/scaling.js';
import { type NameRule, namedPicker, pickNamed } from '../scope.js';

// Naming N things of a list of N one by one grows about linearly through `namedPicker` (`test:scaling`, prism-ui's
// _test/scaling.ts; the superposition's named structures): the list indexed once. `pickNamed` per name is the positive
// control (N × the list). At 500 and 2,000.

interface Item {
  ref: string;
  label: string;
}

const rule: NameRule<Item> = {
  what: 'structure',
  ids: (x) => [x.ref],
  names: (x) => [x.label, x.label.toUpperCase()],
  describe: (x) => `${x.label} (${x.ref})`,
  none: { message: 'Nothing is loaded.' },
};

const listOf = (n: number) => {
  const list = Array.from({ length: n }, (_, i) => ({ ref: `ref-${i}`, label: `Entry ${i}` }));
  // By name (the second tier's lookups as well as the first's), as an operation's caller names them.
  return { list, names: list.map((x, i) => (i % 2 ? x.label : x.ref)) };
};

it('N names picked from N through namedPicker', () => {
  const s = scaling(
    {
      setup: listOf,
      run: ({ list, names }) => {
        const pick = namedPicker(list, rule);
        for (const name of names) pick(name);
      },
    },
    500,
  );
  checkLinear('namedPicker, N names of N', s);
});

it('the control: N names picked from N one by one through pickNamed', () => {
  const s = scaling(
    {
      setup: listOf,
      run: ({ list, names }) => {
        for (const name of names) pickNamed(list, name, rule);
      },
    },
    500,
  );
  checkLinear('pickNamed per name (control)', s, false);
});
