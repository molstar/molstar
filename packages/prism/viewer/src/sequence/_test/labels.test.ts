/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import { labelLines } from '../sequence-control.js';

it("a Seq* hover label is the shell label's lines: the title, then key: value", () => {
  expect(labelLines(null)).toEqual([]);
  expect(
    labelLines({
      title: 'ALA 12',
      lines: [
        { key: 'Chain', value: 'A' },
        { key: 'Note', value: '' },
      ],
    }),
  ).toEqual(['ALA 12', 'Chain: A', 'Note']);
});
