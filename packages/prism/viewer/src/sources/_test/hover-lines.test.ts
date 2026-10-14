/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import { hoverLines } from '../hover-lines.js';

it('Mol* loci labels become plain lines', () => {
  expect(hoverLines([])).toEqual([]);
  expect(hoverLines(['Cytochrome P450 3A4</br><b>A</b> | <i>ALA</i> 12 &amp; more'])).toEqual([
    'Cytochrome P450 3A4',
    'A | ALA 12 & more',
  ]);
  expect(hoverLines(['Sulfate ion <small>|| × 2</small>', { react: true }, ''])).toEqual(['Sulfate ion × 2']);
});
