/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The browser config reaches Chromium with a canvas (the surface and renderer tests of later phases need one).

import { expect, it } from 'vitest';

it('runs in a browser with a 2D canvas', () => {
  expect(document.createElement('canvas').getContext('2d')).not.toBeNull();
});
