/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The DOM UI's default `Hint` never blocks the pointer (seqspec:interaction#accessibility): the cases are
// hint-cases.ts, shared with React's.

import { describe } from 'vitest';
import * as D from '../../src/dom/index.js';
import { hintCases } from './hint-cases.js';

describe('the DOM Hint', () => {
  hintCases(async (host, buttons) => {
    const parts = buttons.map((b, i) =>
      D.Hint(host, {} as never, { label: `Hint ${i}`, children: (x: HTMLElement) => x.append(b) }),
    );
    return () => {
      for (const p of parts) p.dispose?.();
    };
  });
});
