/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

// Every dataset's index page renders, in a sans-serif font, and links each of its views to a story (a view added to a
// dataset without its index entry threw while the page rendered: Genomes' Manhattan, 2026-11-24).

import { describe, expect, it } from 'vitest';
import * as Antibody from '../antibody/index.stories.js';
import * as Complex from '../complex/index.stories.js';
import * as Genomes from '../genomes/index.stories.js';
import * as Globins from '../globins/index.stories.js';
import * as P53 from '../p53/index.stories.js';
import * as Plasmid from '../plasmid/index.stories.js';
import * as Trna from '../trna/index.stories.js';

const Pages = { Antibody, Complex, Genomes, Globins, P53, Plasmid, Trna } as const;

describe('the examples’ index pages', () => {
  for (const [name, page] of Object.entries(Pages)) {
    it(`${name}: renders, sans-serif, every link a story`, () => {
      const render = page.Views.render as () => HTMLElement;
      const root = render();
      const font = root.style.font;
      expect(font, 'a font the browser keeps').not.toBe('');
      expect(font).toMatch(/sans-serif/);
      const links = [...root.querySelectorAll('a')];
      expect(links.length).toBeGreaterThan(0);
      for (const a of links) expect(a.getAttribute('href'), a.textContent ?? '').not.toMatch(/undefined|--$/);
    });
  }
});
