/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { withoutTitle } from '../playback.js';

// A snapshot's description without a first heading that only repeats its title (the card's bar shows the title).

describe('a description without its repeated title', () => {
  it('drops a first heading that is the title, whatever its level and emphasis', () => {
    expect(withoutTitle('# The Birth of a Rogue Kinase\n\nBut in BCR-ABL…', 'The Birth of a Rogue Kinase')).toBe(
      'But in BCR-ABL…',
    );
    expect(withoutTitle('## **Overview**\nText', 'overview')).toBe('Text');
  });

  it('drops a heading the title only adds a part counter or more words to', () => {
    const md = '\n### ATP Binding and Unstoppable Signaling\n\nTo function, every kinase needs ATP.';
    expect(withoutTitle(md, 'ATP Binding and Unstoppable Signaling [1/2]')).toBe(
      'To function, every kinase needs ATP.',
    );
    expect(withoutTitle(md, 'ATP Binding and Unstoppable Signaling (2/2)')).toBe(
      'To function, every kinase needs ATP.',
    );
    expect(withoutTitle('# The End\nBye', 'The End.')).toBe('Bye');
    expect(withoutTitle('# Overview\nText', 'Overview of the active site')).toBe('Text');
  });

  it('keeps a heading that says more, and text with no heading', () => {
    const md = '# The Structural Story of BCR-ABL: A Kinase Out of Control\n\nText';
    expect(withoutTitle(md, 'A Kinase Out of Control')).toBe(md);
    expect(withoutTitle('Just text', 'Just text')).toBe('Just text');
    // The words must be whole: "Arg" is not "Arginine Residues".
    const tbp = '### Arginine Residues\nText';
    expect(withoutTitle(tbp, 'TBP: Arginine [2/5]')).toBe(tbp);
    const conserved = '### TATA-Binding Protein Is Highly Conserved in Eukaryotes\nText';
    expect(withoutTitle(conserved, 'TBP: Highly Conserved in Eukaryotes [1/2]')).toBe(conserved);
  });
});
