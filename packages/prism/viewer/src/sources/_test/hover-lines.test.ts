/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { hoverLines, hoverRich } from '../hover-lines.js';

it('Mol* loci labels become plain lines', () => {
  expect(hoverLines([])).toEqual([]);
  expect(hoverLines(['Cytochrome P450 3A4</br><b>A</b> | <i>ALA</i> 12 &amp; more'])).toEqual([
    'Cytochrome P450 3A4',
    'A | ALA 12 & more',
  ]);
  expect(hoverLines(['Sulfate ion <small>|| × 2</small>', { react: true }, ''])).toEqual(['Sulfate ion × 2']);
});

describe('hoverRich', () => {
  it("keeps Mol*'s emphasis as runs, the first line a title", () => {
    expect(hoverRich(['Cytochrome P450 3A4</br><b>A</b> | <i>ALA</i> 12 &amp; more'])).toEqual([
      'Cytochrome P450 3A4',
      {
        role: 'text',
        runs: [{ text: 'A', bold: true }, { text: ' | ' }, { text: 'ALA', italic: true }, { text: ' 12 & more' }],
      },
    ]);
    expect(hoverRich(['Sulfate ion <small>|| × 2</small>'])).toEqual([
      { role: 'title', runs: [{ text: 'Sulfate ion ' }, { text: '× 2', small: true }] },
    ]);
  });

  it('reads a label with newlines as Markdown', () => {
    expect(
      hoverRich(['# Binding site\nResidues **near** the *ligand*:\n- `HIS 64`\n2. [UniProt](https://www.uniprot.org)']),
    ).toEqual([
      { role: 'title', runs: [{ text: 'Binding site' }] },
      {
        role: 'text',
        runs: [
          { text: 'Residues ' },
          { text: 'near', bold: true },
          { text: ' the ' },
          { text: 'ligand', italic: true },
          { text: ':' },
        ],
      },
      { role: 'item', runs: [{ text: '• ' }, { text: 'HIS 64', code: true }] },
      { role: 'item', runs: [{ text: '2. ' }, { text: 'UniProt', link: true }] },
    ]);
  });

  it('never lets markup through: tags are dropped, their text kept as text', () => {
    const lines = hoverRich([
      '<img src=x onerror=alert(1)>Name<script>alert(2)</script>',
      'Line one\n<iframe src="javascript:alert(3)"></iframe>[x](javascript:alert(4)) `<b>code</b>`',
    ]);
    const text = JSON.stringify(lines);
    expect(text).not.toMatch(/<img|<script|<iframe|onerror/);
    // HTML in Markdown is skipped (v5's skipHtml); a code span keeps its text literally (shown with textContent, never
    // parsed); a javascript: link keeps only its text.
    const runs = lines.flatMap((l) => (typeof l === 'string' ? [] : l.runs));
    expect(runs.find((r) => r.code)?.text).toBe('<b>code</b>');
    expect(runs.find((r) => r.link)?.text).toBe('x');
  });
});
