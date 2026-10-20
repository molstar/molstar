/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { hoverLines, hoverRich, pickHoverLabel } from '../hover-lines.js';

it('Mol* loci labels become plain lines', () => {
  expect(hoverLines([])).toEqual([]);
  expect(hoverLines(['Cytochrome P450 3A4</br><b>A</b> | <i>ALA</i> 12 &amp; more'])).toEqual([
    'Cytochrome P450 3A4',
    'A | ALA 12 & more',
  ]);
  expect(hoverLines(['Sulfate ion <small>|| × 2</small>', { react: true }, ''])).toEqual(['Sulfate ion × 2']);
});

it('decodes only the known entities, never what an object inherits', () => {
  expect(hoverLines(['&constructor; &toString; &lt;b&gt;'])).toEqual(['&constructor; &toString; <b>']);
  expect(hoverRich(['<b>&constructor;</b> &valueOf;'])).toEqual([
    { role: 'title', runs: [{ text: '&constructor;', bold: true }, { text: ' &valueOf;' }] },
  ]);
});

it('decodes numeric references and the common names (a no-break space a plain one)', () => {
  expect(hoverLines(['&#x2B05; back&nbsp;to &mdash; 2&#8201;&Aring;'])).toEqual([
    '\u2B05 back to \u2014 2\u2009\u00c5',
  ]);
  // Markdown (a label with a newline): the same decoding in its text runs.
  const text = hoverRich(['A tooltip\n\n&#x2B05; *Basic&nbsp;workflow*'])
    .map((l) => (typeof l === 'string' ? l : l.runs.map((r) => r.text).join('')))
    .join('|');
  expect(text).toContain('\u2B05 Basic workflow');
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

describe('pickHoverLabel', () => {
  it("prefers Mol*'s label, the sequence's only where its hover does not map", () => {
    const mol = ['CRABP2', '1CBS | A | ARG 59'];
    const seq = ['A R59'];
    expect(pickHoverLabel(mol, seq, true)).toBe(mol);
    expect(pickHoverLabel(mol, seq, false)).toBe(mol);
    expect(pickHoverLabel([], seq, false)).toBe(seq);
    // Mapped, Mol*'s label not in yet: nothing new (the label holds what it had a frame).
    expect(pickHoverLabel([], seq, true)).toEqual([]);
    expect(pickHoverLabel([], [], false)).toEqual([]);
  });
});

describe('labels from documents cannot make the parsing slow', () => {
  const N = 100_000;
  const adversarial = [
    ' '.repeat(N),
    '<'.repeat(N),
    '< '.repeat(N / 2),
    `<${' '.repeat(N)}`,
    `<${' '.repeat(N)}x`,
    '<b'.repeat(N / 2),
    '<'.repeat(N / 2) + '>'.repeat(N / 2),
    '> '.repeat(N / 2),
    '|'.repeat(N),
    `a${' '.repeat(N)}|`,
    '& '.repeat(N / 2),
    `&${'a'.repeat(N)}`,
    `<${' '.repeat(N / 2)}/${' '.repeat(N / 2)}br`,
  ];

  it('reads 100 KB of adversarial text in under 50 ms, each way', () => {
    for (const text of adversarial) {
      for (const read of [hoverLines, hoverRich]) {
        const t0 = performance.now();
        read([text]);
        expect(performance.now() - t0, `${read.name} on ${text.slice(0, 12).trim()}…`).toBeLessThan(50);
      }
    }
  });

  it('reads tags as before: br in its forms, a stray < or >, "||" with its spaces', () => {
    expect(hoverLines(['a< br />b<BR>c</ br>d'])).toEqual(['a', 'b', 'c', 'd']);
    expect(hoverLines(['1 < 2 and 3 > 2'])).toEqual(['1  2']);
    expect(hoverLines(['x <y <br> z'])).toEqual(['x <y', 'z']);
    expect(hoverLines(['a   ||   b || || c'])).toEqual(['a b  c']);
    expect(hoverRich(['3 <  4, <b>bold</b>'])).toEqual([
      { role: 'title', runs: [{ text: '3 < 4, ' }, { text: 'bold', bold: true }] },
    ]);
  });
});
