// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The exported SVG's names (_spec/README.md "Render", Accessibility): the root an image labelled by its title and
// description, a title on every track, item and highlight, checked as axe does (every graphic named, every
// `aria-labelledby` naming something that exists), and the words the viewer's live region shares.

import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { rangesText, renderSvg, spanText } from '../index.js';
import { p53TracksSpec } from './fixtures.js';

const parse = (svg: string) => new DOMParser().parseFromString(svg, 'image/svg+xml');

const render = (emphasis?: Parameters<typeof renderSvg>[1]['emphasis']) =>
  renderSvg(p53TracksSpec(), { size: { width: 800, height: 140 }, ...(emphasis ? { emphasis } : {}) });

describe('the exported SVG’s names', () => {
  it('is an image labelled by its title and description, which exist', async () => {
    const r = await render();
    const doc = parse(r.svg);
    const root = doc.documentElement;
    expect(doc.querySelector('parsererror')).toBeNull();
    expect(root.getAttribute('role')).toBe('img');
    const ids = (root.getAttribute('aria-labelledby') ?? '').split(' ');
    expect(ids).toHaveLength(2);
    const [title, desc] = ids.map((id) => doc.getElementById(id));
    expect(title?.tagName).toBe('title');
    expect(title?.textContent).toBe('p53 tracks');
    expect(desc?.tagName).toBe('desc');
    expect(desc?.textContent).toBe('1 section: protein (p53, 3 tracks)');
  });

  it('puts a title on each track and each item: its label and extent', async () => {
    const r = await render();
    const doc = parse(r.svg);
    const tracks = [...doc.querySelectorAll('[data-seq-track]')];
    expect(tracks.map((t) => t.querySelector(':scope > title')?.textContent)).toEqual([
      'p53 (protein)',
      'Domains: 6 items (protein)',
      'Conservation (protein)',
    ]);
    const items = [...doc.querySelectorAll('[data-seq-item]')].map(
      (i) => i.querySelector(':scope > title')?.textContent,
    );
    expect(items).toContain('Transactivation, p53 1–61');
    expect(items).toContain('DNA-binding, p53 94–292');
    expect(items).toHaveLength(6);
  });

  it('labels every graphic: no box of the accessible layer is without a title, each highlight has one', async () => {
    const r = await render({
      highlights: [{ location: rangesLocation('/data/p53', [10, 20]), color: 0xff0000ff }],
    });
    const doc = parse(r.svg);
    const rects = [...doc.querySelectorAll('[data-seq-accessible] rect')];
    expect(rects.length).toBeGreaterThan(6);
    for (const rect of rects) {
      expect(rect.parentElement!.querySelector(':scope > title')?.textContent?.length ?? 0).toBeGreaterThan(0);
      expect(Number(rect.getAttribute('width'))).toBeGreaterThanOrEqual(0);
    }
    const highlights = [...doc.querySelectorAll('[data-seq-highlight] > title')].map((t) => t.textContent);
    expect(highlights).toEqual(['Highlight 1: 10 positions, p53 11–20']);
  });
});

describe('the words', () => {
  it('say positions 1-based and inclusive', () => {
    expect(spanText(0, 1)).toBe('1');
    expect(spanText(11, 20)).toBe('12–20');
    expect(spanText(1199, 1450)).toBe('1,200–1,450');
    expect(rangesText([0, 3, 5, 6, 9, 10, 20, 30], 2)).toBe('1–3, 6, 2 more');
  });
});
