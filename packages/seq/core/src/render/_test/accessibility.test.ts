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
import { createSeqSpec } from '@molstar/seqx-core/lang';
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
    expect(root.getAttribute('role')).toBe('figure');
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

/** The accessible tree of an SVG: each element with a role, its name (`aria-label`, else the elements `aria-labelledby` names). */
function accessibleTree(doc: Document) {
  const out: { role: string; name: string; depth: number }[] = [];
  const nameOf = (e: Element) => {
    const label = e.getAttribute('aria-label');
    if (label) return label;
    const by = e.getAttribute('aria-labelledby');
    return by
      ? by
          .split(' ')
          .map((id) => doc.getElementById(id)?.textContent ?? '')
          .join('. ')
      : '';
  };
  const walk = (e: Element, depth: number) => {
    const role = e.getAttribute('role');
    if (role) out.push({ role, name: nameOf(e), depth });
    for (const c of Array.from(e.children)) walk(c, role ? depth + 1 : depth);
  };
  walk(doc.documentElement, 0);
  return out;
}

describe('what assistive technology reaches', () => {
  it('exposes the figure, its tracks as groups and its items and highlights as named images', async () => {
    const r = await render({ highlights: [{ location: rangesLocation('/data/p53', [10, 20]), color: 0xff0000ff }] });
    const tree = accessibleTree(parse(r.svg));
    expect(tree[0]).toEqual({ role: 'figure', name: 'p53 tracks. 1 section: protein (p53, 3 tracks)', depth: 0 });
    expect(tree.filter((n) => n.role === 'group').map((n) => n.name)).toEqual([
      'p53 (protein)',
      'Domains: 6 items (protein)',
      'Conservation (protein)',
    ]);
    const images = tree.filter((n) => n.role === 'img');
    expect(images.map((n) => n.name)).toContain('DNA-binding, p53 94–292');
    expect(images.map((n) => n.name)).toContain('Highlight 1: 10 positions, p53 11–20');
    // the items are inside their track's group
    expect(images.find((n) => n.name.startsWith('DNA-binding'))!.depth).toBe(2);
    // every role has a name
    expect(tree.every((n) => n.name.length > 0)).toBe(true);
  });

  it('never takes a pointer from the picture beneath: nothing in the layer is a hit target', async () => {
    const r = await render();
    expect(r.svg).not.toContain('pointer-events="all"');
    const doc = parse(r.svg);
    for (const rect of doc.querySelectorAll('[data-seq-accessible] rect')) {
      expect(rect.getAttribute('fill')).toBe('none');
      expect(rect.getAttribute('pointer-events')).toBe('none');
    }
  });

  it('tells two exports of one model on a page apart, and takes the host’s idPrefix', async () => {
    const a = await render();
    const b = await render({ highlights: [{ location: rangesLocation('/data/p53', [10, 20]), color: 0xff0000ff }] });
    const c = await renderSvg(p53TracksSpec(), { size: { width: 600, height: 140 } });
    const idOf = (svg: string) => parse(svg).querySelector('title')!.id;
    expect(new Set([idOf(a.svg), idOf(b.svg), idOf(c.svg)]).size).toBe(3);
    // equal exports are equal
    expect(idOf((await render()).svg)).toBe(idOf(a.svg));
    const own = await renderSvg(p53TracksSpec(), { size: { width: 800, height: 140 }, idPrefix: 'mine-' });
    expect(idOf(own.svg)).toBe('mine-title');
    expect(parse(own.svg).documentElement.getAttribute('aria-labelledby')).toBe('mine-title mine-desc');
  });
});

describe('the names of a large figure', () => {
  const Tracks = 14;
  const Items = 400;
  function manySpec() {
    const b = createSeqSpec({ title: 'many' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(2500) });
    const d = b.data.ranges('d', {
      target: s,
      rows: {
        start: Int32Array.from({ length: Items }, (_, i) => i * 24),
        end: Int32Array.from({ length: Items }, (_, i) => i * 24 + 12),
        id: Array.from({ length: Items }, (_, i) => `i${i}`),
        label: Array.from({ length: Items }, (_, i) => `Item ${i}`),
      },
    });
    b.views.view('main', (v) => {
      const sec = v.section('sec', { of: s });
      for (let t = 0; t < Tracks; t++)
        sec.track(`t${t}`, { label: `Track ${t}` }).layer('l', { type: 'ranges', data: d });
    });
    return b.build();
  }

  it('names at most a share of items a track and a figure’s worth in all, and says how many more', async () => {
    const r = await renderSvg(manySpec(), { size: { width: 4000, height: 600 } });
    const doc = parse(r.svg);
    const tracks = [...doc.querySelectorAll('[data-seq-track]')];
    expect(tracks).toHaveLength(Tracks);
    const per = tracks.map((t) => t.querySelectorAll('[data-seq-item]').length);
    expect(Math.max(...per)).toBeLessThanOrEqual(200);
    const total = per.reduce((a, b) => a + b, 0);
    expect(total).toBeGreaterThan(0);
    expect(total).toBeLessThanOrEqual(2000);
    expect(tracks[0].querySelector(':scope > title')!.textContent).toMatch(/more not named/);
    expect(doc.querySelector('desc')!.textContent).toMatch(/\d[\d,]* more items not named$/);
    // the growth of the export: the accessible layer against the picture without it
    const at = r.svg.indexOf('<g data-seq-accessible');
    const layer = r.svg.length - at - '</svg>'.length;
    expect(layer).toBeGreaterThan(0);
    // stated in the report: bytes the layer adds, per named item, stays small
    expect(layer / Math.max(total, 1)).toBeLessThan(600);
  });
});
