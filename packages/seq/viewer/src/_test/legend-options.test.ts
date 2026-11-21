// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Legend part's options (_spec/README.md "Legends"): `of` as a list, `omit`, `orient`, the compact form's entries,
// keys and accessible names, and the layer's `legend: none` read by the rows and the part.

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { mountLegend } from '../mount-legend.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
beforeEach(stubCanvasContext);
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

const Domains = '/views/main/p/domains';
const Sites = '/views/main/p/sites';
const Score = '/views/main/p/score';

function spec(sitesLegend?: 'auto' | 'none') {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQ' });
  const d = b.data.ranges('d', {
    target: p,
    rows: { start: [0, 6], end: [5, 12], id: ['A', 'B'], category: ['kinase', 'linker'] },
  });
  const values = b.data.values('v', { target: p, values: Float32Array.from({ length: 16 }, (_, i) => i) });
  const sc = b.views.scale('s', { label: 'Score' });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
    s.track('sites', { label: 'Sites' }).layer('d', {
      type: 'ranges',
      data: d,
      color: theme.field('id'),
      ...(sitesLegend ? { legend: sitesLegend } : {}),
    });
    s.track('score', { label: 'Score' }).layer('v', {
      type: 'values',
      data: values,
      scale: sc,
      color: theme.field('value', { scale: sc }),
    });
  });
  return b.build();
}

async function mounted(options: Parameters<typeof mountLegend>[2], sitesLegend?: 'auto' | 'none') {
  frames = installFrames();
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  await viewer.load(spec(sitesLegend));
  viewer.resize(500);
  await viewer.settle();
  frames.flush();
  const el = document.createElement('div');
  document.body.append(el);
  return { el, dispose: mountLegend(el, viewer, options) };
}

const titles = (el: HTMLElement) => [...el.querySelectorAll('[role="group"]')].map((g) => g.getAttribute('aria-label'));

describe('the legend selection', () => {
  it('`of` is an address or a list; `omit` removes after it; the default is every legend', async () => {
    await mounted({});
    const all = legendRows(viewer!).map((r) => r.title);
    expect(all).toEqual(expect.arrayContaining(['Domains', 'Sites', 'Score']));
    expect(legendRows(viewer!, { of: Domains }).map((r) => r.title)).toEqual(['Domains']);
    expect(legendRows(viewer!, { of: [Domains, Score] }).map((r) => r.title)).toEqual(['Domains', 'Score']);
    expect(legendRows(viewer!, { omit: [Domains, Sites] }).map((r) => r.title)).toEqual(
      all.filter((t) => t !== 'Domains' && t !== 'Sites'),
    );
    expect(legendRows(viewer!, { of: [Domains, Sites], omit: [Sites] }).map((r) => r.title)).toEqual(['Domains']);
    // a layer's address and a scale's name the same legends as their track's
    expect(legendRows(viewer!, { of: `${Domains}/d` }).map((r) => r.title)).toEqual(['Domains']);
    expect(legendRows(viewer!, { of: [] })).toEqual([]);
  });

  it('the part follows `of` as a list and `omit`', async () => {
    const a = await mounted({ of: [Domains, Sites] });
    expect(titles(a.el)).toEqual(['Domains', 'Sites']);
    a.dispose();
    const b = await mounted({ omit: [Domains, Sites] });
    expect(titles(b.el)).not.toContain('Domains');
    expect(titles(b.el)).toContain('Score');
  });

  it('a layer’s `legend: none` is read by the rows and the part; `auto` keeps it', async () => {
    const { el } = await mounted({}, 'none');
    expect(legendRows(viewer!).map((r) => r.title)).not.toContain('Sites');
    expect(titles(el)).not.toContain('Sites');
    expect(titles(el)).toContain('Domains');
    viewer!.retire();
    frames!.restore();
    const auto = await mounted({}, 'auto');
    expect(titles(auto.el)).toContain('Sites');
  });
});

describe('the compact form', () => {
  it('a line a legend: its title then unlabelled swatches, each named, with the same roving focus and gestures', async () => {
    const { el } = await mounted({ form: 'compact', of: Domains });
    const group = el.querySelector<HTMLElement>('[role="group"]')!;
    expect(group.style.flexWrap).toBe('nowrap');
    const options = [...group.querySelectorAll<HTMLElement>('[role="option"]')];
    expect(options.map((o) => o.getAttribute('aria-label'))).toEqual(['Domains: kinase', 'Domains: linker']);
    // the label is on hover, never text beside the swatch
    expect(options.map((o) => o.title)).toEqual(['Domains: kinase', 'Domains: linker']);
    expect(options.map((o) => o.textContent)).toEqual(['', '']);
    options[1].click();
    expect(viewer!.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['B']);
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    box.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight' }));
    expect(options.map((o) => o.getAttribute('aria-selected'))).toEqual(['false', 'true']);
    expect(viewer!.highlight.names.value).toContain('legend');
  });

  it('a ramp is a short bar with its domain’s ends, named by its domain', async () => {
    const { el } = await mounted({ form: 'compact', of: Score });
    const ramp = el.querySelector<HTMLElement>('[role="option"]')!;
    expect(ramp.getAttribute('aria-label')).toMatch(/^Score: 0 to 15/);
    expect(ramp.textContent).toBe('015');
    const bar = ramp.querySelector<HTMLElement>('span')!;
    expect(bar.style.width).toBe('40px');
    expect(bar.style.background).toContain('linear-gradient');
  });

  it('the full form keeps its labels; `orient` sets the direction; unknown forms and orientations throw', async () => {
    const full = await mounted({ of: Domains });
    expect([...full.el.querySelectorAll('[role="option"]')].map((o) => o.textContent)).toEqual(['kinase', 'linker']);
    const box = full.el.querySelector<HTMLElement>('[role="listbox"]')!;
    expect(box.style.flexDirection).toBe('column');
    full.dispose();
    const h = await mounted({ orient: 'horizontal' });
    const row = h.el.querySelector<HTMLElement>('[role="listbox"]')!;
    expect(row.style.flexFlow).toBe('row wrap');
    expect(row.dataset.orient).toBe('horizontal');
    const bad = document.createElement('div');
    expect(() => mountLegend(bad, viewer!, { form: 'tiny' as never })).toThrow(RangeError);
    expect(() => mountLegend(bad, viewer!, { orient: 'diagonal' as never })).toThrow(RangeError);
  });
});

describe('the part’s focus', () => {
  it('items have ids and aria-activedescendant follows the arrows; compact shows the focused name while focused', async () => {
    const { el } = await mounted({ form: 'compact', of: Domains });
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    const options = [...el.querySelectorAll<HTMLElement>('[role="option"]')];
    expect(new Set(options.map((o) => o.id)).size).toBe(2);
    expect(options.every((o) => o.id !== '')).toBe(true);
    expect(box.getAttribute('aria-activedescendant')).toBe(options[0].id);
    const shown = () => el.querySelector('[data-focus-name]')!.textContent;
    expect(shown()).toBe('');
    box.focus();
    expect(shown()).toBe('Domains: kinase');
    box.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight' }));
    expect(box.getAttribute('aria-activedescendant')).toBe(options[1].id);
    expect(shown()).toBe('Domains: linker');
    box.blur();
    expect(shown()).toBe('');
  });

  it('past maxItems the compact group says +n; an `of` naming no layer is an empty part', async () => {
    const { el } = await mounted({ form: 'compact', of: Domains, maxItems: 1 });
    expect([...el.querySelectorAll('[role="group"] > span')].map((s) => s.textContent)).toContain('+1');
    expect(legendRows(viewer!, { of: '/views/main/nowhere' })).toEqual([]);
  });
});
