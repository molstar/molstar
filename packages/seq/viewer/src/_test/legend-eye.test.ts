// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's eye in the Legend part (_spec/README.md "Legends", A layer's eye), a test a gesture: the eye after a
// legend's title hides its layer (`layer.hide`, one patch), the row staying in its place dimmed with the eye-off that
// shows it again; a hidden layer's entries have no action; undo is one step; the eye's tooltip is its name; a legend
// over several layers has one button opening a menu of them; `layers: false` drops it. The header's "Hidden layers"
// and the eye read one state (two consumers in one test).

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { headerMenu } from '../header-entries.js';
import { legendRows } from '../legends.js';
import { loadLegendEyes, mountLegend } from '../mount-legend.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/p';
const Domains = `${Section}/domains`;
const DomainsLayer = `${Domains}/d`;
const Pair = `${Section}/pair`;

/** Domains (one layer, by category) above Pair (two layers of one data and field: one legend over both). */
function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: new Int32Array([0, 6, 14]),
      end: new Int32Array([5, 12, 20]),
      id: ['A', 'B', 'C'],
      category: ['binding', 'linker', 'binding'],
      family: ['x', 'y', 'x'],
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
    const t = s.track('pair', { label: 'Pair' });
    t.layer('a', { type: 'ranges', data: d, color: theme.field('family') });
    t.layer('b', { type: 'ranges', data: d, color: theme.field('family') });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
const disposers: (() => void)[] = [];
afterEach(() => {
  for (const d of disposers.splice(0)) d();
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
  document.body.replaceChildren();
});

async function mounted(options: Parameters<typeof mountLegend>[2] = {}) {
  stubCanvasContext();
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, history: uiHistory() });
  viewer = v;
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  // the eyes' chunk landed first (a lazy import: real I/O, never the frames' clock), and its tooltip installed after
  await loadLegendEyes();
  const el = document.createElement('div');
  document.body.append(el);
  disposers.push(mountLegend(el, v, options));
  await loadLegendEyes();
  return { v, el };
}

/** A command's patch landed and the view settled again (its legends derived). */
async function landed(v: SeqViewer, r?: unknown) {
  await r;
  await v.settle();
  frames!.flush();
}

const group = (el: HTMLElement, title: string) =>
  [...el.querySelectorAll<HTMLElement>('[role="group"]')].find((g) => g.getAttribute('aria-label') === title)!;
const eye = (el: HTMLElement, title: string) =>
  group(el, title).querySelector<HTMLButtonElement>('[data-legend-eye] button')!;
const titles = (el: HTMLElement) =>
  [...el.querySelectorAll<HTMLElement>('[role="group"]')].map((g) => g.getAttribute('aria-label'));
const stack = (v: SeqViewer) => (v.spec.history as UiHistory).stack;
const hiddenParam = (v: SeqViewer, address: string) => {
  let n = v.spec.root.value;
  for (const id of address.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === id)!;
  return n.params?.hidden;
};

describe('a layer’s eye (the agreed gestures)', () => {
  it('Click the eye: layer.hide; the row dims; the button reads "Show Domains"', async () => {
    const { v, el } = await mounted();
    const before = titles(el);
    expect(eye(el, 'Domains').getAttribute('aria-label')).toBe('Hide Domains');
    expect(group(el, 'Domains').dataset.hidden).toBeUndefined();
    eye(el, 'Domains').click();
    await landed(v);
    expect(hiddenParam(v, DomainsLayer)).toBe(true);
    // the row stays in its place, dimmed, the button its eye-off
    expect(titles(el)).toEqual(before);
    const g = group(el, 'Domains');
    expect(g.dataset.hidden).toBe('true');
    expect(g.querySelector<HTMLElement>('span')!.style.opacity).toBe('0.45');
    for (const o of g.querySelectorAll<HTMLElement>('[role="option"]')) expect(o.style.opacity).toBe('0.45');
    expect(eye(el, 'Domains').getAttribute('aria-label')).toBe('Show Domains');
    // two consumers, one state: the track header's "Hidden layers" lists the layer the eye hid
    const hiddenLayers = headerMenu(v, { section: Section, track: Domains })!
      .groups.flat()
      .find((i) => i.name === 'seq:hidden-layers')!;
    expect(hiddenLayers.items!.map((i) => [i.label, i.command, i.args])).toEqual([
      ['Show d', 'layer.show', { address: DomainsLayer }],
    ]);
    expect(legendRows(v).find((r) => r.title === 'Domains')!.layers).toEqual([
      { address: DomainsLayer, name: 'Domains', hidden: true },
    ]);
  });

  it('Click the eye-off: layer.show; the row back as it was', async () => {
    const { v, el } = await mounted();
    const before = group(el, 'Domains').outerHTML;
    eye(el, 'Domains').click();
    await landed(v);
    eye(el, 'Domains').click();
    await landed(v);
    expect(hiddenParam(v, DomainsLayer)).toBeUndefined();
    // the same row, element for element (its ids aside: the part numbers its options on each draw)
    const ids = (html: string) => html.replace(/id="[^"]*"/g, '');
    expect(ids(group(el, 'Domains').outerHTML)).toBe(ids(before));
    expect(eye(el, 'Domains').getAttribute('aria-label')).toBe('Hide Domains');
  });

  it('Tab to the button, Enter or Space: as a click (a native button, in the tab order; its keys never select an entry)', async () => {
    const { v, el } = await mounted();
    const b = eye(el, 'Domains');
    expect(b.tagName).toBe('BUTTON');
    expect(b.type).toBe('button');
    expect(b.tabIndex).toBe(0);
    expect(b.disabled).toBe(false);
    b.focus();
    for (const key of ['Enter', ' ']) b.dispatchEvent(new KeyboardEvent('keydown', { key, bubbles: true }));
    expect(v.selection.value.value.elements.keys).toEqual([]);
    // the browser's activation of a focused button by Enter or Space is its click
    b.click();
    await landed(v);
    expect(hiddenParam(v, DomainsLayer)).toBe(true);
    // the focus stays on the eye across the redraw (its new button)
    expect(document.activeElement).toBe(eye(el, 'Domains'));
  });

  it('Click an entry of a shown layer: unchanged, selects its items (the click rule)', async () => {
    const { v, el } = await mounted();
    const linker = group(el, 'Domains').querySelector<HTMLElement>('[role="option"][data-value="linker"]')!;
    expect(linker.style.cursor).toBe('pointer');
    linker.click();
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['B']);
  });

  it('Click an entry of a hidden layer: nothing; the pointer says so (not-allowed)', async () => {
    const { v, el } = await mounted();
    eye(el, 'Domains').click();
    await landed(v);
    const options = [...group(el, 'Domains').querySelectorAll<HTMLElement>('[role="option"]')];
    expect(options.map((o) => o.getAttribute('aria-label'))).toEqual(['Domains: binding', 'Domains: linker']);
    for (const o of options) {
      expect(o.style.cursor).toBe('not-allowed');
      expect(o.getAttribute('aria-disabled')).toBe('true');
      expect(o.dataset.value).toBeUndefined();
      o.click();
      o.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    }
    expect(v.selection.value.value.elements.keys).toEqual([]);
    expect(v.highlight.names.value).not.toContain('legend');
    // the keyboard too: Enter on a hidden entry selects nothing
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    box.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
    expect(v.selection.value.value.elements.keys).toEqual([]);
  });

  it('Undo: the hide or show undone, one step', async () => {
    const { v, el } = await mounted();
    const steps = stack(v).steps.value.length;
    eye(el, 'Domains').click();
    await landed(v);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    await landed(v, stack(v).undo());
    expect(hiddenParam(v, DomainsLayer)).toBeUndefined();
    expect(group(el, 'Domains').dataset.hidden).toBeUndefined();
    expect(eye(el, 'Domains').getAttribute('aria-label')).toBe('Hide Domains');
    // and the show: one step, undone to hidden
    eye(el, 'Domains').click();
    await landed(v);
    eye(el, 'Domains').click();
    await landed(v);
    await landed(v, stack(v).undo());
    expect(hiddenParam(v, DomainsLayer)).toBe(true);
    expect(group(el, 'Domains').dataset.hidden).toBe('true');
  });

  it('Hover the eye: a tooltip with the button’s name', async () => {
    const { el } = await mounted();
    const b = eye(el, 'Domains');
    expect(b.dataset.tip).toBe('Hide Domains');
    expect(b.getAttribute('aria-label')).toBe('Hide Domains');
    // the part installs the tooltip its buttons' `data-tip` shows
    expect(el.querySelector('[role="tooltip"]')).not.toBeNull();
  });
});

describe('a legend over several layers', () => {
  it('has one button opening a menu of its layers, each checked while shown; a hidden one keeps the row, undimmed', async () => {
    const { v, el } = await mounted();
    const row = legendRows(v).find((r) => r.title === 'Pair')!;
    expect(row.layers.map((l) => l.name)).toEqual(['Pair · a', 'Pair · b']);
    const b = eye(el, 'Pair');
    expect(group(el, 'Pair').querySelectorAll('[data-legend-eye]')).toHaveLength(1);
    expect(b.getAttribute('aria-label')).toBe('Layers of Pair');
    expect(b.getAttribute('aria-haspopup')).toBe('menu');
    b.click();
    const menu = await vi.waitFor(() => {
      const m = document.querySelector<HTMLElement>('[role="menu"]');
      if (!m) throw new Error('no menu yet');
      return m;
    });
    const rows = [...menu.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
    expect(rows.map((r) => [r.textContent?.trim(), r.getAttribute('aria-checked')])).toEqual([
      ['Pair · a', 'true'],
      ['Pair · b', 'true'],
    ]);
    rows[1].click();
    await landed(v);
    expect(hiddenParam(v, `${Pair}/b`)).toBe(true);
    // one layer of two hidden: the row stays, not dimmed, its entries still select
    const after = legendRows(v).find((r) => r.title === 'Pair')!;
    expect(after.hidden).toBe(false);
    expect(after.layers.map((l) => [l.address, l.hidden])).toEqual([
      [`${Pair}/a`, false],
      [`${Pair}/b`, true],
    ]);
    expect(group(el, 'Pair').dataset.hidden).toBeUndefined();
  });
});

describe('the host’s option', () => {
  it('`layers: false` drops the eye; the rows and their entries as with it', async () => {
    const { el } = await mounted({ layers: false });
    expect(el.querySelectorAll('[data-legend-eye]')).toHaveLength(0);
    expect(titles(el)).toEqual(['Domains', 'Pair']);
    expect(el.querySelector('[role="tooltip"]')).toBeNull();
  });
});
