// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Legend` part in a DOM (jsdom): placed by the host for one track (`of`), its entries options with accessible
// names, a click selecting a value's items, the keyboard's roving focus hovering them, a highlight's row clearing it,
// and nothing left behind when it is taken away.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { channelsRow } from '../legends.js';
import { mountLegend } from '../mount-legend.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQ' });
  const d = b.data.ranges('d', {
    target: p,
    rows: { start: [0, 6], end: [5, 12], id: ['A', 'B'], category: ['kinase', 'linker'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
  });
  return b.build();
}

describe('the Legend part', () => {
  it('shows one track’s legend, selects by click, hovers by keyboard, clears a highlight, and goes away whole', async () => {
    frames = installFrames();
    viewer = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
    });
    await viewer.load(spec());
    viewer.resize(500);
    await viewer.settle();
    frames.flush();
    const el = document.createElement('div');
    document.body.append(el);
    const dispose = mountLegend(el, viewer, { of: '/views/main/p/domains' });
    const options = [...el.querySelectorAll('[role="option"]')];
    expect(options.map((o) => o.getAttribute('aria-label'))).toEqual(['Domains: kinase', 'Domains: linker']);
    (options[1] as HTMLElement).click();
    expect(viewer.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['B']);
    const box = el.querySelector('[role="listbox"]') as HTMLElement;
    box.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight' }));
    expect(viewer.highlight.names.value).toContain('legend');
    const clear = [...el.querySelectorAll('button')].find(
      (b) => b.getAttribute('aria-label') === 'Clear the legend highlight',
    );
    clear!.click();
    expect(viewer.highlight.names.value).not.toContain('legend');
    dispose();
    expect(el.children.length).toBe(0);
  });

  async function mounted(options: Parameters<typeof mountLegend>[2]) {
    frames = installFrames();
    viewer = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
    });
    await viewer.load(spec());
    viewer.resize(500);
    await viewer.settle();
    frames.flush();
    const el = document.createElement('div');
    document.body.append(el);
    return { el, dispose: mountLegend(el, viewer, options) };
  }
  const slot = (el: HTMLElement) => el.querySelector<HTMLElement>('[data-channels]');

  it('groups one a line: a column, each title muted, 600, 11 px, 12 px between', async () => {
    const { el } = await mounted({});
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    expect(box.style.flexDirection).toBe('column');
    expect(box.style.gap).toBe('12px');
    const title = el.querySelector<HTMLElement>('[role="group"] > span')!;
    expect(title.style.fontWeight).toBe('600');
    expect(title.style.fontSize).toBe('11px');
    expect(title.style.color).toContain('--msui-text-muted');
  });

  it('`reserve` holds the channels row whether or not a channel is active; `collapse` grows as before', async () => {
    const reserved = await mounted({});
    const empty = slot(reserved.el)!;
    expect(empty.style.minHeight).toBe('24px');
    expect(empty.children.length).toBe(0);
    await viewer!.execute('legend.hover', { address: '/views/main/p/domains/d', entry: 'linker' });
    const lit = slot(reserved.el)!;
    expect(viewer!.highlight.names.value).toContain('legend');
    expect(lit.style.minHeight).toBe('24px');
    expect(lit.querySelector('button[data-clear]')).not.toBeNull();
    reserved.dispose();
    viewer!.retire();
    frames!.restore();
    const collapsed = await mounted({ channels: 'collapse' });
    expect(slot(collapsed.el)).toBeNull();
    await viewer!.execute('legend.hover', { address: '/views/main/p/domains/d', entry: 'linker' });
    expect(slot(collapsed.el)).not.toBeNull();
  });

  it('the row model: reserve is the same with none or some, collapse only with some', () => {
    expect(channelsRow('reserve', 0)).toEqual(channelsRow('reserve', 2));
    expect(channelsRow('collapse', 0)).toEqual({ present: false, minHeight: 0 });
    expect(channelsRow('collapse', 1).present).toBe(true);
  });

  it('refuses an unknown channels mode', async () => {
    frames = installFrames();
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    expect(() => mountLegend(document.createElement('div'), viewer!, { channels: 'hold' as never })).toThrow(
      RangeError,
    );
  });
});
