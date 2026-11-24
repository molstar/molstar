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
import { afterAll, afterEach, beforeEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { channelsRow } from '../legends.js';
import { loadLegendEyes, mountLegend } from '../mount-legend.js';
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

// The legend's eyes are a lazy chunk the first mount in a worker starts: the file waits for it, or it lands after the file ends.
afterAll(() => loadLegendEyes());

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

  it('a hover survives no redraw it should not: the entry under the pointer is replaced by the highlight it brings, and leaving the part still clears it', async () => {
    const { el } = await mounted({ of: '/views/main/p/domains' });
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    const entry = el.querySelector<HTMLElement>('[role="option"][data-value="linker"]')!;
    // as a browser sends them: pointerenter on the entry (no bubbling) and pointerover (bubbling)
    entry.dispatchEvent(new PointerEvent('pointerenter'));
    entry.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(viewer!.highlight.names.value).toContain('legend');
    // the highlight's row redrew the legend: the entry under the pointer is a new element, the old one detached
    expect(entry.isConnected).toBe(false);
    const again = el.querySelector<HTMLElement>('[role="option"][data-value="linker"]')!;
    expect(again).not.toBe(entry);
    // the pointer moves to another entry: that one is hovered
    el.querySelector<HTMLElement>('[role="option"][data-value="kinase"]')!.dispatchEvent(
      new PointerEvent('pointerover', { bubbles: true }),
    );
    expect(viewer!.highlight.names.value).toContain('legend');
    // and out of the part: always cleared, whatever was replaced in between (the removed node fires nothing)
    box.dispatchEvent(new PointerEvent('pointerleave'));
    expect(viewer!.highlight.names.value).not.toContain('legend');
  });

  async function mounted(
    options: Parameters<typeof mountLegend>[2],
    config: NonNullable<Parameters<typeof createSeqViewer>[0]['config']> = [],
  ) {
    frames = installFrames();
    viewer = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000], ...config],
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

  it('an entry’s action follows what is drawn: a pan that leaves its items out of view takes its pointer away, and an entry with no action hovers nothing', async () => {
    // a jump is a cut (no motion to wait out)
    const { el } = await mounted({ of: '/views/main/p/domains' }, [[viewerItem('view.motion'), 'cut']]);
    const cursors = () => [...el.querySelectorAll<HTMLElement>('[role="option"]')].map((o) => o.style.cursor);
    // only what is drawn names items (no address does)
    viewer!.legends.addressable = () => false;
    viewer!.legends.setCvdSafe(true);
    expect(cursors()).toEqual(['pointer', 'pointer']);
    const focus = async (id: string) => {
      expect(viewer!.runCommand('address.focus', { address: `/data/d/(item :id ${id})` })).toBe('Focused');
      // the view rests (its settle timer), and the settled frame is drawn
      for (let k = 0; k < 3; k++) {
        await viewer!.settle();
        frames!.timers.fire();
        frames!.flush();
      }
    };
    // focused on B: A (0–5) out of view, its entry no target
    await focus('B');
    expect(cursors()).toEqual(['default', 'pointer']);
    // then on A: the legends are the same (nothing republished), what is drawn moved, the actions with it
    const fired: string[] = [];
    const off = viewer!.legends.legends.listen(() => fired.push('legends'));
    await focus('A');
    off();
    expect(fired).toEqual([]);
    expect(cursors()).toEqual(['pointer', 'default']);
    // no address names B's items and none is drawn: its entry carries no value and hovers nothing
    const linker = el.querySelectorAll<HTMLElement>('[role="option"]')[1];
    expect(linker.dataset.value).toBeUndefined();
    linker.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(viewer!.highlight.names.value).not.toContain('legend');
  });

  it('an entry’s items are read once a settled frame: a hover again makes nothing, a pan reads them again', async () => {
    await mounted({ of: '/views/main/p/domains' }, [[viewerItem('view.motion'), 'cut']]);
    const legends = viewer!.legends;
    const first = legends.entryItems('/views/main/p/domains/d', 'linker');
    expect(new Set(first?.elements.map((k) => k.key))).toEqual(new Set(['B']));
    expect(legends.entryItems('/views/main/p/domains/d', 'linker')).toBe(first);
    expect(viewer!.runCommand('address.focus', { address: '/data/d/(item :id A)' })).toBe('Focused');
    for (let k = 0; k < 3; k++) {
      await viewer!.settle();
      frames!.timers.fire();
      frames!.flush();
    }
    expect(legends.entryItems('/views/main/p/domains/d', 'linker')).toBeUndefined();
  });

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
    // empty, it holds the line a channel's row takes: a hidden copy of one (not announced, not a tab stop)
    expect(empty.querySelectorAll('[data-highlight]').length).toBe(0);
    const copy = empty.firstElementChild as HTMLElement;
    expect(copy.getAttribute('aria-hidden')).toBe('true');
    expect(copy.style.visibility).toBe('hidden');
    expect(copy.querySelector('button')!.tabIndex).toBe(-1);
    await viewer!.execute('legend.hover', { address: '/views/main/p/domains/d', entry: 'linker' });
    const lit = slot(reserved.el)!;
    expect(viewer!.highlight.names.value).toContain('legend');
    expect(lit.style.flexWrap).toBe('nowrap');
    expect(lit.children.length).toBe(1);
    expect(lit.querySelector('button[data-clear]')).not.toBeNull();
    reserved.dispose();
    viewer!.retire();
    frames!.restore();
    const collapsed = await mounted({ channels: 'collapse' });
    expect(slot(collapsed.el)).toBeNull();
    await viewer!.execute('legend.hover', { address: '/views/main/p/domains/d', entry: 'linker' });
    expect(slot(collapsed.el)).not.toBeNull();
  });

  it('the row is outside the listbox, one line for any number of channels: three at a narrow width show two and a +1', async () => {
    const { el } = await mounted({ channels: 'reserve' });
    el.style.width = '120px';
    const box = el.querySelector<HTMLElement>('[role="listbox"]')!;
    expect(box.querySelector('[data-channels]')).toBeNull();
    // the channels' Clear buttons outside it (the legends' own buttons, an eye and a "+n", are in their groups)
    expect(box.querySelector('button[data-clear]')).toBeNull();
    for (const name of ['a', 'b', 'c']) viewer!.runLink(`!seq-highlight=/data/d&seq-channel=${name}`);
    expect(viewer!.highlight.names.value).toEqual(['a', 'b', 'c']);
    frames!.flush();
    const row = slot(el)!;
    expect(row.contains(box)).toBe(false);
    expect(row.style.flexWrap).toBe('nowrap');
    expect(row.style.overflow).toBe('hidden');
    expect(row.querySelectorAll('[data-highlight]').length).toBe(2);
    expect(row.querySelector('[data-more]')?.textContent).toBe('+1');
  });

  it('the row model: reserve is the same with none or some, collapse only with some', () => {
    expect(channelsRow('reserve', 0)).toEqual(channelsRow('reserve', 2));
    expect(channelsRow('collapse', 0)).toEqual({ present: false });
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
