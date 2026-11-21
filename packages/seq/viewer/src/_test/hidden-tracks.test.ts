// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { railSections } from '@molstar/seqx-core/render';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { mountHiddenTracks } from '../mount-hidden-tracks.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The hidden tracks the application lists (_spec/README.md "Hidden tracks"): `viewer.hiddenTracks` from the spec in its
// order with its change event, `track.show` of one, several or all in one patch, and `mountHiddenTracks`'s button.

const A = '/views/main/one/a';
const B = '/views/main/one/b';
const C = '/views/main/two/c';

function spec() {
  const b = createSeqSpec({ title: 'hidden tracks' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  b.views.view('main', (v) => {
    const one = v.section('one', { of: p });
    one.track('a', { label: 'Alpha' }).layer('l', { type: 'letters' });
    one.track('b').layer('l', { type: 'letters' });
    const two = v.section('two', { of: p });
    two.track('c', { label: 'Gamma' }).layer('l', { type: 'letters' });
    two.track('d', { label: 'Delta' }).layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames;
let viewer: SeqViewer | undefined;
beforeEach(() => {
  frames = installFrames();
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  Object.defineProperty(document.defaultView!, 'matchMedia', {
    configurable: true,
    writable: true,
    value: () => ({ matches: false, addEventListener() {}, removeEventListener() {} }),
  });
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  Reflect.deleteProperty(document.defaultView!, 'matchMedia');
  document.body.innerHTML = '';
});

async function make() {
  const v = (viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    history: uiHistory(),
    config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
  }));
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}
const steps = (v: SeqViewer) => (v.spec.history as UiHistory).stack.steps.value.length;

describe('viewer.hiddenTracks', () => {
  it('lists each section’s hidden tracks in spec order (address, label, section), told on a change only', async () => {
    const v = await make();
    expect(v.hiddenTracks.value).toEqual([]);
    const told: number[] = [];
    v.hiddenTracks.listen((l) => told.push(l.length));
    await v.execute('track.hide', { address: C });
    await v.execute('track.hide', { address: B });
    expect(v.hiddenTracks.value).toEqual([
      { address: B, label: 'b', section: '/views/main/one' },
      { address: C, label: 'Gamma', section: '/views/main/two' },
    ]);
    // an edit that changes no hidden track tells nothing
    await v.execute('track.hide', { address: B });
    await v.execute('layer.hide', { address: `${A}/l` });
    expect(told).toEqual([1, 2]);
  });

  it('a section with every track hidden is listed (it has no track left in the resolved view)', async () => {
    const v = await make();
    await v.execute('track.hide', { address: C });
    await v.execute('track.hide', { address: '/views/main/two/d' });
    await v.settle();
    expect(v.view.resolved.value!.sections[1].tracks).toEqual([]);
    expect(v.hiddenTracks.value.map((t) => t.section)).toEqual(['/views/main/two', '/views/main/two']);
  });
});

describe('a section with every track hidden collapses fully', () => {
  async function threeSections() {
    const v = viewer ?? (await make());
    const b = createSeqSpec({ title: 'three sections' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
    b.views.view('main', (w) => {
      w.section('one', { of: p }).track('a').layer('l', { type: 'letters' });
      const two = w.section('two', { of: p });
      two.track('c').layer('l', { type: 'letters' });
      two.track('d').layer('l', { type: 'letters' });
      w.section('three', { of: p }).track('e').layer('l', { type: 'letters' });
    });
    await v.load(b.build());
    await v.settle();
    frames.flush();
    return v;
  }

  it('no gap, no line, no rail, no cursor stop: the sections around it close up', async () => {
    const v = await threeSections();
    const gapped = v.frame!.layout.sections.map((s) => [s.y, s.height]);
    expect(gapped[1][1]).toBeGreaterThan(0);
    await v.execute('track.hide', { address: '/views/main/two/c' });
    await v.execute('track.hide', { address: '/views/main/two/d' });
    await v.settle();
    frames.flush();
    const [one, two, three] = v.frame!.layout.sections;
    const config = v.frame!.layout.config;
    expect([two.height, two.header]).toEqual([0, undefined]);
    // the gap is the one between the two sections that remain
    expect(three.y).toBe(one.y + one.height + config.sectionGap);
    // the navigator rails skip it
    expect(railSections(v.view.resolved.value!)).toEqual(['/views/main/one', '/views/main/three']);
    // the cursor goes from one to three
    v.selection.setCursor({ section: '/views/main/one', system: '/data/p', position: 5 });
    await v.execute('cursor.down', { section: '/views/main/one' });
    expect(v.selection.cursor.value?.section).toBe('/views/main/three');
    await v.execute('cursor.up', { section: '/views/main/three' });
    expect(v.selection.cursor.value?.section).toBe('/views/main/one');
  });
});

describe('a jump into a hidden track', () => {
  it('says "on a hidden track" and moves nothing; a shown track is focused as before', async () => {
    const v = await make();
    for (const a of [C, '/views/main/two/d']) await v.execute('track.hide', { address: a });
    await v.settle();
    const out = v.addresses.focus(C);
    expect(out).toMatchObject({ done: false, message: 'On a hidden track: show it to see it' });
    expect(v.addresses.reveal('/views/main/two/d/l').message).toBe('On a hidden track: show it to see it');
    expect(v.addresses.focus(A).message).not.toBe('On a hidden track: show it to see it');
    await v.execute('track.show', { all: true });
    await v.settle();
    expect(v.addresses.focus(C).message).not.toBe('On a hidden track: show it to see it');
  });
});

describe('a jump to an item whose key is drawn only on hidden tracks', () => {
  async function withItems() {
    const v = viewer ?? (await make());
    const b = createSeqSpec({ title: 'items' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
    const d = b.data.ranges('d', {
      target: p,
      rows: { start: new Int32Array([0, 12]), end: new Int32Array([10, 18]), id: ['A', 'B'] },
    });
    b.views.view('main', (w) => {
      const s = w.section('one', { of: p });
      s.track('x', { label: 'X' }).layer('l', { type: 'ranges', data: d });
      s.track('y', { label: 'Y' }).layer('l', { type: 'ranges', data: d });
      s.track('z').layer('l', { type: 'letters' });
    });
    await v.load(b.build());
    await v.settle();
    frames.flush();
    return v;
  }
  const Item = '/data/d/(item :id B)';
  const Said = 'On a hidden track: show it to see it';

  it('says so and moves nothing; with a shown track drawing it too, it jumps there', async () => {
    const v = await withItems();
    expect(v.addresses.focus(Item).message).not.toBe(Said);
    await v.execute('track.hide', { address: '/views/main/one/x' });
    await v.settle();
    // still drawn on y: a jump
    expect(v.addresses.focus(Item).message).not.toBe(Said);
    await v.execute('track.hide', { address: '/views/main/one/y' });
    await v.settle();
    const before = v.view.get('/views/main/one');
    expect(v.addresses.focus(Item)).toMatchObject({ done: false, message: Said });
    expect(v.addresses.reveal(Item).message).toBe(Said);
    expect(v.view.get('/views/main/one')).toEqual(before);
    await v.execute('track.show', { address: '/views/main/one/y' });
    await v.settle();
    expect(v.addresses.focus(Item).message).not.toBe(Said);
  });
});

describe('track.show', () => {
  it('takes one, several or all, in one patch and one undo step', async () => {
    const v = await make();
    for (const a of [A, B, C]) await v.execute('track.hide', { address: a });
    await v.settle();
    const hidden = () => v.hiddenTracks.value.map((t) => t.address);
    let n = steps(v);
    await v.execute('track.show', { address: A });
    await v.settle();
    expect(hidden()).toEqual([B, C]);
    expect(steps(v)).toBe(n + 1);
    n = steps(v);
    await v.execute('track.show', { addresses: [B, C] });
    await v.settle();
    expect(hidden()).toEqual([]);
    expect(steps(v)).toBe(n + 1);
    for (const a of [A, B, C]) await v.execute('track.hide', { address: a });
    n = steps(v);
    await v.execute('track.show', { all: true });
    await v.settle();
    expect(hidden()).toEqual([]);
    expect(steps(v)).toBe(n + 1);
    await (v.spec.history as UiHistory).stack.undo();
    await v.settle();
    expect(hidden()).toEqual([A, B, C]);
  });

  it('all with nothing hidden, or at an address that is no section, is refused: no empty undo step', async () => {
    const v = await make();
    const n = steps(v);
    expect(await v.execute('track.show', { all: true })).toBe('No hidden tracks');
    await v.execute('track.hide', { address: A });
    await v.settle();
    const m = steps(v);
    expect(await v.execute('track.show', { all: true, address: '/views/main/two' })).toBe('No hidden tracks');
    expect(await v.execute('track.show', { all: true, address: A })).toBe('No hidden tracks');
    expect(steps(v)).toBe(m);
    expect(m).toBe(n + 1);
  });

  it('all within a section: the other section’s stay hidden', async () => {
    const v = await make();
    for (const a of [A, B, C]) await v.execute('track.hide', { address: a });
    await v.execute('track.show', { all: true, address: '/views/main/one' });
    await v.settle();
    expect(v.hiddenTracks.value.map((t) => t.address)).toEqual([C]);
  });
});

describe('mountHiddenTracks', () => {
  const button = (el: HTMLElement) => el.querySelector<HTMLButtonElement>('[data-toolbar-button]')!;

  it('is one button, always there: disabled with the reason while nothing is hidden, then named by the count', async () => {
    const v = await make();
    const el = document.createElement('div');
    document.body.append(el);
    const off = mountHiddenTracks(v, el);
    expect(button(el).getAttribute('aria-disabled')).toBe('true');
    expect(button(el).getAttribute('aria-label')).toBe('Hidden tracks');
    const badge = () => el.querySelector<HTMLElement>('.seqx-hidden-tracks-badge')!;
    expect([badge().textContent, badge().style.display]).toEqual(['', 'none']);
    await v.execute('track.hide', { address: A });
    await v.execute('track.hide', { address: C });
    expect(button(el).getAttribute('aria-disabled')).not.toBe('true');
    expect(button(el).getAttribute('aria-label')).toBe('Hidden tracks (2): show');
    // the count shows on the button
    expect([badge().textContent, badge().style.display]).toEqual(['2', '']);
    // taken away by what it returns
    off();
    expect(el.childElementCount).toBe(0);
  });

  it('its menu lists each hidden track then Show all; a row shows one, Show all every one in one patch', async () => {
    const v = await make();
    const el = document.createElement('div');
    document.body.append(el);
    const off = mountHiddenTracks(v, el);
    for (const a of [A, B, C]) await v.execute('track.hide', { address: a });
    const rows = () => [...document.querySelectorAll<HTMLElement>('[role="menu"] [role^="menuitem"]')];
    const open = () => {
      button(el).focus();
      button(el).dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }));
    };
    open();
    expect(rows().map((r) => r.textContent)).toEqual(['Alpha', 'b', 'Gamma', 'Show all']);
    rows()[1].click();
    await v.settle();
    expect(v.hiddenTracks.value.map((t) => t.address)).toEqual([A, C]);
    const n = steps(v);
    open();
    expect(rows().map((r) => r.textContent)).toEqual(['Alpha', 'Gamma', 'Show all']);
    rows()[2].click();
    await v.settle();
    expect(v.hiddenTracks.value).toEqual([]);
    expect(steps(v)).toBe(n + 1);
    off();
  });
});
