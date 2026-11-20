// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Accessibility of the viewer (_spec/README.md "Input", Announcements and Forced colors): the live region says the
// keyboard cursor's place, a hover that rests and a selection's change, one text at a time within a rate limit; forced
// colors read the system colors into the tokens; the rows' order is reachable from the keyboard.

import { uiHistory } from '../ui-history.js';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { builtInTokens, colorHex } from '@molstar/seqx-core/render';
import { rangesLocation } from '@molstar/seqx-core/query';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { forcedColorTokens, readCssTokens } from '../canvas/index.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];
let clock = 1000;

beforeEach(() => {
  frames = installFrames();
  clock = 1000;
  vi.spyOn(DomHost, 'now').mockImplementation(() => clock);
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

const Section = '/views/v/s';

function domainsSpec() {
  const b = createSeqSpec({ title: 'announcements' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFS' });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: new Int32Array([2, 4]),
      end: new Int32Array([8, 12]),
      id: ['A', 'B'],
      label: ['Domain A', 'Domain B'],
    },
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: p });
    s.track('letters', { label: 'p' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
  });
  return b.build();
}

async function make() {
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] as never,
  });
  viewers.push(v);
  await v.load(domainsSpec());
  v.resize(600);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

/** The pointer's y on the domains track of the section. */
function domainsY(v: SeqViewer): number {
  const s = v.frame!.layout.sections[0];
  const t = s.tracks[1];
  return t.box.y + t.box.height / 4;
}

describe('the announcer', () => {
  it('says the keyboard cursor’s place: section, system, position, letter and the items there', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0 });
    // the first key places the cursor where the view starts; the next moves it on
    await v.execute('cursor.next', { section: Section });
    expect(a.message.value.text).toBe('s, p 1, letter M');
    clock += 1000;
    await v.execute('cursor.next', { section: Section });
    await v.execute('cursor.next', { section: Section });
    // position 3 (index 2) is inside Domain A only
    expect(a.message.value.text).toBe('s, p 3, letter E, items: Domain A, p 3–8');
    await v.execute('cursor.next', { section: Section });
    await v.execute('cursor.next', { section: Section });
    // index 4 is in both domains
    expect(a.message.value.text).toBe('s, p 5, letter Q, items: Domain A, p 3–8; Domain B, p 5–12');
    a.dispose();
  });

  it('is rate-limited: a held key says where it is each gap, the latest waiting its turn', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 700 });
    const said: string[] = [];
    a.message.listen((m) => said.push(m.text));
    for (let i = 0; i < 6; i++) await v.execute('cursor.next', { section: Section });
    // the first at once, the rest waiting as one (the latest)
    expect(said).toHaveLength(1);
    clock += 700;
    frames.timers.fire();
    expect(said).toHaveLength(2);
    expect(said[1]).toMatch(/^s, p 6, letter S/);
    a.dispose();
  });

  it('says a hover only when it rests, and what is under it', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0, restMs: 500 });
    const layout = v.frame!.layout;
    const axis = layout.sections[0].axis!;
    const x = (axis.x(3) + axis.x(4)) / 2;
    v.input.move(x, domainsY(v));
    frames.flush();
    // moving says nothing; the rest does
    expect(a.message.value.text).toBe('');
    frames.timers.fire();
    expect(a.message.value.text).toBe('Domain A, p 3–8');
    a.dispose();
  });

  it('drops a hover waiting its turn when the pointer leaves, and a rest that never came', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 700, restMs: 500 });
    const said: string[] = [];
    a.message.listen((m) => said.push(m.text));
    const axis = v.frame!.layout.sections[0].axis!;
    const x = (axis.x(3) + axis.x(4)) / 2;
    // a message has just been said: the hover that rests waits for the gap
    await v.execute('cursor.next', { section: Section });
    expect(said).toHaveLength(1);
    v.input.move(x, domainsY(v));
    frames.flush();
    frames.timers.fire();
    expect(said).toHaveLength(1);
    // the pointer leaves: the waiting hover is dropped, nothing is said when the gap is over
    v.input.leave();
    frames.flush();
    clock += 1000;
    frames.timers.fire();
    expect(said).toHaveLength(1);
    // a hover that left before it rested says nothing at all
    v.input.move(x, domainsY(v));
    v.input.leave();
    frames.flush();
    frames.timers.fire();
    expect(said).toHaveLength(1);
    a.dispose();
  });

  it('keeps a waiting selection over a later hover and cursor, which never replace it', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 700, restMs: 500 });
    const said: string[] = [];
    a.message.listen((m) => said.push(m.text));
    await v.execute('cursor.next', { section: Section });
    const loc = rangesLocation('/data/p', [2, 8]);
    v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    // the cursor and the hover come later, within the gap
    a.say('a hover', 'hover');
    a.say('a cursor', 'cursor');
    expect(said).toHaveLength(1);
    clock += 700;
    frames.timers.fire();
    expect(said).toEqual(['s, p 1, letter M', 'Selected 6 positions, p 3–8']);
    a.dispose();
  });

  it('says nothing of the cursor or hover during a drag, and the selection once when it ends', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0, restMs: 500 });
    const said: string[] = [];
    a.message.listen((m) => said.push(m.text));
    const loc = rangesLocation('/data/p', [2, 8]);
    v.selection.beginDrag({ system: '/data/p', position: 2 } as never);
    v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    await v.execute('cursor.next', { section: Section });
    expect(said).toEqual([]);
    clock += 10;
    frames.timers.fire();
    // still dragging: still nothing
    expect(said).toEqual([]);
    v.selection.endDrag();
    frames.timers.fire();
    expect(said).toEqual(['Selected 6 positions, p 3–8']);
    a.dispose();
  });

  it('says a key that leaves the cursor where it was: the same words, a new message', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0 });
    const said: { id: number; text: string }[] = [];
    a.message.listen((m) => said.push(m));
    await v.execute('cursor.prev', { section: Section });
    // ↑ in the top section stays (← at the first position goes onto the track's header now: a key that moves)
    await v.execute('cursor.up', { section: Section });
    expect(said.map((m) => m.text)).toEqual(['s, p 1, letter M', 's, p 1, letter M']);
    expect(said[1].id).toBeGreaterThan(said[0].id);
    a.dispose();
  });

  it('says a selection’s change once, its size and extent', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0 });
    const said: string[] = [];
    a.message.listen((m) => said.push(m.text));
    const loc = rangesLocation('/data/p', [2, 8]);
    v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    expect(said).toEqual(['Selected 6 positions, p 3–8']);
    // the same selection again is no change
    v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    expect(said).toHaveLength(1);
    clock += 10;
    v.selection.clear();
    expect(said.at(-1)).toBe('Selection cleared');
    a.dispose();
  });
});

describe('the mounted live region', () => {
  it('is one polite status the mount owns, saying the cursor after a key', async () => {
    vi.stubGlobal(
      'ResizeObserver',
      class {
        observe() {}
        unobserve() {}
        disconnect() {}
      },
    );
    vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
      return ((this as any).__ctx ??= new RecordingContext());
    });
    const target = document.createElement('div');
    document.body.appendChild(target);
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: { canvas: true as const, history: uiHistory } }, target);
    viewers.push(v);
    await v.load(domainsSpec());
    await v.settle();
    frames.flush();
    const live = target.querySelectorAll('[role="status"][aria-live="polite"]');
    expect(live).toHaveLength(1);
    await v.execute('cursor.next', { section: Section });
    expect(live[0].textContent).toBe('s, p 1, letter M');
    // a key that leaves the cursor in place clears and sets the region again: the same words are heard again
    const writes: string[] = [];
    new MutationObserver(() => {}).observe(live[0], { childList: true });
    const spy = vi
      .spyOn(live[0], 'textContent', 'set')
      .mockImplementation((t: string | null) => void writes.push(t ?? ''));
    clock += 1000;
    // ↑ in the top section stays (← at the first position goes onto the track's header now)
    await v.execute('cursor.up', { section: Section });
    expect(writes).toEqual(['', 's, p 1, letter M']);
    spy.mockRestore();
    v.retire();
    expect(target.querySelector('[role="status"]')).toBeNull();
  });
});

describe('forced colors, whatever the theme source', () => {
  it('maps a configured theme’s tokens to the system colors when the platform forces them', async () => {
    vi.stubGlobal(
      'ResizeObserver',
      class {
        observe() {}
        unobserve() {}
        disconnect() {}
      },
    );
    vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
      return ((this as any).__ctx ??= new RecordingContext());
    });
    const win = document.defaultView!;
    let on = true;
    const listeners = new Set<() => void>();
    Object.defineProperty(win, 'matchMedia', {
      configurable: true,
      writable: true,
      value: (media: string) => ({
        media,
        get matches() {
          return on && media === '(forced-colors: active)';
        },
        addEventListener: (_: string, fn: () => void) => media === '(forced-colors: active)' && listeners.add(fn),
        removeEventListener: (_: string, fn: () => void) => listeners.delete(fn),
      }),
    });
    const target = document.createElement('div');
    document.body.appendChild(target);
    const v = createSeqViewer(
      {
        ...DefaultSeqViewerSpec(),
        ui: { canvas: true as const, history: uiHistory },
        config: [[viewerItem('theme.source'), 'config']] as never,
      },
      target,
    );
    viewers.push(v);
    const plain = builtInTokens('light');
    expect(v.tokens.version).not.toBe(plain.version);
    // the accent, the hover and the window are one color, the system's Highlight; the accent was another
    expect(v.tokens.colors['msui-accent']).not.toBe(plain.colors['msui-accent']);
    expect(colorHex(v.tokens.colors['seq-hover']).slice(0, 7)).toBe(
      colorHex(v.tokens.colors['msui-accent']).slice(0, 7),
    );
    // the platform stops forcing: the configured table is back
    on = false;
    for (const fn of listeners) fn();
    expect(v.tokens.version).toBe(plain.version);
    Reflect.deleteProperty(win, 'matchMedia');
  });
});

describe('forced colors', () => {
  const system = { Canvas: '#101010', CanvasText: '#f0f0f0', Highlight: '#1aebff', GrayText: '#3fe2e2' };

  it('takes the chrome and the emphasis from the system colors and keeps the data’s', () => {
    const base = builtInTokens('light');
    const forced = forcedColorTokens(base, system);
    const hex = (name: keyof typeof base.colors) => colorHex(forced.colors[name]).slice(0, 7);
    expect(hex('msui-bg')).toBe('#101010');
    expect(hex('msui-text')).toBe('#f0f0f0');
    expect(hex('seq-sequence')).toBe('#f0f0f0');
    expect(hex('msui-accent')).toBe('#1aebff');
    expect(hex('seq-hover')).toBe('#1aebff');
    expect(hex('seq-window')).toBe('#1aebff');
    expect(hex('msui-text-muted')).toBe('#3fe2e2');
    expect(hex('seq-cell-edge')).toBe('#f0f0f0');
    // the data's colors are kept: the meaning is in them
    for (const name of ['seq-aa-hydrophobic', 'seq-cat-0', 'seq-ramp-3', 'seq-nt-a'] as const)
      expect(forced.colors[name]).toBe(base.colors[name]);
    expect(forced.theme).toBe(base.theme);
    expect(forced.version).not.toBe(base.version);
  });

  it('is read from the root when the media query matches, and not otherwise', () => {
    const root = document.createElement('div');
    document.body.appendChild(root);
    const win = document.defaultView!;
    const match = (on: boolean) =>
      Object.defineProperty(win, 'matchMedia', {
        configurable: true,
        writable: true,
        value: (media: string) => ({ media, matches: on && media === '(forced-colors: active)' }),
      });
    match(false);
    const plain = readCssTokens(root);
    expect(plain.version).toBe(builtInTokens('light').version);
    match(true);
    const forced = readCssTokens(root);
    expect(forced.version).not.toBe(plain.version);
    // the system colors the page could not resolve (no theme in a test DOM) are the platform's defaults
    expect(colorHex(forced.colors['msui-bg']).slice(0, 7)).toBe('#ffffff');
    expect(colorHex(forced.colors['msui-text']).slice(0, 7)).toBe('#000000');
    Reflect.deleteProperty(win, 'matchMedia');
  });
});

describe('the keyboard to the row columns’ titles', () => {
  async function rowsViewer() {
    const b = createSeqSpec({ title: 'order' });
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: ['NGB', 'CYG', 'MYG'], letters: ['MEEPQ', 'MEEPQSS', 'MEE'] },
    });
    b.views.view('v', (vw) =>
      vw.section('s', { of: aln }, (s) => {
        s.track('rows', { label: '', height: 100 }, (t) =>
          t.layer('r', {
            type: 'rows',
            data: aln,
            columns: [
              { field: 'label', side: 'right' },
              { field: 'length', side: 'right' },
            ],
          } as never),
        );
      }),
    );
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    return v;
  }
  const order = (v: SeqViewer) =>
    v.spec.root.value.children![1].children![0].children![0].children![0].children![0].params?.order;

  it('Shift+O cycles the fields: ascending, descending, then the next column', async () => {
    const v = await rowsViewer();
    const track = '/views/v/s/rows';
    const step = async () => {
      const message = await v.execute('rows.order-next', { track });
      await v.settle();
      frames.flush();
      return message;
    };
    expect(await step()).toBe('Ordered by label, ascending');
    expect(order(v)).toEqual({ field: 'label', direction: 'ascending' });
    expect(await step()).toBe('Ordered by label, descending');
    expect(order(v)).toEqual({ field: 'label', direction: 'descending' });
    expect(await step()).toBe('Ordered by length, ascending');
    expect(order(v)).toEqual({ field: 'length', direction: 'ascending' });
    await step();
    await step();
    expect(order(v)).toEqual({ field: 'label', direction: 'ascending' });
  });

  it('is on a key: Shift+O on the canvas orders the rows', async () => {
    const v = await rowsViewer();
    expect(v.config.view.keys['rows.order-next']).toEqual(['Shift+O']);
    expect(v.input.key({ key: 'O', shiftKey: true })).toBe(true);
    await v.settle();
    expect(order(v)).toEqual({ field: 'label', direction: 'ascending' });
  });
});
