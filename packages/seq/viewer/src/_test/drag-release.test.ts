// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { mountNavigator } from '../mount-navigator.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { proteinOnGene } from './fixture.js';
import { type Frames, installFrames } from './frames.js';

// A drag ends wherever the button is released (_spec/README.md "Input"): the mounts listen on their window while a
// drag runs, so a release outside the element, a lost button, a blurred window or a hidden document end it, and
// nothing is left listening afterwards.

let frames: Frames;
let viewer: SeqViewer | undefined;
let added: string[] = [];
let removed: string[] = [];

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
  added = [];
  removed = [];
  const add = EventTarget.prototype.addEventListener;
  const remove = EventTarget.prototype.removeEventListener;
  vi.spyOn(EventTarget.prototype, 'addEventListener').mockImplementation(function (
    this: EventTarget,
    ...a: [any, any, any]
  ) {
    if (this === window || this === document) added.push(a[0]);
    return add.apply(this, a);
  });
  vi.spyOn(EventTarget.prototype, 'removeEventListener').mockImplementation(function (
    this: EventTarget,
    ...a: [any, any, any]
  ) {
    if (this === window || this === document) removed.push(a[0]);
    return remove.apply(this, a);
  });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  vi.restoreAllMocks();
  document.body.replaceChildren();
});

const ev = (type: string, x: number, init: MouseEventInit = {}) =>
  new MouseEvent(type, { clientX: x, clientY: 18, button: 0, bubbles: true, cancelable: true, ...init });

/** The temporary listeners a drag leaves on the window and the document: added and not yet removed. */
const live = () => {
  const open = [...added];
  for (const t of removed) {
    const i = open.indexOf(t);
    if (i >= 0) open.splice(i, 1);
  }
  return open.filter((t) => ['pointerup', 'pointercancel', 'pointermove', 'blur', 'visibilitychange'].includes(t));
};

interface Subject {
  element: HTMLElement;
  press(): void;
  dragging(): boolean;
  retire(): void;
}

async function canvas(): Promise<Subject> {
  const target = document.createElement('div');
  document.body.append(target);
  viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: { canvas: true } }, target);
  await viewer.load(proteinOnGene());
  await viewer.settle();
  frames.flush();
  const element = target.querySelector<HTMLElement>('[role="application"]')!;
  const v = viewer;
  return {
    element,
    press: () => element.dispatchEvent(ev('pointerdown', 40)),
    dragging: () => v.input.dragging,
    retire: () => v.retire(),
  };
}

async function navigator(): Promise<Subject> {
  viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(30) });
  b.views.view('v', (view) => {
    const s = view.section('protein', { of: p });
    s.track('letters').layer('letters', { type: 'letters' });
    view.navigator('strip', stripNavigator(s));
  });
  await viewer.load(b.build());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  viewer.view.setRange('/views/v/protein', 0, 30);
  viewer.draw();
  const host = document.createElement('div');
  document.body.append(host);
  const dispose = mountNavigator(host, viewer, 'strip');
  frames.flush();
  const element = host.querySelector<HTMLElement>('[data-navigator="strip"]')!;
  vi.spyOn(element, 'getBoundingClientRect').mockReturnValue({ left: 0, top: 0 } as DOMRect);
  const model = viewer.navigatorModel('strip')!;
  return {
    element,
    // the window's body
    press: () => element.dispatchEvent(ev('pointerdown', 20)),
    dragging: () => model.cursor === 'grabbing',
    retire: dispose,
  };
}

describe.each([
  ['the canvas mount', canvas],
  ['the navigator mount', navigator],
])('%s: a drag ends wherever the button is released', (_, make) => {
  const ends: [string, (s: Subject) => void][] = [
    ['a pointerup on the window', () => window.dispatchEvent(ev('pointerup', 900))],
    ['a pointerup on another element', () => document.body.dispatchEvent(ev('pointerup', 900))],
    ['a pointercancel on the window', () => window.dispatchEvent(ev('pointercancel', 900))],
    ['a move with no button held', () => document.body.dispatchEvent(ev('pointermove', 900, { buttons: 0 }))],
    ['the window losing focus', () => window.dispatchEvent(new Event('blur'))],
    [
      'the document turning hidden',
      () => {
        vi.spyOn(document, 'visibilityState', 'get').mockReturnValue('hidden');
        document.dispatchEvent(new Event('visibilitychange'));
      },
    ],
    ['a lost pointer capture', (s) => s.element.dispatchEvent(ev('lostpointercapture', 20))],
  ];
  for (const [what, end] of ends)
    it(`ends on ${what}, and leaves no listener`, async () => {
      const s = await make();
      expect(live()).toEqual([]);
      s.press();
      expect(s.dragging()).toBe(true);
      expect(live().length).toBeGreaterThan(0);
      end(s);
      expect(s.dragging()).toBe(false);
      expect(live()).toEqual([]);
    });

  it('a move with a button held does not end it', async () => {
    const s = await make();
    s.press();
    document.body.dispatchEvent(ev('pointermove', 100, { buttons: 1 }));
    expect(s.dragging()).toBe(true);
    window.dispatchEvent(ev('pointerup', 100));
    expect(s.dragging()).toBe(false);
  });

  it('a second press after a lost release starts a new drag, and the first one’s listeners are not doubled', async () => {
    const s = await make();
    s.press();
    const first = live().length;
    // no release arrived: the next press begins again
    s.press();
    expect(s.dragging()).toBe(true);
    expect(live().length).toBe(first);
    window.dispatchEvent(ev('pointerup', 900));
    expect(s.dragging()).toBe(false);
    expect(live()).toEqual([]);
    s.press();
    expect(s.dragging()).toBe(true);
    s.element.dispatchEvent(ev('pointerup', 20));
    expect(s.dragging()).toBe(false);
  });

  it('an unmount during a drag takes its listeners away', async () => {
    const s = await make();
    s.press();
    expect(live().length).toBeGreaterThan(0);
    s.retire();
    expect(live()).toEqual([]);
  });
});
