// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigator part in a DOM (jsdom): `mountNavigator` into a host's element (its canvases, its slider's value, its
// events), the default placement by `place` around the canvas, the rail toggle; everything taken away on retire.

import { createSeqSpec, overviewNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { mountNavigator, railToggleMount } from '../mount-navigator.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames;
let viewer: SeqViewer | undefined;

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  document.body.replaceChildren();
});

function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(30) });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('letters').layer('letters', { type: 'letters' });
    v.navigator('strip', stripNavigator(s));
    v.navigator('side', { ...overviewNavigator(s), content: 'axis', axes: 'positions', place: 'right' });
  });
  return b.build();
}

describe('the navigator part', () => {
  it('the default placement: a strip above the canvas, a side part beside it; taken away on retire', async () => {
    const target = document.createElement('div');
    document.body.append(target);
    viewer = createSeqViewer(
      { ...DefaultSeqViewerSpec(), ui: { canvas: true, navigators: DefaultSeqViewerSpec().ui!.navigators! } },
      target,
    );
    await viewer.load(spec());
    await viewer.settle();
    frames.flush();
    const parts = [...target.querySelectorAll<HTMLElement>('[data-navigator]')];
    expect(parts.map((p) => p.dataset.navigator)).toEqual(['strip', 'side']);
    const canvas = target.querySelector('[role="application"]')!;
    // the strip before the canvas in the document's order, the side part after it
    expect(parts[0].compareDocumentPosition(canvas) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    expect(parts[1].compareDocumentPosition(canvas) & Node.DOCUMENT_POSITION_PRECEDING).toBeTruthy();
    expect(parts[0].getAttribute('role')).toBe('slider');
    expect(parts[0].getAttribute('aria-valuetext')).toMatch(/^p:1–\d+ of p$/);
    viewer.retire();
    viewer = undefined;
    expect(target.querySelectorAll('[data-navigator]').length).toBe(0);
  });

  it('mountNavigator into a host’s element: its events move the view; its disposer takes it away', async () => {
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    await viewer.load(spec());
    viewer.resize(600);
    await viewer.settle();
    frames.flush();
    viewer.view.setRange('/views/v/protein', 0, 30);
    viewer.draw();
    const el = document.createElement('div');
    document.body.append(el);
    const dispose = mountNavigator(el, viewer, 'strip');
    frames.flush();
    const part = el.querySelector<HTMLElement>('[data-navigator="strip"]')!;
    vi.spyOn(part, 'getBoundingClientRect').mockReturnValue({ left: 0, top: 0 } as DOMRect);
    const ev = (type: string, x: number) =>
      new MouseEvent(type, { clientX: x, clientY: 10, button: 0, bubbles: true }) as unknown as PointerEvent;
    part.dispatchEvent(ev('pointerdown', 450));
    part.dispatchEvent(ev('pointerup', 450));
    const v = viewer.view.get('/views/v/protein')!;
    expect((v.start + v.end) / 2).toBeCloseTo((450 / 600) * 300, 0);
    dispose();
    expect(el.children.length).toBe(0);
  });

  it('the rail toggle hides and shows the rails: the reader’s session, not the spec', async () => {
    const target = document.createElement('div');
    document.body.append(target);
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: { canvas: true, railToggle: railToggleMount } }, target);
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(30) });
    b.views.view('v', (v) => v.section('protein', { of: p }).track('letters').layer('l', { type: 'letters' }));
    await viewer.load(b.build());
    await viewer.settle();
    frames.flush();
    expect(viewer.frame!.layout.sections[0].rail).toBeDefined();
    const button = target.querySelector('button')!;
    expect(button.getAttribute('aria-pressed')).toBe('true');
    button.click();
    frames.flush();
    expect(viewer.railsShown.value).toBe(false);
    expect(viewer.frame!.layout.sections[0].rail).toBeUndefined();
    expect(button.getAttribute('aria-pressed')).toBe('false');
  });
});
