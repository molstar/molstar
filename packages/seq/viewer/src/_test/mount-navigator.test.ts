// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigator part in a DOM (jsdom): `mountNavigator` into a host's element (its canvases, its slider's value set on
// settle and never for a move made elsewhere, its events in the platform's order: a double-click one step), the
// default placement by `place` around the canvas, the rail toggle, the canvas's cursor over a rail; everything taken
// away on retire.

import { createSeqSpec, overviewNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { viewerItem } from '../config.js';
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
    const ev = (type: string, x: number, detail = 0) =>
      new MouseEvent(type, { clientX: x, clientY: 10, button: 0, bubbles: true, detail }) as unknown as PointerEvent;
    part.dispatchEvent(ev('pointerdown', 450));
    part.dispatchEvent(ev('pointerup', 450));
    part.dispatchEvent(ev('click', 450, 1));
    const v = viewer.view.get('/views/v/protein')!;
    expect((v.start + v.end) / 2).toBeCloseTo((450 / 600) * 300, 0);
    expect(viewer.view.history.value.back).toBe(1);
    // a double-click elsewhere, in the platform's order: one step more
    viewer.draw();
    frames.flush();
    for (const detail of [1, 2]) {
      part.dispatchEvent(ev('pointerdown', 100));
      part.dispatchEvent(ev('pointerup', 100));
      part.dispatchEvent(ev('click', 100, detail));
    }
    part.dispatchEvent(ev('dblclick', 100, 2));
    expect(viewer.view.get('/views/v/protein')).toMatchObject({ start: 0, end: 300 });
    expect(viewer.view.history.value.back).toBe(2);
    dispose();
    expect(el.children.length).toBe(0);
  });

  it('its position readout is the chrome’s surface on one line with a space after the colon: shown over a hover, hidden after it', async () => {
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    await viewer.load(spec());
    viewer.resize(600);
    await viewer.settle();
    frames.flush();
    viewer.view.setRange('/views/v/protein', 0, 30);
    viewer.draw();
    const el = document.createElement('div');
    document.body.append(el);
    mountNavigator(el, viewer, 'strip');
    frames.flush();
    const part = el.querySelector<HTMLElement>('[data-navigator="strip"]')!;
    vi.spyOn(part, 'getBoundingClientRect').mockReturnValue({ left: 0, top: 0 } as DOMRect);
    const readout = part.querySelector<HTMLElement>('.seqx-navigator-readout')!;
    // nothing said: no backdrop drawn
    expect(readout.style.display).toBe('none');
    part.dispatchEvent(
      new MouseEvent('pointermove', { clientX: 450, clientY: 2, bubbles: true }) as unknown as PointerEvent,
    );
    viewer.draw();
    frames.flush();
    expect(readout.textContent).toMatch(/^p: [\d,]+$/);
    expect(readout.style.display).toBe('block');
    // the surface's inset and backdrop (`msui-surface`), one line cut by an ellipsis, inside the part's box
    expect(readout.classList.contains('msui-surface')).toBe(true);
    expect(readout.style).toMatchObject({ whiteSpace: 'nowrap', textOverflow: 'ellipsis', overflow: 'hidden' });
    expect(readout.style.maxWidth).toContain('100%');
    part.dispatchEvent(new MouseEvent('pointerleave', { bubbles: true }) as unknown as PointerEvent);
    viewer.draw();
    frames.flush();
    expect(readout.style.display).toBe('none');
  });

  it('its slider’s value is set on settle after its own gesture or key, never for a move made elsewhere', async () => {
    // jumps with no motion: the view rests once its frame is drawn
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
    await viewer.load(spec());
    viewer.resize(600);
    await viewer.settle();
    frames.flush();
    viewer.view.setRange('/views/v/protein', 0, 30);
    viewer.draw();
    const el = document.createElement('div');
    document.body.append(el);
    mountNavigator(el, viewer, 'strip');
    frames.timers.fire();
    frames.flush();
    const part = el.querySelector<HTMLElement>('[data-navigator="strip"]')!;
    const said = part.getAttribute('aria-valuetext');
    expect(said).toBe('p:1–30 of p');
    // focused, the view moved elsewhere: frames drawn, the view rested, nothing said
    part.dispatchEvent(new FocusEvent('focus'));
    viewer.view.setRange('/views/v/protein', 100, 130);
    frames.flush();
    frames.timers.fire();
    frames.flush();
    expect(part.getAttribute('aria-valuetext')).toBe(said);
    // its own key: said once the view rests
    part.dispatchEvent(new KeyboardEvent('keydown', { key: 'Home' }));
    frames.flush();
    expect(part.getAttribute('aria-valuetext')).toBe(said);
    frames.timers.fire();
    frames.flush();
    expect(part.getAttribute('aria-valuetext')).toBe('p:1–30 of p');
    part.dispatchEvent(new KeyboardEvent('keydown', { key: 'End' }));
    frames.flush();
    frames.timers.fire();
    frames.flush();
    expect(part.getAttribute('aria-valuetext')).toBe('p:271–300 of p');
  });

  it('the canvas’s cursor over a rail: grab, grabbing while dragged', async () => {
    const target = document.createElement('div');
    document.body.append(target);
    viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: { canvas: true } }, target);
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(30) });
    b.views.view('v', (v) => {
      const s = v.section('protein', { of: p });
      s.track('ruler').layer('ruler', { type: 'ruler' });
      s.track('letters').layer('l', { type: 'letters' });
    });
    await viewer.load(b.build());
    await viewer.settle();
    frames.flush();
    viewer.view.setRange('/views/v/protein', 100, 150);
    frames.flush();
    const canvas = target.querySelector<HTMLElement>('[role="application"]')!;
    vi.spyOn(canvas, 'getBoundingClientRect').mockReturnValue({ left: 0, top: 0 } as DOMRect);
    const rail = viewer.frame!.layout.sections[0].rail!;
    const x = rail.x + (125 / rail.total) * rail.width;
    const ev = (type: string) =>
      new MouseEvent(type, {
        clientX: x,
        clientY: rail.bottom - 2,
        button: 0,
        bubbles: true,
      }) as unknown as PointerEvent;
    canvas.dispatchEvent(ev('pointermove'));
    expect(canvas.style.cursor).toBe('grab');
    canvas.dispatchEvent(ev('pointerdown'));
    expect(canvas.style.cursor).toBe('grabbing');
    canvas.dispatchEvent(ev('pointerup'));
    expect(canvas.style.cursor).toBe('grab');
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
