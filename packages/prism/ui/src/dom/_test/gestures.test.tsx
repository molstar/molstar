// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { lostReleases, stubPointerCapture } from '../../testing/gestures.js';
import { render } from '../render.js';
import { flushFrame, onCleanup } from '../scope.js';
import { drag } from '../scoped.js';

// Gestures attach once and every end path releases (prism AGENTS.md; `drag`, dom/scoped.ts): pressed, moved and left
// without a release N times, a gesture leaves no listener behind, and the next one works. The positive control (a
// listener added on the press, removed only on the release) is caught.

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

describe('drag (dom/scoped.ts)', () => {
  it('ends on a release, a cancel, the capture lost, a move with no button, and with its scope', () => {
    const ends: string[] = [];
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    const dispose = render(() => {
      el = document.createElement('div');
      drag(el, { start: () => ({}), end: (_s, how, e) => ends.push(`${how}:${e?.type ?? 'scope'}`) });
      return el;
    }, box);
    const down = () => el.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, buttons: 1, pointerId: 1 }));
    down();
    el.dispatchEvent(new PointerEvent('pointerup', { bubbles: true, pointerId: 1 }));
    down();
    el.dispatchEvent(new PointerEvent('pointercancel', { bubbles: true, pointerId: 1 }));
    down();
    el.dispatchEvent(new PointerEvent('lostpointercapture', { bubbles: true, pointerId: 1 }));
    down();
    el.dispatchEvent(new PointerEvent('pointermove', { bubbles: true, buttons: 0, pointerId: 1 }));
    down();
    down();
    dispose();
    expect(ends).toEqual([
      'up:pointerup',
      'cancel:pointercancel',
      'cancel:lostpointercapture',
      'cancel:pointermove',
      'cancel:pointerdown',
      'cancel:scope',
    ]);
  });

  it('a scroll area’s thumb: no listener left by lost releases, and the next drag scrolls', () => {
    const box = document.createElement('div');
    document.body.append(box);
    let viewport!: HTMLElement;
    render(
      () => (
        <ScrollArea
          viewportRef={(el) => {
            viewport = el;
          }}
        >
          <div style="height: 4000px" />
        </ScrollArea>
      ),
      box,
    );
    // jsdom lays nothing out: the bar's sizes as a browser would have them.
    Object.defineProperty(viewport, 'scrollHeight', { value: 4000, configurable: true });
    Object.defineProperty(viewport, 'clientHeight', { value: 400, configurable: true });
    const track = () => box.querySelector<HTMLElement>('[data-orientation=vertical]')!;
    Object.defineProperty(track(), 'clientHeight', { value: 400, configurable: true });
    viewport.dispatchEvent(new Event('scroll'));
    flushFrame();
    const r = lostReleases({
      target: track,
      from: { x: 5, y: 0 },
      to: { x: 5, y: 40 },
      reset: () => (viewport.scrollTop = 0),
      worked: () => viewport.scrollTop > 0 && !track().hasAttribute('data-dragging'),
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });

  it('a listener added on the press and removed only on the release is caught (positive control)', () => {
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    let moved = 0;
    render(() => {
      el = document.createElement('div');
      const move = () => moved++;
      const up = () => {
        document.removeEventListener('pointermove', move);
        document.removeEventListener('pointerup', up);
      };
      el.addEventListener('pointerdown', () => {
        document.addEventListener('pointermove', move);
        document.addEventListener('pointerup', up);
      });
      onCleanup(() => {});
      return el;
    }, box);
    const r = lostReleases({ target: () => el, from: { x: 0, y: 0 }, to: { x: 5, y: 5 }, worked: () => moved > 0 });
    expect(r.grown).toBeGreaterThan(0);
  });
});
