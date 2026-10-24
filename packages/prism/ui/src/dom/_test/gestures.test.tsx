// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
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
    // The capture lost while the button is held: a cancel at once.
    el.dispatchEvent(new PointerEvent('lostpointercapture', { bubbles: true, buttons: 1, pointerId: 1 }));
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

  it("waits a moment after a capture lost with no button held: Chrome's release can come just after it", () => {
    vi.useFakeTimers();
    const ends: string[] = [];
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    const dispose = render(() => {
      el = document.createElement('div');
      el.append(document.createElement('canvas'));
      drag(el, { start: () => ({}), end: (_s, how, e) => ends.push(`${how}:${e?.type ?? 'scope'}`) });
      return el;
    }, box);
    const down = () => el.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, buttons: 1, pointerId: 1 }));
    const lost = () =>
      el.dispatchEvent(new PointerEvent('lostpointercapture', { bubbles: true, buttons: 0, pointerId: 1 }));
    // The capture lost, then the release on a child (the order Chrome sometimes takes): the release ends it.
    down();
    lost();
    expect(ends).toEqual([]);
    el.firstElementChild!.dispatchEvent(new PointerEvent('pointerup', { bubbles: true, pointerId: 1 }));
    expect(ends).toEqual(['up:pointerup']);
    vi.advanceTimersByTime(1000);
    expect(ends).toEqual(['up:pointerup']);
    // No release after it: a cancel, a moment later.
    down();
    lost();
    vi.advanceTimersByTime(1000);
    expect(ends).toEqual(['up:pointerup', 'cancel:lostpointercapture']);
    dispose();
    vi.useRealTimers();
  });

  it("takes a release elsewhere in the page only while a lost capture waits for it, never another gesture's", () => {
    vi.useFakeTimers();
    const ends: string[] = [];
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    const dispose = render(() => {
      el = document.createElement('div');
      drag(el, { start: () => ({}), end: (_s, how, e) => ends.push(`${how}:${e?.type ?? 'scope'}`) });
      return el;
    }, box);
    const send = (at: Element, type: string, pointerId = 1, buttons = 0) =>
      at.dispatchEvent(new PointerEvent(type, { bubbles: true, buttons, pointerId }));
    // The capture lost, then the release elsewhere (another pointer's is not this drag's): the release.
    send(el, 'pointerdown', 1, 1);
    send(el, 'lostpointercapture');
    send(document.body, 'pointerup', 2);
    send(document.body, 'pointerup');
    expect(ends).toEqual(['up:pointerup']);
    // A release lost with the capture held (nothing comes): a click elsewhere later is another gesture, never this
    // drag's release; it ends nothing until this drag's own end (here, a press on it again).
    send(el, 'pointerdown', 1, 1);
    send(document.body, 'pointerdown', 1, 1);
    send(document.body, 'pointerup');
    expect(ends).toEqual(['up:pointerup']);
    // While a lost capture waits, a press anywhere is another gesture: a cancel.
    send(el, 'lostpointercapture');
    send(document.body, 'pointerdown', 1, 1);
    expect(ends).toEqual(['up:pointerup', 'cancel:pointerdown']);
    dispose();
    vi.useRealTimers();
  });

  it('follows only the pointer that started it: a second finger neither moves, ends nor restarts it', () => {
    const log: string[] = [];
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    const dispose = render(() => {
      el = document.createElement('div');
      drag(el, {
        start: (e) => {
          log.push(`start:${e.pointerId}`);
          return {};
        },
        move: (e) => log.push(`move:${e.pointerId}`),
        end: (_s, how, e) => log.push(`${how}:${e?.pointerId ?? 'scope'}`),
      });
      return el;
    }, box);
    const send = (type: string, pointerId: number, buttons = 1) =>
      el.dispatchEvent(new PointerEvent(type, { bubbles: true, buttons, pointerId, pointerType: 'touch' }));
    send('pointerdown', 1);
    send('pointerdown', 2);
    send('pointermove', 2);
    send('pointermove', 1);
    send('pointerup', 2, 0);
    send('pointercancel', 2, 0);
    send('lostpointercapture', 2, 0);
    send('pointerup', 1, 0);
    dispose();
    expect(log).toEqual(['start:1', 'move:1', 'up:1']);
  });

  it('a press by another pointer after the first one ended unseen (its capture gone) starts anew', () => {
    const log: string[] = [];
    const box = document.createElement('div');
    document.body.append(box);
    let el!: HTMLElement;
    const dispose = render(() => {
      el = document.createElement('div');
      drag(el, {
        start: (e) => {
          log.push(`start:${e.pointerId}`);
          return {};
        },
        end: (_s, how, e) => log.push(`${how}:${e?.pointerId ?? 'scope'}`),
      });
      return el;
    }, box);
    el.dispatchEvent(
      new PointerEvent('pointerdown', { bubbles: true, buttons: 1, pointerId: 1, pointerType: 'touch' }),
    );
    // The first touch ended without telling the element: its capture went with it.
    el.releasePointerCapture(1);
    el.dispatchEvent(
      new PointerEvent('pointerdown', { bubbles: true, buttons: 1, pointerId: 2, pointerType: 'touch' }),
    );
    dispose();
    expect(log).toEqual(['start:1', 'cancel:2', 'start:2', 'cancel:scope']);
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
