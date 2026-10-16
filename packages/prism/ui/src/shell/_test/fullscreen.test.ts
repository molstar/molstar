/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { FullscreenModel } from '../fullscreen.js';

describe('FullscreenModel', () => {
  const enabled = Object.getOwnPropertyDescriptor(Document.prototype, 'fullscreenEnabled');
  afterEach(() => {
    if (enabled) Object.defineProperty(Document.prototype, 'fullscreenEnabled', enabled);
    else delete (document as any).fullscreenEnabled;
  });
  const allow = (value: boolean) => Object.defineProperty(document, 'fullscreenEnabled', { value, configurable: true });

  it('on a page: Expand, fullscreen beside it', () => {
    allow(true);
    const model = new FullscreenModel(() => document.body);
    expect(model.embedded).toBe(false);
    expect(model.mode).toBe('expand');
    expect(model.unavailableReason).toBeUndefined();
    model.dispose();
  });

  it('where the page forbids it: Expand, fullscreen explained', () => {
    allow(false);
    const model = new FullscreenModel(() => document.body);
    expect(model.mode).toBe('expand');
    expect(model.unavailableReason).toMatch(/embedded/);
    model.dispose();
  });

  it('follows leaving fullscreen by Escape', () => {
    allow(true);
    const model = new FullscreenModel(() => document.body);
    const states: boolean[] = [];
    const sub = model.active.subscribe((a) => states.push(a));
    Object.defineProperty(document, 'fullscreenElement', { value: document.body, configurable: true });
    document.dispatchEvent(new Event('fullscreenchange'));
    Object.defineProperty(document, 'fullscreenElement', { value: null, configurable: true });
    document.dispatchEvent(new Event('fullscreenchange'));
    expect(states).toEqual([false, true, false]);
    sub.unsubscribe();
    model.dispose();
  });

  it('of two viewers on a page: each follows its own element only', () => {
    allow(true);
    const a = document.createElement('div');
    const b = document.createElement('div');
    document.body.append(a, b);
    const one = new FullscreenModel(() => a);
    const two = new FullscreenModel(() => b);
    const change = (element: Element | null) => {
      Object.defineProperty(document, 'fullscreenElement', { value: element, configurable: true });
      (element ?? a).dispatchEvent(new Event('fullscreenchange', { bubbles: true }));
    };
    change(a);
    expect([one.active.value, two.active.value]).toEqual([true, false]);
    // The browser moves fullscreen to the other: the first is no longer in it.
    change(b);
    expect([one.active.value, two.active.value]).toEqual([false, true]);
    change(null);
    expect([one.active.value, two.active.value]).toEqual([false, false]);
    one.dispose();
    two.dispose();
    a.remove();
    b.remove();
  });

  it('the button of one enters fullscreen while another is in it, and leaves only its own', async () => {
    allow(true);
    const a = document.createElement('div');
    const b = document.createElement('div');
    document.body.append(a, b);
    const requested: Element[] = [];
    let exited = 0;
    (a as any).requestFullscreen = async () => void requested.push(a);
    (b as any).requestFullscreen = async () => void requested.push(b);
    (document as any).exitFullscreen = async () => void exited++;
    const first = new FullscreenModel(() => a);
    const two = new FullscreenModel(() => b);
    Object.defineProperty(document, 'fullscreenElement', { value: a, configurable: true });
    await two.toggle();
    expect(requested).toEqual([b]);
    expect(exited).toBe(0);
    await first.toggle();
    expect(exited).toBe(1);
    Object.defineProperty(document, 'fullscreenElement', { value: null, configurable: true });
    first.dispose();
    two.dispose();
    a.remove();
    b.remove();
  });
});
