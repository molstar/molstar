/** @jest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
});
