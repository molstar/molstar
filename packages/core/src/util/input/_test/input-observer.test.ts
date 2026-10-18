// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { InputObserver } from '../input-observer.js';

function setup() {
  const host = document.createElement('div');
  document.body.appendChild(host);
  const root = host.attachShadow({ mode: 'open' });
  const parent = document.createElement('div');
  const element = document.createElement('div');
  parent.appendChild(element);
  root.appendChild(parent);
  element.getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
  const other = document.createElement('input');
  root.appendChild(other);
  return { host, element, other, observer: InputObserver.fromElement(element) };
}

describe('input-observer in a shadow root', () => {
  afterEach(() => {
    document.body.innerHTML = '';
  });

  it('reports a move from inside the shadow root as on the element', () => {
    const { element, other, observer } = setup();
    const moves: boolean[] = [];
    observer.move.subscribe((m) => moves.push(m.onElement));

    // the window listener sees the event retargeted to the host
    element.dispatchEvent(new MouseEvent('mousemove', { clientX: 10, clientY: 10, bubbles: true, composed: true }));
    other.dispatchEvent(new MouseEvent('mousemove', { clientX: 10, clientY: 10, bubbles: true, composed: true }));
    expect(moves).toEqual([true, false]);
    observer.dispose();
  });

  it('reports keys from inside the shadow root as on the element', () => {
    const { element, other, observer } = setup();
    const keys: string[] = [];
    observer.keyDown.subscribe((k) => keys.push(k.key));

    element.dispatchEvent(new MouseEvent('mousemove', { clientX: 10, clientY: 10, bubbles: true, composed: true }));
    element.dispatchEvent(new KeyboardEvent('keydown', { key: 'a', bubbles: true, composed: true }));
    other.dispatchEvent(new KeyboardEvent('keydown', { key: 'b', bubbles: true, composed: true }));
    expect(keys).toEqual(['a']);
    observer.dispose();
  });

  it('keeps the behavior outside shadow roots', () => {
    const element = document.createElement('div');
    document.body.appendChild(element);
    element.getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    const observer = InputObserver.fromElement(element);
    const moves: boolean[] = [];
    observer.move.subscribe((m) => moves.push(m.onElement));
    element.dispatchEvent(new MouseEvent('mousemove', { clientX: 10, clientY: 10, bubbles: true }));
    document.body.dispatchEvent(new MouseEvent('mousemove', { clientX: 10, clientY: 10, bubbles: true }));
    expect(moves).toEqual([true, false]);
    observer.dispose();
  });
});
