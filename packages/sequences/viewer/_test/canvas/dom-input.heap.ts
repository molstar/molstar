/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The DOM engine's moves during a press allocate nothing (seqspec:interaction#input-events; seqspec:viewer#performance):
// a press subscribes the document's moves, each read in the pressing surface's px from its rect (read at the press and
// on a scroll, never per move) into the one reused record. The same event dispatched again each step (an event is the
// page's allocation, not the engine's). The positive control: the same moves allocating garbage each, sized to the
// run (`leaking`), measured alike, are caught.

import { afterEach, expect, it } from 'vitest';
import { DomInput } from '../../src/canvas/dom-input.js';
import type { InputSink } from '../../src/index.js';
import { controlCaught, heapScaling, isFlat, leaking } from './heap-browser.js';

const hosts: HTMLElement[] = [];
afterEach(() => {
  for (const h of hosts.splice(0)) h.remove();
});

function pressed() {
  const el = document.createElement('div');
  el.style.cssText = 'position:absolute;left:100px;top:50px;width:300px;height:200px';
  document.body.appendChild(el);
  hosts.push(el);
  let moves = 0;
  const sink: InputSink = {
    pointer() {
      moves++;
    },
    wheel: () => false,
    key: () => false,
    focus() {},
    leave() {},
    cancel() {},
  };
  const input = new DomInput(el, sink, () => 0);
  const init = { clientX: 140, clientY: 90, isPrimary: true, pointerId: 1, bubbles: true, cancelable: true };
  el.dispatchEvent(new PointerEvent('pointerdown', { ...init, buttons: 1, button: 0 }));
  // Off the surface, over the page: the document's listener has it.
  const move = new PointerEvent('pointermove', { ...init, clientX: 600, clientY: 500, buttons: 1 });
  return { input, step: () => document.body.dispatchEvent(move), moves: () => moves };
}

it("a press's moves over the page, each read into the one record: flat", () => {
  const { input, step, moves } = pressed();
  const r = heapScaling(20_000, step);
  expect(moves()).toBeGreaterThan(100_000);
  expect(isFlat(r), JSON.stringify(r)).toBe(true);
  input.stop(false);
});

it('the same moves allocating garbage sized to the run each are caught (positive control)', () => {
  const { input, step } = pressed();
  expect(controlCaught(() => heapScaling(20_000, leaking(step, 20_000)))).toBe(true);
  input.stop(false);
});
