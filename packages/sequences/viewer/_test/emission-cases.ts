/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the emission guards (seqspec:viewer#performance; emission.heap.ts) and their quick twin
// (emission.quick.ts), each a step a set, a bump or an emission: state sets with two rxjs subscribers, sets queued
// during delivery, live version bumps, and events with a reused record by priority. The heap guard and the quick twin
// share them.

import { event, live, state } from '../src/state.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 20_000;
const MASK = 0xfffff;

let sink = 0;
/** What the steps summed, for the guards' non-vacuity. */
export const sunk = () => sink;

/** A state with two rxjs subscribers, and its step: one set. */
export function setting(): (k: number) => void {
  const s = state(0);
  s.subscribe((v) => {
    sink = (sink + v) & MASK;
  });
  s.subscribe((v) => {
    sink = (sink ^ v) & MASK;
  });
  return (k) => s.set(k & MASK, 'pointer');
}

/** Sets queued during delivery: an odd value sets the next one from inside the delivery, so each outer set queues one. */
export function queued(): (k: number) => void {
  const s = state(0);
  s.subscribe((v) => {
    if (v & 1) s.set((v + 1) & MASK, 'behavior');
  });
  s.subscribe((v) => {
    sink = (sink + v) & MASK;
  });
  return (k) => s.set((2 * k + 1) & MASK, 'pointer');
}

/** Live version bumps with a subscriber. */
export function bumping(): (k: number) => void {
  const l = live({ x: 0 });
  l.version.subscribe((v) => {
    sink = (sink + v + l.current.x) & MASK;
  });
  return (k) => {
    l.current.x = k & MASK;
    l.bump();
  };
}

/** Events with a reused record, by priority. */
export function emitting(): (k: number) => void {
  const e = event<{ x: number; owner: string }>();
  const rec = { x: 0, owner: '' };
  e.subscribe(
    (r) => {
      sink = (sink + r.x) & MASK;
    },
    { priority: 100 },
  );
  e.subscribe(
    (r) => {
      if (r.x === -1) r.owner = 'never';
    },
    { priority: 50 },
  );
  e.subscribe((r) => {
    sink = (sink ^ r.x) & MASK;
  });
  return (k) => {
    rec.x = k & MASK;
    e.emit(rec);
  };
}
