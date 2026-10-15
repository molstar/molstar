/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of emission.heap.ts (`test:perf`): state sets, sets queued during delivery, live version bumps and
// events by priority, the guard's steps (emission-cases.ts), each with its positive control.

import { expect, it } from 'vitest';
import { bumping, emitting, queued, setting, sunk } from './emission-cases.js';
import { quickCase } from './heap.js';

it('state sets, with rxjs subscribers', ({ task }) => {
  quickCase(task, 'emission: state sets', setting());
});

it('state sets queued during delivery', ({ task }) => {
  quickCase(task, 'emission: state sets queued during delivery', queued());
});

it('live version bumps', ({ task }) => {
  quickCase(task, 'emission: live version bumps', bumping());
});

it('events with a reused record, by priority', ({ task }) => {
  quickCase(task, 'emission: events by priority', emitting());
  expect(sunk()).toBeGreaterThanOrEqual(0);
});
