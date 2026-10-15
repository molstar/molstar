/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of channels.heap.ts (`test:perf`): a hover element change through the hover manager and the channels'
// rebuild, over the gene fixture and over an alignment's 10,000 rows, the guard's steps (channels-cases.ts), each
// with its positive control.

import { expect, it } from 'vitest';
import { alignment, genes, sunk } from './channels-cases.js';
import { quickCase } from './heap.js';

it('a nucleotide or a residue hovered, closed, projected and read', async ({ task }) => {
  const { step, ch } = await genes();
  quickCase(task, 'channels: a nucleotide or a residue hovered', step);
  expect(ch.rebuilds.hover).toBeGreaterThan(1_000);
  expect(sunk()).not.toBe(0);
});

it('element changes over 10,000 distinct rows and positions of an alignment', async ({ task }) => {
  const { step, ch } = await alignment();
  // Every row is first used once (its cached map allocates on a miss): a step reaches all 10,000 rows in 10,000 steps.
  quickCase(task, 'channels: an alignment over 10,000 rows', step, { warmSteps: 2_000 });
  expect(ch.rebuilds.hover).toBeGreaterThan(1_000);
  expect(sunk()).not.toBe(0);
});
