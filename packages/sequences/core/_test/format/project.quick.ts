/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of project.heap.ts (`test:perf`): a hover change projected (project-cases.ts, the guard's cases), each
// with its positive control.

import { beforeAll, describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { type ProjectCases, projectCases } from './project-cases.js';

let cases: ProjectCases;

beforeAll(async () => {
  cases = await projectCases();
}, 60_000);

describe('the projection per move, quick', () => {
  it('project: residues onto their codons, nucleotides onto residues, across the link', async ({ task }) => {
    quickCase(task, 'project: codons, residues and the link to the chromosome', cases.steps.project);
    expect(cases.visited()).toBeGreaterThan(0);
  });

  it('project: a source used twice, and a row through its binding', async ({ task }) => {
    quickCase(task, 'project: a source used twice and a row', cases.steps.twice);
  });

  it('close: by the hull search and by containment, then projected', async ({ task }) => {
    quickCase(task, 'project: close, then projected', cases.steps.close);
    expect(cases.added()).toBeGreaterThan(0);
  });
});
