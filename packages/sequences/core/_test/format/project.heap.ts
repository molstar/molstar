/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the projection core (VW5; seqspec:interaction#channels: after a route's first use nothing allocates),
// on RB0b's stepped `heapScaling` (a step is one hover change; untimed, `heapScaling` warms by count): `project` through
// a codon route split by an intron, its inverse image, links across sections, a source used twice and a row through its
// binding; `close` by the hull search (a nucleotide) and by containment (a residue), then the closed extents projected.
// The cases are project-cases.ts's, shared with the quick twin (project.quick.ts).

import { beforeAll, describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../model/heap.js';
import { type ProjectCases, projectCases } from './project-cases.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 50_000;
let cases: ProjectCases;

beforeAll(async () => {
  cases = await projectCases();
}, 60_000);

describe('the projection allocates nothing per move', () => {
  it('project: residues onto their codons, nucleotides onto residues, across the link to the chromosome', () => {
    const { project } = cases.steps;
    const result = heapScaling(N, project);
    expect(cases.visited()).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(N, leaking(project)))).toBe(true);
  });

  it('project: a source used twice, and a row through its binding to its chain', () => {
    const result = heapScaling(N, cases.steps.twice);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it('close: a nucleotide by the hull search, a residue by containment, then projected', () => {
    const result = heapScaling(N, cases.steps.close);
    expect(cases.added()).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });
});
