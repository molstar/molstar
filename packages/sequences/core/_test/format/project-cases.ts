/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The projection core's cases, for their heap guard and its quick twin (project.heap.ts, project.quick.ts): a step is one
// hover change, built over project-fixtures.ts's views. `project` through a codon route split by an intron, its inverse
// image and links across sections, a source used twice and a row through its binding; `close` by the hull search (a
// nucleotide) and by containment (a residue), then the closed extents projected.

import { Extents } from '@molstar/seq-core/model';
import { constructs, gene, open, rows } from './project-fixtures.js';

/** What the steps' projections visited and what `close` added: a guard over steps that find nothing would pass blind. */
export interface ProjectCases {
  /** The cases: residues onto their codons and across the link, a source used twice and a row, `close`. */
  steps: { project: (i: number) => void; twice: (i: number) => void; close: (i: number) => void };
  visited: () => number;
  added: () => number;
}

let loaded: Promise<ProjectCases> | undefined;

/** The views opened once for every case (a step allocates nothing, so the views are shared). */
export function projectCases(): Promise<ProjectCases> {
  loaded ??= build();
  return loaded;
}

async function build(): Promise<ProjectCases> {
  const [g, c, r] = await Promise.all([open(gene({ chromosome: true })), open(constructs()), open(rows())]);
  const sink = {
    n: 0,
    add(s: number, t: number) {
      this.n += t - s;
    },
  };
  let added = 0;

  const ext = Extents.create(4);
  const [protein, unit] = [g.system('protein'), g.system('gene/unit')];
  const blocks = [g.block('gene'), g.block('protein'), g.block('chr')];
  const project = (i: number) => {
    ext.count = 0;
    Extents.add(ext, protein, i % 20, (i % 20) + 1);
    Extents.add(ext, unit, i % 120, (i % 120) + 1);
    for (let b = 0; b < blocks.length; b++) g.p.project(ext, blocks[b], sink);
  };

  const ext2 = Extents.create(4);
  const [ub, diub] = [c.system('ub'), c.block('diub')];
  const [row, x, a] = [r.system('a/0'), r.block('x'), r.block('a')];
  const rowExt = Extents.create(4);
  const twice = (i: number) => {
    ext2.count = 0;
    Extents.add(ext2, ub, i % 76, (i % 76) + 1);
    c.p.project(ext2, diub, sink);
    rowExt.count = 0;
    Extents.add(rowExt, row, i % 4, (i % 4) + 1);
    r.p.project(rowExt, x, sink);
    r.p.project(rowExt, a, sink);
  };

  const ext3 = Extents.create(8);
  const [chr, gnu, block] = [g.system('chr/unit'), g.system('gene/unit'), g.block('chr')];
  const close = (i: number) => {
    ext3.count = 0;
    Extents.add(ext3, i % 2 ? chr : gnu, 100 * (i % 2) + (i % 120), 100 * (i % 2) + (i % 120) + 1);
    Extents.add(ext3, protein, i % 20, (i % 20) + 1);
    g.p.close(ext3);
    added += ext3.count - 2;
    g.p.project(ext3, block, sink);
  };

  return { steps: { project, twice, close }, visited: () => sink.n, added: () => added };
}
