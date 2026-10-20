/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import { computeCarbohydrates } from '../compute.js';
import type { Carbohydrates } from '../data.js';
import { glycanCif, type GlycanParams } from './glycan.js';

// Pins what `computeCarbohydrates` gives: elements (unit, residue, ring, alt id, component, geometry), links in order,
// terminal links and partial elements, on 7qpd's glycans and synthetic polysaccharides (`glycan.ts`). The values were
// recorded before the 05-model-structure P9 change (.v6/perf-audit/05-model-structure.md), which must keep them:
// link order and the link directions, which later links overwrite, included.

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

async function load(cif: string) {
  const parsed = await CIF.parse(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

const glycan = (p: GlycanParams) => load(glycanCif(p));

function hash(xs: ArrayLike<number>) {
  return (hashFnv32a(xs) >>> 0).toString(16);
}

const round = (v: ArrayLike<number>) => Array.from(v, (x) => Math.round(x * 1e5));

function signature(c: Carbohydrates) {
  const elements = c.elements.map((e) => [e.unit.id, e.residueIndex, e.ringIndex, e.altId.charCodeAt(0) || 0]);
  const components = c.elements.map((e) => e.component.abbr).join(',');
  const geometry = c.elements.flatMap((e) => [
    ...round(e.geometry.center),
    ...round(e.geometry.normal),
    ...round(e.geometry.direction),
  ]);
  const links = c.links.flatMap((l) => [l.carbohydrateIndexA, l.carbohydrateIndexB]);
  const terminal = c.terminalLinks.flatMap((l) => [
    l.carbohydrateIndex,
    l.elementIndex,
    l.elementUnit.id,
    l.fromCarbohydrate ? 1 : 0,
  ]);
  const partial = c.partialElements.map((e) => [e.unit.id, e.residueIndex]);
  return [
    `elements=${elements.length}:${hash(elements.flat())} components=${hash(Array.from(components, (ch) => ch.charCodeAt(0)))}`,
    `geometry=${hash(geometry)} links=${c.links.length}:${hash(links)} terminal=${c.terminalLinks.length}:${hash(terminal)}`,
    `partial=${partial.length}:${hash(partial.flat())}`,
  ].join(' ');
}

/** Links as `a-b` pairs in order. */
function linkList(c: Carbohydrates) {
  return c.links.map((l) => `${l.carbohydrateIndexA}-${l.carbohydrateIndexB}`).join(' ');
}

/**
 * The links between rings of different residues of a unit, in order, as the code before P9 found them: every pair of
 * elements tested with `areVertexSetsConnected` at 3 bonds.
 */
function allPairsLinks(c: Carbohydrates) {
  const ret: string[] = [];
  for (let i = 0; i < c.elements.length; i++) {
    const a = c.elements[i];
    for (let j = i + 1; j < c.elements.length; j++) {
      const b = c.elements[j];
      if (a.unit !== b.unit || a.residueIndex === b.residueIndex) continue;
      const { all } = a.unit.rings;
      if (IntAdjacencyGraph.areVertexSetsConnected(a.unit.bonds, all[a.ringIndex], all[b.ringIndex], 3)) {
        ret.push(`${i}-${j}`, `${j}-${i}`);
      }
    }
  }
  return ret.join(' ');
}

/** `linkList` restricted to links between different residues of a unit. */
function intraUnitLinks(c: Carbohydrates) {
  return c.links
    .filter(({ carbohydrateIndexA: a, carbohydrateIndexB: b }) => {
      const ea = c.elements[a],
        eb = c.elements[b];
      return ea.unit === eb.unit && ea.residueIndex !== eb.residueIndex;
    })
    .map((l) => `${l.carbohydrateIndexA}-${l.carbohydrateIndexB}`)
    .join(' ');
}

describe('computeCarbohydrates', () => {
  it('7qpd', async () => {
    const s = await load(fs.readFileSync(path.join(examples, '7qpd.fw2.cif'), 'utf8'));
    const c = computeCarbohydrates(s);
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=10:a96b4e71 components=5f0db5d3 geometry=f0e861 links=16:6e0e8fbd terminal=4:ab74179 partial=0:811c9dc5"`,
    );
    expect(linkList(c)).toMatchInlineSnapshot(`"0-1 1-0 1-2 2-1 3-4 4-3 4-5 5-4 5-6 6-5 6-7 7-6 7-8 8-7 8-9 9-8"`);
  });

  it('a branched polysaccharide in one chain, with links to other chains', async () => {
    const s = await glycan({ strands: 2, length: 10, branchEvery: 3, otherChains: true });
    const c = computeCarbohydrates(s);
    expect(s.units.map((u) => `${u.id}:${u.elements.length}`).join(' ')).toMatchInlineSnapshot(`"0:225 1:9 2:2"`);
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=26:22626df7 components=1c91b563 geometry=f36dfd4c links=48:a5211c9 terminal=2:2ca66f22 partial=0:811c9dc5"`,
    );
    expect(linkList(c)).toMatchInlineSnapshot(
      `"0-1 1-0 1-2 2-1 2-3 3-2 3-4 4-3 4-5 5-4 5-6 6-5 6-7 7-6 6-8 8-6 8-9 9-8 9-10 10-9 10-11 11-10 12-13 13-12 13-14 14-13 14-15 15-14 15-16 16-15 15-17 17-15 17-18 18-17 18-19 19-18 19-20 20-19 19-21 21-19 21-22 22-21 22-23 23-22 23-24 24-23 3-25 25-3"`,
    );
    expect(
      c.terminalLinks.map((l) => `${l.carbohydrateIndex}:${l.elementUnit.id}/${l.elementIndex}:${l.fromCarbohydrate}`),
    ).toMatchInlineSnapshot(`
      [
        "10:2/0:true",
        "10:2/0:false",
      ]
    `);
    // the lookups by atom
    const unit = c.elements[1].unit;
    const ring = unit.rings.all[c.elements[1].ringIndex];
    const atom = unit.elements[ring[0]];
    expect([
      c.getElementIndices(unit, atom),
      c.getLinkIndices(unit, atom),
      c.getTerminalLinkIndices(unit, atom),
    ]).toMatchInlineSnapshot(`
      [
        [
          1,
        ],
        [
          1,
          2,
        ],
        [],
      ]
    `);
  });

  it('alternate locations', async () => {
    const c = computeCarbohydrates(await glycan({ strands: 2, length: 8, branchEvery: 3, altEvery: 3 }));
    expect(c.elements.map((e) => e.altId || '.').join('')).toMatchInlineSnapshot(`".AB...AB...AB.AB...AB...AB"`);
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=26:7520749a components=5f71111b geometry=14012505 links=56:c240c035 terminal=0:811c9dc5 partial=0:811c9dc5"`,
    );
    expect(linkList(c)).toMatchInlineSnapshot(
      `"0-1 1-0 0-2 2-0 1-3 3-1 2-3 3-2 3-4 4-3 4-5 5-4 4-6 6-4 4-7 7-4 6-8 8-6 7-8 8-7 8-9 9-8 9-10 10-9 9-11 11-9 9-12 12-9 13-14 14-13 13-15 15-13 14-16 16-14 15-16 16-15 16-17 17-16 17-18 18-17 17-19 19-17 17-20 20-17 19-21 21-19 20-21 21-20 21-22 22-21 22-23 23-22 22-24 24-22 22-25 25-22"`,
    );
  });

  it('rings 4 bonds apart are not linked', async () => {
    const c = computeCarbohydrates(await glycan({ strands: 1, length: 8, branchEvery: 3, spacer: true }));
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=10:1cee7033 components=e2a7ae0f geometry=7b304f65 links=14:d3a6b61b terminal=0:811c9dc5 partial=0:811c9dc5"`,
    );
    expect(linkList(c)).toMatchInlineSnapshot(`"0-1 1-0 1-2 2-1 2-3 3-2 3-5 5-3 5-6 6-5 6-7 7-6 7-9 9-7"`);
  });

  it('links between residues of a unit are all pairs at most 3 bonds apart, in order', async () => {
    const s7qpd = await load(fs.readFileSync(path.join(examples, '7qpd.fw2.cif'), 'utf8'));
    const structures = [
      s7qpd,
      await glycan({ strands: 2, length: 10, branchEvery: 3, otherChains: true }),
      await glycan({ strands: 2, length: 8, branchEvery: 3, altEvery: 3 }),
      await glycan({ strands: 2, length: 9, branchEvery: 2, spacer: true }),
    ];
    for (const s of structures) {
      const c = computeCarbohydrates(s);
      expect(intraUnitLinks(c)).toEqual(allPairsLinks(c));
    }
  });

  it('strands without branches, several in a chain', async () => {
    const c = computeCarbohydrates(await glycan({ strands: 3, length: 30, branchEvery: 0 }));
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=90:62fcd813 components=e9d2f147 geometry=61666e2d links=174:44e7d4f terminal=0:811c9dc5 partial=0:811c9dc5"`,
    );
  });

  it('a large polysaccharide in one chain', async () => {
    const c = computeCarbohydrates(await glycan({ strands: 20, length: 100, branchEvery: 3, otherChains: true }));
    expect(signature(c)).toMatchInlineSnapshot(
      `"elements=2660:cec615cf components=364bf603 geometry=3c361c12 links=5280:313fc561 terminal=2:286d66a2 partial=0:811c9dc5"`,
    );
  });

  // scaling, not a pin: every pair of ring elements of a unit was searched (quadratic: 16x for 4x)
  it('the links of a polysaccharide scale near-linearly', { tags: ['heavy'] }, async () => {
    const time = async (strands: number) => {
      const s = await glycan({ strands, length: 100, branchEvery: 3 });
      for (const u of s.units) if (Unit.isAtomic(u)) void u.rings;
      void s.interUnitBonds;
      let best = Infinity;
      for (let i = 0; i < 3; i++) {
        const t0 = performance.now();
        computeCarbohydrates(s);
        best = Math.min(best, performance.now() - t0);
      }
      return best;
    };
    const small = await time(5);
    const large = await time(20);
    expect(large / small).toBeLessThan(10);
  });
});
