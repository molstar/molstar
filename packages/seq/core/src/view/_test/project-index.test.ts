/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A projection's work is in its parts and their links, never in every system the view shows times every part
// (_spec/README.md "Resolved view", Projection): a genome of 3,000 blocks on a ring and a row shows 3,000 systems, and a
// drag across 12 chromosomes projects 12 parts a frame. The view's systems and links are indexed once a resolved view;
// what is carried equals the walk over every system (the oracle below, the rule as written).

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { EmptyLocation, type Location, type LocationPart, location } from '@molstar/seqx-core/query';
import { loadGenome, projectRanges, resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { toyGenome } from '../../render/_test/genome-fixtures.js';
import {
  type ProjectOptions,
  projectOntoView,
  type RegionProjector,
  type ResolvedView,
  resolveView,
} from '../index.js';
import { p53GeneSpec } from '../../render/_test/fixtures.js';
import { moleculesSpec } from '../../render/_test/molecules-fixtures.js';

beforeAll(async () => {
  await loadGenome();
});

/**
 * The rule as written (the walk the index replaced, kept here as the oracle): every system a section shows, every part,
 * the view's link from the part's system, else (no section shows it) the model's static path, else a map read by region.
 */
function oracle(model: StaticModel, view: ResolvedView, loc: Location, options: ProjectOptions = {}): Location {
  if (!loc.parts.length && !loc.rects.length) return EmptyLocation;
  const systems = [
    ...new Set(view.sections.flatMap((s) => (s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system]))),
  ];
  const owned = new Set(loc.parts.map((p) => p.system));
  const flat = loc.rects.length ? location({ parts: [...loc.parts, ...loc.rects.flatMap((r) => [r.a, r.b])] }) : loc;
  const src = options.close ? model.close(flat, systems) : flat;
  // each part's links and whether a section shows its system (what the walk asks of every system, read once)
  const linksOf = src.parts.map((part) => view.links.filter((l) => l.fromSystem === part.system));
  const shownOf = src.parts.map((part) => systems.includes(part.system));
  const parts: { system: string; ranges: ArrayLike<number> }[] = [];
  for (const system of systems) {
    if (!options.own && owned.has(system)) continue;
    for (let k = 0; k < src.parts.length; k++) {
      const part = src.parts[k];
      if (part.system === system) {
        parts.push(part);
        const out: number[] = [];
        if (options.linked !== false && options.regions?.project(model, system, part.ranges, system, out))
          if (out.length) parts.push({ system, ranges: HalfOpenRanges.normalize(out) });
        continue;
      }
      if (options.linked === false) continue;
      const links = linksOf[k];
      const path = links.length
        ? links.find((l) => l.toSystem === system)?.path
        : shownOf[k] || !model.systems.has(part.system) || !model.systems.has(system)
          ? undefined
          : model.path(part.system, system).path;
      if (path) parts.push({ system, ranges: projectRanges(path, part.ranges, []).slice() });
      else if (options.regions) {
        const out: number[] = [];
        if (options.regions.project(model, part.system, part.ranges, system, out) && out.length)
          parts.push({ system, ranges: HalfOpenRanges.normalize(out) });
      }
    }
  }
  return parts.length ? location({ parts }) : EmptyLocation;
}

/** A view whose `sections` and `links` count their reads. */
function counted(view: ResolvedView) {
  const reads = { sections: 0, links: 0 };
  const proxy = new Proxy(view, {
    get(target, key, receiver) {
      if (key === 'sections') reads.sections++;
      if (key === 'links') reads.links++;
      return Reflect.get(target, key, receiver);
    },
  });
  return { view: proxy, reads };
}

describe('a projection over a view of 3,000 blocks', () => {
  it('a drag across 12 chromosomes: indexed once, then read in its parts alone; equal to the walk over every system', () => {
    const b = createSeqSpec();
    const g = b.data.assembly('g', assembly('toy', toyGenome(24, 2976, 7), { order: 'source' } as never));
    b.views.view('v', (v) => {
      genome(v, 'ring', g, { layout: 'rings' });
      genome(v, 'row', g);
    });
    const model = resolveSpec(b.build());
    const resolved = resolveView(model, '/views/v')!;
    expect(resolved.sections[1].blocks!.length).toBe(3000);
    const { view, reads } = counted(resolved);
    // chr11 to chr22, a part a chromosome, as a drag's frames grow it
    const across = (n: number) =>
      location({
        parts: Array.from({ length: n }, (_, k) => ({
          system: resolved.sections[1].blocks![10 + k].axis.system,
          ranges: [0, 1000 + n],
        })),
      });
    projectOntoView(model, view, across(1));
    const indexed = { ...reads };
    for (let n = 2; n <= 12; n++) {
      for (const own of [false, true]) {
        const loc = across(n);
        expect(projectOntoView(model, view, loc, { own })).toEqual(oracle(model, resolved, loc, { own }));
      }
    }
    // the frames after the first read neither the sections nor the links: the index holds them
    expect(reads).toEqual(indexed);
  });

  it('a complex’s blocks and a sequence linked to them: carried as the walk carries them, both ways', () => {
    const model = resolveSpec(moleculesSpec());
    const view = resolveView(model, '/views/v')!;
    for (const s of view.sections)
      for (const system of s.blocks ? s.blocks.map((x) => x.axis.system) : [s.axis.system]) {
        const loc = location({ parts: [{ system, ranges: [0, 2] }] });
        for (const own of [false, true])
          expect(projectOntoView(model, view, loc, { own })).toEqual(oracle(model, view, loc, { own }));
      }
  });

  it('rects, `close`, a map read by region, reversed and collapsed axes: as the walk carries them', () => {
    // p53 on its gene (a codon map, the minus strand), the gene's introns collapsed (`pieces: coding`) or not; and the
    // gene shown reversed beside it
    const reversed = () => {
      const b = createSeqSpec();
      const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
      const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: 'ACGT'.repeat(30) });
      b.data.map('m', { from: p53, to: gene, method: 'codon', exons: { start: [3, 60], end: [40, 105] }, strand: '-' });
      b.views.view('main', (v) => {
        v.section('protein', { of: p53 }).track('t').layer('l', { type: 'letters' });
        v.section('gene', { of: gene, direction: 'reversed' } as never)
          .track('t')
          .layer('l', { type: 'letters' });
      });
      return b.build();
    };
    // a map read by region, faked: p53 carried onto the gene shifted, and each part's other image on itself
    const regions: RegionProjector = {
      version: { value: 0 },
      reads: () => true,
      project(_model, from, ranges, to, out) {
        if (!from.endsWith('p53') && !to.endsWith('p53')) return false;
        for (let i = 0; i < ranges.length; i++) out.push(ranges[i] + (from === to ? 7 : 2));
        return true;
      },
    };
    let checked = 0;
    const kinds = new Set<string>();
    for (const spec of [p53GeneSpec({ coding: true }), p53GeneSpec(), reversed()]) {
      const model = resolveSpec(spec);
      const view = resolveView(model, '/views/main')!;
      const [p, g] = view.sections.map((s) => s.axis.system);
      // what each fixture is for: the gene collapsed to its coding pieces, or shown reversed
      kinds.add(`${view.sections[1].axis.pieces} ${view.sections[1].axis.dir}`);
      const locs = [
        location({ parts: [{ system: p, ranges: [3, 9] }] }),
        location({ parts: [{ system: g, ranges: [10, 20, 40, 41] }] }),
        location({
          parts: [
            { system: p, ranges: [1, 2] },
            { system: g, ranges: [50, 70] },
          ],
        }),
        location({ rects: [{ a: { system: p, ranges: [2, 6] }, b: { system: g, ranges: [12, 30] } }] }),
        location({
          parts: [{ system: g, ranges: [5, 6] }],
          rects: [{ a: { system: p, ranges: [0, 1] }, b: { system: p, ranges: [4, 5] } }],
        }),
      ];
      for (const loc of locs)
        for (const own of [false, true])
          for (const close of [false, true])
            for (const linked of [true, false])
              for (const r of [undefined, regions]) {
                const options: ProjectOptions = { own, close, linked, ...(r ? { regions: r } : {}) };
                expect(
                  projectOntoView(model, view, loc, options),
                  JSON.stringify({ own, close, linked, r: !!r }),
                ).toEqual(oracle(model, view, loc, options));
                checked++;
              }
    }
    expect(checked).toBe(3 * 5 * 16);
    // (the minus strand's follower and the explicit `reversed` both read -1; `coding` the introns collapsed)
    expect([...kinds].sort()).toEqual(['coding -1', 'whole -1']);
  });
});
