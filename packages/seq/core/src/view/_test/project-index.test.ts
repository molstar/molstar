/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A projection's work is in its parts and their links, never in every system the view shows times every part
// (_spec/README.md "Resolved view", Projection): a genome of 3,000 blocks on a ring and a row shows 3,000 systems, and a
// drag across 12 chromosomes projects 12 parts a frame. The view's systems and links are indexed once a resolved view;
// what is carried equals the walk over every system (the oracle below, the rule as written).

import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { EmptyLocation, type Location, type LocationPart, location } from '@molstar/seqx-core/query';
import { loadGenome, projectRanges, resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { toyGenome } from '../../render/_test/genome-fixtures.js';
import { projectOntoView, type ResolvedView, resolveView } from '../index.js';
import { moleculesSpec } from '../../render/_test/molecules-fixtures.js';

beforeAll(async () => {
  await loadGenome();
});

/** The rule as written: every system a section shows, every part, the view's link else the model's static path. */
function oracle(model: StaticModel, view: ResolvedView, loc: Location, own: boolean): Location {
  const systems = [
    ...new Set(view.sections.flatMap((s) => (s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system]))),
  ];
  const owned = new Set(loc.parts.map((p) => p.system));
  const flat = loc.rects.length ? location({ parts: [...loc.parts, ...loc.rects.flatMap((r) => [r.a, r.b])] }) : loc;
  // each part's links and whether a section shows its system (what the walk asks of every system, read once)
  const linksOf = flat.parts.map((part) => view.links.filter((l) => l.fromSystem === part.system));
  const shownOf = flat.parts.map((part) => systems.includes(part.system));
  const parts: LocationPart[] = [];
  for (const system of systems) {
    if (!own && owned.has(system)) continue;
    for (let k = 0; k < flat.parts.length; k++) {
      const part = flat.parts[k];
      if (part.system === system) {
        parts.push(part);
        continue;
      }
      const links = linksOf[k];
      const shown = shownOf[k];
      const path = links.length
        ? links.find((l) => l.toSystem === system)?.path
        : shown || !model.systems.has(part.system) || !model.systems.has(system)
          ? undefined
          : model.path(part.system, system).path;
      if (path) parts.push({ system, ranges: projectRanges(path, part.ranges, []).slice() });
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
        expect(projectOntoView(model, view, loc, { own })).toEqual(oracle(model, resolved, loc, own));
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
          expect(projectOntoView(model, view, loc, { own })).toEqual(oracle(model, view, loc, own));
      }
  });
});
