/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The marks' and the projector's per-move updates, for their heap guard and its quick twin (marks.heap.ts,
// marks.quick.ts): a hover element closed and projected into every section's hover set with its primary and the
// overlay's queries (`hover`); loci and items converted and projected (`loci`); `Coordinates` through a registered
// translator (`coordinates`); 200 runs out of order sorted at close and read (`sorting`); row rectangles and a ligand's
// atoms (`rowRects`). `sink` is what the steps read, for a check that they read something.

import { CoordinateSystemKind as K, type Loci, OrderedSet, type SeqIndex, type UnitId } from '@molstar/seq-core/model';
import { type SystemDescriptor, TranslateMaps, Translator } from '@molstar/seq-core/translate';
import { case7, case12, case14 } from '../../core/_test/format/cases.js';
import { open as core, gene, rows } from '../../core/_test/format/project-fixtures.js';
import { Axis } from '../src/axis.js';
import { MarkBit } from '../src/index.js';
import { MarkSet, RowMarks, SectionMarks } from '../src/marks.js';
import { entries, over, shown } from './projection-fixtures.js';

let sink = 0;

/** 200 runs out of order into `m`, sorted and merged at close, then read by `has`. */
export const sorting = (m: MarkSet) => (i: number) => {
  m.reset(2);
  m.begin(1);
  for (let k = 199; k >= 0; k--) m.add(3 * k + (i % 3), 3 * k + 2);
  m.close(0, 600);
  for (let x = 0; x < 600; x += 7) sink += m.has(1, x) ? 1 : 0;
};

/** The projections opened, and a step per update. */
export async function marksCases() {
  const [a, b] = await Promise.all([core(gene({ chromosome: true })), core(rows())]);
  const g = { ...over(a.view, a.set), system: a.system };
  const r = { ...over(b.view, b.set), system: b.system };
  const [c14, c7, c12] = await Promise.all([shown(case14()), shown(case7(), 'complex-each'), shown(case12())]);

  /** A hover: closed, projected into each section's hover set, its primary, the overlay's queries. */
  const hover = () => {
    const [protein, unit, chr] = [g.system('protein'), g.system('gene/unit'), g.system('chr/unit')];
    const sections = Array.from({ length: g.view.sections.count }, (_, s) => ({
      marks: new SectionMarks(),
      axis: new Axis(g.view, g.set, s, 28, 3),
      blocks: g.view.sections.blockCount[s],
    }));
    const { p, ext } = g;
    return (i: number) => {
      ext.count = 0;
      const s = i % 3 === 0 ? protein : i % 3 === 1 ? unit : chr;
      const x = s === protein ? i % 20 : s === unit ? i % 120 : 100 + (i % 120);
      ext.kind[0] = s.kind;
      ext.id[0] = s.id;
      ext.row[0] = s.row;
      ext.start[0] = x;
      ext.stop[0] = x + 1;
      ext.count = 1;
      p.close(ext);
      for (let k = 0; k < sections.length; k++) {
        const S = sections[k];
        p.project(ext, k, S.marks.hover);
        S.marks.updatePrimary(S.axis, S.blocks);
        sink += S.marks.bitsAt(0, x) + S.marks.primary(0) + (S.marks.any(0, 40, 80, MarkBit.hover) ? 1 : 0);
      }
    };
  };

  /** Loci: a residue's seq loci and the domain's items, converted and projected. */
  const loci = () => {
    const t = c14;
    const [gu, pu] = [t.unit('gene/unit'), t.unit('protein/unit')];
    const list = Array.from({ length: 120 }, (_, i): Loci => {
      const unit = (i % 2 ? pu : gu) as UnitId;
      return {
        kind: 'seq-loci',
        set: t.set,
        elements: [{ unit, indices: OrderedSet.ofSingleton((i % 2 ? i % 20 : i) as SeqIndex) }],
      };
    });
    const domain = t.ctx.refs.toLoci('domain') as Loci;
    const marks = new MarkSet();
    return (i: number) => {
      t.ext.count = 0;
      t.p.fromLoci(list[i % 120], t.ext);
      if (i % 16 === 0) t.p.fromLoci(domain, t.ext);
      t.p.close(t.ext);
      t.p.project(t.ext, 0, marks);
      sink += marks.count;
    };
  };

  /** Coordinates through a registered translator: a chain's residue onto its unit, then each block. */
  const coordinates = () => {
    const t = c7;
    const tr = new Translator();
    const d: SystemDescriptor = {
      key: '3d/A',
      kind: 'external',
      start: 0,
      stop: 200,
      refs: {},
      members: null,
      structure: null,
      numberings: [],
      owner: 'molstar',
    };
    const a = { kind: K.unit, id: t.unit('1brs/A'), row: -1 };
    tr.register(d, (s) => [TranslateMaps.bySegments(s, a, Int32Array.of(5, 0, 105))]);
    t.p.systems = tr;
    const c = { systems: ['3d/A'], runs: new Int32Array(3) };
    const marks = new MarkSet();
    return (i: number) => {
      c.runs[1] = 5 + (i % 100);
      c.runs[2] = 6 + (i % 100);
      t.ext.count = 0;
      t.p.fromLoci(c, t.ext);
      for (let s = 0; s < t.view.sections.count; s++) t.p.project(t.ext, s, marks);
      sink += marks.first(0);
    };
  };

  /** Row rectangles from a row's residues, and a ligand's atoms at its one position. */
  const rowRects = () => {
    const row = r.system('a/0');
    const rects = new RowMarks();
    const one = entries([row.kind, row.id, 0, 0, 1]);
    const heme = c12.unit('dimer.heme');
    const marks = new MarkSet();
    return (i: number) => {
      one.start[0] = i % 4;
      one.stop[0] = (i % 4) + 1;
      rects.clear(MarkBit.hover);
      r.p.rows(one, 0, MarkBit.hover, rects);
      c12.ext.count = 1;
      c12.ext.kind[0] = K.unit;
      c12.ext.id[0] = heme;
      c12.ext.row[0] = -1;
      c12.ext.start[0] = i % 40;
      c12.ext.stop[0] = (i % 40) + 2;
      c12.p.project(c12.ext, 0, marks);
      sink += rects.count + marks.count;
    };
  };

  return { hover, loci, coordinates, rowRects, sink: () => sink };
}
