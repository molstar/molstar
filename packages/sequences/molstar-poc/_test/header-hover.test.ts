/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track header's hover reaches 3D whole (seqspec:interaction#hover, seqspec:libraries#the-mol-adapter): in the p53
// demo, the PTMs header lights in 1TUP exactly the residues its PTMs light one by one, in one highlight; with the link's
// hover `maxItems` below the set's size (read before it resolves), the whole chains those PTMs touch.

import type { SeqSet } from '@molstar/seq-core/model';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  Hit,
  HoverKind,
  manualScheduler,
} from '@molstar/seq-viewer';
import { StructureElement } from '@molstar/model/model/structure';
import { expect, it } from 'vitest';
import { helpers } from '../../viewer/_test/layout-fixtures.js';
import { fixtureMeasurer } from '../../viewer/src/testing/index.js';
import { MolstarSvsEntry, readStructure, StructureLink, structureOf } from '../src/index.js';
import { linkViews } from '../src/viewer.js';
import { bytesOf, demoLoaders, P53 } from './fixtures.js';
import { P53_TRACKS, p53DemoDoc } from './p53-demo.js';

const SECTION = 'main.main';

/** The residues of a Mol* loci, as `<chain> <seq id>`. */
function residues(loci: StructureElement.Loci | null): string[] {
  const out = new Set<string>();
  if (loci)
    StructureElement.Loci.forEachLocation(loci, (l) => {
      const h = l.unit.model.atomicHierarchy;
      const chain = h.chains.label_asym_id.value(h.chainAtomSegments.index[l.element]);
      out.add(`${chain} ${h.residues.auth_seq_id.value(h.residueAtomSegments.index[l.element])}`);
    });
  return [...out];
}

async function demo(hoverOptions?: { maxItems: number }) {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
    { scheduler, loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(p53DemoDoc(P53.sequence));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1200, 2000, 1);
  expect((await done).ok).toBe(true);
  await ctx.graph.update();
  scheduler.advance(50);
  const set = ctx.graph.assembled.value?.set as SeqSet;
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  const link = StructureLink.attach(ctx.systems, '3d', structure, set, '1tup');
  const lit: (StructureElement.Loci | null)[] = [];
  const unlink = linkViews(
    ctx,
    {
      highlight: (l) => lit.push(l),
      select: () => {},
      focus: () => {},
      onHover: () => () => {},
      onSelect: () => () => {},
    },
    link,
    { select: false, focus: false, hover: hoverOptions },
  );
  const t = helpers(ctx);
  const hit = Hit.create();
  const hover = (x: number, y: number) => {
    lit.length = 0;
    ctx.layout.hitTest(SECTION, x, y, hit);
    ctx.hover.setFromHit(hit, 'pointer');
    return [...lit];
  };
  const close = () => {
    unlink();
    link.detach();
    ctx.dispose();
  };
  return { ctx, set, t, hit, hover, close };
}

it("a track header's hover lights in 3D exactly what its items light one by one, as one highlight", async () => {
  const { set, t, hit, hover, close } = await demo();
  const y = t.yOf(SECTION, P53_TRACKS.ptms);
  const ptms = set.annotations.find((a) => a && set.strings.get(a.key) === 'ptms');
  if (!ptms) throw new Error('no ptms');
  // Each PTM by itself, on its item: the union of what they light in 1TUP.
  const each = new Set<string>();
  for (let r = 0; r < ptms.count; r++)
    for (const x of residues(hover(t.xOf(SECTION, 0, ptms.start[r]), y).at(-1) ?? null)) each.add(x);
  expect(each.size).toBeGreaterThan(1);
  const sent = hover(10, y);
  expect(hit.kind).toBe(HoverKind.header);
  expect(sent.length).toBe(1);
  expect(residues(sent[0]).sort()).toEqual([...each].sort());
  close();
}, 60_000);

it("past the link's hover maxItems (the size read before resolving), the whole chains the set touches", async () => {
  const { ctx, t, hover, close } = await demo({ maxItems: 1 });
  const y = t.yOf(SECTION, P53_TRACKS.ptms);
  const sent = hover(10, y);
  expect(ctx.hover.size()?.items).toBeGreaterThan(1);
  const chains = new Set(residues(sent[0]).map((r) => r.split(' ')[0]));
  const atoms = sent[0] ? StructureElement.Loci.size(sent[0]) : 0;
  // Whole chains: every residue of each chain touched, far more than the PTMs' own.
  expect(chains.size).toBeGreaterThan(0);
  expect(residues(sent[0]).length).toBeGreaterThan(100);
  expect(atoms).toBeGreaterThan(1000);
  close();
}, 60_000);

it("a hover on a track's data draws its row band and reaches 3D as its element alone, never as the header's set", async () => {
  const { ctx, set, t, hit, hover, close } = await demo();
  // The band is on (the hoverRow behavior of the demo's viewer spec); it is a drawing, not a header's hover.
  expect(ctx.ui.rowHover.value?.(ctx.layout.sectionIndex(SECTION))?.band).toBe('entry');
  const y = t.yOf(SECTION, P53_TRACKS.ptms);
  const header = residues(hover(10, y)[0] ?? null);
  const ptms = set.annotations.find((a) => a && set.strings.get(a.key) === 'ptms');
  if (!ptms) throw new Error('no ptms');
  // Each PTM by itself: the item, one highlight of its own residues (none for one outside 1TUP), never the header's.
  let lit = 0;
  for (let r = 0; r < ptms.count; r++) {
    const one = hover(t.xOf(SECTION, 0, ptms.start[r]), y);
    const h = ctx.hover.live.current;
    // At most one highlight (none when the element is the one hovered already).
    expect([hit.kind, h.kind, h.header.on, one.length <= 1]).toEqual([HoverKind.item, HoverKind.item, false, true]);
    const own = residues(one[0] ?? null);
    expect(own.length).toBeLessThan(header.length);
    expect(own.every((x) => header.includes(x))).toBe(true);
    if (own.length > 0) lit++;
  }
  expect(lit).toBeGreaterThan(0);
  // A continuous track's position: its residue in each chain, one sequence id.
  const at = hover(t.xOf(SECTION, 0, 120), t.yOf(SECTION, P53_TRACKS.hydropathy));
  expect([hit.kind, ctx.hover.live.current.header.on]).toEqual([HoverKind.seq, false]);
  expect(new Set(residues(at.at(-1) ?? null).map((r) => r.split(' ')[1]))).toEqual(new Set(['121']));
  close();
}, 60_000);
