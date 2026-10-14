/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A dot plot's cell reaches 3D as both its residues (seqspec:interaction#hover, seqspec:libraries#the-mol-adapter): in
// the p53 demo with a dot plot track added over p53, a hover on the cell (R175, R248) lights both hotspots in 1TUP's
// chains, in one highlight.

import type { DocNode } from '@molstar/seq-core/format';
import type { SeqSet } from '@molstar/seq-core/model';
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, Hit, manualScheduler } from '@molstar/seq-viewer';
import { StructureElement } from '@molstar/model/model/structure';
import { expect, it } from 'vitest';
import { helpers } from '../../viewer/_test/layout-fixtures.js';
import { fixtureMeasurer } from '../../viewer/src/testing/index.js';
import { MolstarSvsEntry, readStructure, StructureLink, structureOf } from '../src/index.js';
import { linkViews } from '../src/viewer.js';
import { bytesOf, demoLoaders, P53 } from './fixtures.js';
import { p53DemoDoc } from './p53-demo.js';

const SECTION = 'main.main';

/** The demo with a dot plot of one pair over p53 (R175 against R248), its last track. */
function doc() {
  const d = p53DemoDoc(P53.sequence);
  const kids = (n: DocNode) => (n.children ?? []) as DocNode[];
  const data = kids(d.root as DocNode).find((n) => n.kind === 'data') as DocNode;
  kids(data).push({
    kind: 'annotation',
    ref: 'hotspots-pair',
    params: { geometry: 'pairs', target: 'p53', value: { start: [174], stop: [175], b_start: [247], b_stop: [248] } },
  });
  const view = kids(d.root as DocNode).find((n) => n.ref === 'main') as DocNode;
  const sec = kids(view).find((n) => n.ref === SECTION) as DocNode;
  kids(sec).push({
    kind: 'track',
    ref: 't-dots',
    params: { label: 'Dots' },
    children: [{ kind: 'layer', ref: 'l-dots', params: { type: 'dot_plot', data: 'hotspots-pair' } }],
  });
  return d;
}

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

it("a dot plot's cell lights both its residues in 3D, in one highlight", async () => {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
    { scheduler, loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc());
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1200, 2400, 1);
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
    { select: false, focus: false },
  );
  const t = helpers(ctx);
  // The cell (175, 248): 175's column, and down the plot where 248's x falls across the data area.
  const l = t.layout(SECTION);
  const v = ctx.layout.axis(SECTION);
  if (!v) throw new Error('no axis');
  const e = l.entryCount - 1;
  const top = l.scrolledTop + l.entryTop[e] - v.scrollY;
  const y = top + ((t.xOf(SECTION, 0, 247) - v.left) / (v.right - v.left)) * l.entryHeight[e];
  const hit = Hit.create();
  ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, 174), y, hit);
  expect([hit.index, hit.pairIndex]).toEqual([174, 247]);
  lit.length = 0;
  ctx.hover.setFromHit(hit, 'pointer');
  expect(lit.length).toBe(1);
  const got = residues(lit[0]);
  expect(got).toEqual(expect.arrayContaining(['C 175', 'C 248']));
  expect(new Set(got.map((r) => r.split(' ')[1]))).toEqual(new Set(['175', '248']));
  unlink();
  link.detach();
  ctx.dispose();
}, 60_000);
