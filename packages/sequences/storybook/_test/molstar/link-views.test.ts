/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `linkViews` and the replay rule (seqspec:interaction#replays-and-carries; VS-07): linking sends nothing the 3D view
// did not ask for (a replay clears neither its hover nor its selection), and a carry of an unchanged Seq* selection to a
// new set (a label or a data patch) neither selects again nor moves the 3D camera. A selection made in Seq* is sent and
// focused, as before. A selection made in 3D that Seq* holds only in part (an atom: its residue; an atom of a unit no
// Seq* system covers) is not sent back by a carry either (LB-04): Mol*'s own selection stays.

import type { SeqSet } from '@molstar/seq-core/model';
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from '@molstar/seq-viewer';
import { OrderedSet } from '@molstar/core/data/int';
import { StructureElement, Unit } from '@molstar/model/model/structure';
import { expect, it } from 'vitest';
import { fixtureMeasurer } from '@molstar/seq-viewer/testing';
import {
  MolstarSvsEntry,
  readStructure,
  StructureLink,
  structureOf,
  linkViews,
} from '@molstar/prism-viewer/sequence/molstar';
import { bytesOf, demoLoaders, P53 } from './fixtures.js';
import { P53_STRUCTURES, p53DemoDoc, p53ShownPatch } from '../../src/_stories/support/seqstar-p53-demo.js';

type DocNode = { ref?: string; params?: Record<string, unknown>; children?: DocNode[] };
const find = (n: DocNode | undefined, ref: string): DocNode | undefined =>
  n?.ref === ref ? n : (n?.children ?? []).map((c) => find(c, ref)).find(Boolean);

it('a replay or a carry sends nothing to 3D; a Seq* selection is sent and focused', async () => {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
    { scheduler, loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(p53DemoDoc(P53.sequence));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1200, 2400, 1);
  await done;
  await ctx.graph.update();
  scheduler.advance(50);
  const set = ctx.graph.assembled.value?.set as SeqSet;
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  const link = StructureLink.attach(ctx.systems, '3d', structure, set, '1tup');
  const calls: string[] = [];
  const unlink = linkViews(
    ctx,
    {
      highlight: (l) => calls.push(`highlight(${l === null ? 'null' : 'loci'})`),
      select: (l) => calls.push(`select(${l === null ? 'null' : 'loci'})`),
      focus: () => calls.push('focus'),
      onHover: () => () => {},
      onSelect: () => () => {},
    },
    link,
  );
  // Linked with nothing in Seq*: the 3D view's own hover and selection stay.
  expect(calls).toEqual([]);
  ctx.selection.set(ctx.refs.toLoci('p53', { start: 100, stop: 110 }), 'pointer');
  expect(calls).toEqual(['select(loci)', 'focus']);
  calls.length = 0;
  // A label patch, then a data patch: the selection is carried (origin `graph`), the same residues.
  await ctx.applyPatch(p53ShownPatch(P53_STRUCTURES[0]));
  scheduler.advance(50);
  const hp = find(ctx.document.shown.value?.root as DocNode, 'hydropathy');
  const value = structuredClone(hp?.params?.value) as { value?: number[] } | number[];
  if (Array.isArray(value)) value[0] += 1;
  else if (Array.isArray(value.value)) value.value[0] += 1;
  const r = await ctx.applyPatch({
    format: 'seq_view_spec_patch',
    ops: [{ op: 'set', ref: 'hydropathy', params: { ...hp?.params, value } }],
  });
  expect(r.ok).toBe(true);
  scheduler.advance(50);
  expect(ctx.selection.current.value.origin).toBe('graph');
  expect(calls).toEqual([]);
  unlink();
  link.detach();
  ctx.dispose();
}, 60_000);

it("a 3D selection Seq* holds only in part (an atom, another unit's atom) survives a data edit", async () => {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
    { scheduler, loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(p53DemoDoc(P53.sequence));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1200, 2400, 1);
  await done;
  await ctx.graph.update();
  scheduler.advance(50);
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  const link = StructureLink.attach(ctx.systems, '3d', structure, ctx.graph.assembled.value?.set as SeqSet, '1tup');
  const calls: string[] = [];
  const hooks: { select: ((l: StructureElement.Loci | null) => void) | null } = { select: null };
  const unlink = linkViews(
    ctx,
    {
      highlight: () => {},
      select: (l) => calls.push(`select(${l === null ? 'null' : StructureElement.Loci.size(l)})`),
      focus: () => calls.push('focus'),
      onHover: () => () => {},
      onSelect: (f) => {
        hooks.select = f;
        return () => {};
      },
    },
    link,
  );
  let edits = 0;
  const dataPatch = async () => {
    const hp = find(ctx.document.shown.value?.root as DocNode, 'hydropathy');
    const value = structuredClone(hp?.params?.value) as { value?: number[] } | number[];
    if (Array.isArray(value)) value[0] += ++edits;
    else if (Array.isArray(value.value)) value.value[0] += ++edits;
    const r = await ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'hydropathy', params: { ...hp?.params, value } }],
    });
    scheduler.advance(50);
    return r.ok;
  };
  const u = structure.units.find((x) => Unit.isAtomic(x)) as Unit;
  const last = structure.units[structure.units.length - 1];
  expect(last).not.toBe(u);
  for (const loci of [
    StructureElement.Loci(structure, [{ unit: u, indices: OrderedSet.ofSingleton(0 as never) }]),
    StructureElement.Loci(structure, [
      { unit: u, indices: OrderedSet.ofSingleton(0 as never) },
      { unit: last, indices: OrderedSet.ofSingleton(0 as never) },
    ]),
  ]) {
    calls.length = 0;
    hooks.select?.(loci);
    expect(ctx.selection.current.value.origin).toBe('host:3d');
    expect(await dataPatch()).toBe(true);
    expect(ctx.selection.current.value.origin).toBe('graph');
    expect(calls).toEqual([]);
  }
  unlink();
  link.detach();
  ctx.dispose();
}, 60_000);
