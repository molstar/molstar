/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { NameHit, OpenedSource } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { EmptyKeys, rangesLocation } from '@molstar/seqx-core/query';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { blockPositionAt, blockSpan, placeOfUnit, type ResolvedSection } from '@molstar/seqx-core/view';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { acrossBlocks } from '../behaviors/select.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { BuiltInHeaderEntries, type HeaderFacts, HeaderEntryNames } from '../header-entries.js';
import { goToLocus } from '../locus.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// Genome coordinates in the viewer (_spec/README.md "Behaviors", A section of blocks; "Commands", A block's commands, A
// locus; "ViewModel", History): the cursor crossing blocks (EC-154, amended), a drag across blocks, the selection kept
// over block switches, the view's history by `{ system, position }`, a block's menu, the locus box's fallback to names.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Genome = '/views/v/genome';
const Mixed = '/views/v/mixed';
const Lengths = [400, 300, 200, 250];

/** A small assembly (four sequences, letters inline) as a genome section, and a section of mixed blocks over it. */
function spec() {
  const b = createSeqSpec({ title: 'small genome' });
  const g = b.data.assembly(
    'g',
    assembly(
      'Small',
      Lengths.map((length, i) => ({ name: `chr${i + 1}`, length })),
    ),
  );
  const seq = (n: number) => `/data/g/(sequence :name chr${n})`;
  b.views.view('v', (v) => {
    genome(v, 'genome', g);
    // one reversed, one scaled: display order is not positions
    v.section('mixed', {}, (s) => {
      s.block('m1', { of: seq(1), range: [0, 200] } as never);
      s.block('m2', { of: seq(2), direction: 'reversed' } as never);
      s.block('m3', { of: seq(3), scale: 2 } as never);
      s.block('m4', { of: seq(4), range: [50, 250] } as never);
      s.track('ruler').layer('ruler', { type: 'ruler' } as never);
    });
  });
  return b.build();
}

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(900);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

const sectionOf = (v: SeqViewer, address: string) =>
  v.view.resolved.value!.sections.find((s) => s.address === address)!;
const name = (n: number) => `/data/g/(sequence :name chr${n})`;

/** Every position a set of cells covers, by system: the oracle a drag's selection is checked against. */
function shownBetween(s: ResolvedSection, lo: number, hi: number): Map<string, Set<number>> {
  const out = new Map<string, Set<number>>();
  for (const b of s.blocks!) {
    const o = b.offset ?? 0;
    const k = b.scale ?? 1;
    for (let u = o; u < o + blockSpan(b); u += k) {
      if (u < lo || u >= hi) continue;
      const p = Math.floor(blockPositionAt(b, u + k / 2, 'start')!);
      let set = out.get(b.axis.system);
      if (!set) out.set(b.axis.system, (set = new Set()));
      set.add(p);
    }
  }
  return out;
}

/** The last block starting at or before a unit. */
function lastAtOrBefore(s: ResolvedSection, u: number): number {
  let k = 0;
  while (k + 1 < s.blocks!.length && (s.blocks![k + 1].offset ?? 0) <= u) k++;
  return k;
}

function setOf(loc: { parts: readonly { system: string; ranges: ArrayLike<number> }[] }): Map<string, Set<number>> {
  const out = new Map<string, Set<number>>();
  for (const part of loc.parts) {
    const set = new Set<number>();
    for (let k = 0; k + 1 < part.ranges.length; k += 2)
      for (let p = part.ranges[k]; p < part.ranges[k + 1]; p++) set.add(p);
    out.set(part.system, set);
  }
  return out;
}

describe('a genome in the viewer', () => {
  it('EC-154 (amended): the keyboard cursor crosses blocks in display order, a reversed block’s positions down', async () => {
    const v = await make();
    v.selection.setCursor({ section: Mixed, system: name(1), position: 199 });
    expect(v.runCommand('cursor.next', {})).toBe('Moved, chr2, 300');
    // into the reversed block: its first in display order is its last position
    expect(v.selection.cursor.value).toMatchObject({ system: name(2), position: 299 });
    v.runCommand('cursor.next', {});
    expect(v.selection.cursor.value).toMatchObject({ system: name(2), position: 298 });
    v.runCommand('cursor.prev', {});
    v.runCommand('cursor.prev', {});
    expect(v.selection.cursor.value).toMatchObject({ system: name(1), position: 199 });
    // at the section's last position it stays
    v.selection.setCursor({ section: Mixed, system: name(4), position: 249 });
    v.runCommand('cursor.next', {});
    expect(v.selection.cursor.value).toMatchObject({ system: name(4), position: 249 });
  });

  it('a drag across blocks selects exactly the shown positions between, a part a system (property)', async () => {
    const v = await make();
    const s = sectionOf(v, Mixed);
    const total = s.blocks!.reduce((t, b) => Math.max(t, (b.offset ?? 0) + blockSpan(b)), 0);
    const rng = new PCG(5);
    let checked = 0;
    for (let n = 0; n < 200; n++) {
      const a = Math.floor(((rng.int() % 10_000) / 10_000) * total);
      const c = Math.floor(((rng.int() % 10_000) / 10_000) * total);
      const from = lastAtOrBefore(s, a);
      const loc = acrossBlocks(s, from, a, c);
      const to = lastAtOrBefore(s, c);
      if (to === from) {
        expect(loc).toBeUndefined();
        continue;
      }
      // the two ends' whole cells
      const cell = (u: number, k: number) => {
        const b = s.blocks![k];
        const scale = b.scale ?? 1;
        return (b.offset ?? 0) + Math.floor((u - (b.offset ?? 0)) / scale) * scale;
      };
      const ca = cell(a, from);
      const cc = cell(c, to);
      const lo = Math.min(ca, cc);
      const hi = Math.max(ca, cc) + (s.blocks![to > from ? to : from].scale ?? 1);
      expect(setOf(loc!), `${a} → ${c}`).toEqual(shownBetween(s, lo, hi));
      checked++;
    }
    expect(checked).toBeGreaterThan(100);
  });

  it('switching blocks keeps the selection’s addresses (property: next, previous, reverse, include)', async () => {
    const v = await make();
    const rng = new PCG(9);
    v.selection.set({
      elements: EmptyKeys,
      location: rangesLocation(name(2), [10, 40]),
      shown: rangesLocation(name(2), [10, 40]),
    });
    const before = v.addresses.selectionAddress();
    expect(before).toBeDefined();
    const ops = ['block.next', 'block.previous', 'block.reverse', 'section.include'] as const;
    for (let n = 0; n < 24; n++) {
      const op = ops[rng.int() % ops.length];
      const blocks = sectionOf(v, Genome).blocks!;
      const k = rng.int() % blocks.length;
      if (op === 'block.reverse') await v.execute(op, { address: blocks[k].address });
      else if (op === 'section.include')
        await v.execute(op, { address: Genome, include: rng.int() % 3 ? [`chr${(rng.int() % 4) + 1}`] : null });
      else await v.execute(op, { section: Genome });
      await v.settle();
      frames!.flush();
      expect(v.addresses.selectionAddress(), `${n}: ${op}`).toBe(before);
    }
  });

  it('the view’s history round-trips by { system, position }: a step back over a changed block list lands on the same places (property)', async () => {
    const v = await make();
    const rng = new PCG(13);
    for (let n = 0; n < 12; n++) {
      // a window inside one of the later chromosomes, then chr1 left out (every offset after it moves)
      await v.execute('section.include', { address: Genome, include: null });
      await v.settle();
      frames!.flush();
      const s = sectionOf(v, Genome);
      const b = s.blocks![1 + (rng.int() % 3)];
      const o = b.offset ?? 0;
      const span = blockSpan(b);
      const start = o + Math.floor(((rng.int() % 1000) / 1000) * (span / 2));
      const end = start + 10 + (rng.int() % Math.floor(span / 3));
      v.view.jump(() => v.view.setRange(Genome, start, end), true, 'focus');
      frames!.flush();
      const at = v.view.get(Genome)!;
      const ends = [placeOfUnit(s, at.start, 'start'), placeOfUnit(s, at.end, 'end')];
      v.view.jump(() => v.view.fit(Genome), true, 'fit');
      await v.execute('section.include', { address: Genome, include: ['chr2', 'chr3', 'chr4'] });
      await v.settle();
      frames!.flush();
      expect(v.view.back()).toBe(true);
      frames!.flush();
      const now = v.view.get(Genome)!;
      const s2 = sectionOf(v, Genome);
      const back = [placeOfUnit(s2, now.start, 'start'), placeOfUnit(s2, now.end, 'end')];
      expect(back[0]?.system).toBe(ends[0]?.system);
      expect(back[0]!.position).toBeCloseTo(ends[0]!.position, 3);
      expect(back[1]!.position).toBeCloseTo(ends[1]!.position, 3);
    }
  });

  it('a block’s menu: Reverse, Show only, Show all; a derived block reversed by an override, one undo step', async () => {
    const v = await make();
    const entry = (n: string) => BuiltInHeaderEntries.find((e) => e.entry.name === n)!.entry;
    const facts = {
      on: 'section',
      address: Genome,
      label: 'genome',
      section: Genome,
      system: Genome,
      linear: true,
      tracks: [],
      layers: [],
      hiddenTracks: [],
      block: { address: `${Genome}/(block :name chr2)`, label: 'chr2', name: 'chr2', reversed: false },
      assemblies: true,
    } satisfies HeaderFacts;
    const label = (n: string, f: HeaderFacts) => {
      const l = entry(n).label;
      return typeof l === 'function' ? l(f) : l;
    };
    expect(label(HeaderEntryNames.reverseBlock, facts)).toBe('Reverse chr2');
    expect(entry(HeaderEntryNames.onlyBlock).args!(facts)).toEqual({ address: Genome, include: ['chr2'] });
    expect(label(HeaderEntryNames.onlyBlock, facts)).toBe('Show only chr2');
    expect(entry(HeaderEntryNames.allBlocks).when!(facts)).toBe(false);
    expect(entry(HeaderEntryNames.allBlocks).when!({ ...facts, included: true })).toBe(true);
    const block = sectionOf(v, Genome).blocks![1];
    expect(block.address).toBe(`${Genome}/(block :name chr2)`);
    await v.execute('block.reverse', { address: block.address });
    await v.settle();
    frames!.flush();
    const node = v.spec.root.value.children!.find((c) => c.kind === 'views')!.children![0].children![0];
    const override = node.children!.find((c) => c.kind === 'block')!;
    expect(override.params).toEqual({ of: name(2), direction: 'reversed' });
    expect(sectionOf(v, Genome).blocks![1].axis.dir).toBe(-1);
    // again: turned back (the override's direction unset)
    await v.execute('block.reverse', { address: sectionOf(v, Genome).blocks![1].address });
    await v.settle();
    expect(sectionOf(v, Genome).blocks![1].axis.dir).toBe(1);
  });

  it('a run of small blocks: its hover names the block under the point and the run’s count, not read at this zoom', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    const b = createSeqSpec({ title: 'contigs' });
    const big = [{ name: 'chr1', length: 1_000_000 }];
    const small = Array.from({ length: 400 }, (_, i) => ({ name: `ctg${i}`, length: 100 }));
    const g = b.data.assembly('g', assembly('Draft', [...big, ...small], { order: 'source' } as never));
    b.views.view('v', (v) => genome(v, 'genome', g));
    await viewer.load(b.build());
    viewer.resize(900);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    const p = viewer.frame!.layout.sections[0].placement as import('@molstar/seqx-core/render').BlockPlacement;
    expect(p.runN).toBe(1);
    const x = p.runX[0] + p.runW[0] / 2;
    const h = viewer.hit(x, (p.tracksTop + p.tracksBottom) / 2)!;
    expect(h).toMatchObject({ target: 'run', count: 400 });
    expect(viewer.view.resolved.value!.sections[0].blocks![h.block!].label).toMatch(/^ctg\d+$/);
    viewer.input.move(x, (p.tracksTop + p.tracksBottom) / 2);
    expect(viewer.hover.current).toMatchObject({ target: 'run', count: 400, section: '/views/v/genome' });
    const { Describer } = await import('../describe.js');
    const card = new Describer({
      spec: viewer.spec,
      store: viewer.store,
      addresses: viewer.addresses,
      view: viewer.view,
      layout: () => viewer.frame?.layout,
    }).describe(viewer.hover.current)!;
    expect(card.rows).toContainEqual({ label: 'Run', value: 'one of 400' });
    // a contig a device pixel wide is read (its strip); one narrower is not, and says so
    if (p.strips.atX(x) < 0) expect(card.rows).toContainEqual({ label: 'Data', value: 'not read at this zoom' });
    expect(card.text).toMatch(/^ctg\d+, one of 400/);
    let unread = -1;
    for (let px = Math.ceil(p.runX[0]); px < p.runX[0] + p.runW[0]; px += 0.25)
      if (p.strips.atX(px) < 0) {
        unread = px;
        break;
      }
    if (unread >= 0)
      expect(viewer.hit(unread, (p.tracksTop + p.tracksBottom) / 2)).toMatchObject({ state: 'not read at this zoom' });
  });

  it('the locus box: a locus, an alias by the prefix rule, and a gene name a source finds (EC-228)', async () => {
    const v = await make();
    expect(await goToLocus(v, 'chr2:11-20')).toEqual({ kind: 'locus' });
    frames!.flush();
    const s = sectionOf(v, Genome);
    const at = v.view.get(Genome)!;
    expect(placeOfUnit(s, (at.start + at.end) / 2)?.system).toBe(name(2));
    expect((await goToLocus(v, '3')).kind).toBe('locus');
    // a source that answers names (a GFF3 read whole, synthetic): its hit on chr4, named as the source knows it
    const genes: OpenedSource = {
      coverage: {
        systems: [{ name: 'chr4', length: 250 }],
        shapes: [{ system: 'chr4', shape: 'ranges' }],
        levels: [1],
        names: ['Name'],
      },
      read: () => Promise.reject(new Error('not read here')),
      dispose() {},
      findNames: (text) =>
        Promise.resolve(
          text === 'TP53'
            ? [{ system: 'chr4', start: 100, end: 140, id: 'g1', label: 'TP53', field: 'Name' } satisfies NameHit]
            : [],
        ),
    };
    expect(v.store.store.attach('/data/genes', genes)).toBe(true);
    const r = await goToLocus(v, 'TP53');
    expect(r.kind).toBe('name');
    frames!.flush();
    const now = v.view.get(Genome)!;
    const c = placeOfUnit(s, (now.start + now.end) / 2)!;
    expect(c.system).toBe(name(4));
    expect(c.position).toBeGreaterThanOrEqual(100);
    expect(c.position).toBeLessThan(140);
    expect(await goToLocus(v, 'nothing-here')).toEqual({ kind: 'none', message: 'Nothing is named nothing-here' });
  });
});
