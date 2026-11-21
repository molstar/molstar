/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { NameHit, OpenedSource } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { EmptyKeys } from '@molstar/seqx-core/query';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { blockPositionAt, blockSpan, placeOfUnit, type ResolvedSection } from '@molstar/seqx-core/view';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { acrossBlocks } from '../behaviors/select.js';
import { viewerItem } from '../config.js';
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

  it('a block of a plain sequence is named a Sequence on its card, never a Chain', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    const b = createSeqSpec({ title: 'plain' });
    b.data.sequence('s1', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
    b.data.sequence('s2', { alphabet: 'dna', letters: 'TTGA'.repeat(50) });
    b.views.view('v', (v) => {
      v.section('plain', {}, (s) => {
        s.block('a', { of: '/data/s1' } as never);
        s.block('b', { of: '/data/s2' } as never);
        s.track('ruler').layer('ruler', { type: 'ruler' } as never);
      });
    });
    await viewer.load(b.build());
    viewer.resize(900);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    const p = viewer.frame!.layout.sections[0].placement as import('@molstar/seqx-core/render').BlockPlacement;
    const a = p.strips.axes[1];
    viewer.input.move(a.x0 + a.width / 2, p.strips.bandY[p.strips.bands] + p.strips.bandH[p.strips.bands] / 2);
    const { Describer } = await import('../describe.js');
    const card = new Describer({
      spec: viewer.spec,
      store: viewer.store,
      addresses: viewer.addresses,
      view: viewer.view,
      layout: () => viewer.frame?.layout,
    }).describe(viewer.hover.current)!;
    expect(card.rows.find((r) => r.value === 'b' || r.value === 's2')?.label).toBe('Sequence');
    expect(card.rows.some((r) => r.label === 'Chain')).toBe(false);
  });

  it('the keyboard cursor stays on the copy it is on: two copies of chr1, then chr2', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    const b = createSeqSpec({ title: 'copies' });
    const g = b.data.assembly(
      'g',
      assembly(
        'Small',
        Lengths.map((length, i) => ({ name: `chr${i + 1}`, length })),
      ),
    );
    const seq = (n: number) => `/data/g/(sequence :name chr${n})`;
    b.views.view('v', (v) => {
      v.section('copies', { copies: 'each' } as never, (s) => {
        s.block('a', { of: seq(1), range: [0, 100] } as never);
        s.block('b', { of: seq(1), range: [0, 100], label: 'chr1 again' } as never);
        s.block('c', { of: seq(2), range: [0, 50] } as never);
        s.track('ruler').layer('ruler', { type: 'ruler' } as never);
      });
    });
    void g;
    await viewer.load(b.build());
    viewer.resize(900);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    const Copies = '/views/v/copies';
    expect(sectionOf(viewer, Copies).blocks!.length).toBe(3);
    viewer.selection.setCursor({ section: Copies, system: name(1), position: 99 });
    // onto the second copy: named as another block
    expect(viewer.runCommand('cursor.next', {})).toBe('Moved, chr1, 1');
    expect(viewer.selection.cursor.value).toMatchObject({ system: name(1), position: 0 });
    // along the second copy to its end, then onto chr2 (never back into the first copy)
    for (let k = 0; k < 99; k++) viewer.runCommand('cursor.next', {});
    expect(viewer.selection.cursor.value).toMatchObject({ system: name(1), position: 99 });
    viewer.runCommand('cursor.next', {});
    expect(viewer.selection.cursor.value).toMatchObject({ system: name(2), position: 0 });
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

  it('switching blocks keeps the selection: its positions, what it shows, its address (property: next, previous, reverse, include)', async () => {
    const v = await make();
    const rng = new PCG(9);
    // a drag across chr1's end into chr2 on the genome: exactly the positions shown between
    const s0 = sectionOf(v, Genome);
    const end1 = blockSpan(s0.blocks![0]);
    const loc = acrossBlocks(s0, 0, end1 - 30, (s0.blocks![1].offset ?? 0) + 40)!;
    expect(setOf(loc)).toEqual(
      new Map([
        [name(1), new Set(Array.from({ length: 30 }, (_, i) => end1 - 30 + i))],
        [name(2), new Set(Array.from({ length: 41 }, (_, i) => i))],
      ]),
    );
    v.selection.set({ elements: EmptyKeys, location: loc, shown: loc });
    const before = v.addresses.selectionAddress()!;
    expect(before).toBeDefined();
    const selected = setOf(loc);
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
      const now = v.selection.value.value;
      // its positions and what it shows: the same, whatever blocks show them now
      expect(setOf(now.location), `${n}: ${op}`).toEqual(selected);
      expect(setOf(now.shown), `${n}: ${op}`).toEqual(selected);
      // its address the same, and read back the same positions
      expect(v.addresses.selectionAddress(), `${n}: ${op}`).toBe(before);
      expect(setOf(v.addresses.locate(before).location), `${n}: ${op}`).toEqual(selected);
      // drawn: the cells of the blocks now showing a selected system are selected exactly at the selected positions
      const s = sectionOf(v, Genome);
      for (const b of s.blocks!) {
        const want = selected.get(b.axis.system);
        if (!want) continue;
        const drawn = shownBetween({ ...s, blocks: [b] }, -Infinity, Infinity).get(b.axis.system)!;
        const hit = new Set([...drawn].filter((p) => want.has(p)));
        expect(hit, `${n}: ${op} ${b.label}`).toEqual(new Set([...want].filter((p) => drawn.has(p))));
        expect(hit.size, `${n}: ${op} ${b.label}`).toBe(want.size);
      }
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
      // every other step the window's own block reversed instead (its ends then map mirrored: the window their hull)
      const reverse = n % 2 === 1;
      if (reverse) await v.execute('block.reverse', { address: b.address });
      else await v.execute('section.include', { address: Genome, include: ['chr2', 'chr3', 'chr4'] });
      await v.settle();
      frames!.flush();
      expect(v.view.back()).toBe(true);
      frames!.flush();
      const now = v.view.get(Genome)!;
      const s2 = sectionOf(v, Genome);
      const back = [placeOfUnit(s2, now.start, 'start'), placeOfUnit(s2, now.end, 'end')];
      const sorted = (p: typeof back) => p.map((x) => x!.position).sort((x, y) => x - y);
      expect(back[0]?.system, `${n}`).toBe(ends[0]?.system);
      expect(back[1]?.system, `${n}`).toBe(ends[1]?.system);
      const [a0, a1] = sorted(ends);
      const [b0, b1] = sorted(back);
      expect(b0, `${n} ${reverse ? 'reverse' : 'include'}`).toBeCloseTo(a0, 0);
      expect(b1, `${n} ${reverse ? 'reverse' : 'include'}`).toBeCloseTo(a1, 0);
      // the block turned back for the next step
      if (reverse) {
        await v.execute('block.reverse', {
          address: sectionOf(v, Genome).blocks!.find((x) => x.axis.dir === -1)!.address,
        });
        await v.settle();
        frames!.flush();
      }
    }
  });

  it('the history signature is a few bytes a step whatever the blocks (3,000 contigs)', async () => {
    const { blocksSig } = await import('../models/view-history.js');
    const blocks = Array.from({ length: 3000 }, (_, i) => ({
      axis: { system: `/data/g/(sequence :name ctg${i})`, dir: 1 },
      offset: i * 1000,
    })) as never[];
    const a = blocksSig(blocks);
    expect(a.length).toBeLessThan(32);
    // one block's direction, or two swapped: another signature
    const turned = blocks.map((b: { axis: object }, i) => (i === 7 ? { ...b, axis: { ...b.axis, dir: -1 } } : b));
    expect(blocksSig(turned as never[])).not.toBe(a);
    const swapped = [...blocks];
    [swapped[3], swapped[4]] = [swapped[4], swapped[3]];
    expect(blocksSig(swapped)).not.toBe(a);
    expect(blocksSig([...blocks])).toBe(a);
  });

  it('acrossBlocks stop keeps the chromosome by its system: a block list changed under it never moves the window to another', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('navigate.acrossBlocks'), 'stop']],
    });
    viewers.push(viewer);
    await viewer.load(spec());
    viewer.resize(900);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    const s = sectionOf(viewer, Genome);
    const chr3 = s.blocks![2];
    viewer.view.setRange(Genome, chr3.offset! + 20, chr3.offset! + 60);
    frames.flush();
    // chr2 left out: chr3 is now the second block (the third index is chr4's)
    await viewer.execute('section.include', { address: Genome, include: ['chr1', 'chr3', 'chr4'] });
    await viewer.settle();
    frames.flush();
    const s2 = sectionOf(viewer, Genome);
    const at = viewer.view.get(Genome)!;
    // (no jump since: the stop kept the window on chr3 over the new block list)
    expect(placeOfUnit(s2, (at.start + at.end) / 2)?.system).toBe(name(3));
    // a pan past chr3's end stops there: the window stays on chr3, never chr4
    viewer.view.pan(Genome, 500);
    frames.flush();
    const now = viewer.view.get(Genome)!;
    const placeA = placeOfUnit(s2, now.start, 'start');
    const placeB = placeOfUnit(s2, now.end, 'end');
    expect([placeA?.system, placeB?.system]).toEqual([name(3), name(3)]);
    // a zoom out keeps it on chr3 whole
    viewer.view.zoom(Genome, 1 / 1000);
    frames.flush();
    const z = viewer.view.get(Genome)!;
    expect([z.start, z.end]).toEqual([s2.blocks![1].offset, (s2.blocks![1].offset ?? 0) + blockSpan(s2.blocks![1])]);
    void at;
  });

  it('block.reverse flips the override already written, however spelled (its id the name, no `of`): never a second', async () => {
    // an override written by its id alone (the name), as a writer would
    const written = spec() as { children?: { kind: string; id?: string; children?: unknown[] }[] };
    const views = written.children!.find((c) => c.kind === 'views') as {
      children: { children: { id?: string; children?: unknown[] }[] }[];
    };
    const genomeNode = views.children[0].children.find((c) => c.id === 'genome')!;
    genomeNode.children = [
      ...(genomeNode.children ?? []),
      { kind: 'block', id: 'chr2', params: { direction: 'reversed' } },
    ];
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(written as never);
    v.resize(900);
    await v.settle();
    frames.flush();
    const blocks = () => sectionOf(v, Genome).blocks!;
    const overrides = () =>
      v.spec.root.value
        .children!.find((c) => c.kind === 'views')!
        .children![0].children!.find((c) => c.id === 'genome')!
        .children!.filter((c) => c.kind === 'block');
    if (blocks()[1].axis.dir !== -1) return expect.unreachable('the id override applies');
    expect(overrides().length).toBe(1);
    await v.execute('block.reverse', { address: blocks()[1].address });
    await v.settle();
    expect(blocks()[1].axis.dir).toBe(1);
    expect(overrides().length).toBe(1);
    expect(overrides()[0].params?.direction).toBeUndefined();
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
    // over the ruler at the fit the rail is the target (its gestures), carrying the run's count and state
    const h = viewer.hit(x, (p.tracksTop + p.tracksBottom) / 2)!;
    expect(h).toMatchObject({ target: 'window', part: 'outside', count: 400 });
    expect(viewer.view.resolved.value!.sections[0].blocks![h.block!].label).toMatch(/^ctg\d+$/);
    viewer.input.move(x, (p.tracksTop + p.tracksBottom) / 2);
    expect(viewer.hover.current).toMatchObject({ target: 'window', count: 400, section: '/views/v/genome' });
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

  it('the locus box: an aborted name search settles as cancelled; a hit on a sequence no section shows is said, not focused', async () => {
    // the genome alone (no other section shows chr4)
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    const b = createSeqSpec({ title: 'genome alone' });
    const g = b.data.assembly(
      'g',
      assembly(
        'Small',
        Lengths.map((length, i) => ({ name: `chr${i + 1}`, length })),
      ),
    );
    b.views.view('v', (w) => genome(w, 'genome', g));
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    let release: () => void = () => {};
    const genes: OpenedSource = {
      coverage: {
        systems: [{ name: 'chr4', length: 250 }],
        shapes: [{ system: 'chr4', shape: 'ranges' }],
        levels: [1],
        names: ['Name'],
      },
      read: () => Promise.reject(new Error('not read here')),
      dispose() {},
      findNames: (text, _o, signal) =>
        text === 'SLOW'
          ? new Promise((resolve, reject) => {
              release = () => resolve([]);
              signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')));
            })
          : Promise.resolve(
              text === 'TP53'
                ? [{ system: 'chr4', start: 100, end: 140, id: 'g1', label: 'TP53', field: 'Name' } satisfies NameHit]
                : [],
            ),
    };
    expect(v.store.store.attach('/data/genes', genes)).toBe(true);
    const ac = new AbortController();
    const slow = goToLocus(v, 'SLOW', { signal: ac.signal });
    ac.abort();
    await expect(slow).resolves.toEqual({ kind: 'cancelled' });
    release();
    // chr4 left out of the only section showing it: the hit is said, the view stays
    await v.execute('section.include', { address: Genome, include: ['chr1', 'chr2'] });
    await v.settle();
    frames!.flush();
    const before = v.view.get(Genome)!;
    const r = await goToLocus(v, 'TP53');
    expect(r).toMatchObject({ kind: 'hidden', message: 'TP53 is on chr4, which no section shows' });
    frames!.flush();
    expect(v.view.get(Genome)).toEqual(before);
  });
});
