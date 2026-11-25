/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { rondo, rondoRing } from '@molstar/seqx-ext-rondo';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { ribbonsOf, strongestStem } from '../views.js';

// `.` and `,` in one viewer of a ring and a section of tracks: the feature keys on the tracks (at the last feature too,
// where they say so and stop), the ring's partner walk on a stem (_spec/README.md "Bindings", the feature keys).

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const Marks = [100_000, 300_000, 500_000];

async function open() {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const chroms = Array.from({ length: 4 }, (_, i) => ({ name: `chr${i + 1}`, length: 1_000_000 }));
  const rows = {
    a_name: [] as string[],
    b_name: [] as string[],
    a_start: [] as number[],
    a_end: [] as number[],
    b_start: [] as number[],
    b_end: [] as number[],
  };
  for (let k = 0; k < 60; k++) {
    const a = k % 4;
    const b = (a + 1 + (k % 3)) % 4;
    rows.a_name.push(chroms[a].name);
    rows.b_name.push(chroms[b].name);
    rows.a_start.push(k * 10_000);
    rows.a_end.push(k * 10_000 + 500);
    rows.b_start.push(k * 9_000);
    rows.b_end.push(k * 9_000 + 500);
  }
  const b = createSeqSpec({ title: 'a ring and tracks', extensions: [rondo] });
  const g = b.data.assembly('toy', assembly('Toy', chroms));
  const links = b.data.pairs('links', { a: g, rows } as never);
  const marks = b.data.ranges('marks', {
    target: g,
    rows: { target: Marks.map(() => 'chr1'), start: Marks, end: Marks.map((m) => m + 1000), id: ['a', 'b', 'c'] },
  } as never);
  let lin = '';
  b.views.view('main', (v) => {
    rondoRing(v, { assembly: g, links, square: false });
    lin = genome(v, 'lin', g, { include: ['chr1'] }, (s) =>
      s.track('marks', { label: 'Marks' }).layer('marks', { type: 'ranges', data: marks } as never),
    ).address;
  });
  const base = DefaultSeqViewerSpec();
  viewer = createSeqViewer({ ...base, ui: undefined, extensions: [...(base.extensions ?? []), rondoViewer] });
  void viewer.load(b.build());
  viewer.resize(800);
  await vi.advanceTimersByTimeAsync(3000);
  viewer.draw();
  await vi.advanceTimersByTimeAsync(3000);
  viewer.draw();
  return { v: viewer, lin };
}

describe('the feature keys and the ring’s partner keys in one viewer', () => {
  it('on the tracks: the next feature, and at the last one no partner walk; on a stem: the partners', async () => {
    const { v, lin } = await open();
    const section = v.view.resolved.value!.sections.find((s) => s.address === lin)!;
    const system = section.blocks![0].axis.system;
    const at = (position: number) => v.selection.setCursor({ section: lin, system, position, keyboard: true });
    const here = () => v.selection.cursor.value!;
    at(0);
    expect(v.input.key({ key: '.' })).toBe(true);
    expect(here().position).toBe(Marks[0]);
    v.input.key({ key: '.' });
    v.input.key({ key: '.' });
    expect(here().position).toBe(Marks[2]);
    // the last feature: the key is the feature command's (it stops), never `rondo:partner`'s
    expect(v.input.key({ key: '.' })).toBe(true);
    expect([here().section, here().system, here().position]).toEqual([lin, system, Marks[2]]);
    expect(v.input.key({ key: ',' })).toBe(true);
    expect(here().position).toBe(Marks[1]);
    // the ring: the same key walks the stem's partners
    const l = ribbonsOf(v)!;
    const { stem } = strongestStem(l);
    v.selection.setCursor({ section: l.section, system: l.systems[l.stemSystem[stem]], position: l.stemStart[stem] });
    expect(v.input.key({ key: '.' })).toBe(true);
    expect(here().section).toBe(l.section);
    expect(here().system).not.toBe(l.systems[l.stemSystem[stem]]);
  });
});
