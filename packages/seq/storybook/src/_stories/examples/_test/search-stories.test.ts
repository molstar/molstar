/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search stories, headless on a fake clock, composed as the stories compose them (the example view on the default
// spec; the box's runs as the box sends them, its count in its words: the box's DOM is support/_test's): p53's PXXP on
// the protein, lit on its gene through the codon map; the plasmid's site across the origin, one hit in two runs, on the
// plus strand and, reversed, on the minus.

import { P53Views, PlasmidViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { searchCount } from '../../support/search-box.js';
import { AcrossTheOrigin, AcrossTheOriginMinus } from '../plasmid/search.js';
import { openExample } from './draw-checks.js';

const raf = DomHost.requestAnimationFrame;
const viewers: SeqViewer[] = [];
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  for (const v of viewers.splice(0)) v.retire();
  vi.useRealTimers();
});

/** The clock run while a run's chunks land (their imports are real I/O, never the fake clock's) until it ends. */
async function ran(v: SeqViewer) {
  for (let i = 0; i < 4000; i++) {
    await vi.advanceTimersByTimeAsync(16);
    await new Promise((r) => setImmediate(r));
    const s = v.search!.status.value;
    if (s === 'done' || s === 'failed') return;
  }
}

async function story(view: Parameters<typeof openExample>[0], initial: string) {
  const v = await openExample(view, { fixtures: true });
  viewers.push(v);
  // what the box sends for a pattern (its matcher and strands the configuration's)
  const type = async (pattern: string) => {
    void v.execute('search.run', { pattern });
    await ran(v);
  };
  await type(initial);
  return { v, count: () => searchCount(v), type };
}

describe('p53, search', () => {
  it('PXXP on the protein: its hits tinted and lit on the gene through the codon map; Next steps', async () => {
    const s = await story(P53Views.find((v) => v.id === 'protein-on-gene')!, 'PXXP');
    const search = s.v.search!;
    expect(search.status.value).toBe('done');
    const n = search.count.value;
    // p53's proline-rich region holds several PXXP motifs
    expect(n).toBeGreaterThan(3);
    expect(s.count()).toBe(`${n} hits`);
    const { hits, systems } = search.hits.value;
    const protein = systems[hits.target[0]];
    expect(s.v.spec.model.value.systems.get(protein)?.alphabet).toBe('protein');
    expect(s.v.highlight.get('search')?.location.parts.map((p) => p.system)).toEqual([protein]);
    await s.v.execute('search.next', {});
    expect(search.current.value).toBeGreaterThanOrEqual(0);
    expect(s.count()).toMatch(new RegExp(`^\\d+ of ${n}$`));
    expect(search.describe(search.current.value)).toMatch(/ \d+ to \d+$/);
  });
});

describe('Plasmid, search across the origin', () => {
  it('twelve letters across the origin: one hit in two runs, on the plus strand; reversed, on the minus', async () => {
    const s = await story(PlasmidViews.find((v) => v.id === 'crossing')!, AcrossTheOrigin);
    const search = s.v.search!;
    expect(search.count.value).toBe(1);
    const { hits, lengths } = search.hits.value;
    // held unrolled: its end past the length
    expect(hits.end[0]).toBeGreaterThan(lengths[0]);
    expect(hits.strand[0]).toBe(1);
    expect(search.locationOf(0)?.parts[0].ranges.length).toBe(4);
    expect(search.addressOf(0)).toMatch(/\(range :from (\d+) :to (\d+)\)$/);
    const [, from, to] = /:from (\d+) :to (\d+)/.exec(search.addressOf(0)!)!;
    expect(Number(from)).toBeGreaterThan(Number(to));
    // the reverse complement: the same place, on the minus strand
    await s.type(AcrossTheOriginMinus);
    expect(search.count.value).toBe(1);
    expect(search.hits.value.hits.strand[0]).toBe(-1);
    expect(search.hits.value.hits.start[0]).toBe(hits.start[0]);
    expect(s.count()).toBe('1 hit');
  });
});
