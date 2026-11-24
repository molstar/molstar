/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search in the viewer, headless on a fake clock (_spec/README.md "Search"): a run over what the view draws, its
// hits on two channels of its own, stepping revealed through the address model, a saved set one patch; each edge case
// by name.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { changeKindOf } from '../commands.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

/** A plasmid of 100: GAATTC across its origin, GGATCC at 20, CTAAGG at 56 (CCTTAG on the minus strand). */
export const PlasmidLetters = `TTC${'A'.repeat(17)}GGATCC${'A'.repeat(30)}CTAAGG${'A'.repeat(35)}GAA`;
const Gene = 'ATGGAGGAGCCC'.repeat(3);

function spec(): Node {
  const b = createSeqSpec({ title: 'search' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: PlasmidLetters, cyclic: true });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPMEEPMEEP' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [0], end: [36] }, strand: '+' });
  b.views.view('main', (v) => {
    v.section('plasmid', { of: c }).track('l').layer('l', { type: 'letters' });
    const protein = v.section('protein', { of: p });
    protein.track('l').layer('l', { type: 'letters' });
    v.section('gene', { of: g, follow: protein }).track('l').layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(options: { history?: boolean; config?: [unknown, unknown][] } = {}) {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    ...(options.history ? { history: uiHistory() } : {}),
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('navigator.form'), 'none'],
      ...((options.config ?? []) as never[]),
    ],
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const said = (v: SeqViewer) => {
  const out: string[] = [];
  v.notice.listen(({ text }) => out.push(text));
  return out;
};
const hitsOf = (v: SeqViewer) => {
  const { hits, systems } = v.search!.hits.value;
  return Array.from({ length: hits.length }, (_, i) => [
    systems[hits.target[i]],
    hits.start[i],
    hits.end[i],
    hits.strand[i],
  ]);
};

describe('a search over what the view draws', () => {
  it('runs every system a letters layer draws, whole; the count, the words, the channel', async () => {
    const v = await make();
    const words = said(v);
    const out = await v.execute('search.run', { pattern: 'GGATCC' });
    expect(out).toBe('1 hit for GGATCC on 3 sequences');
    expect(words.at(-1)).toBe(out);
    expect(hitsOf(v)).toEqual([['/data/c', 20, 26, 0]]);
    expect(v.search!.status.value).toBe('done');
    expect(v.search!.complete.value).toBe(true);
    expect(v.highlight.get('search')?.location.parts.map((p) => [p.system, [...p.ranges]])).toEqual([
      ['/data/c', [20, 26]],
    ]);
  });

  it('EC-9, EC-10: a site across the plasmid’s origin is one hit; its address wraps on the cyclic system', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'GAATTC' });
    expect(hitsOf(v)).toEqual([['/data/c', 97, 103, 0]]);
    expect(v.search!.count.value).toBe(1);
    // two ranges in the channel, one hit in the count
    expect(v.highlight.get('search')?.location.parts[0].ranges).toEqual(Int32Array.from([0, 3, 97, 100]));
    expect(v.search!.addressOf(0)).toBe('/data/c/(range :from 98 :to 3)');
    expect(v.search!.describe(0)).toBe('c 98 to 3, both strands');
  });

  it('both strands: a minus-strand site is found and said so', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'CCTTAG' });
    expect(hitsOf(v)).toEqual([['/data/c', 56, 62, -1]]);
    expect(v.search!.describe(0)).toBe('c 57 to 62, minus strand');
    await v.execute('search.run', { pattern: 'CCTTAG', strands: 'forward' });
    expect(hitsOf(v)).toEqual([]);
  });

  it('a compile error is the result’s message and issue: nothing changes', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'GGATCC' });
    const before = v.search!.hits.value;
    expect(await v.execute('search.run', { pattern: '(A)\\1', matcher: 'regex' })).toBe(
      'Backreferences are not searched',
    );
    expect(v.search!.hits.value).toBe(before);
    expect(v.search!.query.value?.pattern).toBe('GGATCC');
  });

  it('a matcher the composition does not register is refused, named, no fallback', async () => {
    const v = await make();
    expect(await v.execute('search.run', { pattern: 'A', matcher: 'app:blast' })).toBe('No search matcher app:blast');
    expect(v.search!.query.value).toBeUndefined();
  });

  it('a protein pattern’s hits light the gene through the codon map (projected, not translated)', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'EP' });
    expect(hitsOf(v)).toEqual([
      ['/data/p', 2, 4, 0],
      ['/data/p', 6, 8, 0],
      ['/data/p', 10, 12, 0],
    ]);
    frames!.flush();
    const lists = recordDrawLists(v.frame!.lists, ['highlight']);
    // the gene section (the third) has highlight rects too: the channel projected
    const gene = v.frame!.layout.sections[2];
    const rects = lists
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .filter((l) => {
        const y = Number(l.split(' ')[2]);
        return y >= gene.y - 4 && y < gene.y + gene.height;
      });
    expect(rects.length).toBeGreaterThan(0);
  });
});

describe('stepping through the hits', () => {
  it('EC-161: "k of n"; the current again clears it', async () => {
    const v = await make();
    const words = said(v);
    await v.execute('search.run', { pattern: 'EP' });
    await v.execute('search.go', { index: 1 });
    expect(v.search!.current.value).toBe(1);
    expect(words.at(-1)).toBe('Hit 2 of 3: p 7 to 8');
    expect(v.highlight.get('search-current')?.style).toBe('outline');
    expect(v.highlight.get('search-current')?.location.parts[0].ranges).toEqual(Int32Array.from([6, 8]));
    await v.execute('search.go', { index: 1 });
    expect(v.search!.current.value).toBe(-1);
    expect(v.highlight.get('search-current')).toBeUndefined();
  });

  it('EC-160: an index outside the hits is refused, never clipped', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'EP' });
    expect(await v.execute('search.go', { index: 3 })).toBe('No hit 4: there are 3');
    expect(await v.execute('search.go', { index: -1 })).toBe('No hit 0: there are 3');
    expect(v.search!.current.value).toBe(-1);
  });

  it('EC-159: next reveals, panning only (the zoom kept), and wraps, said', async () => {
    const v = await make();
    const words = said(v);
    await v.execute('search.run', { pattern: 'AAAA', strands: 'forward' });
    const section = '/views/main/plasmid';
    v.view.setRange(section, 30, 40);
    const zoom = () => {
      const vp = v.view.get(section)!;
      return vp.end - vp.start;
    };
    const width = zoom();
    // from no current hit: the first at or after the window's start
    await v.execute('search.next', {});
    const first = v.search!.current.value;
    expect(v.search!.hits.value.hits.start[first]).toBeGreaterThanOrEqual(30);
    expect(zoom()).toBeCloseTo(width, 6);
    // to the last, then round to the first, said
    const n = v.search!.hits.value.hits.length;
    await v.execute('search.go', { index: n - 1 });
    await v.execute('search.next', {});
    expect(v.search!.current.value).toBe(0);
    expect(words.at(-1)).toMatch(/wrapped to the first hit$/);
    expect(zoom()).toBeCloseTo(width, 6);
    // previous from the first wraps to the last
    await v.execute('search.previous', {});
    expect(v.search!.current.value).toBe(n - 1);
  });

  it('a step is one jump of the view’s history: Back returns', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'GAATTC' });
    const section = '/views/main/plasmid';
    v.view.setRange(section, 30, 40);
    await v.execute('search.next', {});
    const after = v.view.get(section)!.start;
    expect(after).not.toBe(30);
    await v.execute('view.back', { section });
    expect(v.view.get(section)!.start).toBeCloseTo(30, 6);
  });
});

describe('the keys', () => {
  it('EC-155: Mod+g and F3 step on, with Shift back, in the canvas’s scope; Mod+f is never taken', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'EP' });
    expect(v.input.key({ key: 'g', ctrlKey: true, metaKey: true })).toBe(true);
    expect(v.search!.current.value).toBe(0);
    expect(v.input.key({ key: 'F3' })).toBe(true);
    expect(v.search!.current.value).toBe(1);
    expect(v.input.key({ key: 'F3', shiftKey: true })).toBe(true);
    expect(v.search!.current.value).toBe(0);
    expect(v.input.key({ key: 'G', ctrlKey: true, metaKey: true, shiftKey: true })).toBe(true);
    expect(v.search!.current.value).toBe(2);
    // the browser's find stays the host's; a host's box never reaches these (the canvas's scope)
    expect(v.input.key({ key: 'f', ctrlKey: true, metaKey: true })).toBe(false);
    for (const id of ['search.next', 'search.previous'])
      expect(v.commands.find((c) => c.id === id)?.scope).toBe('canvas');
    // no key runs a search: it needs text
    expect(v.commands.find((c) => c.id === 'search.run')?.shortcut).toBeUndefined();
  });
  it('with no hits Mod+g, Mod+Shift+G and F3 are the page’s (the browser’s find next), as Back is with no step', async () => {
    const v = await make();
    expect(v.input.key({ key: 'g', ctrlKey: true, metaKey: true })).toBe(false);
    expect(v.input.key({ key: 'G', ctrlKey: true, metaKey: true, shiftKey: true })).toBe(false);
    expect(v.input.key({ key: 'F3', shiftKey: true })).toBe(false);
    await v.execute('search.run', { pattern: 'EP' });
    expect(v.input.key({ key: 'g', ctrlKey: true, metaKey: true })).toBe(true);
  });
});

describe('what is said', () => {
  it('a run’s end, a step and a wrap in words; no hits; a cancelled run says nothing', async () => {
    const v = await make();
    const words = said(v);
    expect(await v.execute('search.run', { pattern: 'WWWWWW', matcher: 'exact' })).toBe('No hits for WWWWWW');
    await v.execute('search.run', { pattern: 'EP' });
    expect(words.at(-1)).toBe('3 hits for EP on 3 sequences');
    await v.execute('search.go', { index: 2 });
    expect(words.at(-1)).toBe('Hit 3 of 3: p 11 to 12');
    await v.execute('search.next', {});
    expect(words.at(-1)).toBe('Hit 1 of 3: p 3 to 4, wrapped to the first hit');
    const n = words.length;
    const first = v.execute('search.run', { pattern: 'AAAA' }) as Promise<string>;
    await v.execute('search.clear', {});
    expect(await first).toBe('');
    expect(words.length).toBe(n);
    expect(v.search!.status.value).toBe('idle');
  });
});

describe('saving the hits', () => {
  it('EC-162: a saved set is one patch, one undo step: a ref of the hits’ addresses', async () => {
    const v = await make({ history: true });
    const stack = (v.spec.history as UiHistory).stack;
    await v.execute('search.run', { pattern: 'EP' });
    const out = await v.execute('search.save', {});
    expect(out).toBe('Saved 3 hits as Search: EP (3)');
    expect(stack.steps.value.map((x) => x.label)).toEqual(['Save the search']);
    const ref = v.addresses.refs.value.find((r) => r.address === '/refs/search-1');
    expect(ref?.label).toBe('Search: EP (3)');
    await stack.undo();
    expect(v.addresses.refs.value.find((r) => r.address === '/refs/search-1')).toBeUndefined();
  });

  it('is refused, never cut: past search.maxSaved, or incomplete', async () => {
    const v = await make({ config: [[viewerItem('search.maxSaved'), 2]] });
    await v.execute('search.run', { pattern: 'EP' });
    expect(await v.execute('search.save', {})).toBe('3 hits: a saved set holds at most 2');
    const w = await make({ config: [[viewerItem('search.maxHits'), 2]] });
    await w.execute('search.run', { pattern: 'EP' });
    expect(w.search!.count.value).toBe(3);
    expect(w.search!.hits.value.hits.length).toBe(2);
    expect(w.search!.complete.value).toBe(false);
    expect(await w.execute('search.save', {})).toBe('The result is incomplete: not saved');
  });
});

describe('the channels and the legend', () => {
  it('EC-157: removing the search channel clears the search; clearing the search clears both channels', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'EP' });
    await v.execute('search.go', { index: 0 });
    // the legend's channel removal runs `highlight.clear` on it
    await v.execute('highlight.clear', { channel: 'search' });
    expect(v.search!.query.value).toBeUndefined();
    expect(v.highlight.get('search-current')).toBeUndefined();
    await v.execute('search.run', { pattern: 'EP' });
    v.addresses.highlightAddress('/data/p/(range :from 1 :to 2)', 'other');
    await v.execute('search.clear', {});
    expect(v.highlight.names.value).toEqual(['other']);
  });

  it('the commands’ kinds: a run, a step and a clear session state (a step’s reveal a jump of its own), a save a patch', () => {
    expect(['run', 'next', 'previous', 'go', 'clear', 'save'].map((c) => changeKindOf(`search.${c}`))).toEqual([
      'session',
      'session',
      'session',
      'session',
      'session',
      'patch',
    ]);
  });
});

describe('results follow the spec', () => {
  it('new letters for a searched system mark the result stale, still drawn, until the next run', async () => {
    const v = await make();
    await v.execute('search.run', { pattern: 'EP' });
    const words = said(v);
    // the protein's letters edited: another entry for its system
    const next = JSON.parse(JSON.stringify(v.spec.shown.value)) as Node;
    const data = next.children!.find((c) => c.id === 'data')!;
    const p = data.children!.find((c) => c.id === 'p')!;
    (p.params as Record<string, unknown>).letters = 'MEEPMEEPMEEQ';
    await v.load(next);
    await v.settle();
    expect(v.search!.stale.value).toBe(true);
    expect(words).toContain('Search results may be out of date');
    expect(v.highlight.get('search')).toBeDefined();
    await v.execute('search.run', { pattern: 'EP' });
    expect(v.search!.stale.value).toBe(false);
    expect(v.search!.count.value).toBe(2);
  });
});

describe('progress', () => {
  it('is said once, after a second, never per batch', async () => {
    const v = await make({ config: [[viewerItem('search.chunkLetters'), 8]] });
    const words = said(v);
    const run = v.execute('search.run', { pattern: 'A' }) as Promise<string>;
    // the fake clock's second passes while it runs
    frames!.timers.fire();
    await run;
    expect(words.filter((w) => w.startsWith('Searching')).length).toBeLessThanOrEqual(1);
    expect(words.at(-1)).toMatch(/hits for A on 3 sequences$/);
  });
});
