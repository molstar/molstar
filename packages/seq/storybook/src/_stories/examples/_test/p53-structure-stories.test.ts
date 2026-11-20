/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// p53's structure stories (2OCJ), headless on a fake clock, composed as the stories compose them (a viewer of the view
// with `links.hover` listed, and its SVG export): the blocks and chips, every zinc link with both ends in view an arc at
// its data, a zinc's chip lighting its four ligands, the alignment's map taken by `via` (no issue), both themes alike
// in the SVG.

import { type BlockPlacement, DarkTokens, type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { Chains, positionOf, ZincSites } from '@molstar/seqx-examples/p53/structure';
import { LinksHoverChannel, LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const open = async (id: string) => {
  viewer = await openExample(P53Views.find((v) => v.id === id)!, {
    fixtures: true,
    behaviors: [{ provider: LinksHoverProvider, params: undefined }],
  });
  return viewer;
};

const Section = '/views/main/structure';
const runs = (r: ArrayLike<number>) => Array.from(r);
const arcs = (v: SeqViewer) =>
  v.frame!.layout.sections[1].sectionLays.reduce(
    (n, l) => n + (l.type.name === 'pairs' ? (l.lay as PairsLay).n : 0),
    0,
  );

describe('p53’s structure stories, headless', () => {
  it('merged: the chains one block ×4, the zincs one chip ×4; every link an arc on the merged blocks', async () => {
    const v = await open('structure');
    expect(v.layoutIssues.value).toEqual([]);
    const p = v.frame!.layout.sections[1].placement as BlockPlacement;
    expect(p.strips.axes.map((a, i) => [a.system, p.strips.chip[i]])).toEqual([
      ['/data/2ocj/A', 0],
      ['/data/2ocj/E', 1],
    ]);
    expect(p.labels.map((l) => l.replace(/^.*(×\d)$/, '$1'))).toEqual(['×4', '×4']);
    // every unit's links on the block that shows it (merged: the four copies' arcs on one block)
    expect(arcs(v)).toBe(Chains.length * ZincSites.length);
  });

  it('split: a block a chain and a chip a zinc; every zinc’s four ligands an arc, the domains through the alignment’s map', async () => {
    const v = await open('structure-split');
    expect(v.layoutIssues.value).toEqual([]);
    const p = v.frame!.layout.sections[1].placement as BlockPlacement;
    expect(p.strips.axes.map((a) => a.system)).toEqual([
      ...Chains.map((c) => `/data/2ocj/${c.chain}`),
      ...Chains.map((c) => `/data/2ocj/${c.zinc}`),
    ]);
    expect(arcs(v)).toBe(Chains.length * ZincSites.length);
  });

  it('the zinc’s chip lights its four ligands, C176 the zinc', async () => {
    const v = await open('structure-split');
    const owner = {};
    const lit = () => v.highlight.get(LinksHoverChannel)?.location.parts.map((q) => [q.system, runs(q.ranges)]);
    v.hover.set({ section: Section, system: '/data/2ocj/E', position: 0 }, owner);
    expect(lit()).toEqual([['/data/2ocj/A', ZincSites.flatMap((s) => [positionOf(s.auth), positionOf(s.auth) + 1])]]);
    v.hover.set({ section: Section, system: '/data/2ocj/B', position: positionOf(176) }, owner);
    expect(lit()).toEqual([['/data/2ocj/F', [0, 1]]]);
    v.hover.clear(owner);
  });

  for (const id of ['structure', 'structure-split'])
    it(`${id}: the SVG draws the canvas’s lists in both themes, with no issue`, async () => {
      const v = await open(id);
      v.setRailsShown(false);
      for (const tokens of [undefined, DarkTokens]) {
        if (tokens) v.setTokens(tokens);
        v.draw();
        const svg = v.renderSvg({ idPrefix: 's-' });
        await vi.advanceTimersByTimeAsync(400);
        const r = await svg;
        expect(r.issues).toEqual([]);
        expect(recordDrawLists(r.lists, ['data'])).toBe(recordDrawLists(v.frame!.lists, ['data']));
      }
    });
});
