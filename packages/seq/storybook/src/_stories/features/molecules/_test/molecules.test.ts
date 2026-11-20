/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type BlockPlacement, DarkTokens, recordDrawLists } from '@molstar/seqx-core/render';
import { MoleculeViews } from '@molstar/seqx-examples';
import { LinksHoverChannel, LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// The molecules' stories, headless on a fake clock, each composed as its story composes it (a viewer of the view with
// `links.hover` listed, and its SVG export): the canvas and the SVG draw the same lists in both themes, every chip its
// label, and the hover the first story asks for lights what it says.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(id: string): Promise<SeqViewer> {
  viewer = await openExample(MoleculeViews.find((v) => v.id === id)!, {
    delayMs: 0,
    wait: 50,
    behaviors: [{ provider: LinksHoverProvider, params: undefined }],
  });
  return viewer;
}

async function svgOf(v: SeqViewer) {
  const p = v.renderSvg({ idPrefix: 'x-' });
  await vi.advanceTimersByTimeAsync(50);
  return p;
}

describe('the molecules’ stories, headless', () => {
  for (const view of MoleculeViews)
    it(`${view.title}: the canvas and the SVG draw the same lists in both themes; every chip its label`, async () => {
      const v = await open(view.id);
      expect(v.layoutIssues.value).toEqual([]);
      v.setRailsShown(false);
      v.draw();
      const light = await svgOf(v);
      expect(light.issues).toEqual([]);
      for (const pass of ['data', 'labels'] as const)
        expect(recordDrawLists(light.lists, [pass])).toBe(recordDrawLists(v.frame!.lists, [pass]));
      const data = recordDrawLists(v.frame!.lists, ['data']);
      for (const s of v.frame!.layout.sections) {
        const blocks = s.section.blocks;
        if (!blocks) continue;
        const p = s.placement as BlockPlacement;
        for (let i = 0; i < p.strips.n; i++)
          if (p.strips.chip[i]) expect(data).toContain(`text "${blocks[p.strips.block[i]].label}"`);
      }
      v.setTokens(DarkTokens);
      v.draw();
      const dark = await svgOf(v);
      expect(dark.issues).toEqual([]);
      expect(recordDrawLists(dark.lists, ['data'])).toBe(recordDrawLists(v.frame!.lists, ['data']));
    });

  it('Modified residues: `code` shows the code where it fits (16 positions: 2MG), the dot where it does not', async () => {
    const v = await open('molecules-modified');
    v.draw();
    const data = recordDrawLists(v.frame!.lists, ['data']);
    expect(data).toContain('text "2MG"');
    // the four whole-tRNA sections: no code fits, every modified residue its dot
    expect(data.match(/text "PSU"/g)).toBeNull();
  });

  it('A section of a complex: the glycan’s chip lights N47; N47 lights the glycan', async () => {
    const v = await open('molecules-complex');
    const owner = {};
    const section = '/views/main/complex';
    v.hover.set({ section, system: '/data/complex/glycan', position: 0 }, owner);
    const lit = () => v.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    expect(lit()).toEqual([['/data/complex/A', [46, 47]]]);
    v.hover.set({ section, system: '/data/complex/A', position: 46 }, owner);
    expect(lit()).toEqual([['/data/complex/glycan', [0, 5]]]);
    // drawn: the highlight pass outlines the glycan's chip, and nothing else
    v.draw();
    const p = v.frame!.layout.sections[0].placement as BlockPlacement;
    const chip = p.strips.axes.find((x) => x.system === '/data/complex/glycan')!;
    const rects = recordDrawLists(v.frame!.lists, ['highlight'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' ').slice(1, 5).map(Number));
    expect(rects).toHaveLength(1);
    const [x, , w] = rects[0];
    expect(x).toBeGreaterThanOrEqual(Math.floor(chip.x0) - 1);
    expect(x + w).toBeLessThanOrEqual(Math.ceil(chip.x0 + chip.width) + 1);
    expect(w).toBeGreaterThan(chip.width / 2);
    v.hover.clear(owner);
  });
});
