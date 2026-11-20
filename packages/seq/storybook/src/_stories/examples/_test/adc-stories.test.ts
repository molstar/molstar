/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The T-DM1 construct's story (generated), headless on a fake clock, composed as the story composes it (a viewer with
// `links.hover` and its SVG export): the occupancies summing to the DAR, every conjugation link an arc colored by its band of
// occupancy (a more occupied site a different color than a less occupied one), the payload lighting every site.

import { type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { AntibodyViews } from '@molstar/seqx-examples';
import { AdcChains, conjugationSites, Dar, PayloadComponents } from '@molstar/seqx-examples/antibody/adc';
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

const open = async () => {
  viewer = await openExample(AntibodyViews.find((v) => v.id === 'adc')!, {
    gene: '',
    behaviors: [{ provider: LinksHoverProvider, params: undefined }],
  });
  return viewer;
};

describe('T-DM1, headless', () => {
  it('the occupancies sum to the DAR; every link an arc, colored by its band of occupancy; the payload’s ids placeholders', async () => {
    const sites = conjugationSites();
    expect(sites.reduce((n, s) => n + s.occupancy, 0)).toBeCloseTo(Dar, 9);
    expect(sites.every((s) => s.occupancy > 0 && s.occupancy <= 1)).toBe(true);
    expect(PayloadComponents.map((c) => c.code)).toEqual(['UNL', 'UNL']);
    const v = await open();
    expect(v.layoutIssues.value).toEqual([]);
    const lay = v.frame!.layout.sections[0].sectionLays.find((l) => l.type.name === 'pairs')!.lay as PairsLay;
    expect(lay.n).toBe(sites.length);
    // the most and the least occupied sites drawn in different colors
    const strokes = new Set(
      recordDrawLists(v.frame!.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('path'))
        .map((l) => /#[0-9a-f]{8}/.exec(l)?.[0]),
    );
    expect(strokes.size).toBeGreaterThan(1);
    const svg = v.renderSvg({ idPrefix: 'adc-' });
    await vi.advanceTimersByTimeAsync(400);
    expect((await svg).issues).toEqual([]);
  });

  it('the payload lights every conjugation site (`links.hover`)', async () => {
    const v = await open();
    const owner = {};
    v.hover.set({ section: '/views/main/adc', system: '/data/adc/payload', position: 0 }, owner);
    const lit = v.highlight.get(LinksHoverChannel)!.location.parts;
    let count = 0;
    for (const q of lit) for (let k = 0; k < q.ranges.length; k += 2) count += q.ranges[k + 1] - q.ranges[k];
    expect(lit.map((q) => q.system).sort()).toEqual(AdcChains.map((c) => `/data/adc/${c.id}`).sort());
    expect(count).toBe(conjugationSites().length);
    v.hover.clear(owner);
  });
});
