/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// Features › Track groups, headless on a fake clock, composed as its stories compose them: each fold state's bands
// and headers, a header's fold by the command its click runs, Collapse all and back, and the export drawing the folds
// shown. The layout's own tests are `core/src/render/_test/groups.test.ts`.

const Section = '/views/main/protein';
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(id: string): Promise<SeqViewer> {
  viewer = await openExample(FeatureViews.find((v) => v.id === id)!, { delayMs: 0 });
  viewer.resize(900);
  await vi.advanceTimersByTimeAsync(100);
  viewer.draw();
  return viewer;
}

const bands = (v: SeqViewer) =>
  v.frame!.layout.sections[0].tracks.map((t) => t.track.address.slice(Section.length + 1));
const headers = (v: SeqViewer) =>
  v.frame!.layout.sections[0].groups!.map((g) => `${g.folded ? '▸' : '▾'} ${g.label}${g.count ? ` ${g.count}` : ''}`);

describe('Features › Track groups', () => {
  it('Opening: Features open, Variants folded to its authored summary, Signals to three strips', async () => {
    const v = await open('groups');
    expect(v.frame!.layout.issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect(bands(v)).toEqual([
      'letters',
      'features/domains',
      'features/secondary',
      'variants/density',
      'signals#summary',
      'signals#strip0',
      'signals#strip1',
      'signals#strip2',
    ]);
    expect(headers(v)).toEqual(['▾ Features', '▸ Variants 2 tracks', '▸ Signals 3 tracks']);
  });

  it('All open: a header line each, the tracks indented, the authored summary left out', async () => {
    const v = await open('groups-open');
    expect(bands(v)).toEqual([
      'letters',
      'features/domains',
      'features/secondary',
      'variants/clinvar',
      'variants/rare',
      'signals/conservation',
      'signals/plddt',
      'signals/disorder',
    ]);
    expect(headers(v)).toEqual(['▾ Features', '▾ Variants', '▾ Signals']);
    expect(v.frame!.layout.sections[0].tracks.slice(1).every((t) => t.indent === 12)).toBe(true);
  });

  it('All folded: one row a group; a header’s click command opens one, expand all is a jump back undoes', async () => {
    const v = await open('groups-folded');
    expect(bands(v)).toEqual([
      'letters',
      'features#summary',
      'variants/density',
      'signals#summary',
      'signals#strip0',
      'signals#strip1',
      'signals#strip2',
    ]);
    // Features' lane: the domains' and the secondary structure's layers in their summary form
    expect(v.frame!.layout.sections[0].tracks[1].track.layers.map((l) => l.type)).toEqual(['ranges', 'secondary']);
    await v.execute('groups.toggle', { address: `${Section}/features` });
    v.draw();
    expect(bands(v).slice(1, 3)).toEqual(['features/domains', 'features/secondary']);
    await v.execute('groups.expandAll', { section: Section });
    v.draw();
    expect(headers(v)).toEqual(['▾ Features', '▾ Variants', '▾ Signals']);
    expect(v.view.back()).toBe(true);
    v.draw();
    expect(headers(v)).toEqual(['▾ Features', '▸ Variants 2 tracks', '▸ Signals 3 tracks']);
    // the export draws the folds shown
    const p = v.renderSvg({ idPrefix: 'g-' });
    await vi.advanceTimersByTimeAsync(400);
    const { svg } = await p;
    expect(svg).toContain('Signals');
    expect(svg).toContain('3 tracks');
    expect(svg).not.toContain('Conservation');
  });
});
