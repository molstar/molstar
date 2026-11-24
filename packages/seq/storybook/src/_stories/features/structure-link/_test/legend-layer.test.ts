/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag medium

import { ComplexViews } from '@molstar/seqx-examples';
import { legendRows } from '@molstar/seqx-viewer/legends';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { openView } from './open.js';

// The Focus and back story's Sequence row (FX29; _spec/README.md "Themes", Legends): it reads as the layer it keys
// ("Sequence letters", its eye "Hide Sequence letters"), lists everything the layer draws (the ATP rule, the chains'
// scheme groups, the other ligands' otherwise), keeps its title and place when the layer is hidden, and the track's
// header keeps its label whole while the layer is hidden. A whole example view settled: medium.

const Letters = '/views/main/complex/letters/letters';
let viewer: Awaited<ReturnType<typeof openView>>;
beforeAll(async () => {
  viewer = await openView(ComplexViews.find((v) => v.id === 'units')!);
}, 20_000);
afterAll(() => viewer.retire());

const row = () => legendRows(viewer, { maxItems: Number.POSITIVE_INFINITY }).find((r) => r.layer === Letters)!;

describe('the Focus and back story’s Sequence legend row', () => {
  it('is titled by its layer, its eye names the same, and it lists what the layer draws', () => {
    const r = row();
    expect(r.title).toBe('Sequence letters');
    expect(r.layers).toEqual([{ address: Letters, name: 'Sequence letters', hidden: false }]);
    const labels = r.entries.map((e) => e.label);
    expect(labels[0]).toBe('ATP');
    expect(labels).toContain('hydrophobic');
    expect(labels.at(-1)).toBe('otherwise');
  });

  it('keeps its title and its place when hidden, and the track’s header keeps its label whole', async () => {
    const label = () =>
      viewer.frame!.layout.sections[0].tracks.find((x) => x.track.address.endsWith('/letters'))!.labelText;
    const before = label();
    await viewer.execute('layer.hide', { address: Letters });
    await viewer.settle();
    const r = row();
    expect(r).toMatchObject({ title: 'Sequence letters', hidden: true });
    expect(r.layers[0]).toMatchObject({ name: 'Sequence letters', hidden: true });
    expect(label()).toBe(before);
    expect(before.endsWith('…')).toBe(false);
    await viewer.execute('layer.show', { address: Letters });
    await viewer.settle();
    expect(row().hidden).toBe(false);
  });
});
