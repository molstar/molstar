/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The palette control ("Color ▸" in a track's menu): only where the layer reads a palette, one place (the menu), and a
// palette change is a patch on the layer's own color: choosing the original entry again draws what was drawn, an undo too.

import { recordDrawLists } from '@molstar/seqx-core/render';
import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { type HeadersBehavior, HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { openView } from '../../structure-link/_test/open.js';
import type { UiHistory } from '@molstar/seqx-viewer/ui-history';

let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const menuNames = (v: SeqViewer, section: string, track: string) =>
  v
    .behavior<HeadersBehavior>(HeadersProvider)!
    .menuOf({ section, track })!
    .groups.flat()
    .map((i) => i.name);

describe('the palette control', () => {
  it('is in a track whose layer reads a palette, not in one that ignores it: a letters track of no known alphabet', async () => {
    viewer = await openView(P53Views.find((x) => x.id === 'structure')!);
    const sections = viewer.view.resolved.value!.sections;
    const track = (label: string) => {
      for (const s of sections) for (const t of s.tracks) if (t.label === label) return [s.address, t.address] as const;
      throw new Error(label);
    };
    // 2OCJ's letters: its block's alphabet is not one sequence's, so no scheme of categories is offered over letters
    expect(menuNames(viewer, ...track('2OCJ'))).not.toContain('seq:color');
    // P04637's letters (a protein) and a ranges track by category keep theirs
    expect(menuNames(viewer, ...track('P04637'))).toContain('seq:color');
    expect(menuNames(viewer, ...track('Domains'))).toContain('seq:color');
  });

  it('is one control a track: the menu’s "Color", never also a slot in the header', async () => {
    viewer = await openView(P53Views.find((x) => x.id === 'structure')!);
    const h = viewer.behavior<HeadersBehavior>(HeadersProvider)!;
    const s = viewer.view.resolved.value!.sections[0];
    const track = s.tracks.find((t) => t.label === 'Domains')!;
    const items = h.menuOf({ section: s.address, track: track.address })!.groups.flat();
    expect(items.filter((i) => i.name === 'seq:color')).toHaveLength(1);
    expect(items.find((i) => i.name === 'seq:color')!.slot).toBeFalsy();
    const geometry = viewer.headers!.value.find((g) => g.target.track === track.address)!;
    expect(geometry.slots.map((x) => x.name)).not.toContain('seq:color');
  });

  it('a palette change is a patch: change and change back draw what was drawn first, an undo too (the Variants track)', async () => {
    const view = P53Views.find((v) => v.id === 'themes')!;
    viewer = await openExample(view, { history: true });
    const v = viewer;
    const section = v.view.resolved.value!.sections[0];
    const layer = section.tracks.find((t) => t.label === 'Variants')!.layers[0].address!;
    const draw = async () => {
      await vi.advanceTimersByTimeAsync(600);
      v.draw();
      return recordDrawLists(v.frame!.lists, ['data']);
    };
    const first = await draw();
    const run = async (scheme: string) => {
      await v.execute('layer.channel', { address: layer, channel: 'color', scheme });
      return draw();
    };
    const changed = await run('okabe-ito');
    expect(changed).not.toBe(first);
    // not grey: the same marks, other colors (a field it cannot read draws every mark in the missing color)
    const fills = (s: string) => new Set(s.split('\n').flatMap((l) => l.match(/#[0-9a-f]{8}/g) ?? []));
    expect(fills(changed).size).toBeGreaterThan(2);
    // the menu says which entry it is: the palette the layer's form has
    const entries = () =>
      v
        .behavior<HeadersBehavior>(HeadersProvider)!
        .menuOf({ section: section.address, track: section.tracks.find((t) => t.label === 'Variants')!.address })!
        .groups.flat()
        .find((i) => i.name === 'seq:color')!.items!;
    expect(entries().find((i) => i.checked)?.args?.scheme).toBe('okabe-ito');
    // the original entry again: the first drawing exactly
    expect(await run('categories')).toBe(first);
    // and an undo of the change: the first drawing exactly
    await run('okabe-ito');
    const stack = (v.spec.history as UiHistory).stack;
    expect((await stack.undo()).ok).toBe(true);
    expect(await draw()).toBe(first);
  });

  it('a ranges track’s Color change redraws with the new palette', async () => {
    const view = FeatureViews.find((v) => v.id === 'themes-editing') ?? P53Views.find((v) => v.id === 'themes')!;
    viewer = await openExample(view, { history: true });
    const v = viewer;
    const section = v.view.resolved.value!.sections[0];
    const track = section.tracks.find((t) => t.layers.some((l) => l.type === 'ranges'));
    expect(track).toBeDefined();
    const layer = track!.layers.find((l) => l.type === 'ranges')!.address!;
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    const before = recordDrawLists(v.frame!.lists, ['data']);
    await v.execute('layer.channel', { address: layer, channel: 'color', scheme: 'okabe-ito' });
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    expect(recordDrawLists(v.frame!.lists, ['data'])).not.toBe(before);
  });
});
