/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { mergeRegistryEntries } from '../merge.js';

const a = { name: 'a' } as any;
const b = { name: 'b' } as any;
const c = { name: 'c' } as any;

describe('mergeRegistryEntries', () => {
  it('concatenates the lists of the same registry in order', () => {
    const merged = mergeRegistryEntries(
      { structure: { representations: [a] } },
      { structure: { representations: [b], themes: { color: [c] } } },
      { volume: { representations: [c] } },
    );
    expect(merged).toEqual({
      structure: { representations: [a, b], themes: { color: [c] } },
      volume: { representations: [c] },
    });
  });

  it('lists a provider once', () => {
    const merged = mergeRegistryEntries(
      { structure: { representations: [a, b] } },
      { structure: { representations: [b, a, c] } },
    );
    expect(merged.structure!.representations).toEqual([a, b, c]);
  });

  it('merges nested records and leaves the input unchanged', () => {
    const one: PluginRegistryEntry = { structure: { themes: { color: [a] }, presets: { hierarchy: [b] } } };
    const two: PluginRegistryEntry = { structure: { themes: { size: [c] }, presets: { hierarchy: [a] } } };
    const merged = mergeRegistryEntries(one, two);
    expect(merged.structure!.themes).toEqual({ color: [a], size: [c] });
    expect(merged.structure!.presets!.hierarchy).toEqual([b, a]);
    expect(one.structure!.presets!.hierarchy).toEqual([b]);
    expect(two.structure!.themes).toEqual({ size: [c] });
  });

  it('merges nothing into an empty entry', () => {
    expect(mergeRegistryEntries()).toEqual({});
  });
});
