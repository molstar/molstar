/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PartialOptions } from '@molstar/prism-viewer/spec';

// The URL's options (app.md §4, README §3.2): v5's names kept. The load parameters (`pdb`, `url`, …) are read by the
// loaders feature once the viewer is ready; these are the ones read at creation.

export function optionsFromUrl(search: string): PartialOptions {
  const q = new URLSearchParams(search);
  const num = (name: string) => {
    const v = q.get(name);
    if (v === null || v.trim() === '') return undefined;
    const n = Number(v);
    return Number.isFinite(n) ? n : undefined;
  };
  const flag = (name: string) => {
    const v = q.get(name);
    if (v === null) return undefined;
    return v === '' || v === '1' || v.toLowerCase() === 'true';
  };
  const pick = <T extends string>(name: string, allowed: readonly T[]) => {
    const v = q.get(name);
    return v !== null && (allowed as readonly string[]).includes(v) ? (v as T) : undefined;
  };
  return {
    renderer: {
      pixelScale: num('pixel-scale'),
      pickScale: num('pick-scale'),
      transparency: pick('transparency', ['auto', 'blended', 'wboit', 'dpoit'] as const),
      resolutionMode: pick('resolution-mode', ['auto', 'scaled', 'native'] as const),
      powerPreference: pick('power-preference', ['default', 'high-performance', 'low-power'] as const),
      preferWebgl1: flag('prefer-webgl1') === undefined ? undefined : flag('prefer-webgl1') ? 'yes' : 'no',
      allowMajorPerformanceCaveat: flag('allow-major-performance-caveat'),
      disableAntialiasing: flag('disable-antialiasing'),
      illumination: flag('illumination'),
    },
    canvas: {
      background: q.get('bg-color') ?? q.get('viewport-background-color') ?? undefined,
      focusBehavior: pick('viewport-focus-behavior', ['default', 'secondary-zoom', 'disabled'] as const),
    },
    data: {
      pdbProvider: pick('pdb-provider', ['pdbe', 'rcsb', 'pdbj'] as const),
      emdbProvider: pick('emdb-provider', ['pdbe', 'rcsb'] as const),
      volumeStreamingServer: q.get('volume-streaming-server') ?? undefined,
      volumeStreamingDisabled: flag('volume-streaming-disabled'),
      stateServer: q.get('state-server') ?? undefined,
    },
    ui: {
      collapsed: flag('collapse-left-panel') ?? flag('hide-controls'),
    },
  };
}
