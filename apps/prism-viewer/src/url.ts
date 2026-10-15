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
      // Prism's own (README §3.7): `?palette=colorblind`.
      palette: pick('palette', ['default', 'colorblind'] as const),
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

/**
 * Every parameter the URL table knows (app.md §4), besides the load parameters (`LoadUrlParameters`, the loaders
 * feature's): the options' names, and the ones the app or v5 handled that the viewer drops or reads elsewhere.
 */
export const KnownUrlParameters: readonly string[] = [
  'pixel-scale',
  'pick-scale',
  'transparency',
  'resolution-mode',
  'power-preference',
  'prefer-webgl1',
  'allow-major-performance-caveat',
  'disable-antialiasing',
  'illumination',
  'bg-color',
  'viewport-background-color',
  'viewport-focus-behavior',
  'palette',
  'pdb-provider',
  'emdb-provider',
  'map-provider',
  'volume-streaming-server',
  'volume-streaming-disabled',
  'state-server',
  'collapse-left-panel',
  'hide-controls',
  'show-toggle-fullscreen',
  'debug-mode',
  'timing-mode',
  'pick-padding',
];

/** The parameters of a query string that are in neither list (names as typed, once each; matched case-insensitively). */
export function unknownUrlParameters(search: string, loadParameters: readonly string[]): string[] {
  const known = new Set([...KnownUrlParameters, ...loadParameters]);
  const unknown = new Set<string>();
  for (const name of new URLSearchParams(search).keys()) {
    if (!known.has(name.toLowerCase())) unknown.add(name);
  }
  return [...unknown];
}
