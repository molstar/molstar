/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PartialOptions } from '@molstar/prism-viewer/spec';

// The URL's options (app.md §4): v5's names kept. The load parameters (`pdb`, `url`, …) are read by the
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
  // v5's `map-provider` names the volume streaming server; `volume-streaming-server` wins when both are given.
  const mapProvider = pick('map-provider', ['pdbe', 'rcsb'] as const);
  const fullscreen = flag('show-toggle-fullscreen');
  return {
    renderer: {
      pixelScale: num('pixel-scale'),
      pickScale: num('pick-scale'),
      transparency: pick('transparency', ['auto', 'blended', 'wboit', 'dpoit'] as const),
      resolutionMode: pick('resolution-mode', ['auto', 'scaled', 'native'] as const),
      powerPreference: pick('power-preference', ['default', 'high-performance', 'low-power'] as const),
      graphicsBackend: pick('graphics-backend', ['auto', 'webgpu', 'webgl', 'webgl2', 'webgl1'] as const),
      allowMajorPerformanceCaveat: flag('allow-major-performance-caveat'),
      disableAntialiasing: flag('disable-antialiasing'),
      illumination: flag('illumination'),
    },
    canvas: {
      background: q.get('bg-color') ?? q.get('viewport-background-color') ?? undefined,
      focusBehavior: pick('viewport-focus-behavior', ['default', 'secondary-zoom', 'disabled'] as const),
    },
    representation: {
      // The default palette (accessibility.md): `?palette=colorblind` (or `classic`, `prism`).
      defaultPalette: pick('palette', ['prism', 'classic', 'colorblind'] as const),
    },
    data: {
      pdbProvider: pick('pdb-provider', ['pdbe', 'rcsb', 'pdbj'] as const),
      emdbProvider: pick('emdb-provider', ['pdbe', 'rcsb'] as const),
      volumeStreamingServer:
        q.get('volume-streaming-server') ?? (mapProvider ? MapProviderServers[mapProvider] : undefined),
      volumeStreamingDisabled: flag('volume-streaming-disabled'),
      stateServer: q.get('state-server') ?? undefined,
    },
    ui: {
      collapsed: flag('collapse-left-panel') ?? flag('hide-controls'),
      // The empty canvas's splash: `?splash=quantum` (or `light`, `random`, `off`).
      splash: pick('splash', ['random', 'quantum', 'light', 'off'] as const),
      // v5 hides the button unless asked; Prism shows it: the parameter only says so when given. A partial group
      // keeps the other buttons' defaults (`resolveOptions`).
      ...(fullscreen === undefined ? {} : { buttons: { fullscreen } }),
    },
  };
}

/** v5's `map-provider` values: the volume streaming server each names. */
const MapProviderServers = {
  pdbe: 'https://www.ebi.ac.uk/pdbe/densities',
  rcsb: 'https://maps.rcsb.org',
} as const;

/**
 * The page-wide switches (`debug-mode=1`, `timing-mode=1`): Mol*'s `setDebugMode` / `setTimingMode`, which the app
 * calls before it creates the viewer (app.md §4; v5's index.html does the same). Not per viewer.
 */
export function modesFromUrl(search: string): { debug: boolean; timing: boolean } {
  const q = new URLSearchParams(search);
  const on = (name: string) => q.get(name)?.trim() === '1';
  return { debug: on('debug-mode'), timing: on('timing-mode') };
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
  'graphics-backend',
  'allow-major-performance-caveat',
  'disable-antialiasing',
  'illumination',
  'bg-color',
  'viewport-background-color',
  'viewport-focus-behavior',
  'palette',
  // Opens the assistant at start (loads it then).
  'assistant',
  // The embedded layout: the rail hidden behind a button in the canvas's corner.
  'hide-rail',
  // The empty canvas's splash: random, quantum, light or off.
  'splash',
  // This tab's scene comes back after a reload (the journal's light restore).
  'restore',
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
  // How the structure loads build what they load (the loaders feature's `LoadUrlOptions`): read with the loads.
  'all-models',
  'assembly',
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
