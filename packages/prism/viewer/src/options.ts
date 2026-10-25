/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PartialOptions } from './spec.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsBackendChoice } from '@molstar/plugin/util/graphics-backends';

// The viewer's options (viewer-spec.md §2.2): one PD definition with three readers, `createPrismViewer({ options })`,
// the URL (v5's names, app.md §4) and the Settings tab. `sources/compose-plugin.ts` turns a value into config items and
// canvas props. Defaults match the plugin's config (`auto` where the plugin detects the value itself).

export const PrismViewerOptions = {
  renderer: PD.Group({
    pixelScale: PD.Numeric(1, { min: 0.1, max: 4, step: 0.05 }, { description: 'Resolution of the canvas' }),
    pickScale: PD.Numeric(0.25, { min: 0.05, max: 1, step: 0.05 }, { description: 'Resolution of picking' }),
    transparency: PD.Select<'auto' | 'blended' | 'wboit' | 'dpoit'>('auto', [
      ['auto', 'Automatic'],
      ['blended', 'Blended'],
      ['wboit', 'Weighted, blended order-independent'],
      ['dpoit', 'Dual depth peeling'],
    ]),
    resolutionMode: PD.Select<'auto' | 'scaled' | 'native'>('auto', [
      ['auto', 'Automatic'],
      ['scaled', 'Scaled'],
      ['native', 'Native'],
    ]),
    powerPreference: PD.Select<'default' | 'high-performance' | 'low-power'>('high-performance', [
      ['high-performance', 'High performance'],
      ['low-power', 'Low power'],
      ['default', 'Default'],
    ]),
    graphicsBackend: PD.Select<GraphicsBackendChoice>(
      'auto',
      [
        ['auto', 'Automatic'],
        ['webgpu', 'WebGPU'],
        ['webgl', 'WebGL'],
        ['webgl2', 'WebGL 2 only'],
        ['webgl1', 'WebGL 1 only'],
      ],
      { description: 'Automatic: WebGPU where it has the features of WebGL 2, else WebGL; WebGPU falls back to WebGL' },
    ),
    allowMajorPerformanceCaveat: PD.Boolean(false),
    disableAntialiasing: PD.Boolean(false),
    illumination: PD.Boolean(false, { description: 'Path-traced global illumination' }),
  }),
  canvas: PD.Group({
    /** Empty: the viewer's default background. */
    background: PD.Text('', { description: 'A CSS color; empty for the default' }),
    focusBehavior: PD.Select<'default' | 'secondary-zoom' | 'disabled'>('default', [
      ['default', 'Zoom and show interactions'],
      ['secondary-zoom', 'Zoom only, secondary button'],
      ['disabled', 'Off'],
    ]),
  }),
  representation: PD.Group({
    /**
     * accessibility.md, viewer-spec.md §2.8: the colors the categorical themes (chains, entities, …) give by default
     * (`representation.defaultPalette`).
     */
    defaultPalette: PD.Select<'prism' | 'classic' | 'colorblind'>(
      'prism',
      [
        ['prism', 'Prism'],
        ['classic', 'Classic (Mol*)'],
        ['colorblind', 'Colorblind-safe'],
      ],
      { description: 'The default colors of chains, entities and other categories' },
    ),
  }),
  data: PD.Group({
    pdbProvider: PD.Select<'pdbe' | 'rcsb' | 'pdbj'>('pdbe', [
      ['pdbe', 'PDBe'],
      ['rcsb', 'RCSB PDB'],
      ['pdbj', 'PDBj'],
    ]),
    emdbProvider: PD.Select<'pdbe' | 'rcsb'>('pdbe', [
      ['pdbe', 'PDBe'],
      ['rcsb', 'RCSB'],
    ]),
    // PDBe's, as the plugin's default (ds.litemol.org) no longer answers.
    volumeStreamingServer: PD.Text('https://www.ebi.ac.uk/pdbe/densities'),
    volumeStreamingDisabled: PD.Boolean(false),
    stateServer: PD.Text('https://webchem.ncbr.muni.cz/molstar-state'),
    saccharideCompIdMapType: PD.Select<'default' | 'glycam'>('default', [
      ['default', 'Default'],
      ['glycam', 'GLYCAM'],
    ]),
  }),
  ui: PD.Group({
    showSequence: PD.Boolean(true),
    showHistory: PD.Boolean(true),
    showRemoteStates: PD.Boolean(false),
    /** Start with no panel open. */
    collapsed: PD.Boolean(false),
    /**
     * Open Home at start while nothing is loaded, so a person sees where to load from (the app's; off for a viewer a
     * page drives itself).
     */
    homeWhenEmpty: PD.Boolean(false),
    /**
     * The empty canvas's splash, once the viewer has been idle on it a while: light through Prism's mark as particles
     * (`light`), as a wave (`quantum`), either at random each time it shows (the wave two times in three), or none.
     */
    splash: PD.Select<'random' | 'quantum' | 'light' | 'off'>(
      'random',
      [
        ['random', 'Random'],
        ['quantum', 'Waves'],
        ['light', 'Particles'],
        ['off', 'Off'],
      ],
      { description: 'The light behind the empty canvas’s mark' },
    ),
    /**
     * The toolbar's XR button (classic `PluginConfig.Viewport.ShowXR`): `auto` only where WebXR is supported, `always`
     * shows it everywhere (where it is not supported, it says so), `never` leaves it out.
     */
    xr: PD.Select<'auto' | 'always' | 'never'>('auto', [
      ['auto', 'Where supported'],
      ['always', 'Always'],
      ['never', 'Never'],
    ]),
    /**
     * Whether a viewer's canvas pauses out of view, and its WebGL context joins the page's budget (README decision 29):
     * `auto` does, except in an automated browser (`navigator.webdriver`: a capture or test page renders on purpose).
     */
    pauseOutOfView: PD.Select<'auto' | 'always' | 'never'>('auto', [
      ['auto', 'Unless automated'],
      ['always', 'Always'],
      ['never', 'Never'],
    ]),
    /**
     * Whether the canvas' loop stops asking for frames while nothing changes in view (`PluginConfig.General.PauseLoopWhenIdle`;
     * prism AGENTS.md, idle is idle): a host whose code moves the camera or the scene past Mol*'s API without
     * `canvas3d.requestDraw()` turns it off.
     */
    pauseWhenIdle: PD.Boolean(true),
    /**
     * Where the canvas' hover label shows (v5: the corner): `corner`, at the canvas' bottom right above the toolbar, still
     * and out of the way; `cursor`, beside the pointer.
     */
    hoverLabel: PD.Select<'corner' | 'cursor'>('corner', [
      ['corner', 'Bottom right'],
      ['cursor', 'At the cursor'],
    ]),
    /**
     * How many of the latest steps Undo can take back (sequence.md §7.4): only they keep a checkpoint; older steps stay
     * in History, too old to undo.
     */
    undoSteps: PD.Numeric(10, { min: 1, max: 100, step: 1 }, { description: 'How many steps Undo can take back' }),
    /** Records kept in the browser at most (journal.md §5); keeping one more offers to forget the oldest. */
    keptRecordsMax: PD.Numeric(20, { min: 1, max: 200, step: 1 }),
    buttons: PD.Group({
      reset: PD.Boolean(true),
      screenshot: PD.Boolean(true),
      expand: PD.Boolean(true),
      fullscreen: PD.Boolean(true),
      settings: PD.Boolean(true),
      selectionMode: PD.Boolean(true),
      animation: PD.Boolean(true),
      trajectory: PD.Boolean(true),
    }),
  }),
};

export type PrismViewerOptions = typeof PrismViewerOptions;
export type PrismViewerOptionValues = PD.ValuesFor<PrismViewerOptions>;

/** The options' defaults, a new object each call. */
export function defaultOptions(): PrismViewerOptionValues {
  return PD.getDefaultValues(PrismViewerOptions);
}

/** Options with every group filled in from the defaults (a partial group keeps the rest of its defaults). */
export function resolveOptions(partial?: PartialOptions): PrismViewerOptionValues {
  const values = defaultOptions();
  if (!partial) return values;
  for (const group of Object.keys(values) as (keyof PrismViewerOptionValues)[]) {
    const over = partial[group];
    if (!over) continue;
    const merged: Record<string, unknown> = { ...values[group] };
    for (const [k, v] of Object.entries(over)) {
      if (v === undefined) continue;
      const base = merged[k];
      merged[k] =
        base && typeof base === 'object' && !Array.isArray(base) && v && typeof v === 'object'
          ? { ...(base as object), ...(v as object) }
          : v;
    }
    (values as Record<string, unknown>)[group] = merged;
  }
  return values;
}
