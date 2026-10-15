/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

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
    preferWebgl1: PD.Select<'auto' | 'yes' | 'no'>('auto', [
      ['auto', 'Automatic'],
      ['yes', 'Yes'],
      ['no', 'No'],
    ]),
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
    /** README §3.7: the categorical themes' colors (`canvas.palette`). */
    palette: PD.Select<'default' | 'colorblind'>(
      'default',
      [
        ['default', 'Default'],
        ['colorblind', 'Colorblind-safe (Okabe–Ito)'],
      ],
      { description: 'The colors of chains, entities and other categories' },
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
export function resolveOptions(
  partial?: {
    readonly [G in keyof PrismViewerOptionValues]?: Partial<PrismViewerOptionValues[G]>;
  },
): PrismViewerOptionValues {
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
