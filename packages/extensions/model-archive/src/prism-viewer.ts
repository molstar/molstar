/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AnnotationProvider } from '@molstar/prism-viewer/api/annotations';
import { ApiResult } from '@molstar/prism-viewer/api/area';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerApi } from '@molstar/prism-viewer/api/area';
import type { CommandDefinition } from '@molstar/prism-ui';
import type { SequenceTrackProvider, ViewerExtension } from '@molstar/prism-viewer/spec';
import { PLDDTConfidenceColorThemeProvider } from './quality-assessment/color/plddt.js';
import { QmeanScoreColorThemeProvider } from './quality-assessment/color/qmean.js';
import { ModelArchiveApi } from './pae/operations.js';
import { PluginConfig } from '@molstar/plugin/config';
import { QualityAutoPreset } from './quality-assessment/auto-preset.js';
import { QualityAssessment } from './quality-assessment/prop.js';

declare module '@molstar/prism-viewer/api' {
  // The `model-archive` area: the extension's own (api.md §2, extensions add areas).
  interface ExtensionApis {
    'model-archive': ModelArchiveApi;
  }
}

// Model Archive's quality scores in the Prism viewer (README decision 27, parity.md §6): commands first. What
// a predicted model carries is shown through `annotations.apply` (pLDDT, QMEAN as the 3D color), the plugin gets the two
// color themes, and the sequence a pLDDT track under the letters (a track provider, sequence.md §3) whose menu colors
// the 3D view by it. The viewer itself knows nothing of pLDDT. The pairwise metric (PAE) is the `model-archive` area's
// `model-archive.pae` and a card of its own (`pae/`, a chunk loaded with the views).

/** The pLDDT track's bands: AlphaFold DB's four, the official hues softened to the sequence palette's key. */
const Bands = [
  { label: 'Very high (pLDDT > 90)', from: 90, color: '#3f6fc4' },
  { label: 'Confident (70–90)', from: 70, color: '#7ec1e3' },
  { label: 'Low (50–70)', from: 50, color: '#e8cc5a' },
  { label: 'Very low (< 50)', from: 0, color: '#e89468' },
] as const;

/** The track's color theme (Seq*'s `field_palette`): a discrete palette over the value, the bands' lower bounds. */
const PlddtBandsTheme = {
  theme: 'field_palette',
  params: {
    field: 'value',
    palette: { kind: 'discrete', mode: 'absolute', colors: [...Bands].reverse().map((b) => [b.color, b.from]) },
  },
};

export const PlddtTrack = 'model-archive.plddt-track';
const PlddtData = 'model-archive.plddt';

const PlddtTheme = PLDDTConfidenceColorThemeProvider.name;
const QmeanTheme = QmeanScoreColorThemeProvider.name;

export const QualityAnnotations: readonly AnnotationProvider[] = [
  {
    id: 'plddt',
    title: 'pLDDT',
    description: "A predicted model's per-residue confidence, as color",
    isApplicable: (s) => QualityAssessment.isApplicable(s.model, 'pLDDT'),
    apply: (ctx, structure) => ctx.api.representation.color({ target: 'polymer', theme: PlddtTheme, structure }),
  },
  {
    id: 'qmean',
    title: 'QMEAN',
    description: "A model's per-residue QMEAN score, as color",
    isApplicable: (s) => QualityAssessment.isApplicable(s.model, 'qmean'),
    apply: (ctx, structure) => ctx.api.representation.color({ target: 'polymer', theme: QmeanTheme, structure }),
  },
];

/**
 * The pLDDT track (sequence.md §3, a track provider): in by default under the letters of a model that carries pLDDT,
 * the scores as a line in the bands; its menu colors the 3D view by them, its legend has the four bands.
 */
export const PlddtTrackProvider: SequenceTrackProvider = {
  id: 'model-archive.sequence.plddt',
  label: 'pLDDT',
  order: 30,
  track: PlddtTrack,
  applies: (structure) => QualityAssessment.isApplicable(structure.model, 'pLDDT'),
  describe: () => 'Line · confidence bands',
  build: (ctx) => ({
    data: [{ kind: 'computed', ref: PlddtData, params: { provider: 'molstar:plddt', options: { of: ctx.complex } } }],
    track: {
      kind: 'track',
      ref: PlddtTrack,
      params: { label: '', height: 2 },
      children: [
        {
          kind: 'layer',
          ref: `${PlddtData}.values`,
          // The scores as a line in the bands (Seq*'s case 8); Text draws them as a strip.
          params: { type: 'values', data: PlddtData, encoding: 'line', domain: [0, 100], color: PlddtBandsTheme },
        },
      ],
    },
  }),
  menu: [{ id: 'color-plddt', label: 'Color by pLDDT', icon: 'paintbrush', theme: PlddtTheme }],
  legend: [{ title: 'Model confidence', items: Bands.map((b) => ({ label: b.label, color: b.color })) }],
};

export const ModelArchiveExtension: ViewerExtension = {
  id: 'model-archive',
  title: 'Model quality',
  plugin: {
    registry: [
      {
        structure: {
          themes: { color: [PLDDTConfidenceColorThemeProvider, QmeanScoreColorThemeProvider] },
          presets: { representation: [QualityAutoPreset] },
        },
      },
    ],
    behaviors: [],
    // A predicted model opens colored by its confidence (v5's viewer's default preset): the automatic representation,
    // by pLDDT or QMEAN where the model carries one (quality-assessment/auto-preset.ts).
    config: [[PluginConfig.Structure.DefaultRepresentationPreset, QualityAutoPreset.id]],
  },
  api: { 'model-archive': (ctx) => new ModelArchiveApi(ctx) },
  annotations: QualityAnnotations,
  sequence: { tracks: [PlddtTrackProvider] },
  // The color commands run on `annotations.apply`, a feature a minimal viewer may leave out (the host skips the
  // provider there too): without it only the PAE commands are brought.
  requires: ['structure'],
  commands: (api) => [
    ...(api['model-archive'] ? paeCommands(api, api['model-archive']) : []),
    ...(api.area('annotations') ? colorCommands(api) : []),
  ],
  ui: async (ctx) => (await import('./pae/ui.js')).paeUI(ctx),
};

function colorCommands(api: ViewerApi): CommandDefinition<any>[] {
  return [
    operationCommand(api, api.annotations.apply, {
      id: 'model-archive.color-plddt',
      title: 'Color by pLDDT',
      category: 'Structure',
      icon: 'paintbrush',
      alias: ['plddt', 'confidence', 'alphafold'],
      values: { annotation: 'plddt' },
    }),
    operationCommand(api, api.annotations.apply, {
      id: 'model-archive.color-qmean',
      title: 'Color by QMEAN',
      category: 'Structure',
      icon: 'paintbrush',
      alias: ['qmean', 'model quality'],
      values: { annotation: 'qmean' },
    }),
  ];
}

/** The PAE card's id (`pae/ui.tsx`), here so the commands do not load the views. */
const PaeCard = 'model-archive.overlay.pae';

function paeCommands(api: ViewerApi, area: ModelArchiveApi): CommandDefinition<any>[] {
  return [
    operationCommand(api, area.pae, {
      category: 'Structure',
      icon: 'gauge',
      alias: ['pae', 'predicted aligned error', 'alphafold', 'model quality'],
      action: 'Read',
    }),
    operationCommand(api, area.fetchPae, {
      category: 'Structure',
      icon: 'download',
      alias: ['pae', 'predicted aligned error', 'alphafold', 'download pae', 'afdb', 'model archive'],
      action: 'Fetch',
    }),
    {
      id: 'model-archive.show-pae',
      title: 'Show the PAE plot',
      category: 'Structure',
      icon: 'gauge',
      alias: ['pae', 'predicted aligned error', 'pae plot'],
      // Opens the card: a view, not a change.
      async run() {
        const r = await api.view.open({ card: PaeCard });
        if (ApiResult.isFailure(r)) throw new Error(r.message);
        return { message: 'Opened the PAE plot.', ops: `view.open({ card: "${PaeCard}" })` };
      },
    },
  ];
}
