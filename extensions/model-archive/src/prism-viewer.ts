/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AnnotationProvider } from '@molstar/prism-viewer/api/annotations';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { SequenceBehavior, SequenceShown, ViewerExtension } from '@molstar/prism-viewer/spec';
import { PLDDTConfidenceColorThemeProvider } from './quality-assessment/color/plddt.js';
import { QmeanScoreColorThemeProvider } from './quality-assessment/color/qmean.js';
import { QualityAssessment } from './quality-assessment/prop.js';

// Model Archive's quality scores in the Prism viewer (plan 5.1, README decision 27, parity.md §6): commands first. What
// a predicted model carries is shown through `annotations.apply` (pLDDT, QMEAN as the 3D color), the plugin gets the two
// color themes, and the sequence a pLDDT track under the letters (a sequence behavior, sequence.md §3) whose menu colors
// the 3D view by it. The viewer itself knows nothing of pLDDT. The pairwise (PAE) plot waits for a view of its own.

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

const carries = (shown: SequenceShown | undefined, metric: 'pLDDT' | 'qmean') =>
  !!shown && QualityAssessment.isApplicable(shown.structure.model, metric);

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

/** The pLDDT track under the viewer's tracks, for a model that carries pLDDT; its menu and its bands in the legend. */
export const PlddtSequenceBehavior: SequenceBehavior = {
  id: 'model-archive.sequence.plddt',
  title: 'pLDDT track',
  attach(ctx) {
    const sub = ctx.shown.subscribe((shown) => {
      if (!shown || !carries(shown, 'pLDDT')) {
        ctx.offer(undefined);
        return;
      }
      void ctx.patch(shown, [
        {
          op: 'insert',
          parent: shown.data,
          node: {
            kind: 'computed',
            ref: PlddtData,
            params: { provider: 'molstar:plddt', options: { of: shown.complex } },
          },
        },
        {
          op: 'insert',
          parent: shown.section,
          after: shown.after,
          node: {
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
        },
      ]);
      ctx.offer({
        menus: {
          [PlddtTrack]: [{ id: 'color-plddt', label: 'Color by pLDDT', icon: 'paintbrush', theme: PlddtTheme }],
        },
        legend: [{ title: 'Model confidence', items: Bands.map((b) => ({ label: b.label, color: b.color })) }],
      });
    });
    return () => sub.unsubscribe();
  },
};

export const ModelArchiveExtension: ViewerExtension = {
  id: 'model-archive',
  title: 'Model quality',
  plugin: {
    registry: [{ structure: { themes: { color: [PLDDTConfidenceColorThemeProvider, QmeanScoreColorThemeProvider] } } }],
    behaviors: [],
  },
  annotations: QualityAnnotations,
  sequence: [PlddtSequenceBehavior],
  commands: (api) => [
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
  ],
};
