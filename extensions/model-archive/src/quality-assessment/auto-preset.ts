/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { AutoPreset } from '@molstar/plugin/state/builder/structure/representation-presets/auto';
import { StructureRepresentationPresetProvider } from '@molstar/plugin/state/builder/structure/representation-presets/types';
import { PLDDTConfidenceColorThemeProvider } from './color/plddt.js';
import { QmeanScoreColorThemeProvider } from './color/qmean.js';
import { QualityAssessment } from './prop.js';

// The structure's representation out of the box, as v5's viewer did (its "Automatic (w/ Annotation)" preset, the
// viewer's default): Mol*'s automatic representation, colored by the model's confidence when it carries one (pLDDT,
// else QMEAN), so a predicted model opens colored by it rather than being recolored after it loads. Any other structure
// gets the automatic representation as it is.

/** The color theme of the quality score a structure's models carry: pLDDT first, then QMEAN; none for other models. */
export function qualityThemeOf(models: readonly Parameters<typeof QualityAssessment.isApplicable>[0][]) {
  if (models.some((m) => QualityAssessment.isApplicable(m, 'pLDDT'))) return PLDDTConfidenceColorThemeProvider.name;
  if (models.some((m) => QualityAssessment.isApplicable(m, 'qmean'))) return QmeanScoreColorThemeProvider.name;
  return undefined;
}

export const QualityAutoPreset = StructureRepresentationPresetProvider({
  id: 'model-archive.preset.auto',
  display: {
    name: 'Automatic (quality)',
    group: 'Annotation',
    description: "The automatic representation, colored by the model's pLDDT or QMEAN when it carries one.",
  },
  // The default for every structure: those without a score get the automatic representation unchanged.
  isApplicable: () => true,
  params: () => StructureRepresentationPresetProvider.CommonParams,
  async apply(ref, params, plugin) {
    const structure = StateObjectRef.resolveAndCheck(plugin.state.data, ref)?.obj?.data;
    if (!structure) return {};
    const theme = qualityThemeOf(structure.models);
    if (!theme) return await AutoPreset.apply(ref, params, plugin);
    return await AutoPreset.apply(
      ref,
      {
        ...params,
        theme: { ...params.theme, globalName: theme as any, focus: { ...params.theme?.focus, name: theme as any } },
      },
      plugin,
    );
  },
});
