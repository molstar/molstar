/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QualityAssessment } from '@molstar/model-archive-extension/quality-assessment/prop';
import { SbNcbrPartialChargesPropertyProvider } from '@molstar/sb-ncbr-extension/partial-charges/property';
import { StructureRepresentationPresetProvider } from '@molstar/plugin/state/builder/structure/representation-presets/types';
import { AutoPreset } from '@molstar/plugin/state/builder/structure/representation-presets/auto';
import { StateObjectRef } from '@molstar/core/state';

/*
 * The model-archive and SB-NCBR presets belong to behaviors that can be turned off, so they are looked up by id and
 * skipped when they are not registered instead of being imported here.
 */
const QualityAssessmentPLDDTPresetId = 'preset-structure-representation-ma-quality-assessment-plddt';
const QualityAssessmentQmeanPresetId = 'preset-structure-representation-ma-quality-assessment-qmean';
const SbNcbrPartialChargesPresetId = 'sb-ncbr-partial-charges-preset';

export const ViewerAutoPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-viewer-auto',
  display: {
    name: 'Automatic (w/ Annotation)',
    group: 'Annotation',
    description:
      'Show standard automatic representation but colored by quality assessment (if available in the model).',
  },
  isApplicable(a) {
    return (
      !!a.data.models.some((m) => QualityAssessment.isApplicable(m, 'pLDDT')) ||
      !!a.data.models.some((m) => QualityAssessment.isApplicable(m, 'qmean'))
    );
  },
  params: () => StructureRepresentationPresetProvider.CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    const structure = structureCell?.obj?.data;
    if (!structureCell || !structure) return {};

    const builder = plugin.builders.structure.representation;
    const optional: [() => boolean, string][] = [
      [() => structure.models.some((m) => QualityAssessment.isApplicable(m, 'pLDDT')), QualityAssessmentPLDDTPresetId],
      [() => structure.models.some((m) => QualityAssessment.isApplicable(m, 'qmean')), QualityAssessmentQmeanPresetId],
      [
        () => structure.models.some((m) => SbNcbrPartialChargesPropertyProvider.isApplicable(m)),
        SbNcbrPartialChargesPresetId,
      ],
    ];
    for (const [isApplicable, id] of optional) {
      if (!isApplicable()) continue;
      // skipped when the behavior that owns the preset is not active
      const preset = builder.resolveProvider(id);
      if (preset) return await preset.apply(ref, params, plugin);
    }
    return await AutoPreset.apply(ref, params, plugin);
  },
});
