/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  QualityAssessmentPLDDTPreset,
  QualityAssessmentQmeanPreset,
} from '@molstar/model-archive-extension/quality-assessment/behavior';
import { QualityAssessment } from '@molstar/model-archive-extension/quality-assessment/prop';
import { SbNcbrPartialChargesPreset, SbNcbrPartialChargesPropertyProvider } from '@molstar/sb-ncbr-extension';
import { StructureRepresentationPresetProvider } from '@molstar/plugin/state/builder/structure/representation-presets/types';
import { AutoPreset } from '@molstar/plugin/state/builder/structure/representation-presets/auto';
import { StateObjectRef } from '@molstar/core/state';

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

    if (!!structure.models.some((m) => QualityAssessment.isApplicable(m, 'pLDDT'))) {
      return await QualityAssessmentPLDDTPreset.apply(ref, params, plugin);
    } else if (!!structure.models.some((m) => QualityAssessment.isApplicable(m, 'qmean'))) {
      return await QualityAssessmentQmeanPreset.apply(ref, params, plugin);
    } else if (!!structure.models.some((m) => SbNcbrPartialChargesPropertyProvider.isApplicable(m))) {
      return await SbNcbrPartialChargesPreset.apply(ref, params, plugin);
    } else {
      return await AutoPreset.apply(ref, params, plugin);
    }
  },
});
