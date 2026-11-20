/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

/** The structure translation's codes (_spec/README.md "Structures"). */
export const StructureSyncIssues = defineIssues('seq.structure', {
  'seq.structure.complex': {
    severity: 'warning',
    description: 'A binding names no complex of the model: nothing is carried onto it',
    hint: 'Bind a complex the spec or a source holds (its address), or leave `complexes` out for every complex',
  },
  'seq.structure.gated': {
    severity: 'info',
    description:
      'A Location carried onto units is past the size gate (runs or positions): whole units are sent, or nothing',
    hint: 'Raise `limits` if the structure viewer takes it, or `over: skip` to send nothing',
  },
  'seq.structure.unplaced': {
    severity: 'info',
    description:
      'Positions a structure viewer sent are on no unit, unit atoms or molecule of the bound complexes: left out',
  },
  'seq.structure.rows': {
    severity: 'warning',
    description: 'An annotation whose rows are asked is no ranges or pairs annotation of the model: no rows',
  },
  'seq.structure.layer': {
    severity: 'warning',
    description: 'A layer whose colors are asked is no letters, values or ranges layer, or reads no data: no runs',
  },
});
