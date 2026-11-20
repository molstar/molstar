/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

/** Seq io's codes: a structure description read into nodes (_spec/README.md "Seq io"). */
export const StructureIssues = defineIssues('seq.io', {
  'seq.io.structure.columns': {
    severity: 'error',
    description: 'A table of a structure description has columns of unequal lengths: the table is left out',
    hint: 'The adapter fills every column of a table a row each',
  },
  'seq.io.structure.rows': {
    severity: 'warning',
    description:
      'Rows of a structure description name a row of another table that is not there, or a position past its entity: left out',
  },
  'seq.io.structure.id': {
    severity: 'warning',
    description:
      'An entity’s or a unit’s identifier is no id, or is given twice, or is a key the description’s other nodes take: left out',
  },
  'seq.io.structure.letters': {
    severity: 'info',
    description:
      'A polymer’s letters were read from its components, some of which have no standard letter: `X` (protein) or `N` there',
  },
  'seq.io.structure.hydrogen': {
    severity: 'info',
    description: 'A molecule’s atoms listed hydrogens: left out (heavy atoms only, EC-64)',
  },
});
