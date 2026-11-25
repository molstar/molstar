/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

/** The mmCIF runtime's own codes (the rest are `seq.io.structure.*`, said by `structureNodes`, and `seq.source.*`). */
export const MmcifIssues = defineIssues('seq.io.mmcif', {
  'seq.io.mmcif.assembly': {
    severity: 'error',
    description: 'The assembly asked is not in the file: the source is not read',
    hint: 'Name an assembly of `pdbx_struct_assembly_gen`, or none for the deposited model',
  },
  'seq.io.mmcif.operators': {
    severity: 'warning',
    description:
      'An assembly’s operator expression names more operators than the units bound (`limits.maxUnits`): the ids past it are never made, the rest left out',
  },
  'seq.io.mmcif.altloc': {
    severity: 'info',
    description:
      'Atoms at alternate locations: the first `label_alt_id` an atom name has is read, the others dropped (a position is one atom)',
  },
  'seq.io.mmcif.links': {
    severity: 'info',
    description:
      '`struct_conn` rows not read as links: an end on a symmetry mate, or on a residue no unit read has (a water when waters are not read)',
  },
});
