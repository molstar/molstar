/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';

export const proteinEntityTest = MS.core.logic.and([
  MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
  MS.core.str.match([MS.re('(polypeptide|cyclic-pseudo-peptide|peptide-like)', 'i'), MS.ammp('entitySubtype')]),
]);

export const nucleiEntityTest = MS.core.logic.and([
  MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
  MS.core.str.match([MS.re('(nucleotide|peptide nucleic acid)', 'i'), MS.ammp('entitySubtype')]),
]);

/**
 * this is to get non-polymer and peptide terminus components in polymer entities,
 * - non-polymer, e.g. PXZ in 4HIV or generally ACE
 * - carboxy terminus, e.g. FC0 in 4BP9, or ETA in 6DDE
 * - amino terminus, e.g. ARF in 3K4V, or 4MM in 3EGV
 */
export const nonPolymerResidueTest = MS.core.str.match([
  MS.re('non-polymer|(amino|carboxy) terminus|peptide-like', 'i'),
  MS.ammp('chemCompType'),
]);
