/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, it } from 'vitest';

// The real INSDC record L09137.2 (pUC19c) as NCBI serves it in GenBank form: not fetched yet (its download waits for
// the author's approval in the user's own words), so this test is skipped. Once the file is recorded under the
// examples' `fixtures/plasmid/` with its `sources.json` row, it reads it: one cyclic system of 2,686 letters, the
// letters equal to the recorded FASTA's, its `source` feature, no refusal but the header fields and `/translation`.

describe.skip('the real L09137.2 GenBank record (not recorded yet)', () => {
  it('reads as one cyclic sequence equal to the recorded letters, its features each an item', () => {});
});
