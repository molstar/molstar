/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

// The genome runtimes' shared codes (_spec/README.md "Data", Genome runtimes' internals): what the BGZF container says
// whatever the format in it. Not published (`./runtimes/internal/*`).

export const BgzfIssues = defineIssues('seq.source.bgzf', {
  'seq.source.bgzf.eof': {
    severity: 'warning',
    description: 'The BGZF file does not end with its EOF block: it may be cut short; it is read as it is',
    hint: 'Download the file again, or bgzip it again',
  },
});
