/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

/** The source protocol's codes: answers and the `seq` format. */
export const SourceIssues = defineIssues('seq.source', {
  'seq.source.malformed': {
    severity: 'error',
    description: 'An answer, a coverage or a `seq` file is not in the protocol’s shape',
  },
  'seq.source.failed': {
    severity: 'warning',
    description: 'A source could not answer a region (a read failed); retried with backoff, then after a cooldown',
  },
  'seq.source.not_served': {
    severity: 'warning',
    description: 'A source does not serve a region it was asked for; never retried',
  },
  'seq.source.too_large': {
    severity: 'info',
    description: 'A region is too large for a source at its level; a coarser level is asked for',
    hint: 'Zoom in, or let the view use a coarser level',
  },
  'seq.source.no_data': {
    severity: 'error',
    description: 'A `seq` source has neither inline data (`value`) nor a location (`uri`)',
  },
});

/** The data store's codes. */
export const StoreIssues = defineIssues('seq.store', {
  'seq.store.over_budget': {
    severity: 'info',
    description: 'What is in view alone is over the store’s byte bound; a coarser level is asked for the widest layers',
  },
  'seq.store.node_conflict': {
    severity: 'warning',
    description: 'A source sent a node again under its key with other content; the held node is kept',
  },
});
