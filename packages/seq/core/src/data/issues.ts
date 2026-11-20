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
  'seq.source.unreadable': {
    severity: 'error',
    description: 'A source’s file could not be read when it was opened',
  },
  'seq.source.no_data': {
    severity: 'error',
    description: 'A source has neither inline data (`value`) nor a location (`uri`) its format reads',
  },
  'seq.source.rows.field': {
    severity: 'warning',
    description:
      'A field of an alignment’s rows (an inline table’s column) is not one value a row: it is left out, the alignment kept',
    hint: 'Give the column one value a row, `null` where a row has none',
  },
  'seq.source.rows.annotation': {
    severity: 'info',
    description: 'A Stockholm `#=GS` line names a sequence the alignment has no row of: not read',
  },
});

/** The data store's codes. */
export const StoreIssues = defineIssues('seq.store', {
  'seq.store.over_budget': {
    severity: 'info',
    description: 'What is in view alone is over the store’s byte bound; a coarser level is asked for the widest layers',
  },
  'seq.store.too_many_tiles': {
    severity: 'warning',
    description: 'A wanted region covers more tiles than the store asks for at once; it is not requested',
    hint: 'Ask for it at a coarser level',
  },
  'seq.store.limit_exceeded': {
    severity: 'warning',
    description: 'The store reported as many issues as it keeps; later ones are left out',
  },
  'seq.store.node_conflict': {
    severity: 'warning',
    description: 'A source sent a node again under its key with other content; the held node is kept',
  },
});

/** The archive's codes (`archive.ts`; EC-91). */
export const ArchiveIssues = defineIssues('seq.archive', {
  'seq.archive.unreadable': {
    severity: 'error',
    description: 'The bytes are not a zip archive, or it is damaged (cut short, a failed checksum, an unknown method)',
  },
  'seq.archive.bounds': {
    severity: 'error',
    description: 'An archive is past a bound: its entries, their expanded size, one entry’s, or a ratio of expansion',
    hint: 'Raise the bound if the archive is trusted, else do not open it',
  },
  'seq.archive.name': {
    severity: 'error',
    description: 'An entry’s name is no path under the archive’s root (absolute, with `..`, a backslash, a control)',
  },
  'seq.archive.index': {
    severity: 'error',
    description: 'The archive’s index is missing or not in its shape, or lists a file the archive does not hold',
  },
  'seq.archive.spec': {
    severity: 'error',
    description: 'The archive’s spec is missing, past the spec bounds or not a tree',
  },
  'seq.archive.data': {
    severity: 'error',
    description: 'A data file the spec names could not be read to be saved in the archive',
  },
});
