/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printAddress } from '@molstar/core/repr/tree';
import type { SeqSpecBuilder } from '@molstar/seqx-core/lang';

// A UniProt entry as the REST service writes it (`rest.uniprot.org/uniprotkb/<accession>.json`), read by the records
// runtime through the transform below: `sequence.value` the letters, `features[]` the records, `location.start.value`
// and `location.end.value` the ends (from 1, the last position named), `type` the category, `description` the label,
// `featureId` the id; the ends' modifiers (`EXACT`, `OUTSIDE`, `UNSURE`) kept as fields. The exported spec holds the
// transform written out and no service name: a showcase's switch to the service changes only `uri` (files.md).

/** What `uniprotEntry` takes. */
export interface UniprotOptions {
  /** The entry's accession: the name of its own system, and the sequence's `name` (`P04637`). */
  readonly accession: string;
  /**
   * Whether a feature's `location.sequence` (an isoform's accession, `P04637-2`) names the system the row is on: such
   * rows land on that isoform's system, which nothing reads unless a view adds the isoform. Off by default: an entry
   * without isoform features has no such member, and a path found in no record is an error naming it.
   */
  readonly isoforms?: boolean;
}

/** The source node's transform of an entry (the spec's own words, what a reader of the exported spec sees). */
export function uniprotTransform(options: UniprotOptions) {
  return {
    records: 'features',
    letters: 'sequence.value',
    system: options.accession,
    fields: {
      start: 'location.start.value',
      end: 'location.end.value',
      category: 'type',
      label: 'description',
      id: 'featureId',
      start_modifier: 'location.start.modifier',
      end_modifier: 'location.end.modifier',
      ...(options.isoforms ? { target: 'location.sequence' } : {}),
    },
  };
}

/**
 * Writes an entry's nodes to `b`: the source (`entry`, its transform, `uri` the fixture or the service), the sequence
 * it holds (`protein`, by its accession) and one ranges node over all its features (`features`). A track of one
 * feature type is a layer reading `data(type)`: the features with a leading item form, a layer a type from one read.
 */
export function uniprotEntry(b: SeqSpecBuilder, uri: string, options: UniprotOptions) {
  const source = b.data.source('entry', {
    format: 'json',
    uri,
    ...uniprotTransform(options),
    attribution: 'UniProt Consortium, UniProtKB',
    licence: 'CC-BY-4.0',
  });
  const protein = b.data.sequence('protein', { alphabet: 'protein', source, name: options.accession });
  const features = b.data.ranges('features', { target: protein, source });
  return {
    source,
    protein,
    features,
    /** The data of a layer of one feature type (`DNA binding`, `Region`, `Helix`). */
    data: (type: string): string => `${features}${printAddress([{ kind: 'item', params: { category: type } }])}`,
  };
}
