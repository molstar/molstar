/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FilesViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

const view = (id: string) => FilesViews.find((v) => v.id === id)!;
const story = (id: string): StoryObj => ({ render: () => exampleStory(view(id), { files: true }) });

const meta: Meta = { title: 'Features/Files' };
export default meta;

/**
 * p53 from UniProt’s own JSON, the REST entry of P04637 as the service writes it (trimmed, never reshaped): the letters
 * (`sequence.value`) and the features (`features[]`) are read through the source node’s transform, a track a feature
 * type (a layer with a leading item form). The file’s head and the source node are printed under the view.
 */
export const UniprotEntry: StoryObj = story('uniprot');

/** The same six domains from a TSV table: a header row, positions from 1 and the last one named. */
export const TableTsv: StoryObj = story('table-tsv');

/** The table as CSV (RFC 4180): it draws the same as every other form of the file. */
export const TableCsv: StoryObj = story('table-csv');

/** The table as JSON records: an array of objects, a row each. */
export const TableJsonRecords: StoryObj = story('table-json-aos');

/** The table as JSON columns: an object of columns, one array a field. */
export const TableJsonColumns: StoryObj = story('table-json-soa');

/** The table as MessagePack records (its head is the decoded JSON). */
export const TableMessagePackRecords: StoryObj = story('table-msgpack-aos');

/** The table as MessagePack columns. */
export const TableMessagePackColumns: StoryObj = story('table-msgpack-soa');

/**
 * One table for two sequences: a `target` column names each row’s sequence (GFF3’s `seqid` does the same), and a row
 * with an empty target is on the source’s own system. Each section shows only its own rows: p53 its domains (the NLS
 * row among them, with no target), CDK2 its loops and helix.
 */
export const ManySequences: StoryObj = story('many-sequences');

/**
 * The plasmid’s features as GFF3: nine columns, a `region` row saying the construct is circular, and bla and a primer
 * site ending past its length, one row each, as NCBI writes a ring’s features. Each is drawn as one item across the
 * origin, in the ring and in its detail.
 */
export const PlasmidGff3: StoryObj = story('plasmid-gff3');
