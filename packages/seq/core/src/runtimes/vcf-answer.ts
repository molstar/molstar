/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import type {
  MatrixTable,
  RangesField,
  RangesList,
  RangesTable,
  RegionAnswer,
  RegionRequest,
} from '../data/protocol.js';
import { failure } from './internal/answers.js';
import { densityOf } from './internal/levels.js';
import type { matrixAnswer } from './vcf-genotypes.js';
import { column, FixedFields, rowsIn, type VcfRows, type VcfSchema } from './vcf-parse.js';
import type { GenotypeSchema, MatrixBounds } from './vcf-samples.js';

// A VCF's answer of a region (_spec/README.md "Data", VCF), whole and region reads alike: its sites as ranges (a
// density above level 1 by region), its genotypes as a matrix.

/** What a VCF answers besides its sites: its matrix's samples, its bounds, and whether coarse levels are densities. */
export interface VcfAnswerOptions {
  readonly maxRows: number;
  /** The matrix's states and channels; none: no matrix is served (a file without samples). */
  readonly genotypes?: GenotypeSchema;
  readonly matrix: MatrixBounds;
  /** A coarse level's `ranges` as a density (by region: counted from the tile), else the list (read whole). */
  readonly density: boolean;
}

/**
 * The answer of a system's rows over a request (whole and region reads alike): `ranges` a list (or, by region above
 * level 1 and not as items, a density counted from the tile's items), `matrix` the items × the samples held (a summary
 * above level 1).
 */
export async function vcfAnswer(
  rows: VcfRows | undefined,
  schema: VcfSchema,
  request: RegionRequest,
  options: VcfAnswerOptions,
): Promise<RegionAnswer> {
  const ranges = request.shapes.includes('ranges');
  const matrix = request.shapes.includes('matrix') && options.genotypes !== undefined;
  // the matrix's reader loaded on the first matrix request
  const matrixOf: typeof matrixAnswer | undefined = matrix
    ? (await import('./vcf-genotypes.js')).matrixAnswer
    : undefined;
  if (!ranges && !matrix) return failure('not-served', 'A VCF source serves ranges, and a matrix with samples');
  const got = rows ? rowsIn(rows, request.start, request.end) : [];
  if (got.length > options.maxRows) return failure('too-large', `${got.length} rows in one region`);
  const tables: { ranges?: RangesTable; matrix?: MatrixTable } = {};
  if (ranges) {
    const list = sitesList(rows, schema, got);
    tables.ranges = options.density && request.level > 1 && !request.items ? densityOf(list, request, true) : list;
  }
  if (matrix) {
    const m = matrixOf!(rows, got, options.genotypes!, request, options.matrix);
    if (!m.ok) return m;
    tables.matrix = m.tables.matrix;
  }
  return { ok: true, tables };
}

/** The sites list of `got` rows: a row an allele, its id, the fixed and INFO fields. */
function sitesList(rows: VcfRows | undefined, schema: VcfSchema, got: readonly number[]): RangesList {
  const names = [...FixedFields, ...schema.info];
  const fields: RangesField[] = names.map((f, k) => ({
    name: f.name,
    column: rows ? column(rows.values[k], f.type, got) : column([], f.type, []),
  }));
  return {
    form: 'list',
    start: rows ? Int32Array.from(got, (r) => rows.start[r]) : new Int32Array(0),
    end: rows ? Int32Array.from(got, (r) => rows.end[r]) : new Int32Array(0),
    id: categoryColumn(rows ? got.map((r) => rows.key[r]) : []),
    fields,
  };
}
