/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * https://www.ccp4.ac.uk/html/mtzformat.html
 */

/**
 * A single column descriptor from a COL header record.
 *
 * Column types defined by the MTZ format:
 *   H  Miller index h, k, or l
 *   J  intensity
 *   F  structure amplitude
 *   D  anomalous difference
 *   Q  standard deviation of J, F, D or other
 *   G  structure amplitude F(+) or F(-)
 *   L  standard deviation of G
 *   K  intensity I(+) or I(-)
 *   M  standard deviation of K
 *   E  normalised structure factor
 *   P  phase angle in degrees
 *   W  weight
 *   A  phase probability coefficient (Hendrickson-Lattman)
 *   B  batch number
 *   Y  M/ISYM
 *   I  any other integer
 *   R  any other real
 */
export interface MtzColumn {
    /** Column label, up to 30 characters. */
    label: string;
    /** Single-character column type (see above). */
    type: string;
    /** Minimum value stored in this column. */
    min: number;
    /** Maximum value stored in this column. */
    max: number;
    /** Dataset ID this column belongs to. */
    datasetId: number;
    /** 0-based column index within each reflection row. */
    index: number;
}

/** Per-dataset metadata from PROJECT/CRYSTAL/DATASET/DCELL/DWAVEL records. */
export interface MtzDataset {
    /** Integer dataset ID. 0 = HKL_base dataset. */
    id: number;
    /** Project name (from PROJECT record), up to 64 chars. */
    project: string;
    /** Crystal name (from CRYSTAL record), up to 64 chars. */
    crystal: string;
    /** Dataset name (from DATASET record), up to 64 chars. */
    name: string;
    /**
     * Cell dimensions [a, b, c, α, β, γ] in Å and degrees (from DCELL record).
     * Undefined if no DCELL record was present for this dataset.
     */
    cell?: [number, number, number, number, number, number];
    /** Wavelength in Å (from DWAVEL record). Undefined if not present. */
    wavelength?: number;
}

/** All header metadata parsed from the MTZ file. */
export interface MtzHeader {
    /** Version string from the VERS record. */
    version: string;
    /** File title from the TITLE record. */
    title: string;
    /** Number of columns (from NCOL record). */
    nColumns: number;
    /** Number of reflections (from NCOL record). */
    nReflections: number;
    /** Number of batches; > 0 indicates a multi-record (unmerged) file. */
    nBatches: number;
    /**
     * Global cell dimensions [a, b, c, α, β, γ] in Å and degrees (from CELL record).
     * Deprecated in favour of per-dataset DCELL, but used as fallback.
     */
    cell: [number, number, number, number, number, number];
    /** Space-group number from SYMINF record. */
    spaceGroupNumber: number;
    /** Space-group name from SYMINF record. */
    spaceGroupName: string;
    /** Total number of symmetry operations from SYMINF record. */
    nSymOps: number;
    /** Symmetry operations as x,y,z strings from SYMM records. */
    symops: string[];
    /** Column descriptors, one per data column, in order. */
    columns: MtzColumn[];
    /** Dataset descriptors collected from PROJECT/CRYSTAL/DATASET records. */
    datasets: MtzDataset[];
    /**
     * The value used in the data to represent a missing observation (from VALM record).
     * Typically NaN or a very large number.
     */
    missingNumber: number;
}

/**
 * Parsed MTZ file.
 *
 * `data` is laid out in row-major order as `Float32Array` of length
 * `header.nReflections * header.nColumns`.  Element at reflection `r`,
 * column `c` is `data[r * header.nColumns + c]`.
 */
export interface MtzFile {
    /** Source file name. */
    name: string;
    /** Parsed header metadata. */
    header: MtzHeader;
    /**
     * Reflection data as Float32Array, row-major:
     * index = reflectionIndex * nColumns + columnIndex
     */
    data: Float32Array;
}
