/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The recorded cases: the files, the readers over them and the requests the tests (and the recording script, which
// fetches exactly what these read) make. Coordinates are 0-based, half-open, as the stream port's. Recorded 2026-10-12.

import type { StreamRequest } from '@molstar/seq-core/format';
import {
  bigBedProvider,
  bigWigProvider,
  type RangeReader,
  twoBitProvider,
  vcfProvider,
} from '../../../src/stream/gmod/index.js';

export const req = (system: string | null, start = 0, stop = 0, level = -1): StreamRequest => ({
  source: 's',
  system,
  start,
  stop,
  row_start: 0,
  row_stop: 0,
  level,
  labels: false,
  priority: 'view',
});

export const URLS = {
  twobit: 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit',
  bigwig: 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw',
  bigbed: 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV50.bb',
  vcf: 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz',
} as const;
export type Key = keyof typeof URLS | 'tbi';

/** hg38 chromosome sizes the VCF's manifest takes from the host (a ClinVar header has no `##contig` lengths). */
export const LENGTHS = { 1: 248_956_422, MT: 16_569 };

export const makers = {
  twobit: (read: RangeReader) => twoBitProvider(read, { only: (n) => n === 'chr1' || n === 'chrM' }),
  bigwig: (read: RangeReader) => bigWigProvider(read),
  bigbed: (read: RangeReader, maxItems?: number) => bigBedProvider(read, { maxItems }),
  vcf: (read: RangeReader, index: Uint8Array, maxItems?: number) =>
    vcfProvider(read, index, { info: ['CLNSIG'], lengths: LENGTHS, maxItems }),
};

/** [name, request] per source, in the order they are made: the manifest first. */
export const CASES: Record<keyof typeof URLS, [string, StreamRequest][]> = {
  twobit: [
    ['manifest', req(null)],
    ['chr1 start', req('chr1', 9_990, 10_030)], // the last N's of the 10,000-N telomere gap, then the first soft-masked bases
    ['chr1 mid', req('chr1', 100_000_000, 100_004_096)], // deep in the mask list: found by binary search
    ['chrM', req('chrM', 0, 4_096)],
    ['chrM end', req('chrM', 16_000, 17_000)], // past the end of the sequence: cut at 16,569
  ],
  bigwig: [
    ['manifest', req(null)],
    ['chr21 data', req('chr21', 10_000_000, 10_004_096)],
    ['chr1 data', req('chr1', 100_000_000, 100_004_096)],
    ['chrM data', req('chrM', 0, 4_096)],
    ['chr21 level 0', req('chr21', 9_984_000, 10_009_600, 0)], // 256 bins of 100: 390 * 25,600 = 9,984,000
    ['chr1 level 1', req('chr1', 99_942_400, 99_942_400 + 102_400, 1)], // bins of 400
    ['chr1 level 3', req('chr1', 99_942_400, 99_942_400 + 1_638_400, 3)], // bins of 6,400
  ],
  bigbed: [
    ['manifest', req(null)],
    ['chr1 genes', req('chr1', 10_000, 75_536)], // DDX11L1, WASH7P and the first transcripts of the chromosome
    ['chr21 genes', req('chr21', 10_000_000, 10_065_536)],
    ['chr21 level 0', req('chr21', 9_961_472, 10_027_008, 0)], // 256 bins of 256
  ],
  vcf: [
    ['manifest', req(null)],
    ['chr1 variants', req('1', 925_000, 950_000)], // SAMD11
    ['MT variants', req('MT', 5_000, 6_000)],
  ],
};
