/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { cellsSource, delayedSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime, readStockholm } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';

// The Globins dataset: the Pfam PF00042 seed as published (Stockholm, cut to 60 rows by the
// record script), its tree computed by neighbor joining (`tree.ts`) and written inline, UniProt's organisms of its rows
// (a columns file joined by row key), PDB 4HHB as `seq` (`hemoglobin.ts`); and a scale variant of any size whose rows
// are made on demand by index from the seed's per-column profile, lowercase inserts with an `RF`-like column every 25
// (EC-38, EC-41), never held all at once by the source.

/** The seed, as InterPro serves it, cut to 60 rows. */
export const SeedFile = 'fixtures/globins/PF00042-seed.sto';
/** UniProt's organism and class of each row, keyed by the row's label. */
export const ColumnsFile = 'fixtures/globins/family-columns.tsv';
/** 4HHB's chains, helices and heme as `seq`. */
export const HemoglobinFile = 'fixtures/globins/4hhb.seq';

/**
 * The reference row: the seed's row nearest 4HHB's α chain by identity (HBA_HUMAN is not a row of the current seed; a
 * test holds that this one is the nearest). Wild yak's α, 87% identical to human α over its domain.
 */
export const ReferenceRow = 'HBA1_BOSGR/26-136';
/** The row nearest 4HHB's β chain (mouse's embryonic β). */
export const NearestBetaRow = 'HBB0_MOUSE/26-142';

/** The seed's tree, neighbor joining on identity (`tree.ts`), as committed; a test compares it with its generator. */
export const SeedTree =
  '(((GLB1_GLYDI/28-141:0.3075,GLBP1_GLYDI/29-141:0.3222):0.0875,(((LGB1_MEDSA/27-142:0.1748,LGB1_LUPLU/28-146:0.1730):0.0913,(HBP2_CASGL/33-152:0.2232,HBP1_CASGL/27-146:0.2434):0.0399):0.1147,(A0ABY4E6T8_VITST/26-133:0.2608,HMP_ECOLI/26-131:0.2769):0.1225):0.0094):0.0146,(GLB_APLJU/26-139:0.4044,((GLB2_PHRSO/26-129:0.2943,GLB1_PHRSO/26-136:0.3077):0.1214,(HBF1_URECA/28-141:0.3931,(GLBT_CHITH/29-142:0.3491,((GLB73_CHITH/39-157:0.2078,(GLB8_CHITH/25-143:0.2022,GLB2_CHITH/38-156:0.2012):0.0275):0.0361,(GLB1_CHITH/37-154:0.2447,GLB3_CHITH/36-147:0.2689):0.0139):0.0799):0.0407):0.0105):0.0043):0.0025,((GLB2_LUMTE/31-141:0.3302,(GLB4_LUMTE/36-146:0.2658,GLB3_TYLHE/33-143:0.2568):0.0586):0.0759,((GLB1B_ANATR/37-152:0.3872,(GLB_CERRH/29-142:0.2339,GLB_BUSCA/29-147:0.2397):0.1340):0.0188,(((MYG_CYPCA/23-136:0.2783,(MYG_ALLMI/27-143:0.3007,(MYG_GALGA/23-138:0.2391,MYG_HETPO/23-138:0.2350):0.0700):0.0012):0.0825,(HBAM_AQUCT/26-129:0.3375,(((HBA_PANCK/26-137:0.2671,(((HBA1_IGUIG/26-136:0.2201,(HBAT_HORSE/27-137:0.2014,HBA1_BOSGR/26-136:0.1499):0.0412):0.0168,(HBA1_XENBO/27-137:0.2252,HBA1_PLEWA/27-137:0.2523):0.0125):0.0096,((HBAD_ERYML/26-136:0.2219,HBAD_PASMO/26-136:0.2016):0.0303,(HBA3_PLEWA/27-137:0.1972,HBAZ_CAPHI/27-137:0.1992):0.0282):0.0067):0.0071):0.0286,(HBA_LEPPA/26-138:0.2710,(HBA_HETPO/33-143:0.2297,(HBA1_TORMA/26-136:0.2027,HBA_SQUAC/26-136:0.1937):0.0226):0.0398):0.0149):0.0263,(((HBB_AQUCT/19-135:0.2160,HBB1_XENBO/25-141:0.2028):0.0594,((HBB1_CYGMA/26-142:0.2348,HBB_LEPPA/25-142:0.2524):0.0101,(HBB_ALLMI/25-141:0.2391,(HBB0_MOUSE/26-142:0.2100,HBBN_AMMLE/20-136:0.2430):0.0216):0.0070):0.0089):0.0314,(HBB2_XENLA/26-142:0.2481,(HBB_HETPO/25-136:0.2330,HBB_SQUAC/25-137:0.1994):0.0469):0.0222):0.0238):0.0442):0.0283):0.0193,(GLB_PAREP/30-144:0.4174,((GLB_ASCSU/45-162:0.1374,GLB_PSEDC/45-162:0.1423):0.1589,(GLBH_CAEEL/32-147:0.2347,(GLBC_NIPBR/44-161:0.1917,(GLBH_TRICO/45-161:0.1498,A0A3P7CCD7_NIPBR/27-142:0.1346):0.0521):0.0313):0.0699):0.1249):0.0142):0.0055):0.0058):0.0005);';

/** The family's runtimes: its Stockholm file and its joined columns file. */
export function familyRegistry(): readonly SeqRegistryEntry[] {
  return [
    { kind: 'format-runtime', runtime: AlignmentFileRuntime },
    { kind: 'format-runtime', runtime: RecordsRuntime },
  ];
}

const Residues = 'ACDEFGHIKLMNPQRSTVWY';

/** Each column's residue frequencies (a cumulative distribution) and gap rate over an alignment's rows. */
export interface Profile {
  readonly cdf: Float64Array;
  readonly gap: Float64Array;
  readonly columns: number;
}

/** The profile of aligned rows (residues counted case-blind; a letter not of the twenty counted as a gap). */
export function profileOf(letters: readonly string[]): Profile {
  const C = letters[0].length;
  const K = Residues.length;
  const cdf = new Float64Array(C * K);
  const gap = new Float64Array(C);
  for (let j = 0; j < C; j++) {
    let n = 0;
    for (const r of letters) {
      const k = Residues.indexOf(r[j].toUpperCase());
      if (k < 0) gap[j]++;
      else {
        cdf[j * K + k]++;
        n++;
      }
    }
    gap[j] /= letters.length;
    let acc = 0;
    for (let k = 0; k < K; k++) {
      acc += n ? cdf[j * K + k] / n : 1 / K;
      cdf[j * K + k] = acc;
    }
  }
  return { cdf, gap, columns: C };
}

/** The seed's profile from its Stockholm text. */
export const seedProfile = (stockholm: string): Profile =>
  profileOf(readStockholm(stockholm).rows.map((r) => r.letters));

/**
 * The scale variant's row `i` of `columns` (300 by default): the seed's columns first, then columns sampled from its
 * profile in turn; every 25th column an insert (lowercase, gaps mostly); a function of the seed and the index alone.
 */
export function generatedRow(i: number, columns: number, seed: number, profile: Profile, out?: Uint8Array) {
  const K = Residues.length;
  const row = out ?? new Uint8Array(columns);
  const rng = new PCG((seed * 0x9e3779b1 + i * 2654435761) >>> 0);
  for (let j = 0; j < columns; j++) {
    const pj = j % profile.columns;
    const insert = j % 25 === 24;
    if (rng.float() < (insert ? 0.8 : profile.gap[pj])) {
      row[j] = insert ? 46 : 45; // `.` in an insert column, `-` elsewhere
      continue;
    }
    const u = rng.float();
    let k = 0;
    while (k < K - 1 && profile.cdf[pj * K + k] < u) k++;
    const c = Residues.charCodeAt(k);
    row[j] = insert ? c + 32 : c;
  }
  return row;
}

/** A generated row's start number: a function of its index (EC-39). */
export const generatedStart = (i: number) => 1 + (i % 7);

/**
 * A row's label in the scale variant: `GEN_<n>/<start>-<end>` (EC-39), its end the number of its last residue (the
 * row's letters counted, a gap never), so the label's range is the row's own.
 */
export function generatedLabel(i: number, row: Uint8Array): string {
  const start = generatedStart(i);
  let n = 0;
  for (let j = 0; j < row.length; j++) if (row[j] !== 45 && row[j] !== 46) n++;
  return `GEN_${i + 1}/${start}-${start + Math.max(n, 1) - 1}`;
}

/**
 * A balanced tree over `n` rows of the scale variant (for collapsed clades at scale), leaves in data order, each named
 * `GEN_<n>`: a leaf matches its row by the label before its last `/` (EC-43), so the tree needs no row's letters.
 */
export function balancedNewick(n: number): string {
  const parts: string[] = [];
  const grow = (a: number, b: number): void => {
    if (b - a === 1) {
      parts.push(`GEN_${a + 1}:1`);
      return;
    }
    const m = (a + b) >> 1;
    parts.push('(');
    grow(a, m);
    parts.push(',');
    grow(m, b);
    parts.push(`)c${a}_${b}:1`);
  };
  grow(0, n);
  return `${parts.join('')};`;
}

/** The scale variant's format: `value: { rows, columns, seed, profile }`, rows made on demand by index (EC-41). */
export const GeneratedFormat = 'globins-generated';

/** Its runtime: the profile read from `profile` (the seed's file by default); each answer after `delayMs` (0: none). */
export function generatedRuntime(delayMs = 0): FormatRuntime {
  return {
    name: GeneratedFormat,
    formats: [GeneratedFormat],
    async open(source, ctx) {
      const v = (source.value ?? {}) as { rows?: number; columns?: number; seed?: number; profile?: string };
      const rows = v.rows ?? 10_000;
      const columns = v.columns ?? 300;
      const seed = v.seed ?? 1;
      const profile = seedProfile(await ctx.assets.readText(v.profile ?? SeedFile, ctx.signal));
      const inner = cellsSource(rows, columns, (start, end) => {
        const codes = new Uint8Array((end - start) * columns);
        const label: string[] = [];
        const first = new Int32Array(end - start);
        for (let i = start; i < end; i++) {
          const row = codes.subarray((i - start) * columns, (i - start + 1) * columns);
          generatedRow(i, columns, seed, profile, row);
          label.push(generatedLabel(i, row));
          first[i - start] = generatedStart(i);
        }
        return { codes, label, first };
      });
      return delayMs > 0 ? delayedSource(inner, { delayMs }) : inner;
    },
  };
}

/** What a view of the scale variant needs: its runtime. */
export function generatedRegistry(options: { readonly delayMs?: number } = {}): readonly SeqRegistryEntry[] {
  return [{ kind: 'format-runtime', runtime: generatedRuntime(options.delayMs ?? 0) }];
}
