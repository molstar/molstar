/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { parseArgs } from 'node:util';
import type { Bucket, Row, Run } from './rows.ts';

// Baselines and the diff against them (plan §6): a baseline is the counts of a full run, committed, with its rows not
// passing for a dataset (the whole CCD's, thousands of IDs, are kept out of the repository: a local copy under
// `out/baseline/` holds them). Each run is compared with its own: by row when the rows are at hand, by counts otherwise.

/** `baseline/<dataset>/<connectivity>/<hydrogens>/<tier>.json` */
export interface Baseline {
  dataset: string;
  tier: string;
  connectivity: string;
  hydrogens: string;
  components: number;
  counts: Record<string, number>;
  /** every row not passing, with its bucket: committed for a dataset, local only for the whole CCD */
  rows?: Record<string, Bucket>;
}

export interface DiffRow {
  id: string;
  from: Bucket;
  to: Bucket;
  l0?: string[];
}

export interface Diff {
  baseline: string;
  /** the count per bucket, run minus baseline (non-zero ones only; a full run's only) */
  counts: Record<string, number>;
  /** the rows compared, from the committed baseline or the local copy; without them only the counts are */
  byRow: boolean;
  regressions: DiffRow[];
  improvements: DiffRow[];
  /** in or out of `excluded`: neither better nor worse */
  changes: DiffRow[];
}

export const BaselineDir = new URL('../baseline/', import.meta.url).pathname;

export function baselineFile(dataset: string, connectivity: string, hydrogens: string, tier: string) {
  return path.join(BaselineDir, dataset, connectivity, hydrogens, `${tier}.json`);
}

/** The local copy of a baseline with its rows, `<out>/baseline/<dataset>/<connectivity>/<hydrogens>/<tier>.json`. */
export function localBaselineFile(out: string, dataset: string, connectivity: string, hydrogens: string, tier: string) {
  return path.join(out, 'baseline', dataset, connectivity, hydrogens, `${tier}.json`);
}

const read = (file: string): Baseline => JSON.parse(fs.readFileSync(file, 'utf8'));
const sameCounts = (a: Record<string, number>, b: Record<string, number>) =>
  Object.keys({ ...a, ...b }).every((k) => (a[k] ?? 0) === (b[k] ?? 0));

/** worse is higher; `excluded` has no rank (plan §6) */
const Rank: Partial<Record<Bucket, number>> = { pass: 0, known: 1, fail: 2, unsanitizable: 2, error: 3 };

/**
 * The run against the committed baseline, if there is one: its counts, and its rows (all of them, passing included)
 * when the baseline's rows are at hand, committed or in the local copy. A local copy whose counts differ from the
 * committed ones is another baseline's (updated since) and is not used. `counts` is a full run's only (a partial run's
 * do not compare): a partial run without the baseline's rows has no diff.
 */
export function diffAgainstBaseline(
  file: string,
  local: string,
  counts: Record<string, number> | undefined,
  rows: Row[],
): Diff | undefined {
  if (!fs.existsSync(file)) return void 0;
  const baseline = read(file);
  let listed = baseline.rows;
  if (!listed && fs.existsSync(local)) {
    const copy = read(local);
    if (sameCounts(copy.counts, baseline.counts)) listed = copy.rows;
  }
  if (!listed && !counts) return void 0;
  const delta: Record<string, number> = {};
  for (const k of counts ? Object.keys({ ...counts, ...baseline.counts }) : []) {
    const d = (counts![k] ?? 0) - (baseline.counts[k] ?? 0);
    if (d) delta[k] = d;
  }
  const diff: Diff = {
    baseline: path.relative(process.cwd(), file),
    counts: delta,
    byRow: !!listed,
    regressions: [],
    improvements: [],
    changes: [],
  };
  if (!listed) return diff;
  for (const row of rows) {
    const from = listed[row.id] ?? 'pass';
    const to = row.bucket;
    if (from === to) continue;
    const d: DiffRow = { id: row.id, from, to, ...(row.l0 ? { l0: row.l0 } : {}) };
    const a = Rank[from];
    const b = Rank[to];
    if (a === void 0 || b === void 0) diff.changes.push(d);
    else if (b > a) diff.regressions.push(d);
    else if (b < a) diff.improvements.push(d);
  }
  return diff;
}

export function printDiff(diff: Diff, max = 30) {
  const line = (d: DiffRow) => `  ${d.id}: ${d.from} → ${d.to}${d.l0 ? ` (${d.l0.join(', ')})` : ''}`;
  const counts = Object.entries(diff.counts).map(([k, v]) => `${k} ${v > 0 ? '+' : ''}${v}`);
  console.log(`diff against ${diff.baseline}: counts ${counts.length ? counts.join(', ') : 'unchanged'}`);
  if (!diff.byRow) {
    console.log('  by counts only: no local copy of its rows (node src/baseline.ts on a full run writes one)');
    return;
  }
  console.log(
    `  ${diff.regressions.length} regressions, ${diff.improvements.length} improvements, ${diff.changes.length} in or out of excluded`,
  );
  for (const [title, list] of [
    ['regressions', diff.regressions],
    ['improvements', diff.improvements],
    ['changes', diff.changes],
  ] as const) {
    if (!list.length) continue;
    console.log(`${title}:`);
    for (const d of list.slice(0, max)) console.log(line(d));
    if (list.length > max) console.log(`  … ${list.length - max} more in the run's JSON`);
  }
}

/** The baseline of a full run (`<tier>.json` of `run`), with its rows; a `--filter` or `--sample` run is refused. */
export function baselineOf(dataset: string, run: Run): Baseline {
  if (run.selection !== 'all') {
    throw new Error(`tier ${run.tier}: a ${run.selection ?? 'partial'} run, not a full one: no baseline from it`);
  }
  const rows: Record<string, Bucket> = {};
  for (const r of [...run.rows].sort((a, b) => (a.id < b.id ? -1 : a.id > b.id ? 1 : 0))) {
    if (r.bucket !== 'pass') rows[r.id] = r.bucket;
  }
  const { tier, connectivity, hydrogens, components, counts } = run;
  return { dataset, tier, connectivity, hydrogens, components, counts, rows };
}

if (import.meta.main) {
  const { values } = parseArgs({
    options: {
      dataset: { type: 'string', default: 'ccd' },
      tier: { type: 'string' },
      hydrogens: { type: 'string', default: 'stripped' },
      out: { type: 'string', default: 'out' },
    },
  });
  const connectivity = 'ccd';
  const tiers = (values.tier ?? (values.dataset === 'ccd' ? '0,A,B,C' : 'C')).split(',');
  for (const tier of tiers) {
    const runFile = path.join(values.out, values.dataset, connectivity, values.hydrogens, `${tier}.json`);
    if (!fs.existsSync(runFile)) throw new Error(`${runFile}: no run`);
    const baseline = baselineOf(values.dataset, JSON.parse(fs.readFileSync(runFile, 'utf8')));
    const write = (file: string, b: Baseline) => {
      fs.mkdirSync(path.dirname(file), { recursive: true });
      fs.writeFileSync(file, `${JSON.stringify(b, null, 2)}\n`);
    };
    // the whole CCD's rows stay local; a dataset's are committed (the repository keeps it small, plan §6)
    const file = baselineFile(values.dataset, connectivity, values.hydrogens, tier);
    write(file, values.dataset === 'ccd' ? { ...baseline, rows: void 0 } : baseline);
    write(localBaselineFile(values.out, values.dataset, connectivity, values.hydrogens, tier), baseline);
    console.log(`${path.relative(process.cwd(), file)}: ${JSON.stringify(baseline.counts)}`);
  }
}
