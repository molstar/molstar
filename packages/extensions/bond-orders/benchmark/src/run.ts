/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { parseArgs } from 'node:util';
import { baselineFile, diffAgainstBaseline, localBaselineFile, printDiff } from './baseline.ts';
import { type LinkKind, representative } from './candidates.ts';
import { openCcd } from './ccd.ts';
import { comparableInchi, createIdentifiers, createTautomerInchi } from './compare.ts';
import { ccdInchi, checkCcd } from './exclusions.ts';
import { perceiveEntry } from './entry.ts';
import { fetchTierC } from './fetch.ts';
import { type Perceived, perceiveCcd } from './perceive.ts';
import { downloadEntries, readEntry } from './rcsb.ts';
import { DefaultViewer, partialDir, sdfFile, writeReport } from './report.ts';
import type { Row, Run } from './rows.ts';
import { type Hydrogens, type Tier, Tiers, withHydrogens } from './tiers.ts';

const { values } = parseArgs({
  options: {
    tier: { type: 'string' },
    dataset: { type: 'string', default: 'ccd' },
    sample: { type: 'string' },
    filter: { type: 'string' },
    cache: { type: 'string', default: 'cache' },
    out: { type: 'string', default: 'out' },
    hydrogens: { type: 'string', default: 'stripped' },
    viewer: { type: 'string', default: DefaultViewer },
  },
});

/** Cases listed by hand with their reason (plan §5.1): reported as `known`, not counted as failures. */
const manifest: { known: Record<string, string> } = JSON.parse(
  fs.readFileSync(new URL('../manifest.json', import.meta.url), 'utf8'),
);
const ccd = await openCcd(values.cache);
const identifiers = await createIdentifiers();
const tautomerInchi = await createTautomerInchi();
// the only connectivity source so far (plan §3.3)
const connectivity = 'ccd';
const hydrogens = values.hydrogens as Hydrogens;
if (hydrogens !== 'stripped' && hydrogens !== 'kept') throw new Error(`--hydrogens ${hydrogens}: stripped or kept`);
const dataset = values.dataset === 'ccd' ? void 0 : readDataset(values.dataset);
const selection = values.filter ? 'filter' : values.sample ? 'sample' : 'all';
// a `--filter` or `--sample` run in a folder of its own, beside the full run's, which it never overwrites (plan §6)
const runDir = path.join(values.out, values.dataset, connectivity, partialDir(hydrogens, selection));
const tierNames = (values.tier ?? (dataset ? 'C' : '0,A,B')).split(',');
if (dataset && tierNames.some((t) => t !== 'C')) throw new Error(`--dataset: tier C's pipeline only (plan §6)`);

/** What a row perceives: a component alone, in its representative entry (tier C), or a dataset's pair (plan §6). */
interface Item {
  id: string;
  comp: string;
  entry?: string;
  kind?: LinkKind;
  /** tier C: a linked copy only; a dataset: the first complete copy */
  linkedOnly: boolean;
  /** a dataset's pair judged at another tier (`"tier": "B"`) */
  judgedAt?: Tier;
}

for (const name of tierNames) {
  if (!Tiers[name]) throw new Error(`unknown tier ${name} (${Object.keys(Tiers).join(', ')})`);
  const tier = withHydrogens(Tiers[name], hydrogens);
  const dir = path.join(runDir, tier.name);
  // a run's mol blocks only: none left from an earlier run
  fs.rmSync(path.dirname(sdfFile(dir, '')), { recursive: true, force: true });
  fs.mkdirSync(path.dirname(sdfFile(dir, '')), { recursive: true });

  // tier C: the linked components only, the candidates of plan §3.2 (fetched into the cache when missing); a dataset:
  // its pairs, in their own entries
  const t0 = performance.now();
  const items = await itemsOf(tier);

  // perceive: no RDKit; the perceived and reference SDFs of each component, its L0 diff and diagnostics
  const rows: Row[] = [];
  const toCompare: Row[] = [];
  const t1p = performance.now();
  for (const item of items) {
    const { id, comp, judgedAt } = item;
    const block = await ccd.read(comp);
    const check = checkCcd(block, judgedAt ?? tier);
    const row: Row = { id, bucket: 'pass', kind: item.kind };
    if (comp !== id) row.comp = comp;
    if (judgedAt) row.judgedAt = judgedAt.name;
    rows.push(row);
    if (check.excluded) {
      row.bucket = 'excluded';
      row.reason = check.excluded;
      continue;
    }
    if (check.leavingMissing) row.leavingMissing = check.leavingMissing;
    try {
      let result: Perceived;
      if (judgedAt) {
        // the geometry measured is the entry's only when the CCD's model coordinates come from it (plan §6)
        const source = representative(block).entry;
        row.entry = item.entry;
        if (source !== item.entry) {
          throw new Error(`the CCD's model coordinates are from ${source ?? 'no entry'}, not ${item.entry}`);
        }
        result = await perceiveCcd(comp, block, judgedAt);
      } else if (item.entry) {
        const data = readEntry(values.cache, item.entry);
        const entry = data
          ? await perceiveEntry(comp, block, data, tier, item.linkedOnly)
          : { excluded: 'entry not served', linked: false };
        // linked: a cap per link (tier C's copies always are)
        row.linked = 'excluded' in entry ? entry.linked : entry.caps > 0;
        if ('excluded' in entry) {
          row.bucket = 'excluded';
          row.reason = entry.excluded;
          continue;
        }
        row.entry = item.entry;
        row.copy = entry.copy;
        row.caps = entry.caps;
        result = entry;
      } else {
        result = await perceiveCcd(comp, block, tier);
      }
      fs.writeFileSync(sdfFile(dir, id), result.perceived);
      if (result.reference) fs.writeFileSync(sdfFile(dir, id, true), result.reference);
      if (result.l0.bonds.length) {
        row.l0 = result.l0.bonds;
        row.l0Atoms = result.l0.atoms;
      }
      if (result.diagnostics.length) row.diagnostics = result.diagnostics;
      toCompare.push(row);
    } catch (e) {
      row.bucket = 'error';
      row.reason = String(e);
    }
  }
  const t1 = performance.now();

  // compare: tier 0 against the CCD's InChI (the pipeline check); the perception tiers against RDKit's identifiers of
  // the CCD's orders on the same atoms (plan §4.2)
  for (const row of toCompare) {
    const read = (reference: boolean) => fs.readFileSync(sdfFile(dir, row.id, reference), 'utf8');
    // the reference first: the CCD's own graph may not be a structure RDKit accepts (PF6 with a neutral P, an H of the
    // model coordinates within bonding distance of a second atom): nothing to compare with
    const reference = tier.perceive ? identifiers(read(true)) : void 0;
    if (reference && 'error' in reference) {
      row.bucket = 'excluded';
      row.reason = `reference unsanitizable: ${reference.error}`;
      continue;
    }
    const perceived = identifiers(read(false));
    if ('error' in perceived) {
      row.bucket = 'unsanitizable';
      row.reason = perceived.error;
      continue;
    }
    let expectedInchi: string;
    if (reference) {
      if (!row.l0) row.level = 'L0';
      else if (perceived.smiles === reference.smiles) row.level = 'L1';
      if (row.level) continue;
      expectedInchi = reference.inchi;
    } else {
      expectedInchi = ccdInchi(await ccd.read(row.comp ?? row.id))!;
    }
    // charge and protonation are compared in tier 0 only (the pipeline check), never on perceived orders
    const expected = comparableInchi(expectedInchi, !tier.perceive);
    const actual = comparableInchi(perceived.inchi, !tier.perceive);
    if (actual === expected) {
      if (tier.perceive) row.level = 'L2';
    } else if (
      reference &&
      comparableInchi(tautomerInchi(read(false)), false) === comparableInchi(tautomerInchi(read(true)), false)
    ) {
      // keto / enol or a 1,5 shift: a tautomer of the reference, not an order failure (plan §5)
      row.level = 'L3';
    } else {
      row.bucket = 'fail';
      row.expected = expected;
      row.perceived = actual;
    }
  }
  for (const row of toCompare) {
    const reason = manifest.known[row.id];
    if (reason && (row.bucket === 'fail' || row.bucket === 'unsanitizable')) {
      row.bucket = 'known';
      row.reason = reason;
    }
  }
  const t2 = performance.now();

  const run: Run & Record<string, unknown> = {
    dataset: values.dataset,
    selection,
    tier: tier.name,
    connectivity,
    hydrogens,
    components: items.length,
    ...(tier.entry ? { byKind: byKind(rows) } : {}),
    // tier C: the ligands with a linked copy, whatever their bucket (the reports' covalent-ligand flag, plan §4.3)
    ...(tier.entry ? { covalentLigands: rows.filter((r) => r.kind === 'ligand' && r.linked).map((r) => r.id) } : {}),
    counts: tally(rows, (r) => r.bucket),
    levels: tally(
      rows.filter((r) => r.level),
      (r) => r.level!,
    ),
    // a residue without its leaving atoms (tier B) apart: the comparison then covers fewer atoms (plan §5)
    leavingMissing: tally(
      rows.filter((r) => r.leavingMissing),
      (r) => r.bucket,
    ),
    withDiagnostics: tally(
      rows.filter((r) => r.diagnostics),
      (r) => r.bucket,
    ),
    seconds: { fetch: round((t1p - t0) / 1000), perceive: round((t1 - t1p) / 1000), compare: round((t2 - t1) / 1000) },
    maxRssMB: Math.round(process.resourceUsage().maxRSS / 1024),
    // every component not passing, and the passing ones with a diagnostic
    rows: rows.filter((r) => r.bucket !== 'pass' || r.diagnostics),
  };
  // against the baseline, on this run's rows (plan §6)
  const diff = diffAgainstBaseline(
    baselineFile(values.dataset, connectivity, hydrogens, tier.name),
    localBaselineFile(values.out, values.dataset, connectivity, hydrogens, tier.name),
    selection === 'all' ? run.counts : void 0,
    rows,
  );
  if (diff) run.diff = diff;
  fs.writeFileSync(path.join(runDir, `${tier.name}.json`), JSON.stringify(run, null, 2));
  console.log(
    `tier ${tier.name}:`,
    JSON.stringify({ ...run, rows: void 0, diff: void 0, covalentLigands: run.covalentLigands?.length }),
  );
  if (diff) printDiff(diff);
  const report = await writeReport(runDir, tier.name, values.viewer);
  console.log(`report: ${path.join(dir, 'index.html')} (${report.rows} rows, ${report.drawn} drawn)`);
}

/** The tier's items: the CCD's components (tier C: its candidates), or the dataset's pairs, as selected. */
async function itemsOf(tier: Tier): Promise<Item[]> {
  if (dataset) {
    const pairs = select(dataset.pairs, (p) => [`${p.entry}_${p.comp}`, p.entry, p.comp]);
    await downloadEntries(
      pairs.filter((p) => !p.tier).map((p) => p.entry),
      values.cache,
    );
    const items: Item[] = [];
    for (const p of pairs) {
      const judgedAt = p.tier ? withHydrogens(Tiers[p.tier], hydrogens) : void 0;
      if (p.tier && !judgedAt) throw new Error(`${p.entry} ${p.comp}: unknown tier ${p.tier}`);
      const kind = representative(await ccd.read(p.comp)).kind;
      items.push({ id: `${p.entry}_${p.comp}`, comp: p.comp, entry: p.entry, kind, linkedOnly: false, judgedAt });
    }
    return items;
  }
  const ids = select(ccd.ids, (id) => [id]);
  if (!tier.entry) return ids.map((id) => ({ id, comp: id, linkedOnly: true }));
  const { reps, candidates } = await fetchTierC(ids, (id) => ccd.read(id), values.cache);
  return candidates.map((id) => ({
    id,
    comp: id,
    entry: reps.get(id)!.entry,
    kind: reps.get(id)!.kind,
    linkedOnly: true,
  }));
}

/** `--filter` (any of an item's names), or `--sample N` spread evenly, or all. */
function select<T>(all: T[], names: (item: T) => string[]): T[] {
  if (values.filter) {
    const wanted = new Set(values.filter.split(',').map((f) => f.trim().toUpperCase()));
    return all.filter((item) => names(item).some((n) => wanted.has(n.toUpperCase())));
  }
  if (values.sample) {
    const n = Number(values.sample);
    const step = all.length / n;
    return Array.from({ length: Math.min(n, all.length) }, (_, i) => all[Math.floor(i * step)]);
  }
  return all;
}

/** `datasets/<name>.json` (plan §6) */
interface Dataset {
  name: string;
  source: string;
  pairs: { entry: string; comp: string; tier?: string; reason?: string }[];
}

function readDataset(name: string): Dataset {
  const file = new URL(`../datasets/${name}.json`, import.meta.url);
  if (!fs.existsSync(file)) throw new Error(`no dataset ${name} (datasets/${name}.json)`);
  return JSON.parse(fs.readFileSync(file, 'utf8'));
}

/** Tier C: the counts per link kind (plan §5: polymer residues, glycans, covalent ligands). */
function byKind(rows: Row[]) {
  const kinds: Record<string, Record<string, number>> = {};
  for (const kind of new Set(rows.map((r) => r.kind ?? 'ligand'))) {
    kinds[kind] = tally(
      rows.filter((r) => (r.kind ?? 'ligand') === kind),
      (r) => r.bucket,
    );
  }
  return kinds;
}

function tally<T>(items: T[], key: (item: T) => string) {
  const counts: Record<string, number> = {};
  for (const item of items) counts[key(item)] = (counts[key(item)] ?? 0) + 1;
  return counts;
}

function round(v: number) {
  return Math.round(v * 10) / 10;
}
