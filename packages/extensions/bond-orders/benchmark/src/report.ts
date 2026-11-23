/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { parseArgs } from 'node:util';
import { loadRDKit, type RDKit } from './compare.ts';
import type { Run } from './rows.ts';

/** The example's dev server (`pnpm dev:examples -- bond-providers`, plan §7) */
export const DefaultViewer = 'http://localhost:1338/examples/bond-providers/build/';

const Size = { width: 320, height: 240 };
/** the bonds that differ at L0, and their atoms */
const Highlight: [number, number, number, number] = [1, 0.55, 0.55, 1];

/**
 * The HTML report of a tier (plan §4.3): `<tier>/index.html`, static, opened from disk. One row per component not
 * passing, and per passing one with a diagnostic; the expected and perceived 2D drawings side by side
 * (`<tier>/svg/`, loaded lazily: inlined, a tier's few thousand rows would make a page of hundreds of MB).
 */
/** A run's folder below its connectivity: `stripped` for a full run, `stripped-filter` for a `--filter` one (plan §6). */
export function partialDir(hydrogens: string, selection: string) {
  return selection === 'all' ? hydrogens : `${hydrogens}-${selection}`;
}

export async function writeReport(runDir: string, tier: string, viewer = DefaultViewer) {
  const run: Run = JSON.parse(fs.readFileSync(path.join(runDir, `${tier}.json`), 'utf8'));
  const dir = path.join(runDir, tier);
  const svgDir = path.join(dir, 'svg');
  fs.rmSync(svgDir, { recursive: true, force: true });
  fs.mkdirSync(svgDir, { recursive: true });
  const rdkit = await loadRDKit();

  const shown = run.rows.filter((r) => r.bucket !== 'excluded');
  const drawn = new Set<string>();
  for (const row of shown) {
    const perceivedFile = path.join(dir, `${row.id}.sdf`);
    if (!fs.existsSync(perceivedFile)) continue;
    const referenceFile = path.join(dir, `${row.id}.ref.sdf`);
    const reference = fs.existsSync(referenceFile) ? fs.readFileSync(referenceFile, 'utf8') : void 0;
    const svgs = drawPair(rdkit, reference, fs.readFileSync(perceivedFile, 'utf8'), row.l0Atoms ?? []);
    if (!svgs) continue;
    if (svgs.reference) fs.writeFileSync(path.join(svgDir, `${row.id}.ref.svg`), svgs.reference);
    fs.writeFileSync(path.join(svgDir, `${row.id}.svg`), svgs.perceived);
    drawn.add(row.id);
  }

  const covalent = covalentLigands(runDir, run);
  const isCovalent = (id: string) => covalent?.has(id) || void 0;
  // the baseline's bucket of each regression (plan §6), marked on its row
  const regressions = new Map(run.diff?.regressions.map((d) => [d.id, d.from]));
  const data = {
    run: { ...run, rows: void 0, covalentLigands: void 0, diff: void 0 },
    viewer,
    covalent: covalent && covalentCounts(run, covalent),
    diff: run.diff && {
      ...Object.fromEntries(Object.entries(run.diff.counts).map(([k, v]) => [k, v > 0 ? `+${v}` : `${v}`])),
      // by row only when the baseline's rows were at hand (a local copy for the whole CCD, plan §6)
      ...(run.diff.byRow
        ? {
            regressions: run.diff.regressions.length,
            improvements: run.diff.improvements.length,
            changes: run.diff.changes.length,
          }
        : { 'by row': 'no local rows' }),
    },
    rows: shown.map((r) => ({
      ...r,
      l0Atoms: void 0,
      drawn: drawn.has(r.id) || void 0,
      covalent: isCovalent(r.id),
      regression: regressions.get(r.id),
    })),
    excluded: run.rows
      .filter((r) => r.bucket === 'excluded')
      .map((r) => ({ id: r.id, kind: r.kind, reason: r.reason, covalent: isCovalent(r.id) })),
  };
  fs.writeFileSync(path.join(dir, 'index.html'), page(run, data));
  return { rows: shown.length, drawn: drawn.size };
}

/**
 * The covalent ligands: the non-polymer components linked in their entry, which only tier C knows (plan §3.2). Their
 * CCD definition is the unbound molecule, often not the bound one. From this run, or from tier C's of the same
 * connectivity (linking does not depend on H); `undefined` before any tier C run.
 */
function covalentLigands(runDir: string, run: Run): Set<string> | undefined {
  if (run.covalentLigands) return new Set(run.covalentLigands);
  const connectivityDir = path.dirname(runDir);
  const dirs = [runDir, ...fs.readdirSync(connectivityDir).map((d) => path.join(connectivityDir, d))];
  for (const dir of dirs) {
    const file = path.join(dir, 'C.json');
    if (!fs.existsSync(file)) continue;
    const c: Run = JSON.parse(fs.readFileSync(file, 'utf8'));
    if (c.covalentLigands) return new Set(c.covalentLigands);
  }
  return void 0;
}

/** The covalent ligands' buckets in this run: the rows hold every one not passing, the rest pass. */
function covalentCounts(run: Run, covalent: Set<string>) {
  const counts: Record<string, number> = {};
  let notPassing = 0;
  for (const r of run.rows) {
    if (!covalent.has(r.id) || r.bucket === 'pass') continue;
    counts[r.bucket] = (counts[r.bucket] ?? 0) + 1;
    notPassing++;
  }
  return { total: covalent.size, counts: { pass: covalent.size - notPassing, ...counts } };
}

/**
 * The reference with fresh 2D coordinates, and the perceived side on the same coordinates: both mol blocks are the
 * writer's on the same atoms in the same order (plan §4.1), so the drawings differ by their orders only. The perceived
 * side is parsed without sanitizing (an unsanitizable one is drawn too). Tier 0 has no reference: the perceived side
 * gets its own coordinates.
 */
function drawPair(rdkit: RDKit, reference: string | undefined, perceived: string, l0Atoms: [number, number][]) {
  const options = (block: string) => {
    const bonds = bondIndices(block, l0Atoms);
    const atoms = [...new Set(l0Atoms.flat())];
    return { ...Size, atoms, bonds, highlightColour: Highlight, kekulize: false };
  };
  // the file's own orders on both sides: sanitized, RDKit would draw its aromatic form
  const draw = (block: string, coordinates?: string) => {
    const mol = rdkit.getMol(coordinates ? withCoordinates(block, coordinates) : block, {
      removeHs: false,
      sanitize: false,
      throwOnError: true,
    });
    if (!coordinates) mol.setNewCoords(true);
    const svg = mol.getSVG(options(block));
    mol.delete();
    return svg;
  };
  try {
    if (!reference) return { perceived: draw(perceived) };
    const mol = rdkit.getMol(reference, { removeHs: false, throwOnError: true });
    mol.setNewCoords(true);
    const coordinates = mol.getMolblock();
    mol.delete();
    return { reference: draw(reference, coordinates), perceived: draw(perceived, coordinates) };
  } catch {
    return void 0;
  }
}

/** The V2000 counts line and atom block: the first 3 characters are the atom count, coordinates the first 30 of a line. */
function withCoordinates(target: string, source: string) {
  const lines = target.split('\n');
  const from = source.split('\n');
  const n = Number(lines[3].slice(0, 3));
  for (let i = 0; i < n; i++) lines[4 + i] = from[4 + i].slice(0, 30) + lines[4 + i].slice(30);
  return lines.join('\n');
}

/** The mol block's bond indices (0-based, in file order, as RDKit reads them) of the given atom pairs. */
function bondIndices(block: string, pairs: [number, number][]) {
  if (!pairs.length) return [];
  const lines = block.split('\n');
  const atoms = Number(lines[3].slice(0, 3));
  const count = Number(lines[3].slice(3, 6));
  const wanted = new Set(pairs.map(([a, b]) => (a < b ? `${a}-${b}` : `${b}-${a}`)));
  const found: number[] = [];
  for (let t = 0; t < count; t++) {
    const line = lines[4 + atoms + t];
    const a = Number(line.slice(0, 3)) - 1;
    const b = Number(line.slice(3, 6)) - 1;
    if (wanted.has(a < b ? `${a}-${b}` : `${b}-${a}`)) found.push(t);
  }
  return found;
}

function page(run: Run, data: unknown) {
  const of = run.dataset && run.dataset !== 'ccd' ? `dataset ${run.dataset}` : `tier ${run.tier}`;
  const title = `Bond orders, ${of} (${run.connectivity}, H ${run.hydrogens}${run.selection && run.selection !== 'all' ? `, ${run.selection}` : ''})`;
  const json = JSON.stringify(data).replace(/</g, '\\u003c');
  return `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${title}</title>
<style>
:root { --bg: #fff; --fg: #1d1d1f; --muted: #6e6e73; --line: #d9d9de; --card: #f5f5f7; --fail: #c62828; --known: #8a6d00; --pass: #2e7d32; --flag: #5b3fa8; --flag-bg: #ece6fa; }
@media (prefers-color-scheme: dark) { :root { --bg: #161618; --fg: #ececf0; --muted: #a0a0a8; --line: #34343a; --card: #1f1f23; --fail: #ef7b7b; --known: #e0c050; --pass: #7cc47f; --flag: #cbb8ff; --flag-bg: #2f2645; } }
body { margin: 0; padding: 16px 24px; background: var(--bg); color: var(--fg); font: 14px/1.4 system-ui, sans-serif; }
h1 { font-size: 20px; margin: 0 0 8px; }
.summary { display: flex; flex-wrap: wrap; gap: 8px 24px; color: var(--muted); margin-bottom: 12px; }
.summary b { color: var(--fg); }
.controls { position: sticky; top: 0; background: var(--bg); padding: 8px 0; display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: center; border-bottom: 1px solid var(--line); z-index: 1; }
.controls label { color: var(--muted); }
select, input { font: inherit; background: var(--card); color: var(--fg); border: 1px solid var(--line); border-radius: 4px; padding: 2px 6px; }
table { border-collapse: collapse; width: 100%; }
td, th { border-bottom: 1px solid var(--line); padding: 8px; vertical-align: top; text-align: left; }
th { color: var(--muted); font-weight: 500; }
td.id { white-space: nowrap; font-weight: 600; }
img, .nodraw { width: ${Size.width}px; height: ${Size.height}px; background: #fff; border-radius: 4px; display: block; }
.nodraw { display: flex; align-items: center; justify-content: center; color: #888; background: var(--card); }
.bucket-fail, .bucket-unsanitizable, .bucket-error { color: var(--fail); }
.bucket-known { color: var(--known); }
.bucket-pass { color: var(--pass); }
.regression { color: var(--fail); font-weight: 600; }
.flag { display: inline-block; margin-top: 4px; padding: 0 6px; border-radius: 8px; font-size: 12px; font-weight: 500; background: var(--flag-bg); color: var(--flag); }
code { font-size: 12px; word-break: break-all; }
ul { margin: 0; padding-left: 18px; }
.muted { color: var(--muted); }
button { font: inherit; margin: 12px 0; padding: 4px 12px; background: var(--card); color: var(--fg); border: 1px solid var(--line); border-radius: 4px; cursor: pointer; }
details { margin-top: 24px; }
</style>
</head>
<body>
<h1>${title}</h1>
<div class="summary" id="summary"></div>
<div class="controls">
  <label>Bucket <select id="bucket"></select></label>
  <label id="kindLabel">Link <select id="kind"></select></label>
  <label id="covalentLabel">Covalent ligands <select id="covalent"><option value="">included</option><option value="only">only</option><option value="out">left out</option></select></label>
  <label><input type="checkbox" id="diag"> with diagnostics</label>
  <label>Search <input id="search" placeholder="component, bond, reason"></label>
  <span class="muted" id="count"></span>
</div>
<table>
  <thead><tr><th>Component</th><th>Expected (CCD orders)</th><th>Perceived</th><th>Details</th></tr></thead>
  <tbody id="rows"></tbody>
</table>
<button id="more" hidden>Show more</button>
<details id="excludedBox"><summary>Excluded (<span id="excludedCount"></span>)</summary><table><tbody id="excluded"></tbody></table></details>
<script type="application/json" id="data">${json}</script>
<script>
const data = JSON.parse(document.getElementById('data').textContent);
const { run, rows, excluded, viewer } = data;
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const fmt = (o) => Object.entries(o ?? {}).map(([k, v]) => k + ' <b>' + v + '</b>').join(', ');

$('summary').innerHTML = [
  'components <b>' + run.components + '</b>',
  fmt(run.counts),
  'levels: ' + fmt(run.levels),
  ...Object.entries(run.byKind ?? {}).map(([k, v]) => k + ': ' + fmt(v)),
  data.covalent ? '<span class="flag">covalent ligands</span> ' + data.covalent.total + ': ' + fmt(data.covalent.counts) : 'covalent ligands: no tier C run',
  data.diff ? 'against the baseline: ' + fmt(data.diff) : 'no baseline',
].map((s) => '<span>' + s + '</span>').join('');

const buckets = [...new Set(rows.map((r) => r.bucket))];
$('bucket').innerHTML = '<option value="">not passing</option>' + (data.diff ? '<option value="regressions">regressions</option>' : '') + buckets.map((b) => '<option>' + b + '</option>').join('') + '<option value="*">all</option>';
const kinds = [...new Set(rows.map((r) => r.kind).filter(Boolean))];
if (kinds.length) $('kind').innerHTML = '<option value="">all</option>' + kinds.map((k) => '<option>' + k + '</option>').join('');
else $('kindLabel').hidden = true;
if (!data.covalent) $('covalentLabel').hidden = true;

function viewerUrl(r) {
  // the example's state (plan §7): perception in forceCompute, charges removed are its defaults
  // a dataset's pair judged at tier B: the CCD's model coordinates (plan §6)
  const inEntry = r.entry && !r.judgedAt;
  const q = new URLSearchParams(inEntry ? { pdb: r.entry } : { ccd: r.comp ?? r.id });
  if (inEntry) {
    const m = /^(\\S+) \\S+ \\(auth (\\S+)\\)$/.exec(r.copy ?? '');
    if (m) { q.set('asym', m[1]); q.set('seq', m[2]); }
  } else if ((r.judgedAt ?? run.tier) === 'B') q.set('coords', 'model');
  if (run.hydrogens === 'stripped') q.set('removeH', '1');
  if (run.connectivity === 'distance') q.set('connectivity', 'distance');
  return viewer + '?' + q;
}

function picture(r, suffix) {
  if (!r.drawn || (suffix && run.tier === '0')) return '<div class="nodraw">' + (suffix && run.tier === '0' ? 'tier 0: the CCD InChI' : 'not drawn') + '</div>';
  return '<img loading="lazy" alt="' + esc(r.id) + '" src="svg/' + encodeURIComponent(r.id) + suffix + '.svg">';
}

function details(r) {
  const out = [];
  if (r.reason) out.push('<div>' + esc(r.reason) + '</div>');
  if (r.copy) out.push('<div class="muted">' + esc(r.entry) + ' ' + esc(r.copy) + (r.caps ? ', ' + r.caps + ' cap' + (r.caps > 1 ? 's' : '') : '') + '</div>');
  if (r.l0) out.push('<div>L0 <ul>' + r.l0.map((b) => '<li><code>' + esc(b) + '</code></li>').join('') + '</ul></div>');
  if (r.diagnostics) out.push('<div>Diagnostics <ul>' + r.diagnostics.map((d) => '<li>' + esc(d) + '</li>').join('') + '</ul></div>');
  if (r.expected) out.push('<div class="muted">expected <code>' + esc(r.expected) + '</code></div><div class="muted">perceived <code>' + esc(r.perceived) + '</code></div>');
  if (r.leavingMissing) out.push('<div class="muted">' + r.leavingMissing + ' leaving atom(s) missing</div>');
  return out.join('');
}

let filtered = [];
let shown = 0;
const Page = 100;
function render(reset) {
  if (reset) {
    const bucket = $('bucket').value, kind = $('kind').value, diag = $('diag').checked, covalent = $('covalent').value;
    const q = $('search').value.trim().toLowerCase();
    filtered = rows.filter((r) =>
      (bucket === '*' || (bucket === 'regressions' ? r.regression : bucket ? r.bucket === bucket : r.bucket !== 'pass')) &&
      (!kind || r.kind === kind) && (!diag || r.diagnostics) &&
      (!covalent || (covalent === 'only') === !!r.covalent) &&
      (!q || [r.id, r.reason, ...(r.l0 ?? []), ...(r.diagnostics ?? [])].some((s) => s && s.toLowerCase().includes(q))));
    shown = 0;
    $('rows').innerHTML = '';
  }
  const html = filtered.slice(shown, shown + Page).map((r) =>
    '<tr><td class="id"><a href="' + esc(viewerUrl(r)) + '" target="_blank">' + esc(r.id) + '</a>' +
    '<div class="bucket-' + r.bucket + '">' + r.bucket + (r.level ? ' ' + r.level : '') + '</div>' +
    (r.regression ? '<div class="regression" title="the bucket in the baseline">regression (was ' + r.regression + ')</div>' : '') +
    (r.judgedAt ? '<div><span class="flag" title="the dataset judges this pair on the CCD model coordinates, from this entry (plan §6)">tier ' + r.judgedAt + '</span></div>' : '') +
    (r.kind ? '<div class="muted">' + r.kind + '</div>' : '') +
    (r.covalent ? '<div><span class="flag" title="linked in its entry: the CCD has the unbound molecule">covalent ligand</span></div>' : '') + '</td>' +
    '<td>' + picture(r, '.ref') + '</td><td>' + picture(r, '') + '</td><td>' + details(r) + '</td></tr>').join('');
  $('rows').insertAdjacentHTML('beforeend', html);
  shown += Page;
  $('count').textContent = filtered.length + ' rows';
  $('more').hidden = shown >= filtered.length;
}
for (const id of ['bucket', 'kind', 'diag', 'covalent']) $(id).addEventListener('change', () => render(true));
$('search').addEventListener('input', () => render(true));
$('more').addEventListener('click', () => render(false));
render(true);

$('excludedCount').textContent = excluded.length;
$('excludedBox').addEventListener('toggle', () => {
  if (!$('excludedBox').open || $('excluded').childElementCount) return;
  $('excluded').innerHTML = excluded.map((r) => '<tr><td class="id">' + esc(r.id) + '</td><td>' + esc(r.kind ?? '') + (r.covalent ? ' <span class="flag">covalent ligand</span>' : '') + '</td><td>' + esc(r.reason) + '</td></tr>').join('');
});
</script>
</body>
</html>
`;
}

if (import.meta.main) {
  const { values } = parseArgs({
    options: {
      tier: { type: 'string' },
      dataset: { type: 'string', default: 'ccd' },
      selection: { type: 'string', default: 'all' },
      out: { type: 'string', default: 'out' },
      hydrogens: { type: 'string', default: 'stripped' },
      viewer: { type: 'string', default: DefaultViewer },
    },
  });
  // the only connectivity source so far (plan §3.3)
  const runDir = path.join(values.out, values.dataset, 'ccd', partialDir(values.hydrogens, values.selection));
  for (const tier of (values.tier ?? (values.dataset === 'ccd' ? '0,A,B' : 'C')).split(',')) {
    const t0 = performance.now();
    const { rows, drawn } = await writeReport(runDir, tier, values.viewer);
    const file = path.join(runDir, tier, 'index.html');
    console.log(`${file}: ${rows} rows, ${drawn} drawn, ${Math.round((performance.now() - t0) / 1000)} s`);
  }
}
