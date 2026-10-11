/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Visual tests in a browser (Playwright): builds the visual test harness, serves the repository and runs the cases of
// each backend profile in parallel pages (profiles x shards). Each profile is compared with its reference profile by the
// harness (`webgpu` with `webgl2`, the WebGL profiles with themselves); `--save`/`--compare` additionally store or
// compare the rendered images byte for byte, to check that a change keeps a profile identical.
//
// usage: node scripts/visual-tests.mjs [--profiles webgl2,webgpu] [--changed] [--filter 'geometry/**'] [--tags smoke]
//                                      [--shards 8] [--concurrency 32] [--save <dir>] [--compare <dir>]
//                                      [--images <dir>] [--no-build] [--no-debug] [--channel chrome] [--headed]
//
//   --changed      profiles affected by the uncommitted changes: only `webgpu` when only WebGPU code changed, only the
//                  WebGL profiles when only WebGL code changed, all otherwise
//   --save <dir>   stores the images of the WebGL profiles; --compare <dir> reports images that differ from them
//   --images <dir> writes target, reference and diff PNGs of failed cases (of all cases with --all-images)

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2).filter((a) => a !== '--');
const value = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};

const WebGLProfiles = ['webgl2', 'webgl1', 'webgl2-no-float', 'webgl2-minimum', 'webgl1-minimal'];
const AllProfiles = [...WebGLProfiles, 'webgpu'];

/** Profiles that uncommitted changes can affect */
function getChangedProfiles() {
  const git = (a) => spawnSync('git', a, { cwd: root, encoding: 'utf8' }).stdout.split('\n').filter((l) => l);
  const files = [...git(['diff', '--name-only', 'HEAD']), ...git(['ls-files', '--others', '--exclude-standard'])];
  const sources = files.filter((f) => /\.(ts|tsx|mjs|js)$/.test(f) && !f.startsWith('.v6/'));
  if (sources.length === 0) return [];
  const isWebGPU = (f) => f.includes('/gl/webgpu/') || /webgpu|wgsl/.test(path.basename(f));
  const isWebGL = (f) => f.includes('/gl/webgl/') && !isWebGPU(f);
  if (sources.every(isWebGPU)) return ['webgpu'];
  if (sources.every(isWebGL)) return WebGLProfiles;
  return AllProfiles;
}

const profiles = args.includes('--changed') ? getChangedProfiles() : value('--profiles', AllProfiles.join(',')).split(',');
const filter = value('--filter');
const tags = value('--tags');
const shards = Number(value('--shards', '8'));
const concurrency = Number(value('--concurrency', String(Math.max(2, Math.min(32, os.cpus().length)))));
const saveDir = value('--save');
const compareDir = value('--compare');
const imagesDir = value('--images');
const allImages = args.includes('--all-images');
const channel = value('--channel', 'chrome');

const mime = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json' };

function serve() {
  // the pages of all shards load the same (large) bundle: read each file once
  const cache = new Map();
  const server = http.createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      const resolved = path.resolve(root, `.${pathname}`);
      if (!resolved.startsWith(root + path.sep)) return res.writeHead(403).end();
      const stat = await fs.promises.stat(resolved);
      const file = stat.isDirectory() ? path.join(resolved, 'index.html') : resolved;
      let body = cache.get(file);
      if (!body) cache.set(file, (body = await fs.promises.readFile(file)));
      res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(body);
    } catch {
      res.writeHead(404).end('Not found');
    }
  });
  return new Promise((resolve) => server.listen(0, '127.0.0.1', () => resolve(server)));
}

const fileName = (id) => id.replace(/\//g, '.');

/** Runs one shard of a profile in its own page */
async function runShard(browser, origin, profile, shard) {
  const page = await browser.newPage();
  const consoleErrors = [];
  page.on('console', (m) => {
    if (m.type() === 'error') consoleErrors.push(m.text());
  });
  page.on('pageerror', (e) => consoleErrors.push(String(e)));
  const reference = profile === 'webgpu' ? 'webgl2' : profile;
  const q = new URLSearchParams({ target: profile, reference, shard: `${shard}/${shards}` });
  if (!args.includes('--no-debug')) q.set('debug', '1');
  if (filter) q.set('filter', filter);
  if (tags) q.set('tags', tags);
  // many pages load at once
  await page.goto(`${origin}/tests/build/visual/index.html?${q}`, { timeout: 5 * 60 * 1000 });
  await page.waitForFunction(() => window.__visualTestResults?.done, null, { timeout: 60 * 60 * 1000, polling: 500 });
  const results = await page.evaluate(() => window.__visualTestResults);

  const rows = [];
  for (const r of results.results) {
    const row = { id: r.id, status: r.status, errors: r.errors, reason: r.reason, metrics: r.metrics };
    if (r.status !== 'skipped') {
      const image = (kind) => page.evaluate(([id, k]) => window.__visualTestImage(id, k), [r.id, kind]);
      if ((saveDir || compareDir) && profile !== 'webgpu') {
        const target = await image('target');
        const file = path.join(saveDir ?? compareDir, profile, `${fileName(r.id)}.txt`);
        if (saveDir) {
          fs.mkdirSync(path.dirname(file), { recursive: true });
          fs.writeFileSync(file, target);
        } else {
          const saved = fs.existsSync(file) ? fs.readFileSync(file, 'utf8') : undefined;
          row.compare = saved === undefined ? 'new' : saved === target ? 'same' : 'changed';
        }
      }
      if (imagesDir && (allImages || r.status === 'fail' || r.status === 'error')) {
        const dir = path.join(imagesDir, profile);
        fs.mkdirSync(dir, { recursive: true });
        for (const kind of ['target', 'reference', 'diff']) {
          const url = await image(kind);
          if (url) fs.writeFileSync(path.join(dir, `${fileName(r.id)}.${kind}.png`), Buffer.from(url.split(',')[1], 'base64'));
        }
      }
    }
    rows.push(row);
  }
  await page.close();
  return { rows, consoleErrors };
}

/** Runs `tasks` with at most `limit` at a time */
async function pool(tasks, limit) {
  const results = new Array(tasks.length);
  let next = 0;
  const worker = async () => {
    while (next < tasks.length) {
      const i = next++;
      results[i] = await tasks[i]();
    }
  };
  await Promise.all(Array.from({ length: Math.min(limit, tasks.length) }, worker));
  return results;
}

async function main() {
  if (profiles.length === 0) {
    console.log('no source changes, nothing to run');
    return;
  }
  if (!args.includes('--no-build')) {
    console.log('building the visual tests ...');
    const r = spawnSync(process.execPath, ['scripts/esbuild/visual-tests.mjs', '--no-watch'], { cwd: root, stdio: 'inherit' });
    if (r.status !== 0) process.exit(r.status ?? 1);
  }

  const server = await serve();
  const origin = `http://127.0.0.1:${server.address().port}`;
  const browser = await chromium.launch({
    channel,
    headless: !args.includes('--headed'),
    args: [
      '--enable-gpu',
      '--ignore-gpu-blocklist',
      '--enable-unsafe-webgpu',
      ...(process.platform === 'win32' ? ['--use-angle=d3d11'] : []),
    ],
  });

  const start = Date.now();
  const jobs = profiles.flatMap((profile) => Array.from({ length: shards }, (_, shard) => ({ profile, shard })));
  let shardResults;
  try {
    shardResults = await pool(
      jobs.map((j) => () => runShard(browser, origin, j.profile, j.shard)),
      concurrency,
    );
  } finally {
    await browser.close();
    server.close();
  }

  let failed = false;
  for (const profile of profiles) {
    const parts = shardResults.filter((_, i) => jobs[i].profile === profile);
    const rows = parts.flatMap((p) => p.rows);
    const counts = {};
    for (const r of rows) counts[r.status] = (counts[r.status] ?? 0) + 1;
    const bad = rows.filter((r) => r.status === 'fail' || r.status === 'error');
    const changed = rows.filter((r) => r.compare === 'changed');
    const added = rows.filter((r) => r.compare === 'new');
    failed ||= bad.length > 0 || changed.length > 0;
    console.log(`${profile}: ${Object.entries(counts).map(([k, n]) => `${n} ${k}`).join(', ')}`);
    for (const r of bad) {
      const m = r.metrics ? ` (${r.metrics.diffPixels} px differ)` : '';
      console.log(`  ${r.status} ${r.id}${m}${r.errors?.length ? `: ${r.errors.join(' | ').slice(0, 300)}` : ''}`);
    }
    if (compareDir) {
      console.log(`  changed vs ${compareDir}: ${changed.length ? changed.map((r) => r.id).join(', ') : 'none'}`);
      if (added.length) console.log(`  new (not in ${compareDir}): ${added.map((r) => r.id).join(', ')}`);
    }
    const consoleErrors = [...new Set(parts.flatMap((p) => p.consoleErrors))].filter((e) => !/Failed to load resource/.test(e));
    for (const e of consoleErrors.slice(0, 10)) console.log(`  console: ${e.slice(0, 300)}`);
  }
  console.log(`${jobs.length} pages (${profiles.length} profiles x ${shards} shards) in ${((Date.now() - start) / 1000).toFixed(0)} s`);
  process.exit(failed ? 1 : 0);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
