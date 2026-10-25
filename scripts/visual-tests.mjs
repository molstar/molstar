/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Visual tests in a browser (Playwright): builds the visual test harness, serves the repository and runs the cases of
// each backend profile in parallel pages (profiles x shards). Each profile is compared with its reference profile by the
// harness (the WebGPU profiles with `webgl2`, the WebGL profiles with themselves); `--save`/`--compare` additionally
// store and compare the rendered images, to check that a change keeps a profile identical.
//
// GPU results are not bit-exact from run to run: with ANGLE on D3D11 (NVIDIA), pages rendering in parallel
// occasionally differ by one level in single pixels, and in a few dual depth peeling pixels of
// `transparency/dpoit*/solid-interior` (peeling compares depths for exact equality). Sequential runs and software
// rendering (SwiftShader) are deterministic, so this is driver behavior, not a race in the tests. `--compare` therefore
// ignores differences of one level, and reports changes of the cases in `VariableCases` (solid interiors, and
// illumination, which amplifies single pixel differences) separately, without failing.
//
// usage: node scripts/visual-tests.mjs [--profiles webgl2,webgpu] [--all] [--changed] [--filter 'geometry/**'] [--tags smoke]
//                                      [--shards 8] [--concurrency 32] [--save <dir>] [--compare <dir>]
//                                      [--images <dir>] [--no-build] [--no-debug] [--channel chrome] [--headed]
//                                      [--angle d3d11] [--browser chromium] [--firefox <path>] [--headless]
//
//   --changed      what the uncommitted changes can affect: only the WebGPU profiles when only WebGPU code changed,
//                  only the WebGL profiles when only WebGL code changed, all otherwise; and the cases of the changed
//                  paths (`ChangedCases`: a pass, a group's case file), all cases when a path is not narrowed;
//                  nothing for changes the harness doesn't import (docs, other packages)
//   --save <dir>   stores the images of the WebGL profiles; --compare <dir> reports images that differ from them
//   --compare-head compares with the images of the commit checked out (`HEAD`): a run without uncommitted changes the
//                  harness imports stores its WebGL images as that commit's baseline (`tests/build/visual-baselines/`,
//                  the last few commits kept), so comparing a change needs no run of the commit before it
//   --images <dir> writes target, reference and diff PNGs of failed cases (of all cases with --all-images)
//   --skipped      lists the skipped cases with the reason
//   --draw-lists   the target profiles render with the renderer's draw lists (WebGPU render bundles), each case three
//                  times and the replayed frame compared
//   --dry-run      prints the profiles and the case filter of the run (e.g. what --changed selects) and exits
//   --timing       the time per profile (page loads, renders of target and reference, the slowest shard) and the
//                  slowest cases; a short summary is always printed
//   --angle <name> ANGLE backend (`--use-angle`), default `d3d11` on Windows
//   --browser      `chromium` (Playwright, default) or `firefox`: the installed Firefox (Playwright's build has no
//                  WebGPU), driven over WebDriver BiDi, one window per page; `--firefox <path>` (or `FIREFOX`) names
//                  the executable. Firefox runs headed, as headless Firefox provides no WebGPU adapter; `--headless`
//                  runs it headless anyway (WebGL profiles only)
//
// Profiles run by default: `DefaultProfiles` (WebGL2, WebGL1, WebGPU); `--all` runs `AllProfiles`, adding the reduced
// profiles (no float targets, minimum limits, WebGL1 without extensions, WebGPU core limits), which check the fallbacks;
// others only when named with --profiles (`webgl2-no-ubo`: WebGL2 with uniforms set one by one, the `uniformBlocks`
// option the render benchmark compares).
//
// Shards take every n-th case. Assigning them by earlier times (heaviest first, to the least loaded shard) was measured
// slower (250 s against 174 s for the default profiles): every page then starts with a GPU surface case, whose
// synchronous readbacks wait for the other pages' work, so these cases are best spread over the run as the index does.
//
// The WebGPU profiles compare with `webgl2`: where `webgl2` runs too, its pages write their images (shard by shard)
// and the WebGPU pages of the same shard load them instead of rendering the reference again, except where a fallback
// makes their reference differ. Pages share the GPU process: WebGL commands of all pages run in order, so readbacks
// of one page wait for the others' frames; past about 12 pages at once, runs get slower, not faster (`--concurrency`). Benchmarks run separately:
// `pnpm benchmark` (scripts/benchmark.mjs).

import { spawn, spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';
import { PNG } from 'pngjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2).filter((a) => a !== '--');
const value = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};

const WebGLProfiles = ['webgl2', 'webgl1', 'webgl2-no-float', 'webgl2-minimum', 'webgl1-minimal'];
const WebGPUProfiles = ['webgpu', 'webgpu-core'];
const AllProfiles = [...WebGLProfiles, ...WebGPUProfiles];
const DefaultProfiles = ['webgl2', 'webgl1', 'webgpu'];
/** The reference of the WebGPU profiles, see `runShard` */
const SharedReference = 'webgl2';

const git = (a) =>
  spawnSync('git', a, { cwd: root, encoding: 'utf8' })
    .stdout.split('\n')
    .filter((l) => l);

/** The uncommitted source files the visual tests import (the packages of the harness and the harness) */
function getChangedSources() {
  const files = [...git(['diff', '--name-only', 'HEAD']), ...git(['ls-files', '--others', '--exclude-standard'])];
  return files.filter(
    (f) =>
      /\.(ts|tsx|mjs|js)$/.test(f) &&
      !/(^|\/)_(test|stories|spec)\//.test(f) &&
      /^(packages\/(core|io|model|graphics)\/|extensions\/(alpha-orbitals|geo-export)\/|tests\/src\/visual\/)/.test(f),
  );
}

/** Profiles that uncommitted changes can affect */
function getChangedProfiles(sources) {
  if (sources.length === 0) return [];
  const isWebGPU = (f) => f.includes('/gl/webgpu/') || /webgpu|wgsl/.test(path.basename(f));
  const isWebGL = (f) => f.includes('/gl/webgl/') && !isWebGPU(f);
  if (sources.every(isWebGPU)) return WebGPUProfiles;
  if (sources.every(isWebGL)) return WebGLProfiles;
  return AllProfiles;
}

/**
 * The cases changes of a path can affect, where narrower than all: passes with their own cases, the cases files of a
 * group. Conservative: a path not listed (shader chunks, schemas, the renderer, representations) affects all cases.
 */
const ChangedCases = [
  [
    /^tests\/src\/visual\/cases\/(compute|fallback|geometry|gpu-paths|postprocessing|renderer|transparency)\.ts$/,
    (m) => [`${m[1]}/**`],
  ],
  [/^tests\/src\/visual\/cases\/merged\.ts$/, () => ['geometry/merged/**']],
  [
    /\/canvas3d\/passes\/(postprocessing|ssao|outline|shadow|bloom|dof|fxaa|smaa|cas|background|fog)[^/]*\.ts$/,
    () => ['postprocessing/**'],
  ],
  [/\/canvas3d\/passes\/(dpoit|wboit|transparency)[^/]*\.ts$/, () => ['transparency/**', 'fallback/transparency/**']],
  [/\/canvas3d\/passes\/(illumination|tracing)[^/]*\.ts$/, () => ['renderer/illumination/**', 'fallback/illumination']],
  [/\/canvas3d\/passes\/hi-z[^/]*\.ts$/, () => ['renderer/hi-z/**', 'gpu-paths/hi-z-occlusion']],
  [/\/canvas3d\/passes\/marking[^/]*\.ts$/, () => ['renderer/marking/**']],
  [/\/canvas3d\/passes\/multi-sample[^/]*\.ts$/, () => ['renderer/multisample/**']],
  [/\/gl\/(webgpu\/)?compute\//, () => ['compute/**', 'gpu-paths/**', 'fallback/**']],
];

/** The case filters for `sources`, undefined for all cases */
function getChangedFilters(sources) {
  const filters = new Set();
  for (const f of sources) {
    const rule = ChangedCases.find(([r]) => r.test(f));
    if (!rule) return undefined;
    for (const filter of rule[1](f.match(rule[0]))) filters.add(filter);
  }
  return [...filters].join(',');
}

const changed = args.includes('--changed') ? getChangedSources() : undefined;
const profiles = changed
  ? getChangedProfiles(changed).filter((p) => args.includes('--all') || DefaultProfiles.includes(p))
  : value('--profiles', (args.includes('--all') ? AllProfiles : DefaultProfiles).join(',')).split(',');
const filter = value('--filter') ?? (changed && getChangedFilters(changed));
const tags = value('--tags');
const shards = Number(value('--shards', '8'));
const concurrency = Number(value('--concurrency', String(Math.max(2, Math.min(12, os.cpus().length)))));
const BaselineRoot = path.join(root, 'tests', 'build', 'visual-baselines');
const MaxBaselines = 5;
const head = git(['rev-parse', 'HEAD'])[0];
/** Whether the WebGPU pages load the images of `SharedReference` instead of rendering it */
const sharesReference = profiles.includes(SharedReference) && profiles.some((p) => WebGPUProfiles.includes(p));
/** The baseline of `HEAD`: stored by runs without uncommitted changes, read by `--compare-head` */
const headBaseline = head && path.join(BaselineRoot, head);
const storesBaseline =
  !!headBaseline &&
  !value('--save') &&
  !value('--compare') &&
  !args.includes('--compare-head') &&
  getChangedSources().length === 0;
const saveDir = value('--save') ?? (storesBaseline ? headBaseline : undefined);
const compareDir = value('--compare') ?? (args.includes('--compare-head') ? headBaseline : undefined);
const imagesDir = value('--images');
const allImages = args.includes('--all-images');
const listSkipped = args.includes('--skipped');
const timing = args.includes('--timing');
const channel = value('--channel', 'chrome');
const browserName = value('--browser', 'chromium');
const angle = value('--angle', process.platform === 'win32' ? 'd3d11' : undefined);

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

/** Cases that may vary between runs by more than one level, see above */
const VariableCases = [
  /^transparency\/[^/]+\/solid-interior$/,
  /^renderer\/illumination\//,
  /^fallback\/illumination$/,
];
const isVariable = (id) => VariableCases.some((r) => r.test(id));

/** Whether two PNG data URLs differ by more than one level in a channel of a pixel */
function differs(a, b) {
  if (a === b) return false;
  const [pa, pb] = [a, b].map((url) => PNG.sync.read(Buffer.from(url.split(',')[1], 'base64')));
  if (pa.width !== pb.width || pa.height !== pb.height) return true;
  for (let i = 0; i < pa.data.length; ++i) {
    if (Math.abs(pa.data[i] - pb.data[i]) > 1) return true;
  }
  return false;
}

/**
 * A page of the browser: `goto` loads the harness, `results` waits for its results, `evaluate` calls `fn(arg)` in the
 * page (both JSON-serializable), `consoleErrors` collects the console errors and uncaught exceptions
 */

/** A page without a new result for this long has stalled (a case that never finishes) */
const StallMs = 5 * 60 * 1000;

/**
 * The results of a page once its cases are done; a page stalled in a case (`StallMs` without a new result) gives the
 * results so far and the case it stalled in (`stalled`)
 */
async function waitForResults(evaluate) {
  let count = -1;
  let since = Date.now();
  for (;;) {
    const state = await evaluate(() => ({
      done: !!window.__visualTestResults?.done,
      count: window.__visualTestResults?.results.length ?? 0,
      current: window.__visualTestCurrent,
    }));
    if (state.done) return evaluate(() => window.__visualTestResults);
    if (state.count !== count) {
      count = state.count;
      since = Date.now();
    } else if (Date.now() - since > StallMs) {
      const partial = (await evaluate(() => window.__visualTestResults)) ?? { results: [] };
      return { ...partial, stalled: state.current ?? '(loading the page)' };
    }
    await new Promise((resolve) => setTimeout(resolve, 500));
  }
}

/** Pages of Playwright's Chromium */
async function launchChromium() {
  const browser = await chromium.launch({
    channel,
    headless: !args.includes('--headed'),
    args: [
      '--enable-gpu',
      '--ignore-gpu-blocklist',
      '--enable-unsafe-webgpu',
      ...(angle ? [`--use-angle=${angle}`] : []),
    ],
  });
  return {
    newPage: async () => {
      const page = await browser.newPage();
      const consoleErrors = [];
      page.on('console', (m) => {
        if (m.type() === 'error') consoleErrors.push(m.text());
      });
      page.on('pageerror', (e) => consoleErrors.push(String(e)));
      return {
        consoleErrors,
        goto: (url) => page.goto(url, { timeout: 5 * 60 * 1000 }),
        results: () => waitForResults((fn, arg) => page.evaluate(fn, arg)),
        evaluate: (fn, arg) => page.evaluate(fn, arg),
        close: () => page.close(),
      };
    },
    close: () => browser.close(),
  };
}

function getFirefoxPath() {
  const named = value('--firefox', process.env.FIREFOX);
  if (named) return named;
  if (process.platform === 'win32') return 'C:/Program Files/Mozilla Firefox/firefox.exe';
  if (process.platform === 'darwin') return '/Applications/Firefox.app/Contents/MacOS/firefox';
  return 'firefox';
}

/** Pages of the installed Firefox, over WebDriver BiDi; one window per page (background tabs are throttled) */
async function launchFirefox() {
  const profileDir = fs.mkdtempSync(path.join(os.tmpdir(), 'molstar-visual-firefox-'));
  const prefs = {
    'dom.webgpu.enabled': true,
    'gfx.webgpu.ignore-blocklist': true,
    // windows behind others keep rendering
    'widget.windows.window_occlusion_tracking.enabled': false,
    'browser.shell.checkDefaultBrowser': false,
    'datareporting.policy.dataSubmissionEnabled': false,
  };
  fs.writeFileSync(
    path.join(profileDir, 'user.js'),
    Object.entries(prefs)
      .map(([k, v]) => `user_pref(${JSON.stringify(k)}, ${JSON.stringify(v)});`)
      .join('\n'),
  );
  const port = 9300 + Math.floor(Math.random() * 600);
  const proc = spawn(
    getFirefoxPath(),
    [
      ...(args.includes('--headless') ? ['--headless'] : []),
      '--no-remote',
      '--profile',
      profileDir,
      `--remote-debugging-port=${port}`,
      'about:blank',
    ],
    { stdio: 'ignore' },
  );

  let socket;
  // on Windows `firefox.exe` is a launcher that exits once the browser runs
  for (let i = 0; i < 150 && !socket; ++i) {
    await new Promise((resolve) => setTimeout(resolve, 200));
    socket = await new Promise((resolve) => {
      const ws = new WebSocket(`ws://127.0.0.1:${port}/session`);
      ws.onopen = () => resolve(ws);
      ws.onerror = () => resolve(undefined);
    });
  }
  if (!socket) {
    proc.kill();
    throw new Error(`could not connect to Firefox (${getFirefoxPath()}), see --firefox`);
  }

  let nextId = 0;
  const waiting = new Map();
  /** console errors per browsing context */
  const contextErrors = new Map();
  socket.onmessage = (e) => {
    const m = JSON.parse(e.data);
    if (m.id !== undefined && waiting.has(m.id)) {
      waiting.get(m.id)(m);
      waiting.delete(m.id);
    } else if (m.method === 'log.entryAdded' && m.params.level === 'error') {
      contextErrors.get(m.params.source?.context)?.push(m.params.text ?? '');
    }
  };
  const send = (method, params) =>
    new Promise((resolve, reject) => {
      const id = ++nextId;
      waiting.set(id, (m) => (m.error ? reject(new Error(`${method}: ${m.error}: ${m.message}`)) : resolve(m.result)));
      socket.send(JSON.stringify({ id, method, params }));
    });
  await send('session.new', { capabilities: {} });
  await send('session.subscribe', { events: ['log.entryAdded'] });

  return {
    newPage: async () => {
      const { context } = await send('browsingContext.create', { type: 'window' });
      const consoleErrors = [];
      contextErrors.set(context, consoleErrors);
      // results are passed as JSON strings, which BiDi returns as plain values
      const evaluate = async (fn, arg) => {
        const expression = `(async () => JSON.stringify(await (${fn})(${JSON.stringify(arg)})))()`;
        const r = await send('script.evaluate', { expression, target: { context }, awaitPromise: true });
        if (r.type === 'exception') throw new Error(r.exceptionDetails?.text ?? 'exception in the page');
        return r.result.value === undefined ? undefined : JSON.parse(r.result.value);
      };
      return {
        consoleErrors,
        goto: (url) => send('browsingContext.navigate', { context, url, wait: 'complete' }),
        results: () => waitForResults(evaluate),
        evaluate,
        close: async () => {
          contextErrors.delete(context);
          await send('browsingContext.close', { context });
        },
      };
    },
    close: async () => {
      try {
        await send('browser.close', {});
      } catch {
        proc.kill();
      }
      socket.close();
      // the profile is in use until the browser has exited
      for (let i = 0; i < 50; ++i) {
        try {
          fs.rmSync(profileDir, { recursive: true, force: true });
          break;
        } catch {
          await new Promise((resolve) => setTimeout(resolve, 200));
        }
      }
    },
  };
}

/** Where the images of the shared reference are written for a run, served with the repository */
const sharedDir = path.join(root, 'tests', 'build', 'visual-shared', String(process.pid));
const sharedUrl = (origin) => `${origin}/${path.relative(root, sharedDir).replaceAll(path.sep, '/')}/`;

/**
 * Runs one shard of a profile in its own page; the shard of the shared reference writes its images to `sharedDir`, a
 * shard of a profile comparing with it waits for them (`referenceShard`, the same shard of the reference) and loads them
 */
async function runShard(browser, origin, profile, shard, referenceShard) {
  const reference = WebGPUProfiles.includes(profile) ? SharedReference : profile;
  const shared =
    referenceShard &&
    (await referenceShard.then(
      () => true,
      () => false,
    ));
  const page = await browser.newPage();
  const { consoleErrors } = page;
  const q = new URLSearchParams({ target: profile, reference, shard: `${shard}/${shards}` });
  if (shared) q.set('shared', sharedUrl(origin));
  if (!args.includes('--no-debug')) q.set('debug', '1');
  if (filter) q.set('filter', filter);
  if (tags) q.set('tags', tags);
  if (args.includes('--draw-lists')) q.set('drawLists', '1');
  const pageStart = Date.now();
  // many pages load at once
  await page.goto(`${origin}/tests/build/visual/index.html?${q}`);
  const loadMs = Date.now() - pageStart;
  const results = await page.results();
  if (results.stalled) {
    results.results.push({
      id: results.stalled,
      status: 'error',
      errors: [`stalled: no result within ${StallMs / 60000} min, the shard's later cases not run`],
    });
  }

  const rows = [];
  for (const r of results.results) {
    const row = {
      id: r.id,
      status: r.status,
      errors: r.errors,
      reason: r.reason,
      metrics: r.metrics,
      timeMs: r.timeMs,
      warnings: r.warnings,
    };
    if (r.status !== 'skipped') {
      const image = (kind) => page.evaluate(([id, k]) => window.__visualTestImage(id, k), [r.id, kind]);
      // check cases have no image
      const file = (saveDir || compareDir) && path.join(saveDir ?? compareDir, profile, `${fileName(r.id)}.txt`);
      // a baseline keeps the images it has (of the same commit)
      const stored = storesBaseline && file && fs.existsSync(file);
      const writesShared = profile === SharedReference && sharesReference && r.status !== 'error';
      const target =
        (file && !stored && !WebGPUProfiles.includes(profile)) || writesShared ? await image('target') : undefined;
      if (target && writesShared) {
        fs.mkdirSync(sharedDir, { recursive: true });
        fs.writeFileSync(path.join(sharedDir, `${fileName(r.id)}.png`), Buffer.from(target.split(',')[1], 'base64'));
      }
      if (target && file && !stored) {
        if (saveDir) {
          fs.mkdirSync(path.dirname(file), { recursive: true });
          fs.writeFileSync(file, target);
        } else {
          const saved = fs.existsSync(file) ? fs.readFileSync(file, 'utf8') : undefined;
          row.compare = saved === undefined ? 'new' : differs(saved, target) ? 'changed' : 'same';
        }
      }
      if (imagesDir && (allImages || r.status === 'fail' || r.status === 'error')) {
        const dir = path.join(imagesDir, profile);
        fs.mkdirSync(dir, { recursive: true });
        for (const kind of ['target', 'reference', 'diff']) {
          const url = await image(kind);
          if (url)
            fs.writeFileSync(path.join(dir, `${fileName(r.id)}.${kind}.png`), Buffer.from(url.split(',')[1], 'base64'));
        }
      }
    }
    rows.push(row);
  }
  const totalMs = Date.now() - pageStart;
  await page.close();
  return { rows, consoleErrors, loadMs, totalMs };
}

const seconds = (ms) => `${(ms / 1000).toFixed(1)} s`;

/** Where the time of the run went: per profile, and the slowest cases (target and reference renders) */
function printTiming(jobs, shardResults) {
  const cases = [];
  for (const profile of new Set(jobs.map((j) => j.profile))) {
    const parts = shardResults.filter((_, i) => jobs[i].profile === profile);
    const rows = parts.flatMap((p) => p.rows);
    const target = rows.reduce((s, r) => s + (r.timeMs?.target ?? 0), 0);
    const reference = rows.reduce((s, r) => s + (r.timeMs?.reference ?? 0), 0);
    const load = parts.reduce((s, p) => s + p.loadMs, 0);
    const slowest = Math.max(...parts.map((p) => p.totalMs));
    const run = rows.filter((r) => r.status !== 'skipped').length;
    console.log(
      `  ${profile}: ${run} cases, renders ${seconds(target)} target + ${seconds(reference)} reference, page loads ${seconds(load)}, slowest shard ${seconds(slowest)}`,
    );
    for (const r of rows) {
      const ms = (r.timeMs?.target ?? 0) + (r.timeMs?.reference ?? 0);
      if (ms > 0) cases.push({ id: r.id, profile, ms });
    }
  }
  cases.sort((a, b) => b.ms - a.ms);
  console.log('  slowest cases (target + reference):');
  for (const c of cases.slice(0, timing ? 25 : 8)) console.log(`    ${seconds(c.ms)}  ${c.profile}  ${c.id}`);
  if (timing) {
    // the cases' renders summed over the profiles: what a case costs a run
    const byCase = new Map();
    for (const c of cases) byCase.set(c.id, (byCase.get(c.id) ?? 0) + c.ms);
    const total = [...byCase.values()].reduce((s, ms) => s + ms, 0);
    console.log(`  cases summed over the profiles (${seconds(total)} in all):`);
    for (const [id, ms] of [...byCase].sort((a, b) => b[1] - a[1]).slice(0, 25)) {
      console.log(`    ${seconds(ms)}  ${((100 * ms) / total).toFixed(1)}%  ${id}`);
    }
  }
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

/** Keeps the baselines of the last `MaxBaselines` commits run */
function pruneBaselines() {
  if (!fs.existsSync(BaselineRoot)) return;
  const dirs = fs
    .readdirSync(BaselineRoot)
    .map((name) => ({ name, time: fs.statSync(path.join(BaselineRoot, name)).mtimeMs }))
    .sort((a, b) => b.time - a.time);
  for (const d of dirs.slice(MaxBaselines))
    fs.rmSync(path.join(BaselineRoot, d.name), { recursive: true, force: true });
}

async function main() {
  if (args.includes('--dry-run')) {
    console.log(`profiles: ${profiles.join(', ') || 'none'}`);
    console.log(`cases: ${filter || 'all'}${tags ? ` (tags ${tags})` : ''}`);
    if (changed) {
      console.log('changed sources the harness imports:');
      for (const f of changed) console.log(`  ${f}`);
    }
    return;
  }
  if (profiles.length === 0) {
    console.log('no source changes, nothing to run');
    return;
  }
  if (args.includes('--compare-head') && !fs.existsSync(headBaseline)) {
    console.log(`no baseline of HEAD (${head.slice(0, 9)}): run the cases once without uncommitted changes first`);
    process.exit(1);
  }
  if (storesBaseline) {
    fs.mkdirSync(headBaseline, { recursive: true });
    // the commit run last is kept longest
    fs.utimesSync(headBaseline, new Date(), new Date());
    pruneBaselines();
  }
  if (!args.includes('--no-build')) {
    console.log('building the visual tests ...');
    const r = spawnSync(process.execPath, ['scripts/esbuild/visual-tests.mjs', '--no-watch'], {
      cwd: root,
      stdio: 'inherit',
    });
    if (r.status !== 0) process.exit(r.status ?? 1);
  }

  const server = await serve();
  const origin = `http://127.0.0.1:${server.address().port}`;
  if (browserName !== 'chromium' && browserName !== 'firefox') throw new Error(`unknown browser '${browserName}'`);
  const browser = browserName === 'firefox' ? await launchFirefox() : await launchChromium();

  const start = Date.now();
  // the shared reference first, then the other WebGL profiles: the WebGPU shards wait for its images
  const order = (p) => (p === SharedReference ? 0 : WebGPUProfiles.includes(p) ? 2 : 1);
  const ordered = [...profiles].sort((a, b) => order(a) - order(b));
  const jobs = ordered.flatMap((profile) => Array.from({ length: shards }, (_, shard) => ({ profile, shard })));
  const referenceShards = new Map();
  let shardResults;
  try {
    shardResults = await pool(
      jobs.map((j) => () => {
        const referenceShard =
          sharesReference && WebGPUProfiles.includes(j.profile) ? referenceShards.get(j.shard) : undefined;
        const run = runShard(browser, origin, j.profile, j.shard, referenceShard);
        if (sharesReference && j.profile === SharedReference) referenceShards.set(j.shard, run);
        return run;
      }),
      concurrency,
    );
  } finally {
    await browser.close();
    server.close();
    fs.rmSync(sharedDir, { recursive: true, force: true });
  }

  let failed = false;
  for (const profile of ordered) {
    const parts = shardResults.filter((_, i) => jobs[i].profile === profile);
    const rows = parts.flatMap((p) => p.rows);
    const counts = {};
    for (const r of rows) counts[r.status] = (counts[r.status] ?? 0) + 1;
    const bad = rows.filter((r) => r.status === 'fail' || r.status === 'error');
    const changed = rows.filter((r) => r.compare === 'changed' && !isVariable(r.id));
    const varied = rows.filter((r) => r.compare === 'changed' && isVariable(r.id));
    const added = rows.filter((r) => r.compare === 'new');
    failed ||= bad.length > 0 || changed.length > 0;
    console.log(
      `${profile}: ${Object.entries(counts)
        .map(([k, n]) => `${n} ${k}`)
        .join(', ')}`,
    );
    if (listSkipped) {
      for (const r of rows.filter((r) => r.status === 'skipped')) console.log(`  skipped ${r.id}: ${r.reason ?? ''}`);
    }
    if (args.includes('--draw-lists')) {
      const direct = rows.filter((r) => r.warnings?.some((w) => /replayed nothing/.test(w)));
      console.log(
        `  draw lists: ${direct.length} cases drawn without replay${direct.length ? `: ${direct.map((r) => r.id).join(', ')}` : ''}`,
      );
    }
    for (const r of bad) {
      const m = r.metrics ? ` (${r.metrics.diffPixels} px differ)` : '';
      console.log(`  ${r.status} ${r.id}${m}${r.errors?.length ? `: ${r.errors.join(' | ').slice(0, 300)}` : ''}`);
    }
    // only the WebGL profiles' images are stored and compared
    if (compareDir && !WebGPUProfiles.includes(profile)) {
      console.log(
        `  changed vs ${compareDir === headBaseline ? `HEAD (${head.slice(0, 9)})` : compareDir}: ${changed.length ? changed.map((r) => r.id).join(', ') : 'none'}`,
      );
      if (varied.length) console.log(`  changed, may vary between runs: ${varied.map((r) => r.id).join(', ')}`);
      if (added.length) console.log(`  new (not in ${compareDir}): ${added.map((r) => r.id).join(', ')}`);
    }
    const consoleErrors = [...new Set(parts.flatMap((p) => p.consoleErrors))].filter(
      (e) => !/Failed to load resource/.test(e),
    );
    for (const e of consoleErrors.slice(0, 10)) console.log(`  console: ${e.slice(0, 300)}`);
  }
  console.log('timing:');
  printTiming(jobs, shardResults);
  console.log(
    `${jobs.length} pages (${profiles.length} profiles x ${shards} shards) in ${((Date.now() - start) / 1000).toFixed(0)} s`,
  );
  process.exit(failed ? 1 : 0);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
