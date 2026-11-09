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
//                                      [--groups] [--stall <s>] [--browser-arg <flag>]
//
//   --changed      a quick check of the uncommitted changes, not a full run: only the WebGPU profiles when only WebGPU
//                  code changed, only the WebGL profiles when only WebGL code changed, all otherwise; and the cases of
//                  the changed paths (`ChangedCases`: a pass, a geometry, the gizmos, a group's case file), all cases
//                  when a path is not narrowed; a geometry's own files run its cases only, not the others drawing it
//                  (structures as meshes, the clip ghost); nothing for changes the harness doesn't import (docs, other
//                  packages)
//   --groups       one group of cases after another (the --filter patterns, or each top-level group), a line each as
//                  it ends: which group is slow or failing shows early; a little slower in all (pages load per group)
//   --stall <s>    a case without a result for this long has stalled (default 60): reported, and its shard goes on in
//                  a new page from the next case
//   --budget <s>   the run is aborted (exit 2) once it is projected to take longer: three times what the last finished
//                  run with the same profiles, shards, filter, tags and browser took (`tests/build/visual-timing.json`),
//                  at least 120 s; this many seconds where there is none yet (default 600)
//
// The progress is printed every 10 s (cases done, failed, the case running longest), failures as they come.
//   --save <dir>   stores the images of the WebGL profiles; --compare <dir> reports images that differ from them
//   --compare-head compares with the images of the commit checked out (`HEAD`): a run without uncommitted changes the
//                  harness imports stores its WebGL images as that commit's baseline (`tests/build/visual-baselines/`,
//                  the last few commits kept), so comparing a change needs no run of the commit before it
//   --images <dir> writes target, reference and diff PNGs of failed cases (of all cases with --all-images)
//   --skipped      lists the skipped cases with the reason
//   --require-backends  fails where a profile's backend is not available (by default its cases are skipped): CI, so
//                  that a browser without WebGPU does not pass by skipping it
//   --allow-software  runs on a software renderer (SwiftShader, WARP, Microsoft Basic Render Driver, llvmpipe,
//                  lavapipe, a WebGPU fallback adapter): by default the run fails before its cases, as a browser without
//                  the GPU (e.g. over Remote Desktop) is many times slower; CI's software Vulkan passes it
//   --draw-lists   the target profiles render with the renderer's draw lists (WebGPU render bundles), each case three
//                  times and the replayed frame compared
//   --dry-run      prints the profiles and the case filter of the run (e.g. what --changed selects) and exits
//   --timing       the time per profile (page loads, renders of target and reference, the slowest shard) and the
//                  slowest cases; a short summary is always printed
//   --angle <name> ANGLE backend (`--use-angle`), default `d3d11` on Windows
//   --browser-arg <flag>  an extra Chromium flag (repeatable), e.g. a software WebGPU adapter on a machine without a
//                  GPU (`--browser-arg=--use-webgpu-adapter=swiftshader`: WebGPU on SwiftShader, not combinable with
//                  `--angle swiftshader`, which leaves Chromium without WebGPU)
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
 * The cases of a geometry's own files (its params, renderable and shaders): its cases and the gizmos drawing it. Not the
 * other cases drawing it (structures as meshes or spheres): those are left to a full run; nor every shader variant
 * (`shaders/**`, run on demand).
 */
const GeometryCases = {
  mesh: ['geometry/mesh/**', 'geometry/texture-mesh/**', 'overlay/**'],
  'texture-mesh': ['geometry/texture-mesh/**', 'gpu-paths/**', 'fallback/**'],
  lines: ['geometry/lines/**', 'overlay/**'],
  points: ['geometry/points/**', 'overlay/**'],
  spheres: ['geometry/spheres/**', 'renderer/sphere-near**', 'renderer/sphere-raster/**'],
  cylinders: ['geometry/cylinders/**'],
  text: ['geometry/text/**', 'overlay/**', 'screen/**'],
  image: ['geometry/image/**'],
  'direct-volume': ['geometry/direct-volume/**'],
};

/**
 * The cases changes of a path can affect, where narrower than all: passes with their own cases, a geometry's files
 * (`GeometryCases`), the gizmos and helpers, the cases files of a group. Conservative otherwise: a path not listed
 * (shader chunks, schemas, the renderer, representations) affects all cases.
 */
const ChangedCases = [
  [
    /^tests\/src\/visual\/cases\/(compute|fallback|geometry|gpu-paths|postprocessing|renderer|transparency)\.ts$/,
    (m) => [`${m[1]}/**`],
  ],
  [/^tests\/src\/visual\/cases\/merged\.ts$/, () => ['geometry/merged/**']],
  [/^tests\/src\/visual\/cases\/sphere-raster\.ts$/, () => ['renderer/sphere-raster/**']],
  [/^tests\/src\/visual\/cases\/overlay\.ts$/, () => ['overlay/**', 'screen/**']],
  [/^tests\/src\/visual\/cases\/shader-variants\.ts$/, () => ['shaders/**']],
  [
    /^tests\/src\/visual\/cases\/theme-data\.ts$/,
    () => [
      'geometry/palette**',
      'renderer/wiggle/**',
      'renderer/animation/**',
      'renderer/bump/**',
      'renderer/lod-fade/**',
    ],
  ],
  [
    /\/(geo\/geometry|gl\/renderable|gl\/webgl\/shader|gl\/webgpu\/shader)\/(mesh|texture-mesh|lines|points|spheres|cylinders|text|image|direct-volume)[/.]/,
    (m) => GeometryCases[m[2]],
  ],
  [/\/util\/gizmo\/|\/canvas3d\/helper\/overlay-helper\.ts$/, () => ['overlay/**']],
  [/\/canvas3d\/helper\/(camera-helper|screen-helper|scale-line)\.ts$/, () => ['screen/**']],
  [
    /\/canvas3d\/passes\/(postprocessing|ssao|outline|shadow|bloom|dof|fxaa|smaa|cas|background|fog)[^/]*\.ts$/,
    () => ['postprocessing/**'],
  ],
  [/\/canvas3d\/passes\/(dpoit|wboit|transparency)[^/]*\.ts$/, () => ['transparency/**', 'fallback/transparency/**']],
  [/\/canvas3d\/passes\/(illumination|tracing)[^/]*\.ts$/, () => ['renderer/illumination/**', 'fallback/illumination']],
  [/\/canvas3d\/passes\/hi-z[^/]*\.ts$/, () => ['renderer/hi-z/**', 'gpu-paths/hi-z-occlusion']],
  [/\/canvas3d\/passes\/marking[^/]*\.ts$/, () => ['renderer/marking/**']],
  [/\/canvas3d\/passes\/multi-sample[^/]*\.ts$/, () => ['renderer/multisample/**']],
  [/\/gl\/(compute|webgpu\/shader)\/sphere-raster[^/]*\.ts$/, () => ['renderer/sphere-raster/**']],
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
const requireBackends = args.includes('--require-backends');
const allowSoftware = args.includes('--allow-software');
/** Renderer and adapter names of software rendering */
const SoftwareRenderer = /swiftshader|llvmpipe|lavapipe|softpipe|basic render|\bwarp\b/i;
const timing = args.includes('--timing');
const channel = value('--channel', 'chrome');
const browserName = value('--browser', 'chromium');
const angle = value('--angle', process.platform === 'win32' ? 'd3d11' : undefined);
/** `--browser-arg <flag>` or `--browser-arg=<flag>`, repeatable: extra Chromium flags */
const browserArgs = args.flatMap((a, i) =>
  a.startsWith('--browser-arg=') ? [a.slice('--browser-arg='.length)] : args[i - 1] === '--browser-arg' ? [a] : [],
);

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
 * A page of the browser: `goto` loads the harness, `results` waits for its results (`waitForResults`), `evaluate` calls
 * `fn(arg)` in the page (both JSON-serializable), `consoleErrors` collects the console errors and uncaught exceptions
 */

/**
 * A page without a new result for this long has stalled (a case that never finishes): the slowest case takes seconds,
 * so a minute reports a hang early, and the shard goes on past it (`runShard`)
 */
const StallMs = Number(value('--stall', '60')) * 1000;
/** The run's own times, by its key (`runKey`): what `--budget` compares a run with */
const TimingFile = path.join(root, 'tests', 'build', 'visual-timing.json');
/** A run slower than this many times the last one's is aborted; never before `MinBudgetMs` */
const BudgetFactor = 3;
const MinBudgetMs = 120 * 1000;
/** The projection waits this long: the first pages load and compile */
const ProjectAfterMs = 20 * 1000;
/** Set while an over-budget run closes its browser: the pages' calls it cuts off are not the run's errors */
let aborting = false;
/**
 * A page not answering for this long is hung (a loop in the page): its results so far are those it reported. Not
 * shorter: a page in a GPU surface case's synchronous readbacks waits for the other pages' GPU work, measured up to
 * 29 s with the default profiles in 8 shards, and is not hung
 */
const UnresponsiveMs = 90 * 1000;
/** How often the progress is printed */
const ProgressMs = 10 * 1000;

/** `promise`, or undefined once `ms` passed */
const within = (promise, ms) => Promise.race([promise, new Promise((resolve) => setTimeout(resolve, ms).unref())]);

/**
 * The run's progress: the cases done and failed, of `total` (the cases of the pages loaded so far, `pagesLeft` not
 * yet), and the case each page runs since when, printed every `ProgressMs` (`printProgress`), failures at once
 */
const progress = { total: 0, done: 0, failed: 0, pagesLeft: 0, running: new Set(), start: Date.now() };

function printProgress() {
  const now = Date.now();
  let longest;
  for (const r of progress.running) if (r.id && (!longest || r.since < longest.since)) longest = r;
  const of = `${progress.total}${progress.pagesLeft > 0 ? '+' : ''}`;
  const failed = progress.failed ? `, ${progress.failed} failed` : '';
  const running = longest ? `; longest running ${longest.label} ${longest.id} ${seconds(now - longest.since)}` : '';
  console.log(`  [${seconds(now - progress.start)}] ${progress.done}/${of} cases${failed}${running}`);
}

/**
 * The results of a page once its cases are done (`count`: how many it ran, `label` names it in the progress); a page
 * stalled in a case (`StallMs` without a new result) gives the results so far and the case it stalled in (`stalled`;
 * `loading`: before its first case; `counted`: its cases added to the progress), a hung one only the ids and status it reported (`hung`); `counts` adds its cases
 * to the progress' total
 */
async function waitForResults(evaluate, label, counts) {
  /** id and status of the results reported */
  const seen = [];
  let total;
  let since = Date.now();
  const entry = { label, id: undefined, since };
  progress.running.add(entry);
  try {
    for (;;) {
      const state = await within(
        evaluate(
          (n) => ({
            done: !!window.__visualTestResults?.done,
            total: window.__visualTestCount,
            added: (window.__visualTestResults?.results ?? []).slice(n).map((r) => ({ id: r.id, status: r.status })),
            current: window.__visualTestCurrent,
          }),
          seen.length,
        ),
        UnresponsiveMs,
      );
      const now = Date.now();
      if (!state)
        return {
          results: seen,
          count: total ?? seen.length + 1,
          stalled: entry.id ?? '(loading the page)',
          loading: !entry.id && seen.length === 0,
          counted: !!counts && total !== undefined,
          hung: true,
        };
      if (total === undefined && state.total !== undefined) {
        total = state.total;
        if (counts) {
          progress.total += total;
          progress.pagesLeft--;
        }
      }
      for (const r of state.added) {
        seen.push(r);
        progress.done++;
        if (r.status === 'fail' || r.status === 'error') {
          progress.failed++;
          console.log(`  ${r.status} ${label} ${r.id}`);
        }
      }
      if (state.added.length > 0) since = now;
      if (state.current !== entry.id) Object.assign(entry, { id: state.current, since: now });
      if (state.done) return { ...(await evaluate(() => window.__visualTestResults)), count: total };
      if (now - since > StallMs) {
        const partial = (await within(
          evaluate(() => window.__visualTestResults),
          UnresponsiveMs,
        )) ?? { results: seen };
        // the cases before the stalled one: it may have ended since (and the page gone on), which would then be reported
        // twice and the case after it skipped (the shard goes on from the case after the stalled one)
        const at = state.current === undefined ? -1 : partial.results.findIndex((r) => r.id === state.current);
        const results = at < 0 ? partial.results : partial.results.slice(0, at);
        return {
          ...partial,
          results,
          count: total ?? results.length + 1,
          stalled: state.current ?? '(loading the page)',
          loading: !state.current && results.length === 0,
          counted: !!counts && total !== undefined,
        };
      }
      await new Promise((resolve) => setTimeout(resolve, 500));
    }
  } finally {
    progress.running.delete(entry);
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
      // the sphere raster's `atomic-vec2u-min-max` (Chrome 156)
      '--enable-webgpu-developer-features',
      ...(angle ? [`--use-angle=${angle}`] : []),
      ...browserArgs,
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
        results: (label, counts) => waitForResults((fn, arg) => page.evaluate(fn, arg), label, counts),
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
        results: (label, counts) => waitForResults(evaluate, label, counts),
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
 * Runs one shard of a profile (the cases of `caseFilter`) in its own page; the shard of the shared reference writes its
 * images to `sharedDir`, a shard of a profile comparing with it waits for them (`referenceShard`, the same shard of the
 * reference) and loads them. A stalled case ends its page: the shard goes on in a new page from the case after it.
 */
async function runShard(browser, origin, profile, shard, referenceShard, caseFilter) {
  const reference = WebGPUProfiles.includes(profile) ? SharedReference : profile;
  const shared =
    referenceShard &&
    (await referenceShard.then(
      () => true,
      () => false,
    ));
  const rows = [];
  const consoleErrors = [];
  let loadMs = 0;
  const pageStart = Date.now();
  let reloaded = false;
  let counted = false;
  for (let from = 0; ; ) {
    const page = await browser.newPage();
    const q = new URLSearchParams({ target: profile, reference, shard: `${shard}/${shards}` });
    if (from > 0) q.set('from', String(from));
    if (shared) q.set('shared', sharedUrl(origin));
    if (!args.includes('--no-debug')) q.set('debug', '1');
    if (caseFilter) q.set('filter', caseFilter);
    if (tags) q.set('tags', tags);
    if (args.includes('--draw-lists')) q.set('drawLists', '1');
    const loadStart = Date.now();
    // many pages load at once
    await page.goto(`${origin}/tests/build/visual/index.html?${q}`);
    loadMs += Date.now() - loadStart;
    const results = await page.results(profile, from === 0 && !counted);
    const ran = results.results.length;
    // a page stalled while loading (a busy GPU process): loaded again, once; its cases are not known to go on past
    if (results.loading && !reloaded) {
      reloaded = true;
      counted = !!results.counted;
      console.log(`  stalled ${profile} loading the page, loading it again`);
      await within(page.close(), UnresponsiveMs);
      continue;
    }
    if (results.stalled) {
      progress.done++;
      progress.failed++;
      console.log(`  stalled ${profile} ${results.stalled}${results.hung ? ' (the page hung)' : ''}, going on`);
      results.results.push({
        id: results.stalled,
        status: 'error',
        errors: [`stalled: no result within ${seconds(StallMs)}${results.hung ? ', the page hung' : ''}`],
      });
    }
    await collectRows(page, profile, results, rows);
    consoleErrors.push(...page.consoleErrors);
    await within(page.close(), UnresponsiveMs);
    if (!results.stalled || ran + 1 >= results.count) break;
    from += ran + 1;
  }
  return { rows, consoleErrors, loadMs, totalMs: Date.now() - pageStart };
}

/** The rows of a page's `results` added to `rows`, with their images stored, compared or written as the run asks */
async function collectRows(page, profile, results, rows) {
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
    // a hung page answers nothing: its rows have no images
    if (r.status !== 'skipped' && !results.hung) {
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

/**
 * The groups of `--groups`: the patterns of `caseFilter`, or the top-level groups of the cases (`geometry/**`, ...),
 * read from the case files; the on-demand shader variants only when named
 */
function getGroups(caseFilter) {
  if (caseFilter) return caseFilter.split(',');
  const dir = path.join(root, 'tests', 'src', 'visual', 'cases');
  const groups = new Set();
  for (const file of fs.readdirSync(dir)) {
    if (file === 'shader-variants.ts') continue;
    const source = fs.readFileSync(path.join(dir, file), 'utf8');
    for (const m of source.matchAll(/\bid: [`'"]([a-z0-9-]+)\//g)) groups.add(`${m[1]}/**`);
  }
  return [...groups].sort();
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
    if (args.includes('--groups')) console.log(`groups: ${getGroups(filter).join(', ')}`);
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

  if (!(await checkRenderers(browser, origin, profiles))) {
    await browser.close();
    server.close();
    process.exit(3);
  }

  const start = Date.now();
  progress.start = start;
  // the shared reference first, then the other WebGL profiles: the WebGPU shards wait for its images
  const order = (p) => (p === SharedReference ? 0 : WebGPUProfiles.includes(p) ? 2 : 1);
  const ordered = [...profiles].sort((a, b) => order(a) - order(b));
  const groups = args.includes('--groups') ? getGroups(filter) : [filter];
  const jobs = [];
  const shardResults = [];
  // a run far slower than its last (a busy or hung GPU process) is aborted early with what it saw, not after minutes
  const key = runKey(ordered, groups);
  const timings = readTimings();
  const last = timings[key];
  const budgetMs = last ? Math.max(MinBudgetMs, BudgetFactor * last.ms) : Number(value('--budget', '600')) * 1000;
  const checkBudget = async () => {
    const elapsed = Date.now() - start;
    if (elapsed < ProjectAfterMs || progress.done === 0) return;
    const total = progress.pagesLeft > 0 && last ? Math.max(progress.total, last.cases) : progress.total;
    const projected = (elapsed * total) / progress.done;
    if (Math.max(projected, elapsed) <= budgetMs) return;
    clearInterval(timer);
    aborting = true;
    const expected = last ? `, the last such run ${seconds(last.ms)}` : '';
    console.log(
      `aborted: ${progress.done}/${total} cases in ${seconds(elapsed)}, projected ${seconds(projected)}, over the budget of ` +
        `${seconds(budgetMs)}${expected} (${progress.failed} failed so far; --budget <s> where there is no last run)`,
    );
    printProgress();
    await within(browser.close(), UnresponsiveMs);
    server.close();
    fs.rmSync(sharedDir, { recursive: true, force: true });
    process.exit(2);
  };
  const timer = setInterval(() => {
    printProgress();
    void checkBudget();
  }, ProgressMs);
  try {
    for (const group of groups) {
      const groupStart = Date.now();
      const groupJobs = ordered.flatMap((profile) =>
        Array.from({ length: shards }, (_, shard) => ({ profile, shard })),
      );
      progress.pagesLeft += groupJobs.length;
      const referenceShards = new Map();
      const groupResults = await pool(
        groupJobs.map((j) => () => {
          const referenceShard =
            sharesReference && WebGPUProfiles.includes(j.profile) ? referenceShards.get(j.shard) : undefined;
          const run = runShard(browser, origin, j.profile, j.shard, referenceShard, group);
          if (sharesReference && j.profile === SharedReference) referenceShards.set(j.shard, run);
          return run;
        }),
        concurrency,
      );
      if (groups.length > 1) {
        const counts = ordered.map((profile) => {
          const rows = groupResults.filter((_, i) => groupJobs[i].profile === profile).flatMap((p) => p.rows);
          const bad = rows.filter((r) => r.status === 'fail' || r.status === 'error').length;
          return `${profile} ${rows.filter((r) => r.status !== 'skipped').length - bad} pass${bad ? `, ${bad} failed` : ''}`;
        });
        console.log(`${group}: ${counts.join('; ')} (${seconds(Date.now() - groupStart)})`);
      }
      jobs.push(...groupJobs);
      shardResults.push(...groupResults);
    }
  } finally {
    clearInterval(timer);
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
    // `--require-backends`: a profile whose backend is not available fails (CI: a run on software rendering that
    // silently skipped WebGPU would pass)
    const unavailable = rows.filter((r) => r.status === 'skipped' && /Backend not available/.test(r.reason ?? ''));
    if (requireBackends && unavailable.length > 0) {
      failed = true;
      console.log(`  ${profile}: backend not available for ${unavailable.length} cases (--require-backends)`);
    }
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
  writeTiming(key, Date.now() - start, progress.total);
  console.log(
    `${jobs.length} pages (${profiles.length} profiles x ${shards} shards${groups.length > 1 ? ` x ${groups.length} groups` : ''}) in ${((Date.now() - start) / 1000).toFixed(0)} s`,
  );
  process.exit(failed ? 1 : 0);
}

/**
 * The renderers the run's backends get, printed; false where one is a software renderer and `--allow-software` is not
 * given (the run would be many times slower than on the GPU, e.g. a browser over Remote Desktop)
 */
async function checkRenderers(browser, origin, profiles) {
  const page = await browser.newPage();
  let r;
  try {
    // a document of the run's origin (WebGPU needs a secure context: localhost is one)
    await page.goto(`${origin}/package.json`);
    r = await page.evaluate(async () => {
      const out = { webgl: 'none', webgpu: 'none', fallback: false };
      try {
        const canvas = document.createElement('canvas');
        const gl = canvas.getContext('webgl2') ?? canvas.getContext('webgl');
        const ext = gl?.getExtension('WEBGL_debug_renderer_info');
        if (gl) out.webgl = String(gl.getParameter(ext ? ext.UNMASKED_RENDERER_WEBGL : gl.RENDERER));
      } catch (e) {
        out.webgl = `error: ${e}`;
      }
      try {
        const adapter = await navigator.gpu?.requestAdapter();
        if (adapter) {
          const info = adapter.info ?? {};
          out.webgpu = [info.vendor, info.architecture, info.description].filter(Boolean).join(' ') || 'unknown';
          out.fallback = !!(info.isFallbackAdapter ?? adapter.isFallbackAdapter);
        }
      } catch (e) {
        out.webgpu = `error: ${e}`;
      }
      return out;
    });
  } finally {
    await within(page.close(), UnresponsiveMs);
  }
  const usesWebGL = profiles.some((p) => !WebGPUProfiles.includes(p));
  const usesWebGPU = profiles.some((p) => WebGPUProfiles.includes(p));
  const software = [];
  if (usesWebGL && SoftwareRenderer.test(r.webgl)) software.push(`WebGL on ${r.webgl}`);
  if (usesWebGPU && (r.fallback || SoftwareRenderer.test(r.webgpu))) {
    software.push(`WebGPU on ${r.webgpu}${r.fallback ? ' (a fallback adapter)' : ''}`);
  }
  console.log(`renderers: WebGL ${r.webgl}; WebGPU ${r.webgpu}${r.fallback ? ' (fallback adapter)' : ''}`);
  if (software.length === 0) return true;
  const what = `software rendering: ${software.join('; ')}`;
  if (allowSoftware) {
    console.log(`  ${what} (--allow-software)`);
    return true;
  }
  console.log(
    `  ${what}: not the GPU (a browser over Remote Desktop gets none), runs many times slower; ` +
      'run where the GPU is available, or pass --allow-software',
  );
  return false;
}

/** The key of a run's timing: what decides how long it takes */
function runKey(ordered, groups) {
  return JSON.stringify({ profiles: ordered, shards, groups, tags: tags ?? null, browser: browserName });
}

function readTimings() {
  try {
    return JSON.parse(fs.readFileSync(TimingFile, 'utf8'));
  } catch {
    return {};
  }
}

/** Records a finished run's time and case count for `--budget` (a run aborted records nothing) */
function writeTiming(key, ms, cases) {
  const timings = readTimings();
  timings[key] = { ms, cases };
  fs.mkdirSync(path.dirname(TimingFile), { recursive: true });
  fs.writeFileSync(TimingFile, `${JSON.stringify(timings, null, 2)}\n`);
}

main().catch((e) => {
  if (aborting) return;
  console.error(e);
  process.exit(1);
});
