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
// usage: node scripts/visual-tests.mjs [--profiles webgl2,webgpu] [--changed] [--filter 'geometry/**'] [--tags smoke]
//                                      [--shards 8] [--concurrency 32] [--save <dir>] [--compare <dir>]
//                                      [--images <dir>] [--no-build] [--no-debug] [--channel chrome] [--headed]
//                                      [--angle d3d11] [--browser chromium] [--firefox <path>] [--headless]
//
//   --changed      profiles affected by the uncommitted changes: only the WebGPU profiles when only WebGPU code
//                  changed, only the WebGL profiles when only WebGL code changed, all otherwise
//   --save <dir>   stores the images of the WebGL profiles; --compare <dir> reports images that differ from them
//   --images <dir> writes target, reference and diff PNGs of failed cases (of all cases with --all-images)
//   --skipped      lists the skipped cases with the reason
//   --angle <name> ANGLE backend (`--use-angle`), default `d3d11` on Windows
//   --browser      `chromium` (Playwright, default) or `firefox`: the installed Firefox (Playwright's build has no
//                  WebGPU), driven over WebDriver BiDi, one window per page; `--firefox <path>` (or `FIREFOX`) names
//                  the executable. Firefox runs headed, as headless Firefox provides no WebGPU adapter; `--headless`
//                  runs it headless anyway (WebGL profiles only)
//
// Profiles run by default: `AllProfiles`; others only when named with --profiles (`webgl2-no-ubo`: WebGL2 with uniforms
// set one by one, the `uniformBlocks` option the render benchmark compares). Benchmarks run separately:
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

/** Profiles that uncommitted changes can affect */
function getChangedProfiles() {
  const git = (a) =>
    spawnSync('git', a, { cwd: root, encoding: 'utf8' })
      .stdout.split('\n')
      .filter((l) => l);
  const files = [...git(['diff', '--name-only', 'HEAD']), ...git(['ls-files', '--others', '--exclude-standard'])];
  const sources = files.filter((f) => /\.(ts|tsx|mjs|js)$/.test(f) && !f.startsWith('.v6/'));
  if (sources.length === 0) return [];
  const isWebGPU = (f) => f.includes('/gl/webgpu/') || /webgpu|wgsl/.test(path.basename(f));
  const isWebGL = (f) => f.includes('/gl/webgl/') && !isWebGPU(f);
  if (sources.every(isWebGPU)) return WebGPUProfiles;
  if (sources.every(isWebGL)) return WebGLProfiles;
  return AllProfiles;
}

const profiles = args.includes('--changed')
  ? getChangedProfiles()
  : value('--profiles', AllProfiles.join(',')).split(',');
const filter = value('--filter');
const tags = value('--tags');
const shards = Number(value('--shards', '8'));
const concurrency = Number(value('--concurrency', String(Math.max(2, Math.min(32, os.cpus().length)))));
const saveDir = value('--save');
const compareDir = value('--compare');
const imagesDir = value('--images');
const allImages = args.includes('--all-images');
const listSkipped = args.includes('--skipped');
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
        results: async () => {
          await page.waitForFunction(() => window.__visualTestResults?.done, null, {
            timeout: 60 * 60 * 1000,
            polling: 500,
          });
          return page.evaluate(() => window.__visualTestResults);
        },
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
        results: async () => {
          const start = Date.now();
          while (!(await evaluate(() => !!window.__visualTestResults?.done))) {
            if (Date.now() - start > 60 * 60 * 1000) throw new Error('timeout waiting for the results');
            await new Promise((resolve) => setTimeout(resolve, 500));
          }
          return evaluate(() => window.__visualTestResults);
        },
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

/** Runs one shard of a profile in its own page */
async function runShard(browser, origin, profile, shard) {
  const page = await browser.newPage();
  const { consoleErrors } = page;
  const reference = WebGPUProfiles.includes(profile) ? 'webgl2' : profile;
  const q = new URLSearchParams({ target: profile, reference, shard: `${shard}/${shards}` });
  if (!args.includes('--no-debug')) q.set('debug', '1');
  if (filter) q.set('filter', filter);
  if (tags) q.set('tags', tags);
  // many pages load at once
  await page.goto(`${origin}/tests/build/visual/index.html?${q}`);
  const results = await page.results();

  const rows = [];
  for (const r of results.results) {
    const row = { id: r.id, status: r.status, errors: r.errors, reason: r.reason, metrics: r.metrics };
    if (r.status !== 'skipped') {
      const image = (kind) => page.evaluate(([id, k]) => window.__visualTestImage(id, k), [r.id, kind]);
      // check cases have no image
      const target = (saveDir || compareDir) && !WebGPUProfiles.includes(profile) ? await image('target') : undefined;
      if (target) {
        const file = path.join(saveDir ?? compareDir, profile, `${fileName(r.id)}.txt`);
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
    for (const r of bad) {
      const m = r.metrics ? ` (${r.metrics.diffPixels} px differ)` : '';
      console.log(`  ${r.status} ${r.id}${m}${r.errors?.length ? `: ${r.errors.join(' | ').slice(0, 300)}` : ''}`);
    }
    if (compareDir) {
      console.log(`  changed vs ${compareDir}: ${changed.length ? changed.map((r) => r.id).join(', ') : 'none'}`);
      if (varied.length) console.log(`  changed, may vary between runs: ${varied.map((r) => r.id).join(', ')}`);
      if (added.length) console.log(`  new (not in ${compareDir}): ${added.map((r) => r.id).join(', ')}`);
    }
    const consoleErrors = [...new Set(parts.flatMap((p) => p.consoleErrors))].filter(
      (e) => !/Failed to load resource/.test(e),
    );
    for (const e of consoleErrors.slice(0, 10)) console.log(`  console: ${e.slice(0, 300)}`);
  }
  console.log(
    `${jobs.length} pages (${profiles.length} profiles x ${shards} shards) in ${((Date.now() - start) / 1000).toFixed(0)} s`,
  );
  process.exit(failed ? 1 : 0);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
