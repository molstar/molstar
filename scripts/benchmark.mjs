/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Benchmarks in a browser (Playwright): builds tests/src/benchmark, serves the repository and runs the cases of each
// backend in one page, sequentially and one page at a time (parallel pages would distort the times). The `webgl2`
// page also runs the CPU paths. Timings depend on the machine and its load: run on demand on an otherwise idle machine.
//
// usage: node scripts/benchmark.mjs [--backends webgl2,webgl2-no-ubo,webgl1,webgpu] [--paths cpu,webgl,webgpu] [--filter 'mc/volume/*']
//                                   [--runs 7] [--warmup 2] [--out <file.json>] [--compare <file.json>]
//                                   [--no-build] [--channel chrome] [--headed] [--profile <dir>] [--heap-profile <dir>]
//
//   --gpu-timing     GPU time per pass group of the render cases, per pass of the representation cases (timing mode,
//                    timestamp queries), printed per case
//   --strict-math off  WebGPU shader modules without strict math (`strictMath: false` in their descriptors, a Chrome
//                    developer option: launched with `--enable-webgpu-developer-features`), for comparisons with the
//                    default (strict) run
//   --draw-lists, --no-draw-lists  every render case with draw lists or without (`RendererProps.drawLists`, WebGPU
//                    render bundles; on by default)
//   --data-binding <bytes>  WebGPU: the theme layers' storage binding limit (`maxDataBindingSize`), lower than the
//                    device's, so that large items read their layers from several pages
//   --browser-arg <flag>  an extra browser flag (repeatable), e.g. to choose the GPU of a machine with two
//   --images <dir>   writes the canvas after the last run of each render case (<dir>/<backend>/<case id>.png)
//   --profile <dir>  records a V8 CPU profile of each backend's page (<dir>/<backend>.cpuprofile, open in the DevTools
//                    Performance panel); profile a single case (--filter) for readable results
//   --heap-profile <dir>  records a V8 sampling heap profile of each backend's page with the objects garbage collected
//                    since (<dir>/<backend>.heapprofile, open in the DevTools Memory panel): what is allocated, not only
//                    what stays

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2).filter((a) => a !== '--');
const value = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};

const backends = value('--backends', 'webgl2,webgpu').split(',');
const paths = value('--paths');
const filter = value('--filter');
const runs = value('--runs', '7');
const warmup = value('--warmup', '2');
const outFile = value('--out');
const compareFile = value('--compare');
const channel = value('--channel', 'chrome');
const profileDir = value('--profile');
const heapProfileDir = value('--heap-profile');
const imagesDir = value('--images');
const strictMath = value('--strict-math', 'on');
if (strictMath !== 'on' && strictMath !== 'off') throw new Error(`--strict-math takes on or off, not '${strictMath}'`);

const mime = { '.html': 'text/html', '.js': 'text/javascript', '.map': 'application/json' };

function serve() {
  const server = http.createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      const file = path.join(root, pathname);
      if (!file.startsWith(root)) throw new Error('outside of the repository');
      const data = await fs.promises.readFile(file);
      res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(data);
    } catch {
      res.writeHead(404).end();
    }
  });
  return new Promise((resolve) => server.listen(0, '127.0.0.1', () => resolve(server)));
}

const fmt = (ms) => (ms === undefined ? '-' : ms < 10 ? ms.toFixed(2) : ms < 100 ? ms.toFixed(1) : ms.toFixed(0));

/** One row per case and path: first run, median, min, stages, triangles, speedup against the CPU path */
function report(results, previous) {
  const byCase = new Map();
  for (const r of results) {
    if (!byCase.has(r.id)) byCase.set(r.id, []);
    byCase.get(r.id).push(r);
  }
  const old = new Map(
    (previous ?? []).map((r) => [r.backend ? `${r.id}|${r.path}|${r.backend}` : `${r.id}|${r.path}`, r]),
  );
  // rows by backend where a path runs on several (e.g. the WebGL paths of webgl2, webgl2-no-ubo and webgl1)
  const label = (r) => (r.path === 'cpu' || r.backend === undefined ? r.path : r.backend);
  for (const [id, rows] of byCase) {
    console.log(`\n${id}`);
    const cpu = rows.find((r) => r.path === 'cpu' && r.status === 'ok');
    for (const r of rows) {
      if (r.status !== 'ok') {
        console.log(`  ${label(r).padEnd(13)} ${r.status}: ${r.reason}`);
        continue;
      }
      const speedup = cpu && r !== cpu ? `  ×${(cpu.median / r.median).toFixed(1)} vs cpu` : '';
      const before = old.get(`${r.id}|${r.path}|${r.backend}`) ?? old.get(`${r.id}|${r.path}`);
      const change =
        before?.status === 'ok'
          ? `  (${r.median >= before.median ? '+' : ''}${(((r.median - before.median) / before.median) * 100).toFixed(0)}%)`
          : '';
      const stages = Object.entries(r.stages)
        .map(([label, ms]) => `${label} ${fmt(ms)}`)
        .join(', ');
      const triangles =
        r.drawCalls !== undefined
          ? r.drawCalls > 0
            ? `${r.drawCalls} draw calls`
            : 'draw calls counted in timing mode only'
          : cpu && r !== cpu && Math.abs(r.triangles - cpu.triangles) > 0.01 * cpu.triangles
            ? `${r.triangles} triangles (cpu ${cpu.triangles}!)`
            : `${r.triangles} triangles`;
      console.log(
        `  ${label(r).padEnd(13)} median ${fmt(r.median).padStart(6)} ms, min ${fmt(r.min)}, first ${fmt(r.first)}${speedup}${change}`,
      );
      for (const note of r.notes ?? []) console.log(`                ${note}`);
      if (r.gpuPasses) {
        for (const [path, ms] of Object.entries(r.gpuPasses)) console.log(`                gpu ${fmt(ms)} ms  ${path}`);
      }
      console.log(
        `                ${stages}; ${triangles}${r.bytes === undefined ? '' : `, ${(r.bytes / 1024 / 1024).toFixed(1)} MiB`}`,
      );
    }
  }
}

async function main() {
  if (!args.includes('--no-build')) {
    const r = spawnSync(process.execPath, ['scripts/esbuild/visual-tests.mjs', '--suite', 'benchmark', '--no-watch'], {
      cwd: root,
      stdio: 'inherit',
    });
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
      ...(strictMath === 'off' ? ['--enable-webgpu-developer-features'] : []),
      // garbage collection before each run, so earlier cases' garbage doesn't slow later CPU runs
      '--js-flags=--expose-gc',
      ...(process.platform === 'win32' ? ['--use-angle=d3d11'] : []),
      ...args.flatMap((a, i) => (args[i - 1] === '--browser-arg' ? [a] : [])),
    ],
  });

  const results = [];
  const info = {};
  try {
    for (const backend of backends) {
      const page = await browser.newPage();
      if (strictMath === 'off') {
        await page.addInitScript(() => {
          if (typeof GPUDevice === 'undefined') return;
          const create = GPUDevice.prototype.createShaderModule;
          GPUDevice.prototype.createShaderModule = function (descriptor) {
            return create.call(this, { ...descriptor, strictMath: false });
          };
        });
      }
      const q = new URLSearchParams({ backend, runs, warmup });
      if (paths) q.set('paths', paths);
      if (filter) q.set('filter', filter);
      if (args.includes('--gpu-timing')) q.set('gpuTiming', '1');
      if (args.includes('--draw-lists')) q.set('drawLists', '1');
      if (args.includes('--no-draw-lists')) q.set('drawLists', '0');
      if (value('--data-binding')) q.set('dataBinding', value('--data-binding'));
      if (imagesDir) q.set('images', '1');
      page.on('console', (m) => {
        if (m.type() === 'error' || m.type() === 'warning') console.log(`  [${backend}] ${m.text()}`);
      });
      const cdp = profileDir || heapProfileDir ? await page.context().newCDPSession(page) : undefined;
      if (cdp && profileDir) {
        await cdp.send('Profiler.enable');
        await cdp.send('Profiler.setSamplingInterval', { interval: 100 });
        await cdp.send('Profiler.start');
      }
      if (cdp && heapProfileDir) {
        await cdp.send('HeapProfiler.enable');
        await cdp.send('HeapProfiler.startSampling', {
          samplingInterval: 4096,
          includeObjectsCollectedByMajorGC: true,
          includeObjectsCollectedByMinorGC: true,
        });
      }
      await page.goto(`${origin}/tests/build/benchmark/index.html?${q}`);
      await page.waitForFunction(() => window.__benchmark?.done, null, { timeout: 60 * 60 * 1000, polling: 1000 });
      if (cdp && profileDir) {
        const { profile } = await cdp.send('Profiler.stop');
        fs.mkdirSync(profileDir, { recursive: true });
        fs.writeFileSync(path.join(profileDir, `${backend}.cpuprofile`), JSON.stringify(profile));
      }
      if (cdp && heapProfileDir) {
        const { profile } = await cdp.send('HeapProfiler.stopSampling');
        fs.mkdirSync(heapProfileDir, { recursive: true });
        fs.writeFileSync(path.join(heapProfileDir, `${backend}.heapprofile`), JSON.stringify(profile));
      }
      const state = await page.evaluate(() => window.__benchmark);
      if (state.error) console.log(`${backend}: ${state.error}`);
      info[backend] = state.info;
      for (const [id, url] of Object.entries(state.images ?? {})) {
        const file = path.join(imagesDir, backend, `${id.replaceAll('/', '_')}.png`);
        fs.mkdirSync(path.dirname(file), { recursive: true });
        fs.writeFileSync(file, Buffer.from(url.slice(url.indexOf(',') + 1), 'base64'));
      }
      console.log(`${backend}: ${state.info.renderer}`);
      results.push(...state.results.map((r) => ({ ...r, backend })));
      await page.close();
    }
  } finally {
    await browser.close();
    server.close();
  }

  const previous = compareFile ? JSON.parse(fs.readFileSync(compareFile, 'utf8')).results : undefined;
  report(results, previous);
  if (outFile) {
    fs.writeFileSync(outFile, JSON.stringify({ date: new Date().toISOString(), strictMath, info, results }, null, 2));
    console.log(`\nwrote ${outFile}`);
  }
}

await main();
