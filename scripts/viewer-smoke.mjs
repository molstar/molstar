/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Viewer smoke test in a browser (Playwright): builds the Viewer, serves the repository, and per graphics backend loads
// a local structure, waits for a rendered (non-uniform) canvas, checks for console and page errors (WebGPU errors are
// printed in debug mode), and hovers a grid of points; the hover labels of the backends must match.
//
// usage: node scripts/viewer-smoke.mjs [--backends webgl,webgpu] [--no-build] [--require-webgpu] [--channel chrome]
//                                      [--headed]

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';
import { PNG } from 'pngjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const value = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};
const backends = value('--backends', 'webgl,webgpu').split(',');
const channel = value('--channel', 'chrome');
const requireWebGPU = args.includes('--require-webgpu');
const structure = '/examples/1grm_updated.cif';
const grid = 5;

const mime = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.css': 'text/css',
  '.json': 'application/json',
  '.png': 'image/png',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
};

function serve() {
  const server = http.createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      const resolved = path.resolve(root, `.${pathname}`);
      if (!resolved.startsWith(root + path.sep)) return res.writeHead(403).end();
      const stat = await fs.promises.stat(resolved);
      const file = stat.isDirectory() ? path.join(resolved, 'index.html') : resolved;
      const body = await fs.promises.readFile(file);
      res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(body);
    } catch {
      res.writeHead(404).end('Not found');
    }
  });
  return new Promise((resolve) => server.listen(0, '127.0.0.1', () => resolve(server)));
}

function isUniform(png) {
  const { data } = png;
  for (let i = 4; i < data.length; i += 4) {
    if (data[i] !== data[0] || data[i + 1] !== data[1] || data[i + 2] !== data[2]) return false;
  }
  return true;
}

/** Fraction of pixels differing from the top-left pixel (the background) */
function coverage(png) {
  const { data } = png;
  let n = 0;
  for (let i = 0; i < data.length; i += 4) {
    if (Math.abs(data[i] - data[0]) + Math.abs(data[i + 1] - data[1]) + Math.abs(data[i + 2] - data[2]) > 30) n++;
  }
  return n / (data.length / 4);
}

async function check(browser, origin, backend) {
  const failures = [];
  const page = await browser.newPage({ viewport: { width: 800, height: 600 } });
  const errors = [];
  page.on('console', (m) => {
    if (m.type() === 'error') errors.push(m.text().slice(0, 500));
  });
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 500)));
  // the requested WebGPU devices, to simulate a device loss
  await page.addInitScript(() => {
    if (typeof GPUAdapter === 'undefined') return;
    const requestDevice = GPUAdapter.prototype.requestDevice;
    window.__smokeDevices = [];
    GPUAdapter.prototype.requestDevice = async function (...args) {
      const device = await requestDevice.apply(this, args);
      window.__smokeDevices.push(device);
      return device;
    };
  });

  const q = new URLSearchParams({
    'graphics-backend': backend,
    url: `${origin}${structure}`,
    'url-format': 'mmcif',
    'hide-controls': '1',
    'debug-mode': '1',
    'pixel-scale': '1',
  });
  await page.goto(`${origin}/apps/viewer/build/index.html?${q}`);
  const canvas = await page.waitForSelector('canvas', { timeout: 60000 });

  // the backend the Viewer ended up with (`webgpu` falls back to WebGL where unavailable)
  const actual = await page.evaluate(() =>
    document.querySelector('canvas')?.getContext('webgpu') ? 'webgpu' : 'webgl',
  );
  if (actual !== backend) {
    await page.close();
    if (backend === 'webgpu' && !requireWebGPU) return { backend, skipped: 'WebGPU not available' };
    return { backend, failures: [`expected the ${backend} backend, got ${actual}`] };
  }

  // wait for the structure to render
  let covered = 0;
  const start = Date.now();
  while (Date.now() - start < 60000) {
    const png = PNG.sync.read(await canvas.screenshot());
    covered = isUniform(png) ? 0 : coverage(png);
    if (covered > 0.02) break;
    await page.waitForTimeout(500);
  }
  if (covered <= 0.02) failures.push(`canvas not rendered (coverage ${(covered * 100).toFixed(2)}%)`);
  // let progressive passes (multi-sampling) and the picking pass settle
  await page.waitForTimeout(2000);

  const box = await canvas.boundingBox();
  const label = () =>
    page.evaluate(() => (document.querySelector('.msp-highlight-info')?.textContent ?? '').split(' | ').pop() || '-');
  // warm up picking (the first pick compiles the pick programs)
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.waitForTimeout(1500);

  const labels = [];
  for (let iy = 0; iy < grid; ++iy) {
    for (let ix = 0; ix < grid; ++ix) {
      const x = Math.round(box.x + box.width * (0.35 + (0.3 * ix) / (grid - 1)));
      const y = Math.round(box.y + box.height * (0.25 + (0.5 * iy) / (grid - 1)));
      await page.mouse.move(x, y);
      // the label updates once the asynchronous pick data is read: wait until it is stable
      await page.waitForTimeout(300);
      let current = await label();
      for (let i = 0; i < 10; ++i) {
        await page.waitForTimeout(300);
        const next = await label();
        if (next === current) break;
        current = next;
      }
      labels.push(current);
    }
  }
  const hits = labels.filter((l) => l !== '-').length;
  if (hits < 3) failures.push(`hover found ${hits} of ${labels.length} grid points`);

  let deviceLoss;
  if (backend === 'webgpu') {
    // device loss: the context requests a new device, recreates its resources and the Viewer renders again
    const lost = await page.evaluate(() => {
      const devices = window.__smokeDevices ?? [];
      devices.at(-1)?.destroy();
      return devices.length;
    });
    let restored = false;
    let restoredCoverage = 0;
    const lossStart = Date.now();
    while (Date.now() - lossStart < 20000) {
      await page.waitForTimeout(500);
      restored = await page.evaluate((n) => (window.__smokeDevices ?? []).length > n, lost);
      if (!restored) continue;
      // hovering requests a draw
      await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2 + 1);
      const png = PNG.sync.read(await canvas.screenshot());
      restoredCoverage = isUniform(png) ? 0 : coverage(png);
      if (Math.abs(restoredCoverage - covered) < 0.01) break;
    }
    if (!restored) failures.push('device loss: no new device requested');
    else if (Math.abs(restoredCoverage - covered) >= 0.01) {
      failures.push(`device loss: coverage ${(restoredCoverage * 100).toFixed(2)}% after the restore`);
    } else deviceLoss = `restored (coverage ${(restoredCoverage * 100).toFixed(1)}%)`;
    // the loss itself is reported in debug mode
    for (let i = errors.length - 1; i >= 0; --i) if (errors[i].includes('WebGPU device lost')) errors.splice(i, 1);
  }

  if (errors.length) failures.push(...[...new Set(errors)].slice(0, 10).map((e) => `console: ${e}`));
  await page.close();
  return { backend, coverage: covered, hits, labels, deviceLoss, failures };
}

async function main() {
  if (!args.includes('--no-build')) {
    console.log('building the Viewer ...');
    const r = spawnSync(process.execPath, ['scripts/esbuild/app.mjs', 'apps/viewer'], { cwd: root, stdio: 'inherit' });
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

  const results = [];
  try {
    for (const backend of backends) results.push(await check(browser, origin, backend));
  } finally {
    await browser.close();
    server.close();
  }

  // picking must agree between the backends
  const compared = results.filter((r) => r.labels);
  for (let i = 1; i < compared.length; ++i) {
    const [a, b] = [compared[0], compared[i]];
    const differ = a.labels.map((l, j) => (l !== b.labels[j] ? `${j}: ${l} / ${b.labels[j]}` : '')).filter((d) => d);
    if (differ.length) b.failures.push(`hover labels differ from ${a.backend} at ${differ.join(', ')}`);
  }

  let failed = false;
  for (const r of results) {
    if (r.skipped) {
      console.log(`${r.backend}: skipped (${r.skipped})`);
      continue;
    }
    const ok = r.failures.length === 0;
    failed ||= !ok;
    const info = r.labels
      ? ` coverage ${(r.coverage * 100).toFixed(1)}%, hover hits ${r.hits}/${r.labels.length}${r.deviceLoss ? `, device loss ${r.deviceLoss}` : ''}`
      : '';
    console.log(`${r.backend}: ${ok ? 'ok' : 'FAILED'}${info}`);
    for (const f of r.failures) console.log(`  - ${f}`);
  }
  process.exit(failed ? 1 : 0);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
