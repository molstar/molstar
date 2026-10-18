/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Viewer smoke test in a browser (Playwright): builds the Viewer, serves the repository, and per graphics backend loads
// a local structure, waits for a rendered (non-uniform) canvas, checks for console and page errors (WebGPU errors are
// printed in debug mode), and hovers a grid of points; the hover labels of the backends must match, and on WebGPU
// hovering compiles no pipelines (the idle warm-up prepared picking and marking). At each point it also identifies
// (`canvas3d.asyncIdentify(…).wait()`, the read a click waits for): it must find an element where hovering did, and
// the same on all backends; `identify` in the same call must return that element too, or 'pending' (WebGPU).
//
// usage: node scripts/viewer-smoke.mjs [--backends webgl,webgpu] [--no-build] [--require-webgpu] [--channel chrome]
//                                      [--browser chromium] [--variant webgpu] [--headed]
//
//   --variant  checks a Viewer variant with one graphics backend (`molstar.variants`: `webgl`, `webgpu`), built with
//              `--variants` into `build/<variant>/`; backends `default` unless given
//   --browser  chromium (default, with --channel), firefox or webkit (Playwright's builds); where WebGPU is not
//              available, the `webgpu` run checks the WebGL backend the Viewer falls back to
//   --backends the Viewer's `graphics-backend` values (`webgpu`, `webgl`, `webgl2`, `webgl1`, `auto`) and `default` (no
//              `graphics-backend`, the Viewer's default); `auto`, `default` and `webgl` reported with the backend they
//              chose
//   --fallback then loses contexts for good: on a WebGPU page the device is destroyed with device requests failing
//              (falls back to WebGL2), then WebGL2 is lost without a restore (falls back to WebGL after the restore
//              timeout), then WebGL (nothing left: logged once, no further change); on a WebGL page, losses restored each
//              time until more than allowed within the loss window (falls back to WebGL); each fallback renders with the
//              coverage of the first render
//   --switch   then switches the backend at runtime (`plugin.setGraphicsBackends`, through those the Viewer was
//              created with: WebGL2, WebGPU, WebGL, WebGL2 with `auto`): each must be the backend asked for (WebGPU
//              skipped where unavailable), render with the coverage and hover the grid with the labels of the first
//              render; switched to a backend a second time, the context holds as many resources as after the first
//              switch to it (nothing of the previous contexts or the replaced representations left over)

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';
import { PNG } from 'pngjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const value = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};
/** a Viewer variant with one graphics backend (`molstar.variants`), built with `--variants` into `build/<variant>/` */
const variant = value('--variant');
const backends = value('--backends', variant ? 'default' : 'webgl,webgpu').split(',');
const channel = value('--channel', 'chrome');
const browserName = value('--browser', 'chromium');
const requireWebGPU = args.includes('--require-webgpu');
const switchBackends = args.includes('--switch');
const checkFallback = args.includes('--fallback');
const structure = '/data/examples/1grm_updated.cif';
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
  // the Viewer, to identify (`molstar` is a global of the Viewer bundle)
  await page.addInitScript(() => {
    let molstar;
    Object.defineProperty(window, 'molstar', {
      configurable: true,
      get: () => molstar,
      set: (value) => {
        const create = value.Viewer.create.bind(value.Viewer);
        value.Viewer.create = async (...args) => (window.__smokeViewer = await create(...args));
        molstar = value;
      },
    });
  });
  // the requested WebGPU devices, to simulate a device loss
  await page.addInitScript(() => {
    if (typeof GPUAdapter === 'undefined') return;
    const requestDevice = GPUAdapter.prototype.requestDevice;
    window.__smokeDevices = [];
    GPUAdapter.prototype.requestDevice = async function (...args) {
      // to fail the requests of a restore (`--fallback`)
      if (window.__smokeFailDevices) throw new Error('device requests fail (smoke test)');
      const device = await requestDevice.apply(this, args);
      window.__smokeDevices.push(device);
      return device;
    };
    // the render pipelines created, to check that the idle warm-up compiled those of the first hover
    window.__smokePipelines = [];
    for (const name of ['createRenderPipeline', 'createRenderPipelineAsync']) {
      const create = GPUDevice.prototype[name];
      GPUDevice.prototype[name] = function (descriptor) {
        window.__smokePipelines.push(descriptor.label);
        return create.call(this, descriptor);
      };
    }
  });

  const q = new URLSearchParams({
    ...(backend === 'default' ? {} : { 'graphics-backend': backend }),
    url: `${origin}${structure}`,
    'url-format': 'mmcif',
    'hide-controls': '1',
    'debug-mode': '1',
    'pixel-scale': '1',
  });
  await page.goto(`${origin}/apps/viewer/build/${variant ? `${variant}/` : ''}index.html?${q}`);
  const canvas = await page.waitForSelector('canvas', { timeout: 60000 });

  // the backend the Viewer ended up with (`webgpu`, `webgl2` or `webgl`; `webgpu` falls back to WebGL where unavailable)
  await page.waitForFunction(() => window.__smokeViewer?.plugin.canvas3d, null, { timeout: 60000 });
  const actual = await page.evaluate(() => window.__smokeViewer.plugin.canvas3d.graphics.capabilities.backend);
  // the backend asked for: `webgl1` is the backend id `webgl`, `webgl` any WebGL version
  const expected = backend === 'webgl1' ? 'webgl' : backend;
  let fallback;
  if (backend === 'auto' || backend === 'default' || (backend === 'webgl' && actual !== 'webgpu')) {
    // WebGPU where it has the features of WebGL2, else WebGL
    fallback = actual;
  } else if (actual !== expected) {
    if (backend !== 'webgpu' || requireWebGPU) {
      await page.close();
      return { backend, failures: [`expected the ${backend} backend, got ${actual}`] };
    }
    // WebGPU not available: the Viewer must work on the WebGL backend it fell back to
    fallback = actual;
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
  // hovering: on WebGPU, the idle warm-up compiled its pipelines (picking, marking) ahead
  await page.evaluate(() => window.__smokePipelines?.splice(0));
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.waitForTimeout(1500);

  /** the hover label at `x`, `y` once the asynchronous pick data is read and the label stable */
  const hoverLabel = async (x, y) => {
    await page.mouse.move(x, y);
    await page.waitForTimeout(300);
    let current = await label();
    for (let i = 0; i < 10; ++i) {
      await page.waitForTimeout(300);
      const next = await label();
      if (next === current) break;
      current = next;
    }
    return current;
  };
  const gridPoint = (ix, iy) => [
    Math.round(box.x + box.width * (0.35 + (0.3 * ix) / (grid - 1))),
    Math.round(box.y + box.height * (0.25 + (0.5 * iy) / (grid - 1))),
  ];

  const labels = [];
  const identified = [];
  /** grid points where `identify` returned other than the read pick data, without being 'pending' */
  const identifyNow = [];
  /** grid points where the picked position is far from the picked element (e.g. a wrong depth convention) */
  const misplaced = [];
  for (let iy = 0; iy < grid; ++iy) {
    for (let ix = 0; ix < grid; ++ix) {
      const [x, y] = gridPoint(ix, iy);
      labels.push(await hoverLabel(x, y));
      // as a click: the pick data once read, and `identify` in the same call ('pending' without synchronous reads)
      const [now, picked, offset] = await page.evaluate(
        async ([px, py]) => {
          const canvas3d = window.__smokeViewer?.plugin.canvas3d;
          const format = (p) => (p === 'pending' ? p : p ? `${p.id.objectId}/${p.id.instanceId}/${p.id.groupId}` : '-');
          const now = canvas3d?.identify([px, py]);
          const p = await canvas3d?.asyncIdentify([px, py])?.wait();
          // distance of the picked position outside the bounding sphere of the picked element
          let offset;
          if (p?.position) {
            const sphere = window.molstar.lib.loci.Loci.getBoundingSphere(canvas3d.getLoci(p.id).loci);
            if (sphere) {
              const [cx, cy, cz] = sphere.center;
              const [x, y, z] = p.position;
              offset = Math.hypot(x - cx, y - cy, z - cz) - sphere.radius;
            }
          }
          return [format(now), format(p), offset];
        },
        [x - box.x, y - box.y],
      );
      identified.push(picked);
      if (now !== 'pending' && now !== picked) identifyNow.push(`${identified.length - 1}: ${now} / ${picked}`);
      // representations extend a few Å beyond their elements (cartoon, surfaces)
      if (offset !== undefined && offset > 5) misplaced.push(`${identified.length - 1}: ${offset.toFixed(1)} Å`);
    }
  }
  const hits = labels.filter((l) => l !== '-').length;
  if (hits < 3) failures.push(`hover found ${hits} of ${labels.length} grid points`);
  if (actual === 'webgpu') {
    const compiled = await page.evaluate(() => window.__smokePipelines.splice(0));
    if (compiled.length) failures.push(`hovering compiled ${compiled.length} pipelines: ${compiled.join(', ')}`);
  }
  // switching to illumination for the first time: frames wait for its pipelines (`holdFramesUntilPrepared`), none of
  // them shows (almost) only the background meanwhile, i.e. less than half the coverage of the first render (each drawn
  // frame read back in the task that drew it)
  if (actual === 'webgpu') {
    const [frames, backgroundOnly, compiled] = await page.evaluate(async (covered) => {
      const canvas3d = window.__smokeViewer.plugin.canvas3d;
      const { graphics } = canvas3d;
      const illumination = canvas3d.props.illumination;
      const reads = [];
      // `didDraw` replays the last frame on subscribing, outside the task that drew it
      let replayed = false;
      const sub = canvas3d.didDraw.subscribe(() => {
        if (replayed) reads.push(graphics.readDrawingBuffer());
        replayed = true;
      });
      canvas3d.setProps({ illumination: { ...illumination, enabled: true } });
      const start = performance.now();
      while (performance.now() - start < 1500) await new Promise((resolve) => requestAnimationFrame(resolve));
      sub.unsubscribe();
      const images = await Promise.all(reads);
      // as `coverage`
      const coverage = (image) => {
        const d = image.data;
        let n = 0;
        for (let i = 0; i < d.length; i += 4) {
          if (Math.abs(d[i] - d[0]) + Math.abs(d[i + 1] - d[1]) + Math.abs(d[i + 2] - d[2]) > 30) n++;
        }
        return n / (d.length / 4);
      };
      canvas3d.setProps({ illumination });
      const backgroundOnly = images.filter((image) => coverage(image) < covered / 2).length;
      return [images.length, backgroundOnly, window.__smokePipelines.splice(0).length];
    }, covered);
    if (!compiled) failures.push('switching to illumination compiled no pipelines, already compiled before');
    if (!frames) failures.push('switching to illumination drew no frames');
    if (backgroundOnly) {
      failures.push(`switching to illumination showed only the background in ${backgroundOnly} of ${frames} frames`);
    }
    await page.waitForTimeout(500);
  }

  // props prepared ahead (`Canvas3D.prepare`): frames after setting them compile nothing
  if (actual === 'webgpu') {
    const scenarios = await page.evaluate(async () => {
      const canvas3d = window.__smokeViewer.plugin.canvas3d;
      const p = canvas3d.props;
      const frames = async (n) => {
        for (let i = 0; i < n; ++i) {
          canvas3d.requestDraw();
          await new Promise((resolve) => requestAnimationFrame(resolve));
        }
      };
      const results = [];
      for (const [name, props] of [
        [
          'multi-sample, transparent background',
          { multiSample: { ...p.multiSample, mode: 'on' }, transparentBackground: true },
        ],
      ]) {
        const restore = Object.fromEntries(Object.keys(props).map((k) => [k, p[k]]));
        await canvas3d.prepare({ props });
        window.__smokePipelines.splice(0);
        canvas3d.setProps(props);
        await frames(10);
        results.push([name, window.__smokePipelines.splice(0)]);
        canvas3d.setProps(restore);
        await frames(2);
      }
      return results;
    });
    for (const [name, compiled] of scenarios) {
      if (compiled.length)
        failures.push(`${name}, prepared, compiled ${compiled.length} pipelines: ${compiled.join(', ')}`);
    }
  }
  const missed = labels.map((l, i) => (l !== '-' && identified[i] === '-' ? i : -1)).filter((i) => i >= 0);
  if (missed.length) failures.push(`identify found nothing where hovering did, at ${missed.join(', ')}`);
  if (identifyNow.length) failures.push(`identify differs from the read pick data at ${identifyNow.join(', ')}`);
  if (misplaced.length) failures.push(`picked positions away from the picked element at ${misplaced.join(', ')}`);

  // context loss: the backend recovers (WebGPU: a new device; WebGL: the context the browser restores), recreates its
  // resources and the Viewer renders again
  let contextLoss;
  let lost;
  if (actual === 'webgpu') {
    lost = await page.evaluate(() => {
      const devices = window.__smokeDevices ?? [];
      devices.at(-1)?.destroy();
      return devices.length;
    });
  } else {
    lost = await canvas.evaluate((c) => {
      const ext = (c.getContext('webgl2') ?? c.getContext('webgl'))?.getExtension('WEBGL_lose_context');
      if (!ext) return false;
      ext.loseContext();
      setTimeout(() => ext.restoreContext(), 500);
      return true;
    });
  }
  const isRestored = () =>
    actual === 'webgpu'
      ? page.evaluate((n) => (window.__smokeDevices ?? []).length > n, lost)
      : canvas.evaluate((c) => !(c.getContext('webgl2') ?? c.getContext('webgl'))?.isContextLost());
  let restored = false;
  let restoredCoverage = 0;
  const lossStart = Date.now();
  while (lost !== false && Date.now() - lossStart < 20000) {
    await page.waitForTimeout(500);
    restored = await isRestored();
    if (!restored) continue;
    // hovering requests a draw
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2 + 1);
    const png = PNG.sync.read(await canvas.screenshot());
    restoredCoverage = isUniform(png) ? 0 : coverage(png);
    if (Math.abs(restoredCoverage - covered) < 0.01) break;
  }
  if (lost === false) failures.push('context loss: WEBGL_lose_context not available');
  else if (!restored) failures.push('context loss: not restored');
  else if (Math.abs(restoredCoverage - covered) >= 0.01) {
    failures.push(`context loss: coverage ${(restoredCoverage * 100).toFixed(2)}% after the restore`);
  } else contextLoss = `restored (coverage ${(restoredCoverage * 100).toFixed(1)}%)`;
  // the loss itself is reported in debug mode
  for (let i = errors.length - 1; i >= 0; --i) {
    if (/WebGPU device lost|CONTEXT_LOST_WEBGL|context lost/i.test(errors[i])) errors.splice(i, 1);
  }

  // falling back when a context is lost for good, without loops
  const fellBack = [];
  if (checkFallback) {
    await page.evaluate(() => {
      window.__smokeChanges = [];
      // backend changes only, restores on the same backend are expected before a fallback after repeated losses
      window.__smokeViewer.plugin.canvas3d.contextChanged.subscribe(({ reason }) => {
        if (reason === 'restore') return;
        window.__smokeChanges.push(`${reason} ${window.__smokeViewer.plugin.canvas3d.graphics.capabilities.backend}`);
      });
    });
    const backendNow = () => page.evaluate(() => window.__smokeViewer.plugin.canvas3d.graphics.capabilities.backend);
    const changes = () => page.evaluate(() => window.__smokeChanges.slice());
    /** waits up to `ms` for the backend to become `expected`, then for the coverage of the first render */
    const expectFallback = async (expected, ms, what) => {
      const start = Date.now();
      while (Date.now() - start < ms && (await backendNow()) !== expected) await page.waitForTimeout(250);
      const got = await backendNow();
      if (got !== expected) {
        failures.push(`fallback after ${what}: ${got}, expected ${expected}`);
        return false;
      }
      let c = 0;
      const coverageStart = Date.now();
      while (Date.now() - coverageStart < 20000) {
        await page.mouse.move(box.x + 2, box.y + 2 + (fellBack.length % 2));
        const png = PNG.sync.read(await (await page.waitForSelector('canvas')).screenshot());
        c = isUniform(png) ? 0 : coverage(png);
        if (Math.abs(c - covered) < 0.01) break;
        await page.waitForTimeout(500);
      }
      if (Math.abs(c - covered) >= 0.01) failures.push(`fallback to ${expected}: coverage ${(c * 100).toFixed(2)}%`);
      fellBack.push(`${expected} after ${what} (coverage ${(c * 100).toFixed(1)}%)`);
      return true;
    };
    const loseWebGL = (restore) =>
      page.evaluate((restore) => {
        const c = document.querySelector('canvas');
        const ext = (c.getContext('webgl2') ?? c.getContext('webgl'))?.getExtension('WEBGL_lose_context');
        if (!ext) return false;
        ext.loseContext();
        if (restore) setTimeout(() => ext.restoreContext(), 300);
        return true;
      }, restore);

    if (actual === 'webgpu') {
      await page.evaluate(() => {
        window.__smokeFailDevices = true;
        window.__smokeDevices.at(-1)?.destroy();
      });
      if (await expectFallback('webgl2', 10000, 'a failed WebGPU restore')) {
        await loseWebGL(false);
        if (await expectFallback('webgl', 15000, 'a WebGL2 loss without restore')) {
          // nothing left after WebGL: logged, no further change
          await loseWebGL(false);
          await page.waitForTimeout(8000);
          const all = await changes();
          if (all.length !== 2)
            failures.push(`fallback: ${all.length} backend changes, expected 2 (${all.join(', ')})`);
          if (!errors.some((e) => /no backend left to fall back to/.test(e))) {
            failures.push('fallback: no error logged when no backend was left');
          }
          fellBack.push(`nothing after webgl (${all.length} changes in all)`);
        }
      }
    } else if (actual === 'webgl2') {
      // losses restored each time, until more than allowed within the window (the context loss check above lost one)
      for (let i = 0; i < 4 && (await backendNow()) === 'webgl2'; ++i) {
        await loseWebGL(true);
        await page.waitForTimeout(1500);
      }
      await expectFallback('webgl', 10000, 'repeated losses');
      await page.waitForTimeout(4000);
      const all = await changes();
      if (all.length !== 1) failures.push(`fallback: ${all.length} backend changes, expected 1 (${all.join(', ')})`);
    }
    // the losses and fallbacks are reported
    for (let i = errors.length - 1; i >= 0; --i) {
      if (/Graphics context|WebGPU device lost|CONTEXT_LOST_WEBGL|context lost|device requests fail/i.test(errors[i])) {
        errors.splice(i, 1);
      }
    }
  }

  // runtime switching: the state is kept, the representations are created again with the new context
  const switched = [];
  /** resource counts after the first switch to a backend */
  const counts = new Map();
  if (switchBackends) {
    const available = await page.evaluate(() => window.__smokeViewer.plugin.graphicsBackends.map((b) => b.id));
    const order = ['webgl2', 'webgpu', 'webgl', 'webgl2'].filter((id) => available.includes(id));
    if (order.length < 2) failures.push(`switch: the Viewer was created with ${available.join(', ') || 'no'} backends`);
    for (const target of order.length > 1 ? order : []) {
      const result = await page.evaluate(async (id) => {
        const plugin = window.__smokeViewer.plugin;
        const ok = await plugin.setGraphicsBackends(plugin.graphicsBackends.filter((b) => b.id === id));
        return { ok, backend: plugin.canvas3d.graphics.capabilities.backend };
      }, target);
      if (!result.ok) {
        if (target === 'webgpu' && !requireWebGPU) {
          switched.push(`${target} unavailable`);
          continue;
        }
        failures.push(`switch to ${target}: not switched (${result.backend})`);
        continue;
      }
      if (result.backend !== target) failures.push(`switch to ${target}: got ${result.backend}`);
      const switchedCanvas = await page.waitForSelector('canvas');
      let switchedCoverage = 0;
      const switchStart = Date.now();
      while (Date.now() - switchStart < 30000) {
        await page.mouse.move(box.x + 2, box.y + 2 + (switched.length % 2));
        const png = PNG.sync.read(await switchedCanvas.screenshot());
        switchedCoverage = isUniform(png) ? 0 : coverage(png);
        if (Math.abs(switchedCoverage - covered) < 0.01) break;
        await page.waitForTimeout(500);
      }
      if (Math.abs(switchedCoverage - covered) >= 0.01) {
        failures.push(`switch to ${target}: coverage ${(switchedCoverage * 100).toFixed(2)}%`);
      }
      await page.waitForTimeout(1000);
      const differ = [];
      for (let iy = 0; iy < grid; ++iy) {
        for (let ix = 0; ix < grid; ++ix) {
          const i = iy * grid + ix;
          const l = await hoverLabel(...gridPoint(ix, iy));
          if (l !== labels[i]) differ.push(`${i}: ${labels[i]} / ${l}`);
        }
      }
      if (differ.length) failures.push(`switch to ${target}: hover labels differ at ${differ.join(', ')}`);
      const resources = await page.evaluate(() => ({
        ...window.__smokeViewer.plugin.canvas3d.graphics.stats.resourceCounts,
      }));
      const first = counts.get(target);
      if (!first) counts.set(target, resources);
      else {
        const changed = Object.keys(first)
          .filter((k) => first[k] !== resources[k])
          .map((k) => `${k} ${first[k]} -> ${resources[k]}`);
        if (changed.length) failures.push(`switch to ${target} again: resource counts differ: ${changed.join(', ')}`);
      }
      switched.push(`${target} (coverage ${(switchedCoverage * 100).toFixed(1)}%)`);
    }

    // geometry relying on capabilities: a molecular surface colored by chain has groups varying across triangles
    // where the backend has flat interpolation; after a switch to WebGL1 (none) its meshes must be made uniform
    if (available.includes('webgl') && available.length > 1) {
      const varying = await page.evaluate(async () => {
        const plugin = window.__smokeViewer.plugin;
        const structure = plugin.managers.structure.hierarchy.current.structures[0];
        await plugin.builders.structure.representation.addRepresentation(structure.cell, {
          type: 'molecular-surface',
          color: 'chain-id',
        });
        // the representations' render objects (Canvas3D's are committed later)
        const count = () => {
          let n = 0;
          plugin.state.data.cells.forEach((cell) => {
            for (const o of cell.obj?.data?.repr?.renderObjects ?? []) if (o.values.dVaryingGroup?.ref.value) n++;
          });
          return n;
        };
        const other = plugin.graphicsBackends.filter((b) => b.id !== 'webgl');
        await plugin.setGraphicsBackends(other);
        const before = count();
        await plugin.setGraphicsBackends(plugin.graphicsBackends.filter((b) => b.id === 'webgl'));
        for (let i = 0; i < 50 && count() > 0; ++i) await new Promise((r) => setTimeout(r, 100));
        return { before, after: count() };
      });
      if (varying.before === 0) failures.push('switch: the molecular surface has no varying groups before WebGL1');
      if (varying.after > 0) failures.push(`switch to webgl: ${varying.after} meshes keep varying groups`);
      else switched.push(`webgl with a molecular surface (${varying.before} varying meshes made uniform)`);
    }
  }

  if (errors.length) failures.push(...[...new Set(errors)].slice(0, 10).map((e) => `console: ${e}`));
  await page.close();
  return { backend, fallback, coverage: covered, hits, labels, identified, contextLoss, switched, fellBack, failures };
}

async function main() {
  if (!args.includes('--no-build')) {
    console.log('building the Viewer ...');
    const r = spawnSync(
      process.execPath,
      ['scripts/esbuild/app.mjs', 'apps/viewer', ...(variant ? ['--variants'] : [])],
      { cwd: root, stdio: 'inherit' },
    );
    if (r.status !== 0) process.exit(r.status ?? 1);
  }

  const server = await serve();
  const origin = `http://127.0.0.1:${server.address().port}`;
  const browserType = { chromium, firefox, webkit }[browserName];
  if (!browserType) throw new Error(`unknown browser '${browserName}'`);
  const browser = await browserType.launch({
    ...(browserName === 'chromium' ? { channel } : {}),
    headless: !args.includes('--headed'),
    args:
      browserName === 'chromium'
        ? [
            '--enable-gpu',
            '--ignore-gpu-blocklist',
            '--enable-unsafe-webgpu',
            ...(process.platform === 'win32' ? ['--use-angle=d3d11'] : []),
          ]
        : [],
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
    const differId = a.identified
      .map((p, j) => (p !== b.identified[j] ? `${j}: ${p} / ${b.identified[j]}` : ''))
      .filter((d) => d);
    if (differId.length) b.failures.push(`identified elements differ from ${a.backend} at ${differId.join(', ')}`);
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
      ? ` coverage ${(r.coverage * 100).toFixed(1)}%, hover hits ${r.hits}/${r.labels.length}${r.contextLoss ? `, context loss ${r.contextLoss}` : ''}${r.switched?.length ? `, switched to ${r.switched.join(', ')}` : ''}${r.fellBack?.length ? `, fell back to ${r.fellBack.join(', ')}` : ''}`
      : '';
    const backend = r.fallback
      ? `${r.backend} (${['auto', 'default', 'webgl'].includes(r.backend) ? 'chose' : 'fell back to'} ${r.fallback})`
      : r.backend;
    console.log(`${backend}: ${ok ? 'ok' : 'FAILED'}${info}`);
    for (const f of r.failures) console.log(`  - ${f}`);
  }
  process.exit(failed ? 1 : 0);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
