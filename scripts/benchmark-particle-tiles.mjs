/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The CPU of a particle tiles view following the camera, in the Prism viewer (Chrome, WebGPU, Playwright): a pack from
// `.benchmark/` opened as a dropped file is (a file input's, read by range), the camera along a path (still outside, an
// orbit outside, a dolly into the cell through the levels, an orbit inside, still inside), per segment the frames, the
// long tasks, the stages the code measures in timing mode (`UserTiming.span`: the pick, the step, the store, the cull)
// and the uploads (`GPUQueue.writeBuffer` and `writeTexture`: calls, bytes, time; by buffer, its label or usage). Needs the Prism viewer built (`pnpm build`, or
// `pnpm dev:apps -- prism-viewer` running: its build folder is served). Run on demand on an otherwise idle machine.
//
// usage: node scripts/benchmark-particle-tiles.mjs [--pack yeast_cell_atoms.tiles16.levels.mspk] [--segment-ms 4000]
//                                                   [--channel chrome] [--headed] [--out <file.json>]
//                                                   [--profile <segment> --profile-out <file.cpuprofile>] [--vsync]
//                                                   [--phone [--cpu-slowdown 4]] [--step-ms <ms>] [--budget <particles>]
//                                                   [--trace <segment> --trace-out <file.json>] [--load-ms 30000]
//                                                   [--max-particles <particles>]
//                                                   [--shots <dir>] [--timeline]
//
//   --profile <segment>  a V8 CPU profile of that segment (e.g. orbit-inside; open in DevTools' Performance panel),
//                        its functions by self time printed
//   --trace <segment>    a Chrome trace of that segment (every process: the page's, the GPU's; open in DevTools'
//                        Performance panel or Perfetto), its threads' busy time by top-level event printed
//   --timeline           per segment its frames of 33 ms or more in order: when (from the segment's start), how long,
//                        the stages begun in them

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
const pack = path.resolve(root, '.benchmark', value('--pack', 'yeast_cell_atoms.tiles16.levels.mspk'));
const segmentMs = Number(value('--segment-ms', '4000'));
const loadMs = Number(value('--load-ms', '30000'));
const shotsDir = value('--shots');
const timeline = args.includes('--timeline');
const channel = value('--channel', 'chrome');
const outFile = value('--out');
const profileSegment = value('--profile');
const profileOut = value('--profile-out', 'particle-tiles.cpuprofile');
const traceSegment = value('--trace');

/** a page main thread's top-level task by what it ran: the first of these found in it */
const TaskKinds = [
  ['FireAnimationFrame', 'animation frames'],
  ['TimerFire', 'timers'],
  ['SimpleWatcher::OnHandleReady', 'IPC (the reads) and its continuations'],
  ['Receive mojo message', 'IPC (the reads) and its continuations'],
  ['v8.callFunction', 'other script'],
  ['V8.GC_MC_INCREMENTAL', 'GC'],
  ['V8.GC_SCAVENGER_SCAVENGE_PARALLEL_PHASE', 'GC'],
];

/**
 * the page's main thread by frame (the intervals between its animation frames' callbacks): per frame its top-level
 * tasks' time by kind (`TaskKinds`, else other) and its idle time; the frames past `slowMs` and the others averaged
 */
function traceFrames(events, slowMs) {
  const names = new Map();
  for (const e of events) {
    if (e.ph === 'M' && e.name === 'thread_name') names.set(`${e.pid}:${e.tid}`, e.args.name);
  }
  // the renderer main thread with the most animation frames (the page's)
  const counts = new Map();
  for (const e of events) {
    if (e.name !== 'FireAnimationFrame' || e.ph !== 'X') continue;
    const key = `${e.pid}:${e.tid}`;
    if (names.get(key) === 'CrRendererMain') counts.set(key, (counts.get(key) ?? 0) + 1);
  }
  const main = [...counts].sort((a, b) => b[1] - a[1])[0]?.[0];
  if (!main) return undefined;
  const list = events.filter((e) => e.ph === 'X' && e.dur && `${e.pid}:${e.tid}` === main);
  list.sort((a, b) => a.ts - b.ts || b.dur - a.dur);
  const tasks = [];
  let top;
  for (const e of list) {
    if (!top || e.ts >= top.ts + top.dur) {
      top = { ts: e.ts, dur: e.dur, names: new Set([e.name]) };
      tasks.push(top);
    } else top.names.add(e.name);
  }
  const frameStarts = tasks.filter((t) => t.names.has('FireAnimationFrame')).map((t) => t.ts);
  const kindOf = (t) => TaskKinds.find(([name]) => t.names.has(name))?.[1] ?? 'other';
  const sum = { slow: { frames: 0, ms: 0, kinds: {} }, other: { frames: 0, ms: 0, kinds: {} } };
  let k = 0;
  for (let f = 0; f + 1 < frameStarts.length; ++f) {
    const a = frameStarts[f];
    const b = frameStarts[f + 1];
    const ms = (b - a) / 1000;
    const into = ms > slowMs ? sum.slow : sum.other;
    into.frames += 1;
    into.ms += ms;
    let busy = 0;
    while (k < tasks.length && tasks[k].ts < a) ++k;
    for (let j = k; j < tasks.length && tasks[j].ts < b; ++j) {
      const t = tasks[j];
      const d = Math.min(t.ts + t.dur, b) - t.ts;
      busy += d;
      const kind = kindOf(t);
      into.kinds[kind] = (into.kinds[kind] ?? 0) + d / 1000;
    }
    into.kinds.idle = (into.kinds.idle ?? 0) + Math.max(0, b - a - busy) / 1000;
  }
  return sum;
}
const traceOut = value('--trace-out', 'particle-tiles.trace.json');

/**
 * a trace's threads (busiest first) and per thread its top-level events' time (ms), complete events (`X`) not inside
 * another of the same thread: what a thread was busy with, nested events counted once
 */
function traceThreads(events, n) {
  const names = new Map();
  for (const e of events) {
    if (e.ph === 'M' && e.name === 'thread_name') names.set(`${e.pid}:${e.tid}`, e.args.name);
  }
  const byThread = new Map();
  for (const e of events) {
    if (e.ph !== 'X' || !e.dur) continue;
    const key = `${e.pid}:${e.tid}`;
    let list = byThread.get(key);
    if (!list) byThread.set(key, (list = []));
    list.push(e);
  }
  const out = [];
  for (const [key, list] of byThread) {
    list.sort((a, b) => a.ts - b.ts || b.dur - a.dur);
    const byName = new Map();
    let end = -Infinity;
    let total = 0;
    for (const e of list) {
      if (e.ts < end) continue;
      end = e.ts + e.dur;
      total += e.dur / 1000;
      byName.set(e.name, (byName.get(e.name) ?? 0) + e.dur / 1000);
    }
    const top = [...byName].sort((a, b) => b[1] - a[1]).slice(0, n);
    out.push({ thread: names.get(key) ?? key, total, top });
  }
  return out.sort((a, b) => b.total - a.total);
}

/** a CPU profile's functions by self time (ms), the most first */
function topSelf(profile, n) {
  const byNode = new Map(profile.nodes.map((node) => [node.id, node]));
  const self = new Map();
  const dt = profile.timeDeltas;
  for (let i = 0; i < profile.samples.length; ++i) {
    const node = byNode.get(profile.samples[i]);
    const { functionName, url, lineNumber } = node.callFrame;
    const key = `${functionName || '(anonymous)'} ${url.split('/').pop()}:${lineNumber + 1}`;
    self.set(key, (self.get(key) ?? 0) + (dt[i] ?? 0) / 1000);
  }
  return [...self.entries()].sort((a, b) => b[1] - a[1]).slice(0, n);
}

const mime = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json' };

/** the repository's files (the viewer's build) */
function serve() {
  const server = http.createServer((req, res) => {
    const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
    const file = path.join(root, pathname);
    if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      res.writeHead(404).end();
      return;
    }
    res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' });
    fs.createReadStream(file).pipe(res);
  });
  return new Promise((resolve) => server.listen(0, '127.0.0.1', () => resolve(server)));
}

/**
 * in the page, before it runs: `GPUQueue.writeBuffer` counted (calls, bytes, its own time), in all and by buffer (its
 * label, else its usage and size), `writeTexture` too
 */
function countUploads() {
  const uploads = { calls: 0, bytes: 0, ms: 0, byBuffer: {}, bySize: {} };
  globalThis.__uploads = uploads;
  // the shader modules created (when, their label: the timeline's frames those compiled in, a GPU process' compile)
  const modules = [];
  globalThis.__modules = modules;
  const device = globalThis.GPUDevice?.prototype;
  if (device) {
    const create = device.createShaderModule;
    device.createShaderModule = function (descriptor) {
      modules.push({ t: performance.now(), label: descriptor.label ?? '', length: descriptor.code?.length ?? 0 });
      return create.call(this, descriptor);
    };
  }
  const queue = globalThis.GPUQueue?.prototype;
  if (!queue) return;
  const count = (key, bytes, ms) => {
    uploads.ms += ms;
    uploads.calls += 1;
    uploads.bytes += bytes;
    const b = (uploads.byBuffer[key] ??= { calls: 0, bytes: 0, ms: 0 });
    b.calls += 1;
    b.bytes += bytes;
    b.ms += ms;
    // by the write's size: whether a call's cost goes with its bytes or not
    const size = bytes <= 4096 ? '≤ 4 KiB' : bytes <= 65536 ? '≤ 64 KiB' : bytes <= 1048576 ? '≤ 1 MiB' : '> 1 MiB';
    const z = (uploads.bySize[size] ??= { calls: 0, bytes: 0, ms: 0 });
    z.calls += 1;
    z.bytes += bytes;
    z.ms += ms;
  };
  const write = queue.writeBuffer;
  queue.writeBuffer = function (buffer, offset, data, dataOffset, size) {
    const t = performance.now();
    try {
      return write.call(this, buffer, offset, data, dataOffset, size);
    } finally {
      const bytesPer = data.BYTES_PER_ELEMENT ?? 1;
      const bytes = size !== undefined ? size * bytesPer : data.byteLength - (dataOffset ?? 0) * bytesPer;
      const key = buffer.label || `usage ${buffer.usage}, ${(buffer.size / 1024 / 1024).toFixed(1)} MiB`;
      count(key, bytes, performance.now() - t);
    }
  };
  const writeTexture = queue.writeTexture;
  queue.writeTexture = function (destination, data, layout, size) {
    const t = performance.now();
    try {
      return writeTexture.call(this, destination, data, layout, size);
    } finally {
      const tex = destination.texture;
      // the bytes written: the rows the size covers (a part's rows: the array passed whole, its offset in the layout)
      const height = size.height ?? size[1] ?? 1;
      const layers = size.depthOrArrayLayers ?? size[2] ?? 1;
      const bytes = layout.bytesPerRow ? layout.bytesPerRow * height * layers : (data.byteLength ?? 0);
      count(`texture ${tex.label || tex.format}`, Math.min(bytes, data.byteLength ?? bytes), performance.now() - t);
    }
  };
}

/**
 * in the page, once the viewer is up: the GPU timers' results (timing mode, timestamp queries) summed, a frame's
 * `Canvas3D.render` and its passes, not printed (the animation loop prints them every frame: its own cost)
 */
function collectGpuTimes() {
  const timer = globalThis.molstar.viewer.plugin.canvas3d.graphics.timer;
  const gpu = { frames: 0, ms: 0, max: 0, byPass: {}, byRoot: {}, byLeaf: {} };
  globalThis.__gpu = gpu;
  const resolve = timer.resolve.bind(timer);
  // the passes themselves (the results' leaves: compute passes as render passes), by label
  const leaves = (r) => {
    if (r.children.length === 0) {
      const t = r.gpuElapsed / 1e6;
      if (Number.isFinite(t)) gpu.byLeaf[r.label] = (gpu.byLeaf[r.label] ?? 0) + t;
      return;
    }
    for (const c of r.children) leaves(c);
  };
  timer.resolve = () => {
    for (const r of resolve()) {
      const ms = r.gpuElapsed / 1e6;
      if (!Number.isFinite(ms)) continue;
      // every root (a frame's `Canvas3D.render`, and what is timed outside it)
      gpu.byRoot[r.label] = (gpu.byRoot[r.label] ?? 0) + ms;
      leaves(r);
      if (r.label !== 'Canvas3D.render') continue;
      gpu.frames += 1;
      gpu.ms += ms;
      gpu.max = Math.max(gpu.max, ms);
      for (const c of r.children) {
        const t = c.gpuElapsed / 1e6;
        if (Number.isFinite(t)) gpu.byPass[c.label] = (gpu.byPass[c.label] ?? 0) + t;
      }
    }
    return [];
  };
}

/** in the page: the camera along a segment (`kind`), its frames, long tasks, measures, uploads and GPU times */
async function runSegment({ kind, ms }) {
  const plugin = globalThis.molstar.viewer.plugin;
  const canvas = plugin.canvas3d;
  const { center, radius } = canvas.boundingSphere;
  const c = [center[0], center[1], center[2]];
  const far = 2.6 * radius;
  const near = 0.15 * radius;
  const place = (distance, angle) =>
    canvas.camera.setState(
      {
        position: [c[0] + Math.sin(angle) * distance, c[1], c[2] + Math.cos(angle) * distance],
        target: c,
        up: [0, 1, 0],
      },
      0,
    );
  const raf = () => new Promise((r) => requestAnimationFrame(r));
  // the still segments where the moving ones begin and end
  if (kind === 'still-outside') place(far, 0);
  else if (kind === 'still-inside') place(near, 2.0);
  performance.clearMeasures();
  const uploads0 = structuredClone(globalThis.__uploads);
  const gpu0 = structuredClone(globalThis.__gpu);
  let longMs = 0;
  const observer = new PerformanceObserver((list) => {
    for (const e of list.getEntries()) longMs += e.duration;
  });
  observer.observe({ type: 'longtask' });
  const frames = [];
  /** each frame's interval: the stages within a slow one */
  const intervals = [];
  /** the frames' timestamps (their begin, the display's: a frame missed a gap of two) */
  const stamps = [];
  const t0 = performance.now();
  let last = t0;
  for (;;) {
    const now = performance.now();
    const s = Math.min(1, (now - t0) / ms);
    if (kind === 'orbit-outside') place(far, s * 1.0);
    else if (kind === 'dolly-in') place(far * Math.pow(near / far, s), 1.0);
    else if (kind === 'orbit-inside') place(near, 1.0 + s * 1.0);
    else canvas.requestDraw();
    stamps.push(await raf());
    const t = performance.now();
    frames.push(t - last);
    intervals.push([last, t]);
    last = t;
    if (t - t0 >= ms) break;
  }
  await new Promise((r) => setTimeout(r, 50));
  observer.disconnect();
  const stages = {};
  for (const m of performance.getEntriesByType('measure')) {
    const s = (stages[m.name] ??= { count: 0, ms: 0, max: 0 });
    s.count += 1;
    s.ms += m.duration;
    s.max = Math.max(s.max, m.duration);
  }
  // the frames past a display's frame and a fifth (missed: 20 ms at 60 Hz), the stages that began in them summed
  const measures = performance.getEntriesByType('measure');
  const slow = { frames: 0, ms: 0, stages: {} };
  // the frames of two display frames or more, in order: when, how long, their stages
  const long = [];
  for (const [a, b] of intervals) {
    if (b - a <= 20) continue;
    slow.frames += 1;
    slow.ms += b - a;
    const stages = {};
    for (const m of measures) {
      if (m.startTime < a || m.startTime >= b) continue;
      slow.stages[m.name] = (slow.stages[m.name] ?? 0) + m.duration;
      stages[m.name] = (stages[m.name] ?? 0) + m.duration;
    }
    if (b - a >= 33) {
      // the shader modules created in the frame or the one before (a frame's compile in the GPU process the next's wait)
      const modules = (globalThis.__modules ?? []).filter((m) => m.t >= a - 17 && m.t < b).map((m) => m.label);
      long.push({ at: a - t0, ms: b - a, stages, modules });
    }
  }
  // the display's frames missed: by the timestamps' gaps, in display intervals (the median gap: not the least, two
  // frames stamped about together a gap near 0, every other gap then thousands of intervals)
  const gaps = stamps.slice(1).map((t, i) => t - stamps[i]);
  const interval = [...gaps].sort((a, b) => a - b)[gaps.length >> 1] ?? 1;
  let dropped = 0;
  for (const gap of gaps) dropped += Math.max(0, Math.round(gap / interval) - 1);
  const u = globalThis.__uploads;
  const g = globalThis.__gpu;
  const sorted = [...frames].sort((a, b) => a - b);
  const q = (p) => sorted[Math.min(sorted.length - 1, Math.floor(p * sorted.length))] ?? 0;
  const result = {
    kind,
    ms: performance.now() - t0,
    frames: { count: frames.length, p50: q(0.5), p95: q(0.95), max: sorted[sorted.length - 1] ?? 0 },
    longMs,
    dropped,
    slow,
    long,
    stages,
    gpu: {
      frames: g.frames - gpu0.frames,
      ms: g.ms - gpu0.ms,
      max: g.max,
      byPass: Object.fromEntries(Object.entries(g.byPass).map(([k, v]) => [k, v - (gpu0.byPass[k] ?? 0)])),
      byRoot: Object.fromEntries(Object.entries(g.byRoot).map(([k, v]) => [k, v - (gpu0.byRoot[k] ?? 0)])),
      byLeaf: Object.fromEntries(Object.entries(g.byLeaf).map(([k, v]) => [k, v - (gpu0.byLeaf[k] ?? 0)])),
    },
    uploads: {
      calls: u.calls - uploads0.calls,
      bytes: u.bytes - uploads0.bytes,
      ms: u.ms - uploads0.ms,
      byBuffer: Object.fromEntries(
        Object.entries(u.byBuffer).map(([k, b]) => {
          const b0 = uploads0.byBuffer[k] ?? { calls: 0, bytes: 0, ms: 0 };
          return [k, { calls: b.calls - b0.calls, bytes: b.bytes - b0.bytes, ms: b.ms - b0.ms }];
        }),
      ),
      bySize: Object.fromEntries(
        Object.entries(u.bySize).map(([k, b]) => {
          const b0 = uploads0.bySize[k] ?? { calls: 0, bytes: 0, ms: 0 };
          return [k, { calls: b.calls - b0.calls, bytes: b.bytes - b0.bytes, ms: b.ms - b0.ms }];
        }),
      ),
    },
  };
  result.t0 = t0;
  return result;
}

/**
 * in the page, after a segment begun at `t0` (the still segment outside: the first tiles): the load's time until done,
 * its last step's end, the frames drawn on (not counted in the segment, nor profiled or traced) until no step for half
 * a second, at most `limitMs` (`--load-ms`, 30 s)
 */
async function waitForLoad({ t0, limitMs }) {
  const canvas = globalThis.molstar.viewer.plugin.canvas3d;
  const lastStep = () => {
    const steps = performance.getEntriesByName('tiles: step', 'measure');
    const e = steps[steps.length - 1];
    return e ? e.startTime + e.duration : t0;
  };
  for (;;) {
    canvas.requestDraw();
    await new Promise((r) => requestAnimationFrame(r));
    const now = performance.now();
    if (now - lastStep() > 500 || now - t0 > limitMs) break;
  }
  const end = lastStep();
  // what the steps past the segment did: their count and the stages in them (a view that never settles)
  const after = {};
  for (const m of performance.getEntriesByType('measure')) {
    if (m.startTime < t0 + 4000) continue;
    after[m.name] = (after[m.name] ?? 0) + 1;
  }
  return { ms: end - t0, done: performance.now() - end > 500, after };
}

const fmt = (ms) => (ms < 10 ? ms.toFixed(2) : ms < 100 ? ms.toFixed(1) : ms.toFixed(0));

function report(segments) {
  for (const s of segments) {
    const { frames, uploads } = s;
    console.log(
      `\n${s.kind}: ${frames.count} frames, p50 ${fmt(frames.p50)} ms, p95 ${fmt(frames.p95)}, max ${fmt(frames.max)}; long tasks ${fmt(s.longMs)} ms; display frames missed ${s.dropped}`,
    );
    console.log(
      `  uploads: ${uploads.calls} calls, ${(uploads.bytes / 1024 / 1024).toFixed(1)} MiB, ${fmt(uploads.ms)} ms`,
    );
    if (s.particles !== undefined) console.log(`  particles shown: ${s.particles.toLocaleString('en')}`);
    if (s.load) {
      console.log(
        `  load: steps until ${fmt(s.load.ms)} ms${s.load.done ? '' : ` (not done at ${fmt(loadMs / 1000)} s)`}`,
      );
      if (s.load.after) console.log(`    after the segment, stages by count: ${JSON.stringify(s.load.after)}`);
    }
    if (s.slow?.frames > 0) {
      const parts = Object.entries(s.slow.stages)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 6)
        .map(([k, v]) => `${k} ${fmt(v / s.slow.frames)}`)
        .join(', ');
      console.log(
        `  slow frames (> 20 ms): ${s.slow.frames}, ${fmt(s.slow.ms / s.slow.frames)} ms each; their stages a frame: ${parts}`,
      );
    }
    if (timeline) {
      for (const f of s.long ?? []) {
        const parts = Object.entries(f.stages)
          .sort((a, b) => b[1] - a[1])
          .slice(0, 5)
          .map(([k, v]) => `${k} ${fmt(v)}`)
          .join(', ');
        console.log(`    at ${fmt(f.at)} ms: ${fmt(f.ms)} ms (${parts})`);
        if (f.modules?.length) {
          const by = {};
          for (const l of f.modules) by[l] = (by[l] ?? 0) + 1;
          const list = Object.entries(by)
            .map(([l, n]) => (n > 1 ? `${l} ×${n}` : l))
            .join(', ');
          console.log(`      shader modules ${f.modules.length}: ${list}`);
        }
      }
    }
    const { gpu } = s;
    if (gpu.frames > 0) {
      const passes = Object.entries(gpu.byPass)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 4)
        .map(([k, v]) => `${k} ${fmt(v / gpu.frames)}`)
        .join(', ');
      console.log(`  gpu: ${gpu.frames} frames timed, ${fmt(gpu.ms / gpu.frames)} ms a frame (${passes})`);
    }
    // all the GPU's timed work of the segment: by root (what is timed outside a frame's render too) and by pass
    const total = (o) =>
      Object.entries(o ?? {})
        .filter(([, v]) => v > 0.05)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 8)
        .map(([k, v]) => `${k} ${fmt(v)}`)
        .join(', ');
    if (gpu.byRoot) console.log(`    gpu by root (ms in the segment): ${total(gpu.byRoot)}`);
    if (gpu.byLeaf) console.log(`    gpu by pass (ms in the segment): ${total(gpu.byLeaf)}`);
    const sizes = Object.entries(uploads.bySize ?? {})
      .filter(([, b]) => b.calls > 0)
      .map(([k, b]) => `${k}: ${b.calls} calls ${fmt(b.ms)} ms (${fmt((b.ms * 1000) / b.calls)} µs each)`)
      .join('; ');
    if (sizes) console.log(`    by size: ${sizes}`);
    const buffers = Object.entries(uploads.byBuffer ?? {})
      .filter(([, b]) => b.calls > 0)
      .sort((a, b) => b[1].ms - a[1].ms)
      .slice(0, 8);
    for (const [key, b] of buffers) {
      console.log(
        `    ${key.slice(0, 40).padEnd(40)} ${fmt(b.ms).padStart(7)} ms, ${String(b.calls).padStart(6)} calls, ${(b.bytes / 1024 / 1024).toFixed(1).padStart(6)} MiB`,
      );
    }
    const names = Object.keys(s.stages).sort((a, b) => s.stages[b].ms - s.stages[a].ms);
    for (const name of names) {
      const t = s.stages[name];
      console.log(
        `  ${name.padEnd(20)} ${fmt(t.ms).padStart(7)} ms in ${String(t.count).padStart(4)}, ${fmt(t.ms / t.count)} each, max ${fmt(t.max)}`,
      );
    }
  }
}

async function main() {
  if (!fs.existsSync(pack)) throw new Error(`no pack ${pack}`);
  if (!fs.existsSync(path.join(root, 'apps/prism-viewer/build/index.html'))) {
    throw new Error('the Prism viewer is not built (pnpm build, or pnpm dev:apps -- prism-viewer)');
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
      // frames as fast as they go (the CPU's cost of a frame), or paced by the display (`--vsync`: a user's, the GPU
      // idle between frames: what an upload waits on otherwise)
      ...(args.includes('--vsync') ? [] : ['--disable-gpu-vsync', '--disable-frame-rate-limit']),
      ...(process.platform === 'win32' ? ['--use-angle=d3d11'] : []),
    ],
  });
  try {
    // a phone's (`--phone`): its viewport and pixel ratio (a typical Android's), its CPU slowed by the DevTools protocol
    // once the pack is open (the GPU the machine's: a phone's CPU alone)
    const phone = args.includes('--phone');
    const page = await browser.newPage(
      phone
        ? { viewport: { width: 412, height: 915 }, deviceScaleFactor: 2.625, isMobile: true, hasTouch: true }
        : { viewport: { width: 1920, height: 1080 } },
    );
    page.on('console', (m) => {
      if (m.type() === 'error') console.log(`  [page] ${m.text()}`);
    });
    await page.addInitScript(countUploads);
    const q = new URLSearchParams({ 'timing-mode': '1', 'graphics-backend': 'webgpu', splash: 'off' });
    await page.goto(`${origin}/apps/prism-viewer/build/index.html?${q}`);
    await page.waitForFunction(() => globalThis.molstar?.viewer?.plugin?.canvas3d, null, { timeout: 120_000 });

    // the pack as a dropped file is: a file input's (backed by the disk, read by range), opened as the drop opens it
    await page.evaluate(() => {
      const input = document.createElement('input');
      input.type = 'file';
      input.id = 'benchmark-pack';
      input.style.display = 'none';
      document.body.appendChild(input);
    });
    await page.setInputFiles('#benchmark-pack', pack);
    const t = Date.now();
    await page.evaluate(async () => {
      const file = document.getElementById('benchmark-pack').files[0];
      await globalThis.molstar.api.load.files({ files: [file] });
    });
    await page.waitForFunction(() => globalThis.molstar.viewer.plugin.canvas3d.boundingSphere.radius > 0, null, {
      timeout: 600_000,
    });
    console.log(`${path.basename(pack)} opened in ${((Date.now() - t) / 1000).toFixed(1)} s`);
    await page.evaluate(collectGpuTimes);

    // the views' steps by time (`stepMs` of the following views: their particles by how long the last steps took)
    // and their particles a step at most (`budget`)
    const stepMs = value('--step-ms');
    const budget = value('--budget');
    const maxParticles = value('--max-particles');
    if (stepMs !== undefined || budget !== undefined || maxParticles !== undefined) {
      const set = {
        ...(stepMs !== undefined && { stepMs: Number(stepMs) }),
        ...(budget !== undefined && { budget: Number(budget) }),
        ...(maxParticles !== undefined && { maxParticles: Number(maxParticles) }),
      };
      await page.evaluate(async (set) => {
        const plugin = globalThis.molstar.viewer.plugin;
        const state = plugin.state.data;
        const update = state.build();
        for (const cell of state.cells.values()) {
          const p = cell.transform.params;
          if (p?.tiles?.name !== 'camera') continue;
          update.to(cell).update({ ...p, tiles: { ...p.tiles, params: { ...p.tiles.params, ...set } } });
        }
        await plugin.runTask(state.updateTree(update));
      }, set);
      console.log(`steps: ${JSON.stringify(set)}`);
    }
    if (phone) {
      const rate = Number(value('--cpu-slowdown', '4'));
      const session = await page.context().newCDPSession(page);
      await session.send('Emulation.setCPUThrottlingRate', { rate });
      console.log(`phone: 412 × 915 at 2.625, the CPU ${rate} times slower`);
    }
    const segments = [];
    const cdp = profileSegment ? await page.context().newCDPSession(page) : undefined;
    for (const kind of ['still-outside', 'orbit-outside', 'dolly-in', 'orbit-inside', 'still-inside']) {
      const profiled = cdp && kind === profileSegment;
      const traced = kind === traceSegment;
      if (traced) {
        await browser.startTracing(page, {
          path: traceOut,
          screenshots: false,
          categories: [
            'devtools.timeline',
            'disabled-by-default-devtools.timeline',
            'disabled-by-default-devtools.timeline.frame',
            'toplevel',
            'v8',
            'blink',
            'gpu',
            'gpu.dawn',
            'disabled-by-default-gpu.dawn',
            'disabled-by-default-gpu.service',
            'viz',
            'cc',
            'ipc',
            'mojom',
          ],
        });
      }
      if (profiled) {
        await cdp.send('Profiler.enable');
        await cdp.send('Profiler.setSamplingInterval', { interval: 100 });
        await cdp.send('Profiler.start');
      }
      segments.push(await page.evaluate(runSegment, { kind, ms: segmentMs }));
      if (traced) {
        const buffer = await browser.stopTracing();
        const { traceEvents } = JSON.parse(buffer.toString());
        console.log(`
${kind}: Chrome trace in ${traceOut}, the threads' busy time (top-level events, ms):`);
        for (const t of traceThreads(traceEvents, 6).slice(0, 8)) {
          console.log(`  ${fmt(t.total).padStart(7)} ms  ${t.thread}`);
          for (const [name, ms] of t.top) console.log(`      ${fmt(ms).padStart(7)} ms  ${name}`);
        }
        // the page's main thread a frame: what a slow frame spends its time on, against the others
        const byFrame = traceFrames(traceEvents, 20);
        if (byFrame) {
          for (const [label, f] of [
            ['slow (> 20 ms)', byFrame.slow],
            ['the others', byFrame.other],
          ]) {
            if (f.frames === 0) continue;
            const parts = Object.entries(f.kinds)
              .sort((a, b) => b[1] - a[1])
              .map(([kind, ms]) => `${kind} ${fmt(ms / f.frames)}`)
              .join(', ');
            console.log(
              `  frames ${label}: ${f.frames}, ${fmt(f.ms / f.frames)} ms each; a frame's main thread: ${parts}`,
            );
          }
        }
      }
      if (profiled) {
        const { profile } = await cdp.send('Profiler.stop');
        fs.writeFileSync(profileOut, JSON.stringify(profile));
        console.log(`
${kind}: CPU profile in ${profileOut}, by self time:`);
        for (const [name, ms] of topSelf(profile, 25)) console.log(`  ${fmt(ms).padStart(7)} ms  ${name}`);
      }
      const segment = segments[segments.length - 1];
      if (kind === 'still-outside')
        segment.load = await page.evaluate(waitForLoad, { t0: segment.t0, limitMs: loadMs });
      // where the camera rests (outside after the load, inside): the views' particles, the canvas saved (`--shots`)
      if (kind === 'still-outside' || kind === 'still-inside') {
        segment.particles = await page.evaluate(() => {
          let n = 0;
          for (const cell of globalThis.molstar.viewer.plugin.state.data.cells.values()) {
            if (cell.transform.transformer.definition.name !== 'particle-tiles-view') continue;
            const list = cell.obj?.data;
            if (list) n += list.pages ? list.pages.liveCount : list.count;
          }
          return n;
        });
        if (shotsDir) {
          fs.mkdirSync(shotsDir, { recursive: true });
          await page
            .locator('canvas')
            .first()
            .screenshot({ path: path.join(shotsDir, `${kind}.png`) });
        }
      }
    }
    report(segments);
    // the ready tiles' positions and rotations: the bytes their arrays view and the buffers they keep (a pack read's
    // parts views of a merged read: a read kept while any of them is)
    const held = await page.evaluate(() => {
      const buffers = new Set();
      let views = 0;
      for (const cell of globalThis.molstar.viewer.plugin.state.data.cells.values()) {
        const loader = cell.obj?.data?.loader;
        if (!loader?.getReady) continue;
        for (const { id } of loader.getReady(0)) {
          const list = loader.getTile(id, 0);
          for (const a of [list?.coordinates?.data, list?.rotations]) {
            if (!a?.buffer) continue;
            views += a.byteLength;
            buffers.add(a.buffer);
          }
        }
      }
      let kept = 0;
      for (const b of buffers) kept += b.byteLength;
      // the loaded targets (their sites' bytes the catalog's estimates) and those no ready or loading tile uses
      let targets = { count: 0, bytes: 0, unused: 0, catalog: 0, targetBytes: 0, unusedBytes: 0 };
      for (const cell of globalThis.molstar.viewer.plugin.state.data.cells.values()) {
        const loader = cell.obj?.data?.loader;
        if (!loader?.getLoaded) continue;
        const loaded = loader.getLoaded();
        const unused = loader.getUnusedTargets();
        const catalog = cell.obj.data.set?.catalog?.targets;
        const bytesOf = (ids) => ids.reduce((sum, id) => sum + (catalog?.get(id)?.byteCount ?? 0), 0);
        targets = {
          count: loaded.targets,
          bytes: loaded.bytes,
          unused: unused.length,
          catalog: catalog?.size ?? 0,
          targetBytes: bytesOf([...loader.targets.keys()]),
          unusedBytes: bytesOf(unused),
        };
      }
      return { views, kept, buffers: buffers.size, targets };
    });
    console.log(
      `
tiles held: their positions and rotations ${(held.views / 2 ** 20).toFixed(1)} MiB, the buffers kept ${(held.kept / 2 ** 20).toFixed(1)} MiB (${held.buffers})`,
    );
    const tg = held.targets;
    console.log(
      `targets loaded: ${tg.count} of ${tg.catalog}, ${tg.unused} used by no tile; the loader's bytes (tiles and targets) ${(tg.bytes / 2 ** 20).toFixed(1)} MiB, the targets' ${(tg.targetBytes / 2 ** 20).toFixed(1)} MiB, the unused ${(tg.unusedBytes / 2 ** 20).toFixed(1)} MiB`,
    );
    // the page's JS heap after a collection: what the view keeps at the path's end (the tiles, the representations)
    const heapSession = await page.context().newCDPSession(page);
    await heapSession.send('HeapProfiler.collectGarbage');
    const { usedSize } = await heapSession.send('Runtime.getHeapUsage');
    console.log(`JS heap after a collection: ${(usedSize / 2 ** 20).toFixed(1)} MiB`);
    if (outFile) fs.writeFileSync(outFile, JSON.stringify({ pack: path.basename(pack), segments }, null, 2));
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
