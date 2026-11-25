#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 *
 * Memory profile of a running ModelServer (molstar#1816). Starts the built server (`npx tsc -b servers/model`) with
 * a remote source, sends a workload, waits past the cache TTL and reports the heap after a full GC at each step.
 * A healthy server returns to its baseline after the idle step; a leaking one keeps the structures. The idle step is
 * also read before the GC is forced (`no GC`): what a container's memory stats show while the server is left alone.
 * `--trace-idle` prints that reading every 10 s of the idle wait. `heap` is V8's `heapUsed`; `ext` is Node's `external`,
 * the memory outside the heap (array buffers, the structures' columns among them, are part of it: not added again).
 *
 *   node servers/model/scripts/profile-memory.mjs [--cycles 2] [--concurrency 8] [--ttl 5000] [--cache 512]
 *        [--heap 2048] [--encoding cif] [--format bcif] [--ids 1cbs,1tqn,...] [--snapshot out-dir] [--log]
 *        [--local data-dir] [--idle ms] [--trace-idle]
 *
 * By default the server fetches each entry from models.rcsb.org (bcif) or files.rcsb.org (cif), which exercises the
 * fetch path. `--local <dir>` downloads the entries once into `<dir>` and serves them from disk (`--sourceMap`): no
 * network noise in the numbers. Either way it needs the network: it is a profiling tool, not a test.
 */

import { spawn } from 'node:child_process';
import * as fs from 'node:fs';
import * as path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { parseArgs } from 'node:util';

const here = path.dirname(fileURLToPath(import.meta.url));
const serverBin = path.join(here, '..', 'bin', 'model-server.mjs');
const agent = pathToFileURL(path.join(here, 'profile-agent.mjs')).href;

// The issue's two entries (3j3q: ~2.4M atoms) and mid-size entries, some of them requested twice per pass so the
// cache sees hits as well as misses.
const DefaultIds = [
  '3j3q',
  '7arv',
  '1cbs',
  '1tqn',
  '4hhb',
  '1ubq',
  '1crn',
  '3pqr',
  '1tup',
  '2ocj',
  '1c26',
  '6lu7',
  '6vxx',
  '6m0j',
  '7bv2',
  '1bna',
  '1grm',
  '1kx5',
  '1jj2',
  '1ffk',
  '1aon',
  '5ire',
  '2rh1',
  '3sn6',
  '1gfl',
  '1ema',
  '1mbn',
  '1hrc',
  '2lyz',
  '1stp',
  '4ins',
  '1a3n',
  '1hho',
  '2hhb',
  '5t2c',
];

const { values: args } = parseArgs({
  options: {
    cycles: { type: 'string', default: '2' },
    concurrency: { type: 'string', default: '8' },
    ttl: { type: 'string', default: '5000' },
    idle: { type: 'string' },
    cache: { type: 'string', default: '512' },
    heap: { type: 'string', default: '2048' },
    encoding: { type: 'string', default: 'cif' },
    format: { type: 'string', default: 'bcif' },
    ids: { type: 'string' },
    port: { type: 'string', default: '13370' },
    snapshot: { type: 'string' },
    log: { type: 'boolean', default: false },
    local: { type: 'string' },
    'trace-idle': { type: 'boolean', default: false },
    'server-arg': { type: 'string', multiple: true, default: [] },
  },
});

const ids = args.ids ? args.ids.split(',').map((s) => s.trim()) : DefaultIds;
const ttl = +args.ttl;
const idleMs = args.idle ? +args.idle : ttl + 2000;
const port = +args.port;
const ext = args.format === 'cif' ? 'cif' : 'bcif';
const remoteUrl = (id) =>
  ext === 'cif' ? `https://files.rcsb.org/download/${id}.cif` : `https://models.rcsb.org/${id}.bcif`;
const source = args.local
  ? ['--sourceMap', 'pdb', path.join(path.resolve(args.local), `\${id.toLowerCase()}.${ext}`)]
  : ['--sourceMapUrl', 'pdb', remoteUrl('${id.toLowerCase()}'), ext];

/** Downloads the entries missing from `dir`. */
async function ensureLocalFiles(dir) {
  fs.mkdirSync(dir, { recursive: true });
  for (const id of new Set(ids.map((i) => i.toLowerCase()))) {
    const fn = path.join(dir, `${id}.${ext}`);
    if (fs.existsSync(fn)) continue;
    const res = await fetch(remoteUrl(id));
    if (!res.ok) {
      console.log(`  could not download ${id}: ${res.status}`);
      continue;
    }
    fs.writeFileSync(fn, Buffer.from(await res.arrayBuffer()));
  }
}

const MB = 1024 * 1024;
const fmt = (b) => (b / MB).toFixed(0).padStart(6);

function startServer() {
  const nodeArgs = ['--expose-gc', `--max-old-space-size=${args.heap}`, '--import', agent, serverBin];
  const serverArgs = [
    '--defaultPort',
    '' + port,
    '--cacheMaxSizeInBytes',
    '' + +args.cache * MB,
    '--cacheEntryTimeoutMs',
    '' + ttl,
    '--shutdownTimeoutMinutes',
    '0',
    ...source,
    '--defaultSource',
    'pdb',
    ...args['server-arg'],
  ];
  const child = spawn(process.execPath, [...nodeArgs, ...serverArgs], { stdio: ['ignore', 'pipe', 'pipe', 'ipc'] });
  const counters = { added: 0, expired: 0, evicted: 0, accessed: 0, errors: 0 };
  let stopping = false;
  let crashed = false;
  let resolveReady;
  const ready = new Promise((resolve) => (resolveReady = resolve));
  let buffered = '';
  const onLine = (line) => {
    if (args.log) console.log(`  [server] ${line}`);
    if (line.includes('is running on port')) resolveReady();
    const cache = /\[Cache\] \S+ (added|expired|evicted|accessed)\b/.exec(line);
    if (cache) counters[cache[1]]++;
    if (/error/i.test(line)) counters.errors++;
  };
  const onData = (d) => {
    buffered += d.toString();
    const lines = buffered.split('\n');
    buffered = lines.pop();
    lines.forEach(onLine);
  };
  child.stdout.on('data', onData);
  child.stderr.on('data', onData);
  child.on('exit', (code, signal) => {
    if (!stopping) {
      console.error(`Server exited unexpectedly (code ${code}, signal ${signal}).`);
      process.exitCode = 1;
      crashed = true;
      resolveReady();
    }
  });

  let nextId = 0;
  const pending = new Map();
  child.on('message', (msg) => {
    if (!msg || msg.__profile !== true) return;
    pending.get(msg.id)?.(msg);
    pending.delete(msg.id);
  });
  const request = (cmd, extra = {}) =>
    new Promise((resolve, reject) => {
      if (crashed) return reject(new Error('server crashed'));
      const id = nextId++;
      pending.set(id, resolve);
      child.send({ __profile: true, id, cmd, ...extra });
    });

  return {
    ready,
    counters,
    get crashed() {
      return crashed;
    },
    mem: (gc) => request('mem', { gc }),
    snapshot: (p) => request('snapshot', { path: p }),
    stop() {
      stopping = true;
      child.kill();
    },
  };
}

async function query(id) {
  const t = performance.now();
  try {
    const res = await fetch(`http://localhost:${port}/ModelServer/v1/${id}/full?encoding=${args.encoding}`);
    let bytes = 0;
    for await (const chunk of res.body) bytes += chunk.length;
    return { id, status: res.status, bytes, ms: performance.now() - t };
  } catch (e) {
    return { id, status: 'ERR ' + (e.cause?.code ?? e.message), bytes: 0, ms: performance.now() - t };
  }
}

/** Runs `ids` with at most `concurrency` requests in flight, sampling the server's (un-collected) heap meanwhile. */
async function pass(server, concurrency) {
  let peakHeap = 0;
  let peakExternal = 0;
  let peakRss = 0;
  let sampling = true;
  const sampler = (async () => {
    while (sampling && !server.crashed) {
      try {
        const { mem } = await server.mem(false);
        peakHeap = Math.max(peakHeap, mem.heapUsed);
        peakExternal = Math.max(peakExternal, mem.external);
        peakRss = Math.max(peakRss, mem.rss);
      } catch {
        return;
      }
      await new Promise((r) => setTimeout(r, 100));
    }
  })();

  const queue = [...ids];
  const results = [];
  const t = performance.now();
  await Promise.all(
    Array.from({ length: concurrency }, async () => {
      while (queue.length && !server.crashed) results.push(await query(queue.shift()));
    }),
  );
  sampling = false;
  await sampler;
  const failed = results.filter((r) => r.status !== 200);
  return {
    ms: performance.now() - t,
    peakHeap,
    peakExternal,
    peakRss,
    failed,
    bytes: results.reduce((n, r) => n + r.bytes, 0),
  };
}

/** Waits `idleMs` without requests; with `--trace-idle`, prints the un-collected memory every 10 s. */
async function idle(server) {
  const end = performance.now() + idleMs;
  const start = performance.now();
  while (performance.now() < end && !server.crashed) {
    await new Promise((r) => setTimeout(r, Math.min(10_000, end - performance.now())));
    if (!args['trace-idle'] || server.crashed) continue;
    const { mem } = await server.mem(false);
    console.log(
      `  idle ${((performance.now() - start) / 1000).toFixed(0).padStart(4)} s  heap ${fmt(mem.heapUsed)} MB  ` +
        `ext ${fmt(mem.external)} MB  rss ${fmt(mem.rss)} MB`,
    );
  }
}

async function main() {
  console.log(
    `ids: ${ids.length}, cycles: ${args.cycles}, concurrency: ${args.concurrency}, cache: ${args.cache} MB, ` +
      `ttl: ${ttl} ms, idle: ${idleMs} ms, heap limit: ${args.heap} MB, source: ${args.local ? 'local ' : ''}${ext}, ` +
      `output: ${args.encoding}`,
  );
  if (args.local) await ensureLocalFiles(path.resolve(args.local));
  const server = startServer();
  await server.ready;
  if (server.crashed) return;

  const rows = [];
  const record = async (step, peak, gc = true) => {
    const { mem } = await server.mem(gc);
    const c = server.counters;
    rows.push({
      step,
      heap: mem.heapUsed,
      external: mem.external,
      rss: mem.rss,
      peakHeap: peak?.peakHeap,
      peakExternal: peak?.peakExternal,
      peakRss: peak?.peakRss,
      cache: `+${c.added} -${c.expired + c.evicted} hit ${c.accessed}`,
      note: peak
        ? `${(peak.ms / 1000).toFixed(1)} s, ${fmt(peak.bytes).trim()} MB out, ${peak.failed.length} failed`
        : '',
    });
    const r = rows[rows.length - 1];
    console.log(
      `${r.step.padEnd(22)} heap ${fmt(r.heap)} MB  ext ${fmt(r.external)} MB  rss ${fmt(r.rss)} MB` +
        (peak ? `  | peak heap ${fmt(r.peakHeap)} MB  ext ${fmt(r.peakExternal)} MB  rss ${fmt(r.peakRss)} MB` : '') +
        `  | cache ${r.cache}  ${r.note}`,
    );
    for (const f of peak?.failed ?? []) console.log(`    failed: ${f.id} -> ${f.status}`);
  };

  try {
    await record('baseline');
    for (let cycle = 1; cycle <= +args.cycles && !server.crashed; cycle++) {
      await record(`c${cycle} sequential`, await pass(server, 1));
      if (server.crashed) break;
      await record(`c${cycle} concurrent x${args.concurrency}`, await pass(server, +args.concurrency));
      if (server.crashed) break;
      await idle(server);
      await record(`c${cycle} idle, no GC`, undefined, false);
      await record(`c${cycle} after idle`);
    }
    if (args.snapshot && !server.crashed) {
      const { file } = await server.snapshot(path.resolve(args.snapshot, `model-server-${Date.now()}.heapsnapshot`));
      console.log(`heap snapshot: ${file}`);
    }
  } catch (e) {
    console.error('' + e);
  } finally {
    server.stop();
  }

  const base = rows[0]?.heap ?? 0;
  const last = rows[rows.length - 1];
  if (last) {
    console.log(
      `\nheap after the last idle step: ${fmt(last.heap).trim()} MB (baseline ${fmt(base).trim()} MB, ` +
        `+${fmt(last.heap - base).trim()} MB retained)${server.crashed ? ' -- SERVER CRASHED' : ''}`,
    );
  }
}

main();
