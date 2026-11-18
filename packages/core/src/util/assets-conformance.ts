/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The asset manager's conformance suite (`_spec/assets.md` §9): `checkAssetManager` says how a manager breaks the
// interface's guarantees (`assets.ts`). Data and a function, no test framework: a manager's own test runs it and
// expects no failure, so a manager outside this package (an app's, derived from `DefaultAssetManager` or not) is held
// to the same suite. Downloads go through a fake `fetch` that the suite installs on `globalThis` while it runs and
// restores after: nothing reaches the network. A manager that downloads by other means than `fetch` is out of the
// suite's reach (it cannot be observed), and fails the cases that read the fake's requests.

import { Asset, AssetManager } from './assets.js';

const enc = new TextEncoder();
const dec = new TextDecoder();

/** `gzip('hello conformance hello conformance')`, written out so the suite needs no compressor. */
const GzipHello = new Uint8Array([
  31, 139, 8, 0, 0, 0, 0, 0, 0, 19, 203, 72, 205, 201, 201, 87, 72, 206, 207, 75, 203, 47, 202, 77, 204, 75, 78, 85,
  200, 64, 23, 1, 0, 85, 147, 125, 86, 35, 0, 0, 0,
]);
const HelloText = 'hello conformance hello conformance';

interface Route {
  body?: Uint8Array;
  status?: number;
  /** Where the route redirects to (`redirect: 'error'` fails the request). */
  redirectTo?: string;
  /** The body is sent in chunks of this many bytes, one per read of the stream (the whole body when absent). */
  chunk?: number;
  /** `Content-Length` is declared (default true). */
  declare?: boolean;
  /** The response waits for this before it is given. */
  gate?: Promise<void>;
}

interface Seen {
  url: string;
  init: RequestInit;
  /** Chunks of the body the reader pulled. */
  pulled: number;
  /** The request's signal aborted, or its body cancelled. */
  stopped: boolean;
}

/** A fake `fetch` over routes, recording what it was asked and what its readers took. */
class FakeFetch {
  readonly routes = new Map<string, Route>();
  readonly seen: Seen[] = [];

  serve = (url: string, route: Route | string) => {
    this.routes.set(url, typeof route === 'string' ? { body: enc.encode(route) } : route);
  };

  fetch = async (input: RequestInfo | URL, init: RequestInit = {}): Promise<Response> => {
    const url = String(input);
    const seen: Seen = { url, init, pulled: 0, stopped: false };
    this.seen.push(seen);
    const signal = init.signal ?? undefined;
    const aborted = () => {
      seen.stopped = true;
      return signal!.reason ?? new DOMException('aborted', 'AbortError');
    };
    if (signal?.aborted) throw aborted();

    let route = this.routes.get(url);
    if (route?.redirectTo !== undefined) {
      if (init.redirect === 'error') throw new TypeError('redirect mode is set to error');
      route = this.routes.get(route.redirectTo);
    }
    if (route?.gate) {
      await new Promise<void>((resolve, reject) => {
        signal?.addEventListener('abort', () => reject(aborted()), { once: true });
        route!.gate!.then(resolve, reject);
      });
    }
    if (!route) return new Response('not found', { status: 404 });
    if (route.status !== undefined && route.status >= 400) return new Response('', { status: route.status });

    const body = route.body ?? new Uint8Array(0);
    const size = route.chunk ?? Math.max(body.length, 1);
    let offset = 0;
    const stream = new ReadableStream<Uint8Array>(
      {
        pull(controller) {
          if (signal?.aborted) return controller.error(aborted());
          if (offset >= body.length) return controller.close();
          seen.pulled++;
          controller.enqueue(body.slice(offset, offset + size));
          offset += size;
        },
        cancel() {
          seen.stopped = true;
        },
      },
      { highWaterMark: 0 },
    );
    signal?.addEventListener('abort', () => {
      seen.stopped = true;
    });
    const headers: Record<string, string> = route.declare === false ? {} : { 'Content-Length': String(body.length) };
    return new Response(stream, { status: route.status ?? 200, headers });
  };
}

/** How long a read may take before it counts as not ending (`checkAssetManager`'s `timeoutMs`, for the run). */
let timeoutMs = 2000;

/** The reason a read failed, or `undefined` when it gave data (released). A read that does not end in time fails. */
async function failure(
  read: Promise<Asset.Wrapper<any>>,
  ms = timeoutMs,
): Promise<{ error?: unknown; data?: unknown }> {
  let timer: ReturnType<typeof setTimeout> | undefined;
  const timeout = new Promise<'timeout'>((resolve) => {
    timer = setTimeout(() => resolve('timeout'), ms);
  });
  try {
    const r = await Promise.race([read, timeout]);
    if (r === 'timeout') {
      read.then(
        (w) => w.dispose(),
        () => {},
      );
      return { error: new Error('the read did not end in time') };
    }
    const data = r.data;
    r.dispose();
    return { data };
  } catch (e) {
    return { error: e };
  } finally {
    if (timer !== undefined) clearTimeout(timer);
  }
}

const textOf = (data: unknown) => (typeof data === 'string' ? data : String(data));

type Case = (create: () => AssetManager, fake: FakeFetch, fail: (message: string) => void) => Promise<void>;

const cases: Record<string, Case> = {
  'a read gives the data': async (create, fake, fail) => {
    fake.serve('https://x.org/a', 'abc');
    const m = create();
    const r = await failure(m.resolve(Asset.Url('https://x.org/a'), 'string').run());
    if (r.error || textOf(r.data) !== 'abc') fail(`a text read gave ${r.error ?? JSON.stringify(r.data)}`);
    m.dispose();
  },

  'references are counted and released': async (create, fake, fail) => {
    fake.serve('https://x.org/a', 'abc');
    const m = create();
    const asset = Asset.Url('https://x.org/a');
    const first = await m.resolve(asset, 'string').run();
    const second = await m.resolve(asset, 'string').run();
    if (!m.has(asset)) fail('a resolved asset is not held');
    first.dispose();
    if (!m.has(asset)) fail('an asset with a reference left was let go');
    second.dispose();
    if (m.has(asset)) fail('an asset without a reference is still held');
    m.dispose();
  },

  'a failed read gives its reference back': async (create, fake, fail) => {
    fake.serve('https://x.org/gone', { status: 500 });
    const m = create();
    const asset = Asset.Url('https://x.org/gone');
    const r = await failure(m.resolve(asset, 'string').run());
    if (!r.error) fail('a download that failed gave data');
    if (m.has(asset)) fail('a failed read left the asset held');
    m.dispose();
  },

  'concurrent reads of a URL share one download': async (create, fake, fail) => {
    let open!: () => void;
    const gate = new Promise<void>((resolve) => (open = resolve));
    fake.serve('https://x.org/a', { body: enc.encode('abc'), gate });
    const m = create();
    const a = m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    const b = m.resolve(Asset.Url('https://x.org/a'), 'binary').run();
    await new Promise((resolve) => setTimeout(resolve, 20));
    open();
    const [ra, rb] = await Promise.all([failure(a), failure(b)]);
    if (ra.error || rb.error) fail(`concurrent reads failed: ${ra.error ?? rb.error}`);
    if (fake.seen.length !== 1) fail(`two concurrent reads of a URL made ${fake.seen.length} downloads`);
    m.dispose();
  },

  'a failed download is not cached': async (create, fake, fail) => {
    const m = create();
    fake.serve('https://x.org/a', { status: 503 });
    const first = await failure(m.resolve(Asset.Url('https://x.org/a'), 'string').run());
    if (!first.error) fail('a download that failed gave data');
    fake.serve('https://x.org/a', 'abc');
    const second = await failure(m.resolve(Asset.Url('https://x.org/a'), 'string').run());
    if (second.error || textOf(second.data) !== 'abc') fail(`the read after a failed download gave ${second.error}`);
    m.dispose();
  },

  'gzip bytes are inflated on request': async (create, fake, fail) => {
    fake.serve('https://x.org/a.gz', { body: GzipHello });
    const m = create();
    const inflated = await failure(
      AssetManager.resolveData(m, Asset.Url('https://x.org/a.gz'), 'string', { decompress: true }).run(),
    );
    if (inflated.error || textOf(inflated.data) !== HelloText)
      fail(`a decompressed read gave ${inflated.error ?? JSON.stringify(inflated.data)}`);
    const raw = await failure(m.resolve(Asset.Url('https://x.org/a.gz'), 'binary').run());
    if (
      raw.error ||
      !(raw.data instanceof Uint8Array) ||
      dec.decode(raw.data) === HelloText ||
      raw.data.length !== GzipHello.length
    )
      fail('a read without decompress did not give the bytes as received');
    m.dispose();
  },

  'credentials are honored or the read refused': async (create, fake, fail) => {
    for (const credentials of ['omit', 'include', 'same-origin'] as const) {
      const url = `https://x.org/${credentials}`;
      fake.serve(url, 'abc');
      const m = create();
      const r = await failure(m.resolve(Asset.Url(url), 'string', { credentials }).run());
      const sent = fake.seen.filter((s) => s.url === url);
      if (!r.error && !sent.some((s) => s.init.credentials === credentials))
        fail(`a read with credentials '${credentials}' was made with '${sent[0]?.init.credentials ?? 'the default'}'`);
      m.dispose();
    }
  },

  'redirects are honored or the read refused': async (create, fake, fail) => {
    fake.serve('https://x.org/moved', { redirectTo: 'https://x.org/target' });
    fake.serve('https://x.org/target', 'target');
    const m = create();
    const refused = await failure(m.resolve(Asset.Url('https://x.org/moved'), 'string', { redirect: 'error' }).run());
    if (!refused.error)
      fail(`a read with redirect 'error' followed a redirect and gave ${JSON.stringify(refused.data)}`);
    const followed = await failure(m.resolve(Asset.Url('https://x.org/moved'), 'string', { redirect: 'follow' }).run());
    if (!followed.error && textOf(followed.data) !== 'target')
      fail(`a read with redirect 'follow' gave ${followed.data}`);
    m.dispose();
  },

  'a signal cancels the read': async (create, fake, fail) => {
    const m = create();
    fake.serve('https://x.org/pre', 'abc');
    const done = new AbortController();
    done.abort();
    const before = await failure(m.resolve(Asset.Url('https://x.org/pre'), 'string', { signal: done.signal }).run());
    if (!before.error) fail('a read with a signal already aborted gave data');

    // aborted while the server has not answered: the read ends and the download is stopped
    const gate = new Promise<void>(() => {});
    fake.serve('https://x.org/slow', { body: enc.encode('abc'), gate });
    const controller = new AbortController();
    const slow = failure(m.resolve(Asset.Url('https://x.org/slow'), 'string', { signal: controller.signal }).run());
    await new Promise((resolve) => setTimeout(resolve, 20));
    controller.abort();
    const r = await slow;
    if (!r.error) fail('a read aborted by its signal gave data');
    else if (r.error instanceof Error && r.error.message === 'the read did not end in time')
      fail('a read aborted by its signal did not end');
    const sent = fake.seen.filter((s) => s.url === 'https://x.org/slow');
    if (sent.length && !sent.every((s) => s.stopped)) fail('a read aborted by its signal left its download running');

    // a waiter on a shared download with its own signal leaves alone: the download goes on for the other
    let open!: () => void;
    const shared = new Promise<void>((resolve) => (open = resolve));
    fake.serve('https://x.org/shared', { body: enc.encode('abc'), gate: shared });
    const owner = failure(m.resolve(Asset.Url('https://x.org/shared'), 'string').run());
    await new Promise((resolve) => setTimeout(resolve, 20));
    const waiterSignal = new AbortController();
    const waiter = failure(
      m.resolve(Asset.Url('https://x.org/shared'), 'string', { signal: waiterSignal.signal }).run(),
    );
    await new Promise((resolve) => setTimeout(resolve, 20));
    waiterSignal.abort();
    const w = await waiter;
    if (!w.error) fail('a waiting read aborted by its signal gave data');
    open();
    const o = await owner;
    if (o.error || textOf(o.data) !== 'abc')
      fail(`the read that was not aborted gave ${o.error ?? JSON.stringify(o.data)}`);
    m.dispose();
  },

  'maxBytes bounds a body while it streams': async (create, fake, fail) => {
    const big = new Uint8Array(100).fill(97);
    const m = create();
    // a body of undeclared length: the reader stops after the bound, it does not take the whole body first
    fake.serve('https://x.org/chunked', { body: big, chunk: 10, declare: false });
    const stopped = await failure(m.resolve(Asset.Url('https://x.org/chunked'), 'binary', { maxBytes: 50 }).run());
    if (!stopped.error) fail('a body over maxBytes was read');
    const chunked = fake.seen.find((s) => s.url === 'https://x.org/chunked');
    if (chunked && chunked.pulled > 7) fail(`a body over maxBytes was read on to ${chunked.pulled} of 10 chunks`);
    // a declared length over the bound: refused before the body
    fake.serve('https://x.org/declared', { body: big, chunk: 10 });
    const declared = await failure(m.resolve(Asset.Url('https://x.org/declared'), 'binary', { maxBytes: 50 }).run());
    if (!declared.error) fail('a declared length over maxBytes was read');
    const seen = fake.seen.find((s) => s.url === 'https://x.org/declared');
    if (seen && seen.pulled > 0) fail('a declared length over maxBytes was read before it was refused');
    // within the bound, and exactly at it
    for (const limit of [100, 1000]) {
      const ok = await failure(
        m.resolve(Asset.Url('https://x.org/chunked'), 'binary', { maxBytes: limit, store: false }).run(),
      );
      if (ok.error || !(ok.data instanceof Uint8Array) || ok.data.length !== 100)
        fail(`a read within maxBytes ${limit} gave ${ok.error ?? 'other data'}`);
    }
    // a bound is the read's: the same URL without one is read whole
    const whole = await failure(m.resolve(Asset.Url('https://x.org/chunked'), 'binary').run());
    if (whole.error) fail(`a read without maxBytes failed after one that was bound: ${whole.error}`);
    m.dispose();
  },

  'held data is not served to a read with other credentials': async (create, fake, fail) => {
    fake.serve('https://x.org/a', 'abc');
    const m = create();
    const url = 'https://x.org/a';
    const first = await m.resolve(Asset.Url(url), 'string', { credentials: 'include' }).run();
    // the asset is found as an app finds it, and read again without credentials, the first read still held
    const again = Asset.getUrlAsset(m, url);
    const r = await failure(m.resolve(again, 'string', { credentials: 'omit' }).run());
    const sent = fake.seen.filter((s) => s.url === url);
    if (!r.error && !(sent.length > 1 && sent[sent.length - 1].init.credentials === 'omit'))
      fail('a read with credentials omitted was served what a read with credentials included got');
    first.dispose();
    m.dispose();
  },

  'held data over maxBytes is refused': async (create, fake, fail) => {
    const big = new Uint8Array(100).fill(97);
    fake.serve('https://x.org/a', { body: big });
    const m = create();
    const asset = Asset.Url('https://x.org/a');
    const held = await m.resolve(asset, 'binary').run();
    const hit = await failure(m.resolve(asset, 'binary', { maxBytes: 50 }).run());
    if (!hit.error) fail('held data over maxBytes was read by the asset it is held for');
    // not held in memory any more (a manager that lets go keeps what it reads again from): the bound still holds
    m.letGoOfBytes?.(asset);
    const letGo = await failure(m.resolve(asset, 'binary', { maxBytes: 50 }).run());
    if (!letGo.error) fail('data a manager let go of was read over maxBytes');
    held.dispose();

    const file = Asset.File(new File([big], 'big.bin'));
    const fileRead = await m.resolve(file, 'binary').run();
    const fileHit = await failure(m.resolve(file, 'binary', { maxBytes: 50 }).run());
    if (!fileHit.error) fail('a held file over maxBytes was read');
    fileRead.dispose();
    m.dispose();
  },

  'maxBytes that is not a size refuses the read': async (create, fake, fail) => {
    fake.serve('https://x.org/a', 'abc');
    const m = create();
    for (const maxBytes of [Number.NaN, -1]) {
      const r = await failure(m.resolve(Asset.Url('https://x.org/a'), 'string', { maxBytes }).run());
      if (!r.error) fail(`a read with maxBytes ${maxBytes} gave data`);
    }
    m.dispose();
  },
};

/**
 * How a manager breaks the guarantees of `AssetManager` (none: it conforms): each failure names the case and what
 * went wrong. `create` makes a fresh manager for every case; `timeoutMs` is how long a read may take before it counts
 * as not ending (2000 by default). Runs on a fake `fetch`, installed on `globalThis` for the
 * run (do not run two at once, or other downloads in the same process meanwhile).
 */
export async function checkAssetManager(
  create: () => AssetManager,
  options: { timeoutMs?: number } = {},
): Promise<string[]> {
  const failures: string[] = [];
  const original = globalThis.fetch;
  timeoutMs = options.timeoutMs ?? 2000;
  for (const [name, run] of Object.entries(cases)) {
    const fake = new FakeFetch();
    globalThis.fetch = fake.fetch as typeof fetch;
    // every manager a case made is disposed after it, with what a read that did not end still holds
    const made: AssetManager[] = [];
    try {
      await run(
        () => {
          const m = create();
          made.push(m);
          return m;
        },
        fake,
        (message) => failures.push(`${name}: ${message}`),
      );
    } catch (e) {
      failures.push(`${name}: threw (${String(e)})`);
    } finally {
      globalThis.fetch = original;
      for (const m of made) {
        try {
          m.dispose();
        } catch {
          // a manager that fails to dispose is not what the case tests
        }
      }
    }
  }
  timeoutMs = 2000;
  return failures;
}

/** The names of the conformance cases. */
export const AssetManagerCases: readonly string[] = Object.keys(cases);
