/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { gunzipSync, gzipSync } from 'node:zlib';

/** No personal detail in requests: a generic client name. */
const Headers = { 'User-Agent': 'molstar-bond-orders-benchmark' };
const GraphQL = 'https://data.rcsb.org/graphql';
/** The data API's cap on `entries(entry_ids:)` */
const EntriesPerRequest = 1000;

/**
 * The components each entry lists as bound (`rcsb_entry_info.nonpolymer_bound_components`, plan §3.2), metal
 * coordination included; `null` for an entry the data API does not serve. Cached in `cache/rcsb/bound.json`.
 */
export async function boundComponents(entries: string[], cacheDir: string): Promise<Map<string, string[] | null>> {
  const file = path.join(cacheDir, 'rcsb', 'bound.json');
  const cached: Record<string, string[] | null> = fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) : {};
  const missing = [...new Set(entries)].filter((e) => !(e in cached));
  for (let k = 0; k < missing.length; k += EntriesPerRequest) {
    const ids = missing.slice(k, k + EntriesPerRequest);
    const query = `{ entries(entry_ids: ${JSON.stringify(ids)}) { rcsb_id rcsb_entry_info { nonpolymer_bound_components } } }`;
    const response = await retry(() =>
      fetch(GraphQL, {
        method: 'POST',
        headers: { ...Headers, 'Content-Type': 'application/json' },
        body: JSON.stringify({ query }),
      }),
    );
    const data = (await response.json()) as {
      data?: { entries: { rcsb_id: string; rcsb_entry_info?: { nonpolymer_bound_components?: string[] } }[] };
      errors?: unknown;
    };
    if (!data.data) throw new Error(`${GraphQL}: ${JSON.stringify(data.errors)}`);
    for (const id of ids) cached[id] = null;
    for (const e of data.data.entries) cached[e.rcsb_id] = e.rcsb_entry_info?.nonpolymer_bound_components ?? [];
  }
  if (missing.length) {
    fs.mkdirSync(path.dirname(file), { recursive: true });
    fs.writeFileSync(file, JSON.stringify(cached));
  }
  return new Map(entries.map((e) => [e, cached[e] ?? null]));
}

/** Cached gzipped: BinaryCIF compresses to about a fifth (plan §3.2). */
export function entryFile(cacheDir: string, id: string) {
  return path.join(cacheDir, 'entries', `${id.toUpperCase()}.bcif.gz`);
}

/** The entry's static BinaryCIF, or `undefined` when it is not served (obsolete, withdrawn). */
export function readEntry(cacheDir: string, id: string): Uint8Array | undefined {
  const file = entryFile(cacheDir, id);
  return fs.existsSync(file) ? gunzipSync(fs.readFileSync(file)) : void 0;
}

/**
 * Downloads the entries not cached yet (`models.rcsb.org/<id>.bcif`, whole files, gzipped, plan §3.2), `concurrency`
 * at a time.
 * A 404 or 410 is recorded (`<id>.missing`) so it is not asked again.
 */
export async function downloadEntries(ids: string[], cacheDir: string, concurrency = 4) {
  const dir = path.join(cacheDir, 'entries');
  fs.mkdirSync(dir, { recursive: true });
  const todo = [...new Set(ids.map((i) => i.toUpperCase()))].filter(
    (id) => !fs.existsSync(entryFile(cacheDir, id)) && !fs.existsSync(path.join(dir, `${id}.missing`)),
  );
  let done = 0;
  const worker = async () => {
    for (let id = todo.pop(); id !== void 0; id = todo.pop()) {
      const url = `https://models.rcsb.org/${id.toLowerCase()}.bcif`;
      const response = await retry(() => fetch(url, { headers: Headers }), [404, 410]);
      if (response.status === 404 || response.status === 410) {
        fs.writeFileSync(path.join(dir, `${id}.missing`), String(response.status));
      } else {
        const partial = `${entryFile(cacheDir, id)}.partial`;
        fs.writeFileSync(partial, gzipSync(new Uint8Array(await response.arrayBuffer()), { level: 9 }));
        fs.renameSync(partial, entryFile(cacheDir, id));
      }
      if (++done % 500 === 0) console.log(`downloaded ${done} entries, ${todo.length} left`);
    }
  };
  await Promise.all(Array.from({ length: concurrency }, worker));
}

/** Retries on a network error, 429 and 5xx, backing off (1 s, 2 s, 4 s...); other statuses but `accept` throw. */
async function retry(request: () => Promise<Response>, accept: number[] = []): Promise<Response> {
  for (let attempt = 0; ; attempt++) {
    let response: Response | undefined;
    try {
      response = await request();
    } catch (e) {
      if (attempt >= 5) throw e;
    }
    if (response?.ok || (response && accept.includes(response.status))) return response;
    if (response && response.status !== 429 && response.status < 500) {
      throw new Error(`${response.url}: ${response.status}`);
    }
    if (attempt >= 5) throw new Error(`${response?.url}: ${response?.status}, after ${attempt} retries`);
    await new Promise((r) => setTimeout(r, 1000 * 2 ** attempt));
  }
}
