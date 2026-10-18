/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 2bit: a genome's letters (README, "2bit"). A sequence's size is read for the manifest, the rest of its record lazily,
// and its block lists by binary search, so a region of chr1 does not read its 200,000 mask blocks.

import { type DataPatch, StreamNotServed, type StreamProvider, type StreamRequest } from '@molstar/seq-core/format';
import { join, live, once, paged, type RangeReader, type Shared, share, view } from './ranges.js';

const SIGNATURE = 0x1a412743;
const BASES = 'TCAG'; // codes 0..3
const [UPPER_N, LOWER] = [78, 0x20];

interface Sequence {
  readonly name: string;
  readonly offset: number;
  readonly size: number;
}
/** A record's block lists (the file offsets of their start arrays) and its packed bases. */
interface Record2 {
  readonly nCount: number;
  readonly nAt: number;
  readonly maskCount: number;
  readonly maskAt: number;
  readonly dnaAt: number;
}

/**
 * At most this many sequences kept (after `only`: a draft assembly has 10^5 to 10^7 scaffolds, Picea abies' about 10M,
 * and a host asks for the few it shows; GI-12, LA-33), unless the host passes its own `maxSequences`.
 */
const MAX_SEQUENCES = 1_000_000;
/** The size reads in flight at once, one per sequence kept. */
const SIZE_READS = 16;

/**
 * A 2bit file as a `letters` provider: one coordinate system per sequence, `only` keeping the named ones, at most
 * `maxSequences` of them (a million by default).
 */
export function twoBitProvider(
  read: RangeReader,
  o: { only?: (name: string) => boolean; maxSequences?: number } = {},
): StreamProvider {
  const pages = paged(read);
  const records = new Map<string, Shared<Record2>>();
  const opened = once(async (signal: AbortSignal) => {
    const head = view(await read(0, 16, signal));
    const le = head.getUint32(0, true) === SIGNATURE;
    if (!le && head.getUint32(0, false) !== SIGNATURE) throw new Error('not a 2bit file: no 0x1A412743 signature');
    if (head.getUint32(4, le) !== 0) throw new Error('2bit version 1 (64-bit offsets) is not read');
    // The index: per sequence its name's length, name and offset; its size is the first word of its record.
    const count = head.getUint32(8, le);
    // The index is read as far as the file holds it: a count past it (a corrupt one, up to 2^32) ends at its end.
    const cut = (i: number) => new Error(`2bit: the index ends at sequence ${i} of the ${count} its header counts`);
    const kept: { name: string; offset: number }[] = [];
    const most = o.maxSequences ?? MAX_SEQUENCES;
    for (let at = 16, i = 0; i < count; i++) {
      const [n] = await pages(at, 1, signal);
      if (n === undefined) throw cut(i);
      const name = String.fromCharCode(...(await pages(at + 1, n, signal)));
      const offset = await pages(at + 1 + n, 4, signal);
      if (name.length < n || offset.length < 4) throw cut(i);
      at += 5 + n;
      if (!(o.only?.(name) ?? true)) continue;
      kept.push({ name, offset: view(offset).getUint32(0, le) });
      if (kept.length > most)
        throw new Error(`2bit: over ${most} sequences kept, the bound: name the ones to read (\`only\`)`);
    }
    // Each kept sequence's size, `SIZE_READS` reads at a time (GI-12: 100,000 at once from one manifest request),
    // set in the index's order.
    const sizes = new Array<number>(kept.length);
    let next = 0;
    const worker = async () => {
      for (let k = next++; k < kept.length; k = next++) {
        const s = kept[k] as { offset: number };
        sizes[k] = view(await read(s.offset, 4, signal)).getUint32(0, le);
      }
    };
    await Promise.all(Array.from({ length: Math.min(SIZE_READS, kept.length) }, worker));
    const sequences = new Map<string, Sequence>();
    for (let k = 0; k < kept.length; k++) {
      const { name, offset } = kept[k] as { name: string; offset: number };
      sequences.set(name, { name, offset, size: sizes[k] as number });
    }
    return { le, sequences };
  });
  /** A sequence's record, read once and shared by every request on it (`share`: CB-01). */
  const record = (s: Sequence, le: boolean, signal: AbortSignal) => {
    if (signal.aborted) return Promise.reject(signal.reason);
    let r = records.get(s.name);
    if (!live(r)) {
      const mine = share(async (own) => {
        const nCount = view(await read(s.offset + 4, 4, own)).getUint32(0, le);
        const maskAt = s.offset + 12 + 8 * nCount;
        const maskCount = view(await read(maskAt - 4, 4, own)).getUint32(0, le);
        return { nCount, nAt: s.offset + 8, maskCount, maskAt, dnaAt: maskAt + 8 * maskCount + 4 };
      });
      r = mine;
      records.set(s.name, mine);
      mine.promise.catch(() => records.get(s.name) === mine && records.delete(s.name));
    }
    return join(r, signal);
  };
  /** The blocks of a list (`count` starts, then `count` sizes) overlapping [a, b): a binary search over the sorted starts. */
  async function overlapping(le: boolean, at: number, count: number, a: number, b: number, signal: AbortSignal) {
    const word = async (i: number) => view(await pages(at + 4 * i, 4, signal)).getUint32(0, le);
    let [lo, hi] = [0, count];
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if ((await word(mid)) <= a) lo = mid + 1;
      else hi = mid;
    }
    const out: [number, number][] = [];
    for (let i = Math.max(0, lo - 1); i < count; ) {
      const n = Math.min(64, count - i);
      const starts = view(await pages(at + 4 * i, 4 * n, signal));
      const sizes = view(await pages(at + 4 * (count + i), 4 * n, signal));
      for (let j = 0; j < n; j++) {
        const s = starts.getUint32(4 * j, le);
        if (s >= b) return out;
        const e = s + sizes.getUint32(4 * j, le);
        if (e > a) out.push([s, e]);
      }
      i += n;
    }
    return out;
  }
  return {
    async fetch(r: StreamRequest, signal: AbortSignal): Promise<DataPatch> {
      const { le, sequences } = await opened(signal);
      const at = {
        source: r.source,
        system: r.system,
        start: r.start,
        stop: r.stop,
        row_start: 0,
        row_stop: 0,
        level: r.level,
      };
      if (r.system === null) {
        const systems = [...sequences.values()].map((s) => ({
          key: s.name,
          extent: [0, s.size] as [number, number],
          numbering: null,
          label: s.name,
        }));
        return { ...at, kind: 'manifest', manifest: { systems, levels: [], alphabet: 'dna' } };
      }
      if (r.level !== -1) throw new StreamNotServed(`a 2bit file serves letters (level -1), not level ${r.level}`);
      const s = sequences.get(r.system);
      if (!s) throw new Error(`2bit: no sequence ${r.system}`);
      const [start, stop] = [r.start, Math.min(r.stop, s.size)];
      if (start < 0 || start > stop) throw new Error(`2bit: ${s.name} [${r.start}, ${r.stop}) is outside 0..${s.size}`);
      const rec = await record(s, le, signal);
      const packed = await pages(rec.dnaAt + (start >> 2), ((stop + 3) >> 2) - (start >> 2), signal);
      const letters = new Uint8Array(stop - start);
      for (let i = 0; i < letters.length; i++) {
        const p = start + i;
        letters[i] = BASES.charCodeAt((packed[(p >> 2) - (start >> 2)] >> (6 - 2 * (p & 3))) & 3);
      }
      // N blocks first, then the soft mask (lower case, an N too), as `twoBitToFa` writes them.
      const [ns, mask] = await Promise.all([
        overlapping(le, rec.nAt, rec.nCount, start, stop, signal),
        overlapping(le, rec.maskAt, rec.maskCount, start, stop, signal),
      ]);
      for (const [a, b] of ns) letters.fill(UPPER_N, Math.max(a, start) - start, Math.min(b, stop) - start);
      for (const [a, b] of mask)
        for (let p = Math.max(a, start); p < Math.min(b, stop); p++) letters[p - start] |= LOWER;
      return { ...at, kind: 'letters', letters: new TextDecoder('latin1').decode(letters) };
    },
  };
}
