/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Particle, ParticleList, ParticleRotations } from './particle-list.js';

/** a deferred list's page (`Particle.DeferredParticles`): its particles read together, the GPU cull's page */
export const DeferredPageBits = 13;
export const DeferredPageSize = 1 << DeferredPageBits;

/** the bytes a cache keeps of the pages read (`ParticlePageCache`), least recently used let go first */
export const DeferredPageCacheBytes = 256 * 1024 * 1024;

/** A page's particles as read: its positions' words (on its list's grid), its rotations */
export interface DeferredPage {
  readonly positions: Uint32Array;
  readonly rotations?: ParticleRotations;
}

/** reads particles `[first, first + count)` of a deferred list */
export type DeferredRead = (first: number, count: number, signal?: AbortSignal) => Promise<DeferredPage>;

function pageBytes(page: DeferredPage) {
  return page.positions.byteLength + (page.rotations?.byteLength ?? 0);
}

/**
 * The pages read of deferred lists (a pack opened by range: its templates'), each read once while in flight (its
 * readers await it), kept within its bytes, the least recently used let go first (read again when asked); a failed read
 * not kept. One a pack, not a module's: the readers of its lists (page stores, a loci's sphere) share one copy
 */
export class ParticlePageCache {
  /** the bytes of the pages kept */
  bytes = 0;
  /** the bytes read (every read: again where let go) */
  readBytes = 0;
  /** per list and page its page, by use (the least recently used first) */
  private readonly kept = new Map<number, DeferredPage>();
  private readonly reading = new Map<number, Promise<DeferredPage>>();
  private readonly ids = new WeakMap<ParticleList, number>();
  private nextId = 0;

  constructor(readonly maxBytes = DeferredPageCacheBytes) {}

  /** `list`'s deferred particles (`Particle.setDeferred`): read by `read`, by page through this cache */
  deferred(list: ParticleList, read: DeferredRead): Particle.DeferredParticles {
    return {
      read,
      getPage: (page) => this.get(list, page),
      readPage: (page) => this.read(list, page, read),
    };
  }

  /** page `page` of `list` where kept (made the most recently used), else undefined: no read begun */
  get(list: ParticleList, page: number): DeferredPage | undefined {
    const key = this.key(list, page);
    const kept = this.kept.get(key);
    if (!kept) return undefined;
    this.kept.delete(key);
    this.kept.set(key, kept);
    return kept;
  }

  /** page `page` of `list`: kept, in flight, or read by `read` (its particles' range, the list's last page shorter) */
  read(list: ParticleList, page: number, read: DeferredRead): Promise<DeferredPage> {
    const kept = this.get(list, page);
    if (kept) return Promise.resolve(kept);
    const key = this.key(list, page);
    const inFlight = this.reading.get(key);
    if (inFlight) return inFlight;
    const first = page * DeferredPageSize;
    const count = Math.min(DeferredPageSize, list.count - first);
    if (count <= 0) return Promise.reject(new Error(`particle list '${list.label}': no page ${page}`));
    const reading = read(first, count).then(
      (values) => {
        this.reading.delete(key);
        const bytes = pageBytes(values);
        this.readBytes += bytes;
        this.kept.set(key, values);
        this.bytes += bytes;
        this.trim();
        return values;
      },
      (error) => {
        this.reading.delete(key);
        throw error;
      },
    );
    this.reading.set(key, reading);
    return reading;
  }

  /** the least recently used pages let go past the cache's bytes */
  private trim() {
    for (const [key, page] of this.kept) {
      if (this.bytes <= this.maxBytes) break;
      this.kept.delete(key);
      this.bytes -= pageBytes(page);
    }
  }

  /** a list's page as one number: an id a list (2^21 of them) above its page (2^32 of them) */
  private key(list: ParticleList, page: number) {
    let id = this.ids.get(list);
    if (id === undefined) this.ids.set(list, (id = this.nextId++));
    return id * 2 ** 32 + page;
  }
}
