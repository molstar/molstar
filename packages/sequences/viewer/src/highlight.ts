/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Highlight channels (seqspec:interaction#channels), named by their owner (`links:<link>`, `set:<section>/<n>`,
// `results:<name>`, a host's own): the overlay draws their union, and clearing one never clears another. A target is
// loci or `Coordinates`, closed by its channel's `closure` (default `element`) whenever the channels light it
// (channels.ts); loci are carried to every new set (lineage.ts), `Coordinates` are held by key. Anything else is refused
// at `set` (several loci: their union, or a channel each); a carry that throws drops its channel alone.

import type { Closure } from '@molstar/seq-core/format';
import type { SeqSet } from '@molstar/seq-core/model';
import { carry, type SeqLociAny } from './lineage.js';
import { type Coordinates, isCoordinates } from './projection.js';
import { type SeqState, state } from './state.js';

/** `band`: a section's location set drawn full height in the main pass. */
export type HighlightStyle = 'highlight' | 'accent' | 'band';
export interface HighlightEntry {
  readonly target: SeqLociAny | Coordinates;
  readonly style: HighlightStyle;
  readonly closure: Closure;
  /** Lit in this section only (a section's `highlights`, a bookmark's step); null: every section. */
  readonly section: string | null;
}
export interface HighlightOptions {
  style?: HighlightStyle;
  closure?: Closure;
  section?: string | null;
}
export interface HighlightManager {
  readonly channels: SeqState<ReadonlyMap<string, HighlightEntry>>;
  set(channel: string, target: SeqLociAny | Coordinates | null, opts?: HighlightOptions): void;
  /** That channel only. */
  clear(channel: string): void;
  readonly version: SeqState<number>;
}

/** Loci (a kind and the set they are of), as a host passes them; arrays and other objects are not. */
const isLoci = (t: object) => !Array.isArray(t) && typeof (t as { kind?: unknown }).kind === 'string' && 'set' in t;

export class Highlight implements HighlightManager {
  readonly channels = state<ReadonlyMap<string, HighlightEntry>>(new Map());
  readonly version = state(0);
  private readonly failed: (message: string) => void;

  /** `failed`: a channel's carry threw (its channel is dropped; the others are carried). */
  constructor(failed: (message: string) => void = () => {}) {
    this.failed = failed;
  }

  set(channel: string, target: SeqLociAny | Coordinates | null, opts: HighlightOptions = {}) {
    if (!target) return this.clear(channel);
    if (!isCoordinates(target) && !isLoci(target))
      throw new TypeError(
        `highlight.set("${channel}"): the target is loci or Coordinates; for several loci, set their union (Loci.union) or one channel each`,
      );
    const entry = {
      target,
      style: opts.style ?? 'highlight',
      closure: opts.closure ?? 'element',
      section: opts.section ?? null,
    } as const;
    this.publish(new Map(this.channels.value).set(channel, entry));
  }

  clear(channel: string) {
    if (!this.channels.value.has(channel)) return;
    const next = new Map(this.channels.value);
    next.delete(channel);
    this.publish(next);
  }

  /** On every assembled set: each loci carried to it, a channel left with nothing dropped; null changes nothing. */
  update(set: SeqSet | null) {
    let next: Map<string, HighlightEntry> | null = null;
    for (const [channel, e] of set ? this.channels.value : []) {
      let target: HighlightEntry['target'] | null;
      try {
        target = isCoordinates(e.target) ? e.target : carry(e.target, set as SeqSet);
      } catch (err) {
        this.failed(
          `highlight channel "${channel}" could not be carried to the new set and is dropped: ${String(err)}`,
        );
        target = null;
      }
      if (target === e.target) continue;
      next ??= new Map(this.channels.value);
      if (target) next.set(channel, { ...e, target });
      else next.delete(channel);
    }
    if (next) this.publish(next);
  }

  private publish(next: ReadonlyMap<string, HighlightEntry>) {
    this.channels.set(next);
    this.version.set(this.version.value + 1);
  }

  dispose() {
    this.channels.complete();
    this.version.complete();
  }
}
