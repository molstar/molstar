/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A runner's retention, for the conformance probe (P10): a test-only contract, nothing more. */
export interface RetentionSnapshot {
  readonly entries: number;
  readonly segments: number;
  readonly bytes: number;
  /** Items evicted to the bounded spill cache (`cacheBytes`). */
  readonly spillBytes: number;
}

export interface Inspectable {
  retention(): RetentionSnapshot;
}
