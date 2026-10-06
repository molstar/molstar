/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';

/**
 * Combines entries into one entry, for modules whose entry includes the entries of what they use. Lists are
 * concatenated in order and a provider listed by several entries appears once, so the result registers each provider a
 * single time.
 */
export function mergeRegistryEntries(...entries: readonly PluginRegistryEntry[]): PluginRegistryEntry {
  return mergeRecords(entries) as PluginRegistryEntry;
}

function mergeRecords(records: readonly any[]): any {
  const keys = new Set<string>();
  for (const r of records) for (const k of Object.keys(r)) keys.add(k);

  const out: any = {};
  for (const key of keys) {
    const values = records.map((r) => r[key]).filter((v) => v !== undefined);
    out[key] = values.every(Array.isArray) ? mergeLists(values) : mergeRecords(values);
  }
  return out;
}

function mergeLists(lists: readonly any[][]): any[] {
  const out: any[] = [];
  const seen = new Set<any>();
  for (const list of lists) {
    for (const item of list) {
      if (seen.has(item)) continue;
      seen.add(item);
      out.push(item);
    }
  }
  return out;
}
