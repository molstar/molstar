/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const anchors: unknown[] = [];

/**
 * Keeps `value` alive for the life of the program and returns it. For V8: optimized code holds the hidden classes
 * (maps) it was specialised on only weakly, and is thrown away when a major GC finds no object of such a map. The
 * state objects of a builder (chunked arrays, the mesh builder's and marching cubes' states) live only while a build
 * runs, so after a GC between builds the next one runs unoptimized again (and boxes its doubles) until it re-tiers.
 * One small instance of each state kept here keeps its map (perf-audit-3 08/P12, 08/P15). Make it the way the real
 * ones are made, with the same kinds of values (a double where they hold doubles), or it keeps a map they do not use.
 */
export function anchorMap<T>(value: T): T {
  anchors.push(value);
  return value;
}
