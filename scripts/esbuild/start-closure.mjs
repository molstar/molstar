/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';

/**
 * Bytes an ESM entry loads at start: the entry and the chunks it imports statically, transitively (esbuild's metafile;
 * `entry` is the output's path relative to the build's working directory). Dynamic imports are lazy chunks: not counted.
 */
export function startBytes(meta, entry) {
  const stack = [entry.split(path.sep).join('/')];
  const seen = new Set();
  let bytes = 0;
  while (stack.length) {
    const file = stack.pop();
    const out = meta.outputs[file];
    if (seen.has(file) || !out) continue;
    seen.add(file);
    bytes += out.bytes;
    for (const i of out.imports) if (i.kind === 'import-statement' && !i.external) stack.push(i.path);
  }
  return bytes;
}
