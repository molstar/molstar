/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Relative URIs inside an `.svsx` archive (seqspec:sources#archives): what names a file of it. Apart from svsx.ts so
// the validator can use it without loading the reader.

const SCHEME = /^[A-Za-z][A-Za-z0-9+.-]*:/;

/** A URI that names a file of the archive, not an absolute one (a scheme, a leading `/`, `//host`). */
export const isRelativeUri = (uri: string) => uri !== '' && !SCHEME.test(uri) && !/^[\\/]/.test(uri);

/**
 * The archive path a relative `uri` names from `base`, the directory of the main file (`""`: the root); `.` and empty
 * segments dropped. Null: not a relative URI, a backslash, `..` past the archive's root, or no path left.
 */
export function resolveArchivePath(uri: string, base: string): string | null {
  if (!isRelativeUri(uri) || uri.includes('\\')) return null;
  const path = base === '' ? [] : base.split('/');
  for (const seg of uri.split('/')) {
    if (seg === '' || seg === '.') continue;
    if (seg !== '..') path.push(seg);
    else if (path.pop() === undefined) return null;
  }
  return path.length === 0 ? null : path.join('/');
}
