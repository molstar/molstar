/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `fs` module a Node.js host gives (`setFSModule` of `util/data-source`), read by the downloads of `file://` URLs
// and the range reads of them. Not published (`./util/internal/*` is not exported): `setFSModule` is the public setter.

let _fs: typeof import('fs') | undefined;

export function getFS() {
  if (!_fs) {
    throw new Error(
      "When running in Node.js and reading from files, call mol-util/data-source's setFSModule function first.",
    );
  }
  return _fs;
}

export function setFS(fs: typeof import('fs')) {
  _fs = fs;
}
