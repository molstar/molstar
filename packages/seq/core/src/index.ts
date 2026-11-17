/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The core of Seq* (next): the SeqViewSpec language, the coordinate model, the source protocol and the pure passes
// from a spec to draw lists. What it holds: _spec/README.md. Each folder becomes an entry of its
// own as it is built (`@molstar/seqx-core/<folder>`).

/** The version of the SeqViewSpec language this package reads and writes (_spec/README.md, Versions). */
export { SeqLanguageVersion as LanguageVersion } from './lang/language.js';
