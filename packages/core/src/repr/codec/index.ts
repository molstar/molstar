/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/repr/codec: the codec interface and its conformance suite; no codec (_spec/README.md "Codecs").

export { type Codec, CodecIssues, type DecodeOptions, type Decoded } from './codec.js';
export { ConformanceLanguage, checkCodec, conformanceTrees } from './conformance.js';
