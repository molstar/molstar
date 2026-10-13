/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { NarrationAudio, NarrationEncoder } from './encoder.js';
import type { NarrationProvider, NarrationRequest } from './provider.js';

/** The whole flow: the provider voices the text, the encoder makes the file. */
export async function narrate(
  provider: NarrationProvider,
  encoder: NarrationEncoder,
  request: NarrationRequest,
): Promise<NarrationAudio> {
  const audio = await provider.synthesize(request);
  return encoder.encode(audio, request.signal);
}
