/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export { mp3Duration, pcmDuration, pcmFromBytes, pcmToWav, wavSamples } from './audio.js';
export {
  extensionOf,
  type NarrationAudio,
  type NarrationEncoder,
  passThrough,
  wavEncoder,
} from './encoder.js';
export { geminiNarration } from './gemini.js';
export { narrate } from './narrate.js';
export { openaiNarration } from './openai.js';
export {
  type Fetch,
  NarrationError,
  type NarrationModel,
  type NarrationProvider,
  type NarrationRequest,
  type NarrationVoice,
  type ProviderAudio,
  refusal,
} from './provider.js';
export { narrationProviders } from './providers.js';
