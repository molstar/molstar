# Narration (`@molstar/ai-narration`)

Text to speech for the story mode: providers behind one port, audio encoding behind another. Independent of the rest of
`packages/ai` (no `@molstar/ai-*` imports). No package uses it yet; it is meant to be loaded lazily by the app.

## Ports (`narration/src/`)

- **`NarrationProvider`** (`provider.ts`): `{ id, name, models, voices, defaultVoice, synthesize(request) }`.
  `NarrationRequest { text, voice, model?, style?, signal? }` (`style` is a delivery note). It returns `ProviderAudio`:
  `{ kind: 'pcm', pcm: Int16Array, sampleRate }` (mono 16-bit) or `{ kind: 'encoded', bytes, mimeType, durationMs }`.
  Failures are `NarrationError`s worded for the person (`await refusal(provider, response)`).
- **`NarrationEncoder`** (`encoder.ts`):
  `encode(audio, signal?) → NarrationAudio { bytes, mimeType, extension, durationMs }`; encoded audio passes through.
  `wavEncoder()` (no dependency, Node too); `passThrough(encoded)` turns encoded provider audio into `NarrationAudio`.
- `narrate(provider, encoder, request)` = `encoder.encode(await provider.synthesize(request), request.signal)`; no
  retries, caching or chunking.

## Providers

| Factory                               | id       | Models                                                                          | Output                        |
| ------------------------------------- | -------- | ------------------------------------------------------------------------------- | ----------------------------- |
| `geminiNarration({ apiKey, fetch? })` | `google` | `gemini-3.8-flash-tts` (default), `…-flash-lite-tts`; 30 voices, `Kore` default | PCM (WAV or L16 from the API) |
| `openaiNarration({ apiKey, fetch? })` | `openai` | `gpt-4o-mini-tts`; 13 voices, `sage` default                                    | MP3, kept as it arrives       |

`narrationProviders({ google?, openai? }, fetch?)` gives those with a key, Gemini first.

## Encoding

`mediabunnyEncoder({ prefer = ['aac', 'opus'], bitrate = 64_000, canEncode? })` at `@molstar/ai-narration/mediabunny`
(not re-exported from the main entry, so Mediabunny stays out of a bundle that does not ask for it): AAC in MP4 (`.m4a`)
or Opus in WebM, by the browser's WebCodecs support, falling back to the next codec, then WAV.

**No MP3 encoding, ever**: MP3 encoders are LGPL (LAME or shine). A provider's MP3 is stored as it arrives; only its
duration is read (`mp3Duration`). Audio helpers in `audio.ts`: `wavSamples`, `pcmToWav`, `pcmDuration`, `mp3Duration`.

Tests (`_test/`) use injected `fetch`; nothing reaches a provider.
