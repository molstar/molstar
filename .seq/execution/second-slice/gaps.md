# Gaps (§gaps-queue)

## Built but not agreed

1. **Lane M**: `SizeLimitExceeded` exported from `@molstar/core/util/data-source` (beyond the four agreed options): the
   policy tells a size refusal from other failures by it. Recommendation: keep (a typed error, not a match on a name).
2. **Lane M**: Seq\*'s policy fixes `redirect: 'error'` (no option). Recommendation: keep until a composition needs
   `'follow'`; then a policy option.
3. **Lane L** (wire format): MessagePack extension codes beyond the agreed `0x11`–`0x18`: `0x10` special numbers (NaN,
   ±Inf, −0), `0x19` category (u32 dictionary size, each entry a u32 length and UTF-8, then int32 codes), `0x1a`
   bigint. The review found no collision (BinaryCIF uses no extensions). Recommendation: confirm.
4. **Lane L**: a one-line bug fix inside io's `common/msgpack` (`str()` reads with `utf8ReadValue`, so a leading U+FEFF
   survives), against Gate 0's "no change to io's MessagePack itself": a bug fix with no API change, which BinaryCIF
   shares. Allowed by the orchestrator, for the author to see.

5. **Lane K**: SeqLanguage has no visibility param, so "hide a layer" is `set :where custom :key hidden` in its cases.
   Recommendation: add `hidden: boolean` (default false) to `layer` (a language change: the author's call).
   **Decided** (2026-11-18): added; lane K builds it.
6. **Lane K**: `ChildIndexThreshold` 8 (unmeasured, from immutable.md's "about 8"); `AnyPatchLanguage` id `patch`,
   version `1.0.0`; a place laid out as `:after id`, `:first true` or `:last true` (patches.md to say so).

7. **Lane R**: `FlatThreshold` measured at 16 (the design said about 8); `Vector.indexOf` uses `Object.is` (finds NaN,
   tells −0 from 0, unlike `Array.indexOf`); collision nodes at the level two keys meet (documented).

## Agreed but not built

## Built but not adopted
