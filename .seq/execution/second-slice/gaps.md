# Gaps (§gaps-queue)

## Built but not agreed

1. **Lane M**: `SizeLimitExceeded` exported from `@molstar/core/util/data-source` (beyond the four agreed options): the
   policy tells a size refusal from other failures by it. Recommendation: keep (a typed error, not a match on a name).
2. **Lane M**: Seq\*'s policy fixes `redirect: 'error'` (no option). Recommendation: keep until a composition needs
   `'follow'`; then a policy option.
3. **Lane L** (wire format): MessagePack extension codes beyond the agreed `0x11`–`0x18`: `0x10` special numbers (NaN,
   ±Inf, −0), `0x19` category (u32 dictionary size, each entry a u32 length and UTF-8, then int32 codes), `0x1a`
   bigint. Recommendation: confirm, pending the review.

## Agreed but not built

## Built but not adopted
