## Core: assets

- `@molstar/core`: `AssetManager.ResolveOptions` and `fetchData` take `credentials`, `redirect` and `maxBytes` (enforced
  while the body streams, failing with the new `SizeLimitExceeded`, also exported from `@molstar/core/util/assets`), and
  `AssetManager.resolve` takes a `signal`. `DefaultAssetManager` honors them; with none given it behaves as before. It
  refuses a held asset to a read with other credentials or a smaller `maxBytes`; `fetchData` honors `maxBytes` and the
  signal for `file://` in Node.js and refuses a `maxBytes` that is not a number of bytes. New
  `@molstar/core/util/assets-conformance` exports `checkAssetManager`, a conformance suite any `AssetManager` passes.
