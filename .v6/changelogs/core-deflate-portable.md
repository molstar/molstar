## Core: a deflate with the same bytes on every runtime

- `deflate` (`@molstar/core/util/zip/zip`) takes `native: false`: core's own deflate instead of the platform's
  `CompressionStream`, whose bytes are the same on every runtime (the platform's change with its zlib, e.g. between
  Node versions); the default is unchanged (#2080)
