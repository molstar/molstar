## Seq* (next): extensions

- Seq*: extensions. A descriptor (`defineExtension`, `defineViewerExtension`) whose id is the namespace of everything it
  ships; a layer type named `<id>:<name>` with a fallback, its own params and keys on built-in nodes; schemes; a reader
  without the extension draws the fallback and says which extension (and version) it needs; `viewer.extend` registers
  one at run time and undoes it.
- Language 0.5.0: `spec.extensions`, `layer.fallback`, `layer.type` a word, `<id>:<name>` for `layout`, `heights` and
  `content`.
- core's repr: `KindSpec.extensionParam` and `KindSpec.fallbackParam` keep an extension's params and a newer type's
  params given a fallback.
