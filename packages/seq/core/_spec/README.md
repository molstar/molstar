# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.0.0` until the first slice runs). Why: one language, one
  description, one number for readers to check.
