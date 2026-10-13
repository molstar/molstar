---
title: Mol* v6
description: An introduction to Mol* v6 and what to expect from its documentation.
---

Mol\* is a library and a set of apps for viewing and analyzing molecular structures and volumes. Version 6 is a rewrite
of how the library is packaged and composed, with a new user interface and an optional AI assistant.

:::caution[Work in progress] These are the v6 pages written so far. For everything else, see the
[v5 documentation](/docs/v5/), which describes the current release. :::

## What changes in v6

- **A workspace of ESM packages.** The library is split into layers (`core`, `io`, `model`, `graphics`, `plugin`), each
  with its own dependencies and exports, instead of one `molstar` package with a `lib/` namespace.
- **Explicit plugin composition.** A plugin starts with empty registries, and you choose the representations, themes,
  file formats and behaviors it includes. This keeps what you do not use out of your bundle.
- **Per-app builds.** The viewer, the docking viewer, the mesoscale explorer and the other apps each own their
  dependencies and their build.
- **Prism, a new UI,** and an optional **AI assistant** built on an agent harness.
- **MolViewSpec** has a standalone builder package, separate from the Mol\* runtime.

## Where to find more

- The [Mol\* repository](https://github.com/molstar/molstar) holds the v6 sources, with design notes and plans in its
  `.v6` folder.
- The [v5 documentation](/docs/v5/) covers the plugin, the data access tools and the extensions as of the current
  release.
- [All documentation versions](/docs/versions/).
