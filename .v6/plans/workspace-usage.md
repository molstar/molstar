# Working with the v6 workspace prototype

The repository is a pnpm workspace. Use Node 22 or newer and the pnpm version
recorded in the root `packageManager` field. Native headless backends may require
an older supported Node ABI and platform build tools.

```sh
pnpm install
pnpm build
pnpm check:workspace
pnpm test
pnpm --dir smoke smoke
```

`pnpm build:lib` compiles the TypeScript solution and copies package assets.
`pnpm build:apps` builds apps and browser examples directly from source using the
`molstar-src` export condition; it does not need a preceding library build.
`pnpm build:distribution` assembles the Viewer/MVS Stories classic bundles and
browser ESM modules in `distributions/molstar/build`. `pnpm pack:workspace` creates
local public package tarballs in `build/packages` and checks their version ranges and published entry points/assets.

Use `pnpm dev:apps -- viewer` for a source-based Viewer dev server. Each browser app/example
also has its own `build` and `dev` scripts. Select examples from the root with
`pnpm dev:examples -- basic-wrapper`. Both selectors accept `--port 1340` and
`--help` to list available browser targets. Omit names to watch every browser
app (`pnpm dev:apps`) or example (`pnpm dev:examples`). Supply multiple names to
watch a subset, such as `pnpm dev:apps -- viewer mvs-stories`. All selected targets
share one server, which prints each page URL. Node examples compile with TypeScript.
The root `scripts/workspace/inventory.json` records package ownership; explicit
exports and dependencies live in each package's manifest. When adding/moving a
source module, update those declarations and TypeScript references. The
`scripts/refresh-workspace.mjs` helper regenerates the explicit inventories from
production imports; review its result before committing. It preserves external
and optional dependency choices and is not a runtime resolver.

## Package families

| Location | Public package | Ownership |
| --- | --- | --- |
| `packages/core` | `@molstar/core` | Utility, data, math, tasks, state |
| `packages/io` | `@molstar/io` | Parsers and low-level writers |
| `packages/model` | `@molstar/model` | Models, format adapters, properties, queries |
| `packages/graphics` | `@molstar/graphics` | Existing GL, geometries, themes, representations, Canvas3D |
| `packages/plugin/core` | `@molstar/plugin` | Plugin context, behavior, state, default catalogs |
| `packages/plugin/ui` | `@molstar/plugin-ui` | React UI and styles |
| `packages/plugin/headless` | `@molstar/plugin-headless` | Node screenshots with injected native backends |
| `packages/mvs/builder` | `@molstar/mvs-builder` | Standalone data/schema/builder, MVSX, validation CLI |
| `packages/mvs/runtime` | `@molstar/mvs` | Plugin loading and semantic query validation |
| `distributions/molstar` | `molstar` | Browser distribution artifacts |

Extensions, apps, servers and CLI tools own their manifests in their corresponding
top-level directories. The filename migration inventory is in
[`migration-map.json`](migration-map.json). Cross-package imports use public
package subpaths rather than relative paths into another package's source.

Data-only shapes belong to model; rendering shape factories/helpers belong to
`@molstar/graphics/geo/shape/shape`. Model `Shape.create` takes an explicit group
count; the graphics factory derives it from geometry. Model-aware molecular/volume
writers live in model. GPU density helpers live in graphics. MP4 support is an
explicit extension: `@molstar/mp4-export-extension/headless` provides
`Mp4HeadlessPluginContext`, while the base headless package has no encoder import.

`CollapsableControls`, `CollapsableProps`, and `CollapsableState` now live at
`@molstar/plugin-ui/controls/collapsable`; base React components stay at
`@molstar/plugin-ui/base`. This split prevents initialization cycles between
base components and controls. MVS registration keeps its existing transformer ID
through a leaf `@molstar/mvs/behavior-id` module.

## ESM consumers

Node imports compiled installed JavaScript, with no source condition or loader:

```js
import { Task } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
import { MVSData } from '@molstar/mvs-builder';
```

The browser distribution contains `build/esm/viewer.js`, `viewer.css`, supported
library entry points, shared chunks and `import-map.json`. Serve that `build`
directory at `/build/`, then load the Viewer using native browser modules:

```html
<link rel="stylesheet" href="/build/esm/viewer.css">
<div id="viewer" style="width:800px;height:600px"></div>
<script type="module">
  import { Viewer } from '/build/esm/viewer.js';
  const viewer = await Viewer.create('viewer', { layoutIsExpanded: false });
  await viewer.loadStructureFromUrl('/data/example.pdb', 'pdb');
</script>
```

The import map has exact entries for the supported browser module inventory; it
is not an arbitrary deep-import resolver. Pass `--base-url /your/esm/path/` (or an
absolute HTTP URL) to `scripts/distribution.mjs` to generate map values for a
different deployment. Embed the generated map before any module scripts.
Classic assets remain under `molstar/build/viewer` and
`molstar/build/mvs-stories`. Consumers need no build step; Mol* maintainers build
the published artifacts once.

## Versions and checks

All public packages release together. Change `version.json`, run
`pnpm version:sync`, and update the lockfile. Workspace ranges are `workspace:*`;
`pnpm pack` replaces them with the exact release version. `pnpm version:check`
detects drift. No package is published by these commands.

`smoke/` tests isolated tarball consumers, Node exports, emitted types, source
build resolution, MVS CLI execution, and native browser ESM/classic rendering.
Browser checks serve files extracted from the distribution tarball. Use the
focused `pnpm --dir smoke smoke:*` commands after building, or `pnpm --dir smoke smoke` to build and pack
first. `smoke:headless` is a separate native capture check. Set
`MOLSTAR_SMOKE_BROWSER` to a Chromium executable when needed.

The full plugin composition proposal, slim bundles, renderer redesign, broad
barrel removal and publishing automation remain outside this structural prototype;
see [`workspace-prototype.md`](workspace-prototype.md) for the scope and follow-up.
