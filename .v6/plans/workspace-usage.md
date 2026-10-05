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

`pnpm build:lib` compiles the TypeScript solution with the native TypeScript 7.0.2
`tsc -b` command and copies package assets.
`pnpm build:apps` builds apps and browser examples directly from source using the
`molstar-src` export condition; it does not need a preceding library build.
`pnpm build:distribution` assembles the Viewer/MVS Stories classic bundles and
browser ESM modules in `distributions/molstar/build`. `pnpm pack:workspace` creates
local public package tarballs in `build/packages` and checks their version ranges and published entry points/assets.

Use `pnpm dev` to watch all browser apps and examples on one server, or
`pnpm dev:apps -- viewer` for a source-based Viewer dev server. Each browser app/example
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

Colocated unit tests use `_test/<name>.test.ts`. Jest discovers `.test.ts` files;
package builds and tarballs exclude `_test/` directories, including their helpers.

All public packages release together. Change `version.json`, run
`pnpm version:sync`, and update the lockfile. Workspace ranges are `workspace:*`;
`pnpm pack` replaces them with the exact release version. `pnpm version:check`
detects drift. No package is published by these commands.

Normal builds use `skipLibCheck: true`: each project checks its own source and its
usage of dependency types, without repeatedly checking imported declarations.
Run `pnpm check:types:full` to force a rebuild with declaration checking enabled
in every project using the native TypeScript 7 compiler. The check generates
short-lived configs beside the originals, preserving source paths, package/type
resolution, and project references while overriding `skipLibCheck`. These configs
are removed on success or failure. Separate `*.full.tsbuildinfo` caches keep the
normal incremental caches intact and are not published. CI runs this check after
the library build, and the publish gate also requires it.

Before publishing, run `pnpm check:publish`. This required release gate lints, tests, builds
libraries/apps/distributions, runs the full declaration check, verifies workspace
and version consistency, packs all public packages, and tests packed Node, types,
source-condition, and CLI consumers. It does not publish or install native GL.
Browser smokes remain local checks as described above. Publishing automation is
still deferred; wire this gate into it when that automation is added.

CLI packages expose checked-in `bin/*.mjs` launchers, so a clean install can link
commands before compiling. Run `pnpm build:lib` before using those commands in a
source checkout; installed tarballs already contain their compiled implementations.

Native `gl` and `canvas` are optional peers. Workspace peer auto-installation is
disabled; all ordinary dependencies remain declared explicitly. Normal installation and
`pnpm test` do not install them. Opt in explicitly:

```sh
pnpm native:install                    # gl only, for native tests/capture
pnpm test:native                       # requires the installed backend
pnpm native:install -- --canvas        # gl + canvas, for the rendering CLI
pnpm native:run -- node cli/mvs-render/lib/mvs-render.js --help
pnpm native:run -- pnpm --dir smoke smoke:headless
```

Native modules and their separate npm lockfile live in ignored `.cache/native/`.
The setup command leaves workspace manifests and `pnpm-lock.yaml` unchanged.
`native:run` exposes them to Node commands through `NODE_PATH`. Applications using
installed packages can instead install the optional peers themselves and continue
to inject native modules into the headless context. CI explicitly installs `gl`
and runs native unit tests under Xvfb. Consumer smoke checks run locally.

`smoke/` tests isolated tarball consumers, Node exports, emitted types, source
build resolution, MVS CLI execution, and native browser ESM/classic rendering.
Browser checks serve files extracted from the distribution tarball. Use the
focused `pnpm --dir smoke smoke:*` commands after building, or `pnpm --dir smoke smoke` to build and pack
first. `smoke:headless` is a separate native capture check. Set
`MOLSTAR_SMOKE_BROWSER` to a Chromium executable when needed.

The full plugin composition proposal, slim bundles, renderer redesign, broad
barrel removal and publishing automation remain outside this structural prototype;
see [`workspace-prototype.md`](workspace-prototype.md) for the scope and follow-up.

The workspace uses TypeScript 7.0.2 for compilation and Biome 2.5.15 for linting.
The separate `@typescript/typescript6` 6.0.2 development dependency provides only
AST/config parsing for the workspace inventory, dependency checker, and full-check
config generator. It does not compile or check source. TypeScript 7.0 has no stable
replacement for that programmatic API; migrate these helpers when one is available.
See [Microsoft's compatibility guidance](https://devblogs.microsoft.com/typescript/announcing-typescript-7-0/#running-side-by-side-with-typescript-60).

VS Code recommends the TypeScript 7 extension (`TypeScriptTeam.native-preview`)
and selects the local `node_modules/typescript` package through `js/ts.tsdk.path`.
Other editors should select their native TypeScript language server. Normal build,
project-reference, and packed type-consumer commands still use `tsc`, which now
launches the platform-specific native executable. Retain optional dependencies
when installing: TypeScript distributes that executable through platform packages.
Migration evidence is recorded in [the implementation plan](workspace-prototype.md#10-next-step-typescript-7-and-biome).

### Linting and formatting

```sh
pnpm lint          # Check lint rules only
pnpm lint:fix      # Apply safe lint fixes only
pnpm format:check  # Check Biome and Prettier formatting; required by CI/publish
pnpm format        # Reformat source, configuration, documentation, and styles
```

Biome owns JavaScript/TypeScript, JSON/JSONC, HTML, CSS, and GraphQL. Prettier 3.9.9
owns Markdown/MDX, YAML, and SCSS; its ignore file restricts it to those formats.
For a focused format, use `pnpm exec biome format --write path/to/file.ts` or
`pnpm exec prettier --write path/to/file.scss`. VS Code recommends both extensions
and enables format on save with the matching formatter for each language.

All formats use two spaces, a 120-column line width, and LF endings with a final
newline. MkDocs pages under `docs/docs/` retain four-space nested Markdown blocks
because Python-Markdown requires that indentation to preserve their rendering. JavaScript uses single quotes, JSX/HTML use double quotes, semicolons are
always present, multiline structures get trailing commas everywhere allowed,
and arrow parameters always have parentheses. Object braces have spaces; property
names are quoted only when required. Objects/arrays wrap automatically and JSX
closing brackets occupy their own line in multiline elements. Markdown prose wraps
at 120 columns. Prettier leaves embedded code unchanged to avoid a competing style
for fenced JavaScript/TypeScript examples. `.editorconfig` supplies the shared
indentation/newline defaults to other editors and retains Markdown hard line breaks.

Assist actions, including import organization, are disabled. Generated `lib/` and
`build/` trees, dependencies, `deploy/`, `docs/site/`, `build.mjs`, Git-ignored files,
and lockfiles are excluded. Linting covers JavaScript/TypeScript files (including
tests, scripts, and smoke fixtures), matching the previous ESLint language scope.
The 5 MiB Biome file limit includes the large alpha-orbitals example data.

The settings are committed before the repository-wide formatting pass so changes
to tooling can be reviewed separately. Formatting checks become green after the
bulk pass, which is recorded separately in `.git-blame-ignore-revs`.

Only explicitly selected lint rules are enabled; Biome's recommended preset is
not enabled. The mapping from the previous ESLint rules is:

| Previous rule | Biome rule / treatment |
| --- | --- |
| `eqeqeq: smart` | `suspicious/noDoubleEquals`, allowing null comparisons; other `==` comparisons are rejected, including same-type comparisons ESLint allowed |
| `no-eval` | `security/noGlobalEval` |
| `no-new-func` | `nursery/noImpliedEval`; also rejects string arguments to timers |
| `no-extend-native` | `nursery/noExtendNative` (warning); polyfills have an explicit suppression |
| `no-unsafe-finally` | `correctness/noUnsafeFinally` (warning) |
| `no-self-compare` | `suspicious/noSelfCompare` (warning); the NaN polyfill keeps its explicit suppression |
| `no-var` | `suspicious/noVar` |
| No default export declarations | `style/noDefaultExport`; also covers default re-exports, with a suppression for the image-loader declaration |
| `no-throw-literal` | `style/useThrowOnlyError` |
| `prefer-const` | `style/useConst`; Biome has no matching destructuring/read-before-assignment options |
| `no-constant-binary-expression` | `suspicious/noConstantBinaryExpressions` |
| `@typescript-eslint/prefer-namespace-keyword` | `suspicious/useNamespaceKeyword` (warning) |
| Quotes, semicolons, braces, whitespace, and spacing rules | Formatter; enforced by the separate `pnpm format:check` gate |
| `spaced-comment`, `no-new-wrappers` | No equivalent enabled; comment spacing and wrapper construction are not enforced |

The two nursery rules are explicitly enabled and the Biome version is pinned.
These rules may differ from ESLint and may change when Biome is upgraded.

Package exports use wildcard mappings for regular modules, with explicit root and
directory aliases, TS/TSX exceptions, Sass/CSS patterns, and blocked test paths.
The workspace refresh script preserves this compact form. Workspace and tarball
checks expand the patterns and validate every concrete conditional target.

CIF field-name filters have one source of truth in root `data/cif-field-names/`.
Workspace CLI runs read those files directly. Library asset staging copies them
into ignored `cli/cifschema/lib/data/` for installed consumers, and tarball checks
verify each staged filter against its root source.
