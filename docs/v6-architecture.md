# Mol\* 6.0: packages, ESM, and plugin composition

Proposal against the `molstar@5.11.0` tree. The APIs and paths below describe the target, not features available in 5.x. See the [short summary](v6-summary.md) for the main decisions.

## 1. Scope

Mol* 6.0 moves to a pnpm workspace of ESM packages grouped by layer. Parsers, representations, and themes become explicit plugin features; `DefaultPluginSpec()` remains the full built-in composition. Apps keep esbuild, with dependencies and build configuration owned by each app.

Ship `@molstar/migrate-6` with the release to handle mechanical import changes and report manual work. The unscoped `molstar` package retains the CDN viewer bundle only: no legacy `lib/mol-*` exports, compatibility re-exports, or CommonJS build.

The release also includes a standalone MolViewSpec builder, dependency-cycle removal, maintainer skills, updated developer docs, and workspace CI.

Try JSR source publication with `--allow-slow-types`, beginning with the MVS builder, alongside native npm packages with compiled ESM from the same release commit. Defer fast-type migration and `isolatedDeclarations` to consideration for v7. The [fast-types and distribution analysis](v6-fasttypes.md) preserves the audit for that decision; its annotation work and agent estimates are outside the v6 scope.

Establish a rendering-backend boundary in `@molstar/graphics` for future WebGPU and other targets, retaining WebGL as the working implementation. The [rendering-backend design](v6-webgpu.md) covers the blast radius, minimal contracts, migration, and validation; a production WebGPU renderer and feature parity are later work.

The same geometry/readback boundary should support portable scene extraction for a future [Blender offline-rendering extension](v6-webgpu.md#71-offline-rendering-with-blender). Offline rendering uses scene snapshots and asynchronous jobs, separately from the interactive view contract. Keep Blender dependencies and integration in an optional extension; implementing it is outside the 6.0 scope.

Out of scope:

- Independent package versions or a package per parser/representation.
- Native Node execution of TypeScript source and a repository-wide erasable-syntax rewrite.
- Changes to transformer identifiers or snapshot JSON, except where registration must become explicit.
- Moving Python `molviewspec` into this repository.
- Publishing a stories library or merging the MolViewStories webapp.

## 2. Starting point

Today one package owns library code, apps, servers, and their npm dependencies. Consumers use deep imports such as `molstar/lib/mol-plugin-ui`; the package has no `exports` map.

| Output | Current build |
| --- | --- |
| `lib/` | `tsc`, then `tsc-alias` to complete relative JS paths; includes declarations |
| `lib/commonjs/` | Separate `tsc` CommonJS emit; used by CLI bins and servers |
| `build/<app>/` | esbuild from `src/`, with SCSS, copied assets, version injection, watch, and serve |

Tests use Jest and `esbuild-jest-transform`, with colocated `_spec/*.spec.ts` files. Asset copying and version patching run from root scripts.

`PluginSpec` already selects actions, behaviors, animations, and additional `customFormats`. Registry constructors still preload all built-in formats, representations, and themes. Large transform modules and default-spec imports connect even a small plugin to the full catalog.

The existing source graph is **not a layer DAG**. Besides plugin/state initialization cycles, lower folders contain runtime dependencies on higher layers. Resolving those dependencies is part of the architecture work, not a consequence of moving files.

## 3. Packages and dependency boundaries

### 3.1 Package ownership

Keep recognizable subpaths, removing the `mol-` prefix. The table is the target ownership after the relocations in §3.2.

| Package | Owns |
| --- | --- |
| `@molstar/core` | `util/`, `task/`, `data/`, `math/`, `state/`; no molecular or rendering dependencies |
| `@molstar/io` | Readers and writers; depends on core |
| `@molstar/model` | Structures, volumes, particles, domain formats, properties, and script/query code |
| `@molstar/graphics` | `gl/`, `geo/`, `theme/`, `repr/`, `canvas3d/`, plus rendering code relocated from model/math |
| `@molstar/plugin` | Plugin and plugin-state together: runtime, spec helpers, and explicit default-spec/catalog entry points |
| `@molstar/plugin-ui` | React UI, UI-spec types, and an explicit default UI-spec entry point |
| `@molstar/mvs-builder` | Standalone MolViewSpec schema, builder, serialization, and validation |
| `@molstar/mvs` | MolViewSpec runtime, feature, loader, and rendering CLI |
| `@molstar/<name>-extension` | Individual extensions and their dependencies |
| `@molstar/viewer` | Published Viewer API and app |
| `@molstar/<name>-server` | Model, volume, and plugin-state servers |
| `@molstar/migrate-6` | Migration CLI |
| `molstar` | CDN viewer assets at `build/viewer/` |

The library dependency direction is:

```text
plugin → graphics → model → io → core
```

An arrow means “depends on.” Packages may also depend directly on any lower layer they import. UI, extensions, apps, and servers sit above the layers they use. `DefaultPluginUISpec` in plugin-ui extends `DefaultPluginSpec` in plugin; plugin never depends on UI. Default compositions live behind explicit subpath entry points in those packages, keeping React out of the non-UI default.

`mvs-builder` has no `@molstar/*` dependency. `mvs` depends on the builder and the runtime layers it imports. Servers that only need IO and model should not acquire plugin or graphics through those packages.

### 3.2 Required relocations before packaging

Audit the graph by proposed package ownership before assigning project references. These are known runtime edges that cannot be removed with `import type`:

| Current source | Conflicting edge | Required work |
| --- | --- | --- |
| `mol-data/db/column.ts` | core → IO number parser | Move shared numeric parsing primitives into core and make IO consume them |
| `mol-util/data-source.ts` | core → IO string/UTF-8 helpers | Move general byte/string helpers into core; keep format-specific parsing in IO |
| `mol-math/geometry/gaussian-density/gpu.ts` | core → GL | Keep CPU/domain-independent math in core; move GPU computation and its GL-facing API to graphics |
| `mol-model/shape/shape.ts` | model → geometry, themes, GL | Separate shape data from rendering helpers; place graphics-dependent APIs in graphics |
| `mol-model-formats/shape/*` | model → mesh builders | Move mesh conversion into graphics, retaining raw file readers in IO |
| `mol-model-props/**/representations`, `**/themes`, and label helpers | model → representations/themes | Keep property computation in model; move visuals/themes to graphics and remove rendering dependencies from computation |
| Model particle/unit transform utilities | model → geometry transform helpers | Extract shared matrix/transform primitives into core or move the graphics-specific callers |

This is a starting list, not the full audit. Resolve declaration dependencies too: a type-only edge avoids runtime evaluation, but still requires a package dependency and can create a cycle in `tsc -b` references. Assign exact replacement subpaths during this work and record exceptions in the migration map.

Do not preserve an old path through a lower-layer re-export that recreates the dependency. General leaf paths remain recognizable; APIs that cross these boundaries need explicit migration entries.

### 3.3 Plugin initialization cycles

Within the plugin layer:

- Use `import type` for parser result types in `PluginStateObject`, and for `PluginContext` wherever only its type is needed. Extract a context interface where it helps isolate construction; its dependencies must remain type-only.
- Separate the `PluginBehavior` contract from the classes extending `PluginStateObject`. Behavior detection uses `typeClass === 'Behavior'` without importing those classes.
- Move default specs out of base spec/context modules into `@molstar/plugin/default-spec` and `@molstar/plugin-ui/default-spec`. Base construction takes an explicit spec.
- Split transform modules into leaves. Leaves, builders, and managers import individual transformers, never the `StateTransforms` facade.
- Rebuild that facade as a plain object in `@molstar/plugin/state/transforms` after removing cycles. Delete the lazy getters.

Enforce an acyclic package graph and no value-import cycles within packages. Type-only cycles inside a package are allowed; they must not conceal a package cycle.

## 4. Workspace, imports, and dependencies

### 4.1 Layout

```text
packages/
  core/                       # @molstar/core
  io/
  model/
  graphics/
  plugin/
  plugin-ui/
  mvs/
    builder/                  # @molstar/mvs-builder
    runtime/                  # @molstar/mvs
extensions/<name>/            # @molstar/<name>-extension
apps/
  viewer/                     # published @molstar/viewer
  docking-viewer/             # private workspace apps
  mesoscale-explorer/
  mvs-stories/
examples/<name>/              # private workspace packages
servers/<name>/               # published server packages
cli/<name>/                   # general CLI packages, including migrate-6
molstar/                      # CDN-only package
scripts/                      # shared build/release tooling
.agents/                      # maintainer skills
```

Each library package has `src/`, `lib/`, `package.json`, and a composite `tsconfig.json`. CLI tools tied to a product can live in that product package: MVS validation belongs to the builder, MVS rendering to the runtime. Give other existing bins, including `cif2bcif` and `cifschema`, explicit CLI package owners.

```yaml
# pnpm-workspace.yaml
packages:
  - 'packages/*'
  - 'packages/mvs/*'
  - 'extensions/*'
  - 'apps/*'
  - 'examples/*'
  - 'servers/*'
  - 'cli/*'
  - 'molstar'
```

Use physical moves into the workspace. Keep relocation commits separate from import rewrites where practical; keep each completed PR buildable.

### 4.2 One import convention

Use the same public package specifiers inside and outside the repository:

```ts
// packages/graphics/src/canvas3d/passes/illumination.ts
import { ValueCell } from '@molstar/core/util';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
```

Rules:

1. Crossing a layer-root folder (`util`, `task`, `gl`, `canvas3d`, plugin `state`, etc.) uses an extensionless package subpath, even within the same package.
2. Relative source imports stay within that folder and use the emitted `.js` path: `./names.js`, `../hex.js`, or `./controls.js` for a TSX module compiled with `react-jsx`. TypeScript and esbuild resolve these to source files while building; emitted JS retains the specifiers.
3. Package self-imports resolve through `exports`. Any development `paths` mappings must describe the same public subpaths and source files.
4. Apps, extensions, examples, and servers never reach into another package through relative filesystem paths.

Lint these boundaries. Bare specifiers such as `@molstar/core/util/color` and relative `.js` specifiers remain unchanged in emitted JS. Resolve directory imports explicitly to an index file or an exported subpath.

### 4.3 Dependency declarations

Publish all public Mol* workspace packages at one version. Declare internal dependencies with `workspace:*`, which becomes an exact version on publication. A consumer should not have to install each lower layer manually. Mixing separate Mol* release versions in one application can still duplicate state/transformer identities; lockstep publishing does not prevent consumers from requesting conflicting versions.

Each package declares every external package it directly imports. “Owned by core” does not make a transitive dependency available to another package. Shared versions live in the pnpm catalog; root `package.json` is private tooling.

| Dependency | Ownership |
| --- | --- |
| `rxjs`, `immutable`, `mutative` | Every direct importer, starting with core |
| `tslib` | Each package whose emit uses imported helpers |
| `argparse` | CLI/server packages that parse arguments |
| `express`, `compression`, `cors`, `swagger-ui-dist` | Server packages |
| `h264-mp4-encoder` | MP4 export extension |
| `io-ts` | MVS builder, plus any remaining direct importer until migrated |
| `react-markdown`, `remark-gfm` | Plugin UI |

Use peers for React/React DOM on UI packages, and optional peers for headless dependencies (`gl`, `canvas`, `pngjs`, `jpeg-js`) and optional cloud storage where used. Do not add React to plugin or core.

Types needed only to build a package belong in its dev dependencies. Types referenced by published declarations must be available to consumers through dependencies or declared peers. In particular, React does not install `@types/react`: UI packages need an explicit consumer-facing type dependency/peer policy, verified with a clean TypeScript consumer.

## 5. ESM and TypeScript source

### 5.1 Supported execution modes

Publish ESM JavaScript and declarations in `lib/`, plus source in `src/` for bundlers and debugging. Set `"type": "module"` on published packages and retain the existing Node **22 or later** baseline; raise it only if runtime or tooling requirements demand it.

| Consumer | Resolution | Executes |
| --- | --- | --- |
| Installed library, CLI, or server | Default `import` condition | `lib/*.js` |
| TypeScript consumer | `types` condition | Checks `lib/*.d.ts` |
| In-repo esbuild | `molstar-src` condition | Compiles source, including TSX |

Node executes compiled JavaScript. Keep all published bins on `lib/*.js`; workspace tooling runs as JavaScript or is compiled before execution. `molstar-src` selects source for bundlers and does not promise native Node execution. Publishing source does not require erasable syntax.

Validated JSR packages expose TypeScript source for Deno or compatible tooling. Deno consumers of native npm packages use the compiled exports above. Publishing to either registry does not make browser, Node, or optional native APIs available in every runtime; see the [distribution matrix](v6-fasttypes.md#6-typescript-distribution-through-npm).

### 5.2 Compiler settings

Merge these into the existing strictness settings:

```json
{
  "compilerOptions": {
    "target": "ES2022",
    "module": "NodeNext",
    "moduleResolution": "NodeNext",
    "verbatimModuleSyntax": true,
    "isolatedModules": true,
    "jsx": "react-jsx",
    "declaration": true,
    "declarationMap": true
  }
}
```

Use `.js` relative specifiers in TypeScript source as described in §4.2. No extension-rewriting compiler option is needed. Keep `verbatimModuleSyntax` and explicit `import type` for predictable module dependencies.

During dependency and module refactoring, keep the existing ESM build settings. If `verbatimModuleSyntax` is enabled in the shared config while CommonJS still exists, explicitly set it to `false` in `tsconfig.commonjs.json`. Otherwise ESM syntax fails in that build with TS1287/TS1295. Switch to `NodeNext` and remove the override with the CommonJS build in the ESM phase. See [TypeScript’s module-syntax behavior](https://www.typescriptlang.org/tsconfig/verbatimModuleSyntax.html).

### 5.3 TypeScript syntax and performance

Do not enable `erasableSyntaxOnly`. Both library and app builds compile TypeScript, so existing namespaces, enums, parameter properties, and other compiler-supported syntax can remain. Refactor individual constructs only where module boundaries, composition, or ESM compatibility require it. There is no blanket namespace/enum conversion or syntax-driven CIF schema regeneration phase.

Keep hot `const enum`s under the existing `isolatedModules` constraints. If a necessary refactor changes their use or emit, inspect the generated code and benchmark the affected parse/render paths. Do not assume that replacing enum uses with object properties or module constants preserves performance.

Do not require `isolatedDeclarations` or a broad annotation migration in v6. Preserve existing inferred API precision, including parameter/schema keys, literal unions, overloads, and factory constructor types. Consider fast types for v7 using the [audit](v6-fasttypes.md#3-measured-blast-radius); do not redesign public contracts solely to satisfy JSR fast types in this release.

### 5.4 Convert runtime CommonJS assumptions

Changing compiler flags does not fix `require`, `__dirname`, or `__filename`. Audit source and root scripts before setting `"type": "module"`:

- Convert module loads to ESM imports; use `createRequire(import.meta.url)` only where CommonJS loading is needed, such as an optional native dependency.
- Replace path globals with `import.meta.url`-based paths and account for moved source, emitted files, and packaged assets.
- Update conditional loads such as `servers/model/preprocess.ts`, native loading in `cli/mvs/mvs-render.ts`, cloud storage loading, and `cifschema` data paths.
- Convert `scripts/clean.js` and `scripts/deploy.js` to ESM, or explicitly retain tooling as `.cjs` and update callers. Keep the published library ESM-only.

Smoke-test emitted CLI/server entry points and verify packaged data assets. Typechecking alone will not catch these runtime failures.

### 5.5 Exports and package contents

A representative core export map:

```json
{
  "name": "@molstar/core",
  "version": "6.0.0",
  "type": "module",
  "engines": { "node": ">=22.0.0" },
  "exports": {
    "./task": {
      "types": "./lib/task/index.d.ts",
      "molstar-src": "./src/task/index.ts",
      "import": "./lib/task/index.js"
    },
    "./*": {
      "types": "./lib/*.d.ts",
      "molstar-src": "./src/*.ts",
      "import": "./lib/*.js"
    }
  },
  "files": ["lib", "src"]
}
```

Add conditional entries for other directory indexes, root entry points, and `.tsx` sources; the generic `*.ts` pattern does not cover TSX. Public code uses `@molstar/core/util/color`, without `.js`. Appending `.js` would make this wildcard target `color.js.js`.

Publish only intended source/assets and generated output. Exclude `_test/` and fixtures with a pack step or appropriate nested ignore files; inspect the actual tarball. Export UI skins and built CSS explicitly. Mark modules `sideEffects: false` only after auditing initialization behavior, and retain CSS/asset side effects where needed.

Use project references for package builds. Source conditions drive esbuild; normal consumer types resolve to generated declarations. Prove both from a clean checkout and from packed packages, rather than depending on stale `lib/` output or root hoisting.

`lib/` is generated and ignored by Git, but its JavaScript, declarations, and required assets ship in the npm tarball and remain in the installed package. It is not merely a temporary input that packing removes. JSR source artifacts exclude it. See the [build-output distinction](v6-fasttypes.md#8-is-lib-only-temporary).

### 5.6 npm and JSR publication

Keep one source tree and derive npm/JSR manifests from one package/export inventory. npm retains compiled conditional exports, peers, bins, and source for opt-in bundlers. JSR uses explicit source exports and resolved dependency mappings. Its publisher supports `package.json` projects with `.js`-to-TypeScript import resolution; validate self-subpaths, TSX, assets, and exact dependencies before relying on that path. JSR's npm compatibility layer does not replace native publication to npmjs.com.

Try publication with `deno publish --dry-run --allow-slow-types`, then use the same allowance when publishing validated packages to JSR. Slow types can degrade JSR-generated documentation and npm-compatibility declarations and make consumer checking slower. Keep native npm declarations generated by `tsc`, and verify JSR source consumers without promising equivalent generated documentation/types. The allowance does not skip normal typechecking or other publication requirements.

Publish validated packages to both registries at the same version from the same commit, starting with the MVS builder. Extend JSR coverage in dependency order without implying universal runtime support. Keep per-registry completion records and retry partial releases from unchanged artifacts; the two registries cannot publish atomically. The [dual-publication design](v6-fasttypes.md#7-publishing-to-npm-and-jsr-together) covers normalization, validation, package coverage, and the `deno pack` alternative. Retain pnpm/`tsc -b` and its complete declarations for native npm packaging.

## 6. Plugin composition

### 6.1 Empty registries and explicit specs

`PluginContext` creates empty registries for formats, structure/volume/particle representations, and themes. `PluginSpec` supplies their contents along with existing actions, behaviors, animations, config, canvas, and layout settings.

Keep `customFormats` as a deprecated append-only alias for `formats` during 6.x. Initialize themes, representations, formats, actions, then behaviors/animations. Register themes in the relevant structure, volume, or particle scope.

`@molstar/plugin/spec` and `@molstar/plugin-ui/spec` contain types and composition helpers. They must not import default catalogs. This also applies to `createPluginUI`: its base entry point takes an explicit spec instead of importing `DefaultPluginUISpec` for an omitted argument. The lean package root entry points must not re-export default specs, full catalogs, or the `StateTransforms` facade.

Defaults and catalogs share packages with the runtime but remain separate modules. Slim consumers install those files without importing them; isolation comes from the module graph. Enforce these boundaries with lint/import-graph checks and bundle validation.

For the full built-in composition:

```ts
import { createPluginUI } from '@molstar/plugin-ui';
import { renderReact18 } from '@molstar/plugin-ui/react18';
import { DefaultPluginUISpec } from '@molstar/plugin-ui/default-spec';

const plugin = await createPluginUI({
    target: document.getElementById('app')!,
    spec: DefaultPluginUISpec(),
    render: renderReact18,
});
```

`DefaultPluginSpec` is exported from `@molstar/plugin/default-spec`; `DefaultPluginUISpec` extends it with UI defaults. Both preserve the existing built-in composition. Compositions that import extension packages remain in the Viewer or other apps, so plugin defaults do not acquire dependencies back on extensions.

### 6.2 Features

A `PluginFeature` groups declarative contributions: formats, transformers/actions, representations, themes, behaviors, and any presets they require. Importing a feature does not mutate registries. Format features live at plugin level; representation features may live in graphics and must not depend on plugin at runtime or through published types. Keep their contribution contract in graphics or adapt it structurally at the plugin boundary.

The proposed call shape remains:

```ts
PluginSpec.fromFeatures(Core, Sdf, BallAndStick, {
    canvas3d: { /* custom overrides */ },
});
```

The last argument is optional spec options. Give feature values a discriminator so the helper can distinguish a trailing feature from options; type the API with overloads or tuple rest parameters. Do not declare a parameter after a rest parameter.

Composition must define and test these rules:

- Process features in order; deduplicate identical contributions by registry key/transformer identity.
- Report conflicting providers for the same key instead of silently choosing one.
- Append extra actions/behaviors from options. Apply explicit settings such as canvas/layout overrides without losing feature contributions.
- Preserve existing behavior-based registration, such as interactions, with a defined initialization order.

SDF contributes `SdfProvider`, `TrajectoryFromSDF`, and its action. Ball-and-stick contributes its representation and default element-symbol/uniform color and physical/uniform size themes. `Core` contributes shared data/model/structure transforms, structure representation infrastructure, camera/highlight behaviors, uniform themes, and the minimal hierarchy/preset path used below.

### 6.3 Remove implicit catalog imports

Split `transforms/data.ts` and `transforms/model.ts` into per-operation/per-format modules. Apply the same rule to representation/theme registries, builders, UI controls, and preset helpers: a registry class or base entry point must not value-import all providers.

Use one path convention throughout:

- SDF reader: `@molstar/io/reader/sdf/parser`
- SDF model conversion: `@molstar/model/formats/structure/sdf`
- SDF feature/provider: `@molstar/plugin/state/formats/trajectory/sdf`
- Ball-and-stick feature/provider: `@molstar/graphics/repr/structure/representation/ball-and-stick`
- Full transform facade: `@molstar/plugin/state/transforms`

Keep full catalogs in explicit modules in the packages that own their providers. Default-spec modules assemble them; runtime leaves and base entry points never import the default specs or full catalogs. Internal code uses individual transformer/provider imports; the full facade is a consumer convenience entry point.

Transformer name strings remain unchanged. Importing a transform leaf must retain whatever registration is necessary for snapshots and actions; account for that behavior when auditing side-effect annotations.

### 6.4 Built-in names and preset fallback

Retain `BuiltInTrajectoryFormat` and sibling types for completion. Use `import type` at call sites. Keep names/type metadata in the owning package and verify the full catalogs against it. Deriving a type from a local catalog is also valid if the import is type-only and creates no upward package dependency.

A built-in name can typecheck while being absent from a particular plugin. Runtime lookup should report a clear missing-format/provider error. Preserve extension-provider and string overloads where supported today.

Presets must consult available representations and themes. When a preferred representation is missing, select an applicable registered provider and its registered defaults. Today's `registry.default` is simply the first entry; it does not check applicability. Define the no-applicable-provider case explicitly, and avoid building duplicate fallback visuals for the same component.

### 6.5 Slim-plugin acceptance example

```ts
import { PluginSpec } from '@molstar/plugin/spec';
import { Core } from '@molstar/plugin/features/core';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { BallAndStick } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { createPluginUI } from '@molstar/plugin-ui';
import { renderReact18 } from '@molstar/plugin-ui/react18';
import type { BuiltInTrajectoryFormat } from '@molstar/plugin/state/formats/trajectory/types';

const spec = PluginSpec.fromFeatures(Core, Sdf, BallAndStick, {
    canvas3d: { /* custom overrides */ },
});
const plugin = await createPluginUI({
    target: document.getElementById('app')!,
    spec,
    render: renderReact18,
});
const format: BuiltInTrajectoryFormat = 'sdf';
const data = await plugin.builders.data.download(
    { url: '/ligand.sdf' }, { state: { isGhost: true } }
);
const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
```

This is a release acceptance target: render an SDF ligand while excluding unrelated mmCIF/CCP4 parsers, cartoon/volume representations, and MP4 export from the runtime import graph and app bundle. Verify with an import-graph check, esbuild metadata, and a rendering smoke test. Empty registries alone do not prove the result.

## 7. MolViewSpec, extensions, and apps

Split `extensions/mvs` into:

| Package | API and responsibilities |
| --- | --- |
| `@molstar/mvs-builder` | Schema, `createMVSBuilder`, `MVSData`, MVSJ/MVSX serialization/validation, `mvs-validate`, schema-printing CLI |
| `@molstar/mvs` | `loadMVS`/`loadMVSData`, plugin feature, annotations, cameras, runtime representations, `mvs-render` |

The builder replaces [molviewspec-ts](https://github.com/molstar/mol-view-spec/tree/master/molviewspec-ts) and the JSR `@molstar/molviewspec` distribution after API/parity checks. Publish the replacement to npm and JSR. Own `io-ts` and any archive dependencies in the builder; inline or extract its small general helpers without introducing a Mol* runtime dependency. An optional core peer would not make unconditional core imports optional.

The runtime imports the builder as a dependency; do not bundle a second copy into it. Python stays in the [mol-view-spec repository](https://github.com/molstar/mol-view-spec).

Other extensions become `@molstar/<name>-extension`. Each owns its direct dependencies, exports features/behaviors, and imports UI only when needed. Viewer dependencies and imports define its extension set; a smaller app declares and imports only its selected extensions.

`@molstar/viewer` is published. Docking viewer, mesoscale explorer, MVS stories, and examples remain private workspace apps. Moving parts of [MolViewStories](https://github.com/molstar/mol-view-stories) into a future `packages/mvs/stories` library needs a separate plan.

## 8. Builds and maintenance

### 8.1 Library and app builds

Use `tsc -b` with a root solution config and per-package references, `rootDir: src`, and `outDir: lib`. Cover all published TypeScript packages, including extensions, CLI, servers, and Viewer. Copy assets per package and generate plugin version source once for both library and app builds.

Apps use `scripts/esbuild/app.mjs` for shared plugins, source conditions, watch, and serve. Each app/example owns its entry, output paths, dependencies, themes, and scripts. Resolve imports from the importing workspace package; do not depend on a root `Apps` list or root dependency hoisting.

```text
pnpm --filter @molstar/viewer dev
pnpm -r --filter "./apps/**" --filter "./examples/**" build
```

Keep SCSS, copied HTML/images/icons, version injection, IIFE globals, and existing CDN output names. Shaders remain `.glsl.ts` strings. The deploy script consumes the app outputs after its ESM/tooling conversion. Stage Viewer JS/CSS and supporting assets into `molstar/build/viewer/` when packing the CDN package.

### 8.2 Tests, skills, and docs

Rename `_spec/*.spec.ts` to `_test/*.test.ts`; keep tests colocated. Reserve “spec” for future agent/behavior specifications. Choose a workspace runner with working ESM and export-condition resolution, remove Jest's `moduleDirectories: ["lib"]`, and retain browser/WebGL coverage. Install headless dependencies on the packages that run those tests.

Add `.agents/` skills for `add-extension`, `add-example-app`, `add-app`, `add-format`, `add-representation`, `add-server`, and `update-dependencies`. Root `AGENTS.md` points to them. Each skill covers package ownership, imports, feature registration, tests, and build wiring. The dependency skill covers catalog/version updates, lockfile refresh, validation, and advisory checks. Update skills when architecture changes.

Rewrite mkdocs installation, plugin, examples, formats, extensions, MVS, and clone/build instructions. Add package boundaries, composition, import conventions, migration, app setup, and adding-code guidance. Update source links and branch references. Docs and skills ship with the implementation.

### 8.3 CI and release checks

- Use pnpm with a frozen lockfile, a store cache, supported checkout/setup actions, and the declared minimum Node version plus the release LTS used for validation.
- Run typechecking, lint, package-cycle/value-cycle checks, unit tests, and app/example builds. Enforce module boundaries between lean entry points/runtime leaves and defaults/full catalogs, even within one package. Include extensions, servers, and CLI in boundary checks.
- Verify source-based esbuild app builds and compiled JS consumption. Install tarballs in clean consumers to check exports, declarations, direct dependencies, CSS/assets, and CLI bins.
- For each JSR package, run publication dry runs with `--allow-slow-types` and test its source/dependency graph with the pinned Deno version. Preserve public type precision through normal TypeScript checks and native npm declaration/consumer checks; fast-type compliance is not a v6 release gate.
- Run the slim-plugin acceptance example and the full Viewer; test snapshots with their required features registered.
- Check advisories with dependency review plus `pnpm audit --prod` or OSV; fail high/critical production findings. Track any justified exceptions explicitly.
- Build mkdocs for documentation changes. Before stable release, run the migrator and smoke-test `pdbe-molstar` and `rcsb-molstar`.

Use a root release script or configured Changesets workflow to version public packages together. On npm, publish `6.0.0-dev.N` under `dev`, then stable `6.0.0` under `latest`; scoped packages use public access. Publish matching versions of the validated JSR packages through the coordinated workflow in §5.6. Verify package-name availability and access on both registries before the first prerelease.

## 9. Migration from 5.x

### 9.1 Import map

These are default prefix mappings; the dependency-relocation audit supplies explicit exceptions. Root/index imports need export-map entries as well as prefix rewrites.

| 5.x prefix or API | 6.0 target |
| --- | --- |
| `molstar/lib/mol-util` | `@molstar/core/util` |
| `molstar/lib/mol-task` | `@molstar/core/task` |
| `molstar/lib/mol-data` | `@molstar/core/data` |
| `molstar/lib/mol-math` | `@molstar/core/math` |
| `molstar/lib/mol-state` | `@molstar/core/state` |
| `molstar/lib/mol-io` | `@molstar/io` |
| `molstar/lib/mol-model` | `@molstar/model` |
| `molstar/lib/mol-model-formats` | `@molstar/model/formats` |
| `molstar/lib/mol-model-props` | `@molstar/model/props` |
| `molstar/lib/mol-script` | `@molstar/model/script` |
| `molstar/lib/mol-gl`, `mol-geo`, `mol-theme`, `mol-repr`, `mol-canvas3d` | Corresponding `@molstar/graphics/gl`, `geo`, `theme`, `repr`, `canvas3d` |
| `molstar/lib/mol-plugin` | `@molstar/plugin` |
| `molstar/lib/mol-plugin-state` | `@molstar/plugin/state` |
| `molstar/lib/mol-plugin-ui` | `@molstar/plugin-ui` |
| `DefaultPluginSpec` | `@molstar/plugin/default-spec` |
| `DefaultPluginUISpec` | `@molstar/plugin-ui/default-spec` |
| `StateTransforms` | Leaf transformer imports, or `@molstar/plugin/state/transforms` |
| `molstar/lib/extensions/mvs` | `@molstar/mvs` for runtime; `@molstar/mvs-builder` for builder/schema APIs |
| `molstar/lib/extensions/<name>` | `@molstar/<name>-extension` |
| `molstar/lib/apps/viewer/app` | `@molstar/viewer` |

### 9.2 Migration tool

Ship `@molstar/migrate-6` with dry-run output and a report of unresolved/manual changes. It should:

1. Rewrite imports/re-exports using resolved paths and symbol-aware exceptions, including mixed imports of base spec helpers and default specs.
2. Replace in-repo cross-layer relatives with public package subpaths. Normalize `.js` suffixes on legacy package imports to the new extensionless exports.
3. In repository mode, resolve relative imports to their emitted `.js` paths, including directory indexes, and apply the new compiler settings. For downstream projects, respect their compiler/bundler convention; do not blindly apply the repository's relative-extension policy to every consumer.
4. Convert imports used solely as types to `import type` where analysis can establish that safely.
5. Update Mol* package dependencies and flag unsupported CommonJS usage, implicit default specs, relocated APIs, and full-catalog imports that prevent a slim bundle.

No promise of a fully automatic upgrade. Validate idempotence and representative transformations, then run the tool against `pdbe-molstar` and `rcsb-molstar` and smoke-test the resulting apps.

### 9.3 Compatibility contract

Keep transformer identifiers, snapshot JSON, feature/provider name strings, `PluginSpec.Action`/`Behavior`, and the `createPluginUI`/`Viewer.create` entry-point names. Snapshots still require their referenced features to be loaded. Explicit specs and registration replace implicit catalog loading in base entry points.

Remove CommonJS and legacy deep import paths. Existing applications must migrate imports and dependencies; installing `molstar@6` alone does not upgrade library consumers. Script-tag users continue using the CDN viewer asset paths and update the version. There is no compatibility-shim package or later shim sunset.

## 10. Implementation order

### 10.1 Release gate

Before 6.0 work lands on the default branch:

1. Land the planned 5.x work: particles, MVS changes, and bond-order perception.
2. Release the final 5.x feature minor.
3. Cut `v5` for any continuing fixes; freeze new features on that line.
4. Rename `master` to `main`, updating CI branch filters, docs, and clone instructions in the same change.

Implement 6.0 on `main`, publish `dev` prereleases, then release stable after the acceptance checks. Subsequent 5.x fixes stay on `v5`.

### 10.2 Technical phases

Keep the major workstreams in separate PRs. Each phase ends with a working build; validate compiled packages and source-based app bundles throughout.

Coordinate the [rendering-backend workstream](v6-webgpu.md#6-minimal-implementation-sequence) with dependency cleanup and packaging, and validate its contracts before freezing the v6 graphics API.

Keep the [fast-types workstream](v6-fasttypes.md#5-effort-and-adoption) deferred for possible v7 adoption. For v6, dual-registry release checks belong with packaging and CI, with slow types allowed on JSR and no associated annotation or generator migration phase.

| Phase | Work and exit condition |
| --- | --- |
| 0 | Verify package names; introduce pnpm, the workspace skeleton, and matching CI |
| 1 | Audit and relocate reverse dependencies in §3.2; break plugin cycles and add explicit type imports. Record final API locations and verify the proposed package graph, including declaration edges. Retain existing module settings and override `verbatimModuleSyntax: false` for the temporary CJS build if needed |
| 2 | Split transform/provider modules and default-spec entry points; remove internal facade imports and lazy getters; rename tests |
| 3 | Introduce features, empty registries, explicit base specs, and registry-aware presets; prove the slim example and full default composition |
| 4 | Drop CJS; set `type: module`, `NodeNext`, and `verbatimModuleSyntax`; use `.js` relative specifiers in source. Convert CommonJS globals/tooling and smoke-test emitted bins |
| 5 | Move into grouped workspace packages; add exports, project references, direct dependencies, and per-app esbuild. Verify clean builds and packed consumers; stage the CDN-only `molstar` package |
| 6 | Finish standalone MVS builder/runtime, extension and server/CLI packaging; verify npm/JSR builder parity and replacement instructions |
| 7 | Ship the migrator, skills, mkdocs updates, advisory checks, and downstream smoke tests; publish `6-dev`, then stable when ready |

Maintain the migration map, docs, skills, and relevant checks as each phase lands; phase 7 closes remaining release work. No legacy import shims are introduced at any phase.
