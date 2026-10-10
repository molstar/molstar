# Mol\* 6.0: packages, ESM, and plugin composition

Proposal against the `molstar@5.11.0` tree. The APIs and paths below describe the target, not features available in 5.x.
See the [short summary](summary.md) for the main decisions.

Current status (2026-10-10): the workspace/ESM split, TypeScript 7/Biome tooling, formatting, and plugin composition are
implemented. The [checklist](../plans/checklist.md) tracks remaining work. Alex owns rendering-backend implementation;
full WebGPU/Blender integration and fast-type migration remain deferred beyond v6. MolQL builder exposure, standalone
syntax/field-name validation, and translation of MolScript, PyMOL, VMD, and Jmol expressions are
[implemented](#72-molql-builder-and-validation-design).

## 1. Scope

Mol* 6.0 moves to a pnpm workspace of ESM packages grouped by layer. Parsers, representations, and themes become
explicit plugin registry entries; `DefaultPluginSpec()` remains the full built-in composition. Apps keep esbuild, with
dependencies and build configuration owned by each app.

Ship `@molstar/migrate-6-cli` with the release to handle mechanical import changes and report manual work. The unscoped
`molstar` package retains both Viewer and MVS Stories CDN apps, including their existing script-tag APIs. Existing MVS
HTML viewers loading `molstar@latest` must continue working without edits. Legacy `lib/mol-*` exports, library
compatibility re-exports, and the CommonJS build are removed.

The release also includes a standalone MolViewSpec builder, dependency-cycle removal, maintainer skills, updated
developer docs, and workspace CI.

Try JSR source publication with `--allow-slow-types`, beginning with the MVS builder, alongside native npm packages with
compiled ESM from the same release commit. Defer fast-type migration and `isolatedDeclarations` to consideration for v7.
The [fast-types and distribution analysis](fasttypes.md) preserves the audit for that decision; its annotation work and
agent estimates are outside the v6 scope.

Establish a rendering-backend boundary in `@molstar/graphics` for future WebGPU and other targets, retaining WebGL as
the working implementation. The [rendering-backend design](webgpu.md) covers the blast radius, minimal contracts,
migration, and validation; a production WebGPU renderer and feature parity are later work.

The same geometry/readback boundary should support portable scene extraction for a future
[Blender offline-rendering extension](webgpu.md#71-offline-rendering-with-blender). Offline rendering uses scene
snapshots and asynchronous jobs, separately from the interactive view contract. Keep Blender dependencies and
integration in an optional extension; implementing it is outside the 6.0 scope.

Out of scope:

- Independent package versions or a package per parser/representation.
- Native Node execution of TypeScript source and a repository-wide erasable-syntax rewrite.
- Changes to transformer identifiers or snapshot JSON, except where registration must become explicit.
- Moving Python `molviewspec` into this repository.
- Publishing a new stories library or merging the MolViewStories webapp. Retaining the existing MVS Stories CDN app and
  API is in scope.

## 2. Starting point

At the 5.11.0 proposal baseline, one package owned library code, apps, servers, and their npm dependencies. Consumers
used deep imports such as `molstar/lib/mol-plugin-ui`; the package had no `exports` map.

| Output          | Baseline build                                                                     |
| --------------- | ---------------------------------------------------------------------------------- |
| `lib/`          | `tsc`, then `tsc-alias` to complete relative JS paths; includes declarations       |
| `lib/commonjs/` | Separate `tsc` CommonJS emit; used by CLI bins and servers                         |
| `build/<app>/`  | esbuild from `src/`, with SCSS, copied assets, version injection, watch, and serve |

Tests used Jest and `esbuild-jest-transform`, with colocated `_spec/*.spec.ts` files. Asset copying and version patching
ran from root scripts.

`PluginSpec` already selected actions, behaviors, animations, and additional `customFormats`. Registry constructors
preloaded all built-in formats, representations, and themes. Large transform modules and default-spec imports connected
even a small plugin to the full catalog.

The baseline source graph was **not a layer DAG**. Besides plugin/state initialization cycles, lower folders contained
runtime dependencies on higher layers. The workspace refactor resolved those package edges explicitly, rather than
relying on moving files to remove them.

## 3. Packages and dependency boundaries

### 3.1 Package ownership

Keep recognizable subpaths, removing the `mol-` prefix. The table is the target ownership after the relocations in §3.2.

| Package                     | Owns                                                                                                         |
| --------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `@molstar/core`             | `util/`, `task/`, `data/`, `math/`, `state/`; no molecular or rendering dependencies                         |
| `@molstar/io`               | Readers and writers; depends on core                                                                         |
| `@molstar/query-language`   | MolQL expressions, builders, symbol tables, validation, and MolScript/PyMOL/VMD/Jmol text translation        |
| `@molstar/model`            | Structures, volumes, particles, domain formats, properties, and executable molecular queries                 |
| `@molstar/graphics`         | `gl/`, `geo/`, `theme/`, `repr/`, `canvas3d/`, plus rendering code relocated from model/math                 |
| `@molstar/plugin`           | Plugin and plugin-state together: runtime, spec helpers, and explicit default-spec/catalog entry points      |
| `@molstar/plugin-ui`        | React UI, UI-spec types, and an explicit default UI-spec entry point                                         |
| `@molstar/plugin-headless`  | Reusable Node headless plugin context, screenshot helpers, and output handling                               |
| `@molstar/mvs-builder`      | Standalone MolViewSpec schema, builder, serialization, and validation                                        |
| `@molstar/mvs`              | MolViewSpec runtime, plugin behavior and registry entries, and loader                                        |
| `@molstar/<name>-extension` | Individual extensions and their dependencies                                                                 |
| `@molstar/viewer`           | Published Viewer API and app                                                                                 |
| `@molstar/<name>-server`    | Model, volume, and plugin-state servers                                                                      |
| `@molstar/<name>-cli`       | Command-focused packages: `mvs-render-cli`, `cif2bcif-cli`, `cifschema-cli`, and `migrate-6-cli`             |
| `molstar`                   | Viewer and MVS Stories CDN apps at `build/viewer/` and `build/mvs-stories/`, with compatible script-tag APIs |

The library dependency direction is:

```text
plugin → graphics → model → io → core
                       └→ query-language → core
```

An arrow means “depends on.” Packages may also depend directly on any lower layer they import. UI, headless support,
extensions, apps, CLIs, and servers sit above the layers they use. `DefaultPluginUISpec` in plugin-ui extends
`DefaultPluginSpec` in plugin; plugin never depends on UI or headless support. Default compositions live behind explicit
subpath entry points in those packages, keeping React out of the non-UI default.

`mvs-builder` depends on the standalone `@molstar/query-language` language package, which depends only on core. Its
validation and authoring entry points do not load the molecular query runtime, plugin, or rendering code. `mvs` depends
on the builder and the runtime layers it imports. Servers that only need IO and model should not acquire plugin or
graphics through those packages.

### 3.2 Required relocations before packaging

The workspace refactor resolved the following baseline runtime edges, which could not be removed with `import type`.
Keep auditing source and declaration dependencies when changing package ownership; record replacement paths in the
migration map.

| Current source                                                       | Conflicting edge                   | Required work                                                                                                          |
| -------------------------------------------------------------------- | ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `mol-data/db/column.ts`                                              | core → IO number parser            | Move shared numeric parsing primitives into core and make IO consume them                                              |
| `mol-util/data-source.ts`                                            | core → IO string/UTF-8 helpers     | Move general byte/string helpers into core; keep format-specific parsing in IO                                         |
| `mol-math/geometry/gaussian-density/gpu.ts`                          | core → GL                          | Keep CPU/domain-independent math in core; move GPU computation and its GL-facing API to graphics                       |
| `mol-model/shape/shape.ts`                                           | model → geometry, themes, GL       | Separate shape data from rendering helpers; place graphics-dependent APIs in graphics                                  |
| `mol-model-formats/shape/*`                                          | model → mesh builders              | Move mesh conversion into graphics, retaining raw file readers in IO                                                   |
| `mol-model-props/**/representations`, `**/themes`, and label helpers | model → representations/themes     | Keep property computation in model; move visuals/themes to graphics and remove rendering dependencies from computation |
| Model particle/unit transform utilities                              | model → geometry transform helpers | Extract shared matrix/transform primitives into core or move the graphics-specific callers                             |

This is a starting list, not the full audit. Resolve declaration dependencies too: a type-only edge avoids runtime
evaluation, but still requires a package dependency and can create a cycle in `tsc -b` references. Assign exact
replacement subpaths during this work and record exceptions in the migration map.

Do not preserve an old path through a lower-layer re-export that recreates the dependency. General leaf paths remain
recognizable; APIs that cross these boundaries need explicit migration entries.

### 3.3 Plugin initialization cycles

Within the plugin layer:

- Use `import type` for parser result types in `PluginStateObject`, and for `PluginContext` wherever only its type is
  needed. Extract a context interface where it helps isolate construction; its dependencies must remain type-only.
- Separate the `PluginBehavior` contract from the classes extending `PluginStateObject`. Behavior detection uses
  `typeClass === 'Behavior'` without importing those classes.
- Move default specs out of base spec/context modules into `@molstar/plugin/default-spec` and
  `@molstar/plugin-ui/default-spec`. Base construction takes an explicit spec.
- Split transform modules into leaves. All transformer consumers import individual defining modules.
- Remove the `StateTransforms` convenience facade and its lazy getters; do not replace it with another aggregate object.

Enforce an acyclic package graph and no value-import cycles within packages. Type-only cycles inside a package are
allowed; they must not conceal a package cycle.

### 3.4 Headless plugin support

Move `HeadlessPluginContext` and `HeadlessScreenshotHelper` from `mol-plugin` into `@molstar/plugin-headless`, under
`packages/plugin/headless/`. This is a reusable Node library, separate from plugin-ui and command-line wrappers. It owns
Node filesystem/output handling and headless setup, depending on plugin, graphics, and the lower layers it imports.
Plugin, browser MVS runtime, and graphics must not depend back on it.

Retain native-module injection for embedding applications. The headless package supplies the Node environment and
external modules to the graphics backend; shared device/resource/capture contracts remain in graphics. Adapt screenshots
to the [backend capture/readback design](webgpu.md) as that boundary is extracted, rather than adding another renderer
abstraction here.

Keep MP4 integration in an explicit module of the MP4 extension, composed by the rendering CLI or embedding application.
The base headless context must not import the encoder or register the extension automatically. Record migration of
existing `getAnimation`/`saveAnimation` calls to the explicit integration. Image rendering and snapshot output should
work without loading video support.

## 4. Workspace, imports, and dependencies

### 4.1 Layout

The packaging-first [prototype plan](../plans/workspace-prototype.md) implemented this layout. Plugin composition was
implemented in a subsequent workstream; rendering-backend changes are owned separately by Alex.

```text
packages/
  core/                       # @molstar/core
  io/
  model/
  graphics/
  plugin/
    core/                     # @molstar/plugin
    ui/                       # @molstar/plugin-ui
    headless/                 # @molstar/plugin-headless
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
cli/<name>/                   # @molstar/<name>-cli (mvs-render, cif2bcif, cifschema, migrate-6)
distributions/
  molstar/                    # CDN-only package
scripts/                      # shared build/release tooling
.agents/                      # maintainer skills
```

Each library package has `src/`, `lib/`, `package.json`, and a composite `tsconfig.json`. Command-focused packages use
the `-cli` suffix; their executable names need not. Libraries that also offer commands, such as the MVS builder, retain
their domain names, as do server packages. Lightweight MVS validation/schema commands stay in the builder; the native
rendering stack lives in a separate `cli/mvs-render/` package. See §7.1 for the complete command map.

```yaml
# pnpm-workspace.yaml
packages:
  - 'packages/*'
  - 'packages/plugin/*'
  - 'packages/mvs/*'
  - 'extensions/*'
  - 'apps/*'
  - 'examples/*'
  - 'servers/*'
  - 'cli/*'
  - 'distributions/*'
  - 'smoke'
```

Use physical moves into the workspace. Keep relocation commits separate from import rewrites where practical; keep each
completed PR buildable.

### 4.2 One import convention

Use the same public package specifiers inside and outside the repository:

```ts
// packages/graphics/src/canvas3d/passes/illumination.ts
import { ValueCell } from '@molstar/core/util/value-cell';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
```

Rules:

1. Crossing a layer-root folder (`util`, `task`, `gl`, `canvas3d`, plugin `state`, etc.) uses an extensionless package
   subpath, even within the same package.
2. Relative source imports stay within that folder and use the emitted `.js` path: `./names.js`, `../hex.js`, or
   `./controls.js` for a TSX module compiled with `react-jsx`. TypeScript and esbuild resolve these to source files
   while building; emitted JS retains the specifiers.
3. Package self-imports resolve through `exports`. Any development `paths` mappings must describe the same public
   subpaths and source files.
4. Apps, extensions, examples, and servers never reach into another package through relative filesystem paths.

Lint these boundaries. Bare specifiers such as `@molstar/core/util/color` and relative `.js` specifiers remain unchanged
in emitted JS. Resolve directory imports explicitly to an index file or an exported subpath.

### 4.3 Dependency declarations

Publish all public Mol* workspace packages at one version. Declare internal dependencies with `workspace:*`, which
becomes an exact version on publication. A consumer should not have to install each lower layer manually. Mixing
separate Mol* release versions in one application can still duplicate state/transformer identities; lockstep publishing
does not prevent consumers from requesting conflicting versions.

Each package declares every external package it directly imports. “Owned by core” does not make a transitive dependency
available to another package. Shared versions live in the pnpm catalog; root `package.json` is private tooling.

| Dependency                                          | Ownership                                                      |
| --------------------------------------------------- | -------------------------------------------------------------- |
| `rxjs`, `immutable`, `mutative`                     | Every direct importer, starting with core                      |
| `tslib`                                             | Each package whose emit uses imported helpers                  |
| `argparse`                                          | CLI/server packages that parse arguments                       |
| `express`, `compression`, `cors`, `swagger-ui-dist` | Server packages                                                |
| `h264-mp4-encoder`                                  | MP4 export extension                                           |
| `io-ts`                                             | MVS builder, plus any remaining direct importer until migrated |
| `react-markdown`, `remark-gfm`                      | Plugin UI                                                      |

Use peers for React/React DOM on UI packages. Keep injected headless modules (`gl`, `canvas`, `pngjs`, `jpeg-js`)
optional for consumers of `@molstar/plugin-headless`, declaring optional peers where needed. The ready-to-run
`@molstar/mvs-render-cli` declares its native rendering modules and codecs as direct dependencies and supplies them
explicitly; browser plugin/MVS packages must not acquire them. Declare optional cloud storage where used. Do not add
React to plugin or core.

Types needed only to build a package belong in its dev dependencies. Types referenced by published declarations must be
available to consumers through dependencies or declared peers. In particular, React does not install `@types/react`: UI
packages need an explicit consumer-facing type dependency/peer policy, verified with a clean TypeScript consumer.

### 4.4 No barrel files

Do not ship convenience barrel files. Consumers, examples, and internal code import defining modules through supported
subpaths. Package export maps provide those paths directly; no root or directory aggregation module is needed. This
removes a source of accidental broad imports rather than relying on every caller to choose correctly.

Barrels do not inherently prevent production tree-shaking: ESM bundlers can remove unused exports when usage and side
effects permit it. They can still enlarge the graph processed by tooling and retain initialization effects.
[Vite's warning](https://vite.dev/guide/performance#avoid-barrel-files) concerns extra development-time fetching and
transformation; [esbuild's documentation](https://esbuild.github.io/api/#tree-shaking) describes conditional removal of
unused code. The policy avoids the aggregation risk without assuming every barrel defeats every bundler.

- **Remove aggregation modules**, including named re-export collections, `export *` chains, and type-only convenience
  barrels. Use `import type` from the defining module for types. Exporting locally defined symbols remains normal module
  structure.
- **Preserve cohesive implementation entry points.** A file named `index.ts` is not inherently a barrel. Split mixed
  implementation/re-export modules such as `mol-util/index.ts` into suitable defining modules and record migrated paths.
  An export-map alias may point directly to an implementation file without adding a wrapper barrel.
- **Do not substitute aggregate convenience objects.** Remove `StateTransforms` and migrate consumers to transformer
  leaves. Preserve transformer identifiers and required registration behavior, not the aggregate access syntax.
- **Keep deliberate composition explicit.** Default specs and registration catalogs assemble complete provider sets for
  app composition. They are not general symbol-import entry points. Base runtime modules and individual providers must
  not depend on them; enforce that boundary even within a package.

Enforce the policy through lint/import-graph checks and the public export inventory. Validate that a leaf import cannot
reach unrelated providers, default catalogs, or optional backends through re-exports. Reuse the slim-plugin fixture to
check both processed modules and production output, with a downstream Vite development smoke/profile case alongside
esbuild. Audit real registration and asset side effects before adding purity annotations or `sideEffects: false`; no
consumer barrel-rewriting plugin is required.

The existing CDN app globals are compatibility APIs at the final bundle boundary. Preserve their exported names,
including existing re-exported values, as required by §9.3. This does not introduce library barrels: library code must
not import the app entry points or depend on their global objects.

## 5. ESM and TypeScript source

### 5.1 Supported execution modes

Publish ESM JavaScript and declarations in `lib/`, plus source in `src/` for bundlers and debugging. Set
`"type": "module"` on published packages and retain the existing Node **22 or later** baseline; raise it only if runtime
or tooling requirements demand it.

| Consumer                          | Resolution                 | Executes                       |
| --------------------------------- | -------------------------- | ------------------------------ |
| Installed library, CLI, or server | Default `import` condition | `lib/*.js`                     |
| TypeScript consumer               | `types` condition          | Checks `lib/*.d.ts`            |
| In-repo esbuild                   | `molstar-src` condition    | Compiles source, including TSX |

Node executes compiled JavaScript. Keep all published bins on `lib/*.js`; workspace tooling runs as JavaScript or is
compiled before execution. `molstar-src` selects source for bundlers and does not promise native Node execution.
Publishing source does not require erasable syntax.

Validated JSR packages expose TypeScript source for Deno or compatible tooling. Deno consumers of native npm packages
use the compiled exports above. Publishing to either registry does not make browser, Node, or optional native APIs
available in every runtime; see the [distribution matrix](fasttypes.md#6-typescript-distribution-through-npm).

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

Use `.js` relative specifiers in TypeScript source as described in §4.2. No extension-rewriting compiler option is
needed. Keep `verbatimModuleSyntax` and explicit `import type` for predictable module dependencies.

During dependency and module refactoring, keep the existing ESM build settings. If `verbatimModuleSyntax` is enabled in
the shared config while CommonJS still exists, explicitly set it to `false` in `tsconfig.commonjs.json`. Otherwise ESM
syntax fails in that build with TS1287/TS1295. Switch to `NodeNext` and remove the override with the CommonJS build in
the ESM phase. See
[TypeScript’s module-syntax behavior](https://www.typescriptlang.org/tsconfig/verbatimModuleSyntax.html).

### 5.3 TypeScript syntax and performance

Do not enable `erasableSyntaxOnly`. Both library and app builds compile TypeScript, so existing namespaces, enums,
parameter properties, and other compiler-supported syntax can remain. Refactor individual constructs only where module
boundaries, composition, or ESM compatibility require it. There is no blanket namespace/enum conversion or syntax-driven
CIF schema regeneration phase.

Keep hot `const enum`s under the existing `isolatedModules` constraints. If a necessary refactor changes their use or
emit, inspect the generated code and benchmark the affected parse/render paths. Do not assume that replacing enum uses
with object properties or module constants preserves performance.

Do not require `isolatedDeclarations` or a broad annotation migration in v6. Preserve existing inferred API precision,
including parameter/schema keys, literal unions, overloads, and factory constructor types. Consider fast types for v7
using the [audit](fasttypes.md#3-measured-blast-radius); do not redesign public contracts solely to satisfy JSR fast
types in this release.

### 5.4 Convert runtime CommonJS assumptions

Changing compiler flags does not fix `require`, `__dirname`, or `__filename`. Audit source and root scripts before
setting `"type": "module"`:

- Convert module loads to ESM imports; use `createRequire(import.meta.url)` only where CommonJS loading is needed, such
  as an optional native dependency.
- Replace path globals with `import.meta.url`-based paths and account for moved source, emitted files, and packaged
  assets.
- Update conditional loads such as `servers/model/preprocess.ts`, native loading in `cli/mvs/mvs-render.ts`, cloud
  storage loading, and `cifschema` data paths.
- Convert `scripts/clean.js` and `scripts/deploy.js` to ESM, or explicitly retain tooling as `.cjs` and update callers.
  Keep the published library ESM-only.

Smoke-test emitted CLI/server entry points and verify packaged data assets. Typechecking alone will not catch these
runtime failures.

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
      "types": "./lib/task/task.d.ts",
      "molstar-src": "./src/task/task.ts",
      "import": "./lib/task/task.js"
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

Add conditional entries for cohesive implementation entry points and `.tsx` sources; the generic `*.ts` pattern does not
cover TSX. Root or directory aliases must point directly to defining modules under the
[no-barrel policy](#44-no-barrel-files). Public code uses `@molstar/core/util/color`, without `.js`. Appending `.js`
would make this wildcard target `color.js.js`.

Publish only intended source/assets and generated output. Exclude `_test/` and fixtures with a pack step or appropriate
nested ignore files; inspect the actual tarball. Export UI skins and built CSS explicitly. Mark modules
`sideEffects: false` only after auditing initialization behavior, and retain CSS/asset side effects where needed.

Use project references for package builds. Source conditions drive esbuild; normal consumer types resolve to generated
declarations. Prove both from a clean checkout and from packed packages, rather than depending on stale `lib/` output or
root hoisting.

`lib/` is generated and ignored by Git, but its JavaScript, declarations, and required assets ship in the npm tarball
and remain in the installed package. It is not merely a temporary input that packing removes. JSR source artifacts
exclude it. See the [build-output distinction](fasttypes.md#8-is-lib-only-temporary).

### 5.6 npm and JSR publication

Keep one source tree and derive npm/JSR manifests from one package/export inventory. npm retains compiled conditional
exports, peers, bins, and source for opt-in bundlers. JSR uses explicit source exports and resolved dependency mappings.
Its publisher supports `package.json` projects with `.js`-to-TypeScript import resolution; validate self-subpaths, TSX,
assets, and exact dependencies before relying on that path. JSR's npm compatibility layer does not replace native
publication to npmjs.com.

Try publication with `deno publish --dry-run --allow-slow-types`, then use the same allowance when publishing validated
packages to JSR. Slow types can degrade JSR-generated documentation and npm-compatibility declarations and make consumer
checking slower. Keep native npm declarations generated by `tsc`, and verify JSR source consumers without promising
equivalent generated documentation/types. The allowance does not skip normal typechecking or other publication
requirements.

Publish validated packages to both registries at the same version from the same commit, starting with the MVS builder.
Extend JSR coverage in dependency order without implying universal runtime support. Keep per-registry completion records
and retry partial releases from unchanged artifacts; the two registries cannot publish atomically. The
[dual-publication design](fasttypes.md#7-publishing-to-npm-and-jsr-together) covers normalization, validation, package
coverage, and the `deno pack` alternative. Retain pnpm/`tsc -b` and its complete declarations for native npm packaging.

## 6. Plugin composition

The detailed design is in [plugin-composition.md](plugin-composition.md), and the ordered steps, call-site inventories,
and acceptance checklist are in the [implementation plan](../plans/plugin-composition.md); this section summarizes them.

- `PluginContext` starts with empty registries. A new `PluginSpec.registry` field lists declarative
  `PluginRegistryEntry` records (formats, representations and themes per scope, presets, selection queries, loci labels,
  markdown extensions, drag-and-drop handlers, actions, animations). `init()` registers them in order without resolving
  dependencies, so creating a plugin from a spec stays cheap.
- The spec keeps `behaviors`, `config`, `canvas3d`, and `layout`. Spec-level `actions`, `animations`, and
  `customFormats` move into registry entries, and `DataFormatProvider` gains a required `name`.
- What is included is what the spec imports and lists. Entries have no `requires`; each lists the providers it directly
  needs, and registries count references by provider identity so overlaps are safe. Composition uses only static
  imports, so single-file builds tree-shake like split builds.
- Behaviors are unchanged. They can reuse the entry format through `plugin.register(entry)`, which returns a function
  that undoes the registration.
- Presets and other code import the providers they run and include them in their entries. Policy choices, such as which
  preset a format applies after loading, go through config and the registries. Presets resolve by `id` or an optional
  `alias`; built-in presets keep their short keys as aliases.
- Unregistered representation and theme names keep today's lenient substitution of the registry default, with a warning
  where the plugin sees the name; an empty registry is an error.
- PyMOL, VMD, and Jmol script support is enabled by importing its transpiler module, like transformer registration.
- Transformers stay in the global registry when their module is imported. The `StateTransforms` facade is removed and
  transform modules are split by functionality, with format-specific transformers next to their format providers.
  Snapshot loading checks all transformer ids before changing any state.
- `BuiltInTrajectoryFormat` and sibling name types remain, derived in catalog modules and used elsewhere only through
  `import type`.
- `DefaultPluginSpec` (`@molstar/plugin/default-spec`) and `DefaultPluginUISpec` (`@molstar/plugin-ui/default-spec`)
  assemble the full built-in registry. Base entry points, including `createPluginUI`, take an explicit spec and never
  import defaults or catalogs as values.

The slim-plugin acceptance target (render an SDF ligand without unrelated parsers, representations, presets, or MP4
export in the import graph and bundle) is defined in [plugin-composition §12](plugin-composition.md#12-acceptance).

## 7. MolViewSpec, extensions, and apps

Split `extensions/mvs` into:

| Package                | API and responsibilities                                                                                       |
| ---------------------- | -------------------------------------------------------------------------------------------------------------- |
| `@molstar/mvs-builder` | Schema, `createMVSBuilder`, `MVSData`, MVSJ/MVSX serialization/validation, `mvs-validate`, schema-printing CLI |
| `@molstar/mvs`         | `loadMVS`/`loadMVSData`, plugin behavior and registry entries, annotations, cameras, runtime representations   |

The builder replaces [molviewspec-ts](https://github.com/molstar/mol-view-spec/tree/master/molviewspec-ts) and the JSR
`@molstar/molviewspec` distribution after API/parity checks. Publish the replacement to npm and JSR. Own `io-ts` and any
archive dependencies in the builder; inline or extract its small general helpers without introducing a Mol* runtime
dependency. An optional core peer would not make unconditional core imports optional.

The runtime imports the builder as a dependency; do not bundle a second copy into it. Python stays in the
[mol-view-spec repository](https://github.com/molstar/mol-view-spec).

Other extensions become `@molstar/<name>-extension`. Each owns its direct dependencies, exports registry
entries/behaviors, and imports UI only when needed. Viewer dependencies and imports define its extension set; a smaller
app declares and imports only its selected extensions.

`@molstar/viewer` is published. Docking viewer, mesoscale explorer, MVS Stories, and examples remain private workspace
app packages; the built MVS Stories app is nevertheless distributed through the root `molstar` package. Preserve that
app's existing browser API. Moving parts of [MolViewStories](https://github.com/molstar/mol-view-stories) into a future
`packages/mvs/stories` library needs a separate plan.

### 7.1 CLI packages and executable names

Use `@molstar/<name>-cli` when the package's public interface is a command. Keep existing executable names unchanged.
Packages that primarily provide a library or server retain their domain names even when they expose bins; do not create
a separate package for every executable or one umbrella CLI that installs all tools' dependencies.

| Package                   | Workspace location      | Executables                                                     |
| ------------------------- | ----------------------- | --------------------------------------------------------------- |
| `@molstar/mvs-render-cli` | `cli/mvs-render/`       | `mvs-render`                                                    |
| `@molstar/cif2bcif-cli`   | `cli/cif2bcif/`         | `cif2bcif`                                                      |
| `@molstar/cifschema-cli`  | `cli/cifschema/`        | `cifschema`                                                     |
| `@molstar/migrate-6-cli`  | `cli/migrate-6/`        | `molstar-migrate-6` (new)                                       |
| `@molstar/mvs-builder`    | `packages/mvs/builder/` | `mvs-validate`, `mvs-print-schema`                              |
| `@molstar/model-server`   | `servers/model/`        | `model-server`, `model-server-query`, `model-server-preprocess` |
| `@molstar/volume-server`  | `servers/volume/`       | `volume-server`, `volume-server-query`, `volume-server-pack`    |

`@molstar/mvs-render-cli` composes `@molstar/mvs`, `@molstar/plugin-headless`, the MP4 extension, and its native
modules/codecs. It owns argument parsing, file processing, and output selection. `@molstar/mvs` retains the reusable
browser-capable runtime and never imports the rendering CLI or headless package. Validation/schema commands must remain
usable without native rendering dependencies.

Publish bins as compiled ESM JavaScript with Node shebangs and explicit `package.json` `bin` entries, for example
`"mvs-render": "./lib/index.js"`. CLI directories use the same source/build layout and dependency rules as other
workspace packages. Existing development-only generators and diagnostics may remain private workspace tools until
deliberately promoted to supported commands.

Library and CLI installation moves away from the CDN-only root `molstar` package. Document the replacement owner of
every existing bin. For example, after publication, `npm exec --package @molstar/mvs-render-cli -- mvs-render --help`
selects the rendering package explicitly. Preserve current flags and supported output formats unless a separate
migration entry records a deliberate change.

### 7.2 MolQL builder and validation design

Status: implemented (2026-10-10). Validation covers syntax and bad field names against symbol and argument-definition
tables without the molecular query runtime. Argument type checking, required-argument checks, and selection-result
checking remain outside this pass.

#### Package boundary

`@molstar/query-language` owns the expression language, builder, symbol/argument definitions, text parsers/transpilers
for MolScript, PyMOL, VMD, and Jmol, and standalone validation. Structure-dependent compilation and evaluation stay in
`@molstar/model`. The builder's duplicate expression representation has been removed in favor of the language package's
canonical type. Moved source paths and split symbols are recorded in the migration inventories.

```text
@molstar/mvs-builder → @molstar/query-language
@molstar/model       → @molstar/query-language
@molstar/query-language  → @molstar/core
@molstar/mvs         → @molstar/mvs-builder + @molstar/model + its existing runtime dependencies
```

The split avoids a package cycle: structure bundles, schemas, and loci use the expression builder, while loci and
computed properties also use the executable query compiler. The existing compiler constructs molecular query contexts
and evaluates constant operations, so it is unsuitable for standalone validation. The pre-extraction source graph had
282 processed modules for executable compilation versus nine for JSON authoring; these are module counts, not timing
measurements. Compilation also accepted unknown argument names, and 21 of the standard table's 163 symbols had no
default runtime implementation. Language validity remains separate from execution support.

Text parsers/transpilers have explicit entry points. JSON builder and validation imports do not reach them. The
executable compiler and custom-property runtime registration remain model-owned.

#### Text expressions to MolQL

`compileScript(language, source)` in `@molstar/query-language/compile` returns a serializable MolQL `Expression` for any
of the four supported languages. It parses/translates text, then validates callable and argument names. It does not
compile executable queries or evaluate constants. Existing supported syntax and unsupported-feature errors are
preserved; this does not promise full parity with the original PyMOL, VMD, or Jmol engines.

The text-only `Script` API lives at `@molstar/query-language/script`. Its `toExpression` dispatch preserves explicit
registration: import `@molstar/query-language/transpilers/<lang>` or `transpilers/all` to enable non-MolScript languages
for plugin use. Direct `transpilers/<lang>/parser` imports work without registration. `compileScript` imports all four
translators directly, so MVS authoring never depends on prior plugin initialization and does not modify the language
registry. Model's `Script.toQuery`, `toLoci`, and `getStructureSelection` keep their molecular behavior.

#### Validation pass

`expressionValidationIssues(expression, options?)` in `@molstar/query-language/language/validation` returns
path-qualified issues or `undefined`. It walks expressions without creating a structure context or executing operations.

- Accept finite JSON literals, symbol objects, and applications with positional arrays or argument maps.
- Reject unexpected expression fields and malformed heads; an application head must be a symbol.
- Resolve callable names and nested calls against the standard symbol table.
- Reject dictionary keys outside `MSymbol.args.map`, including misspelled names and undeclared positional indices.
  Positional array indices correspond to dictionary keys such as `'0'` and `'1'`. Variadic lists accept nonnegative
  positional indices, including numeric map keys, and reject named keys.
- Preserve bare symbol references as string values, matching the existing language semantics.
- Reject cycles, sparse arrays, and array fields that would be lost during JSON serialization.

`options.getSymbol(name)` supplies an explicit replacement vocabulary for custom definitions; callers can fall back to
`SymbolMap[name]`. `options.syntaxOnly` checks shapes without resolving names. The language package does not read the
runtime's global table.

Fixed MVS schema codecs enforce expression syntax and the existing application-root rule. `MVSData.validationIssues` and
`isValid` additionally validate names throughout scene and animation trees, multiple snapshots, and primitive positions
with `structure_ref`. Supply custom definitions through `options.getMolQLSymbol`. The default CLI knows only the
standard vocabulary and reports custom callable names as unknown. CLI failures include the file and expression path,
return a nonzero status, and do not prevent later files from being checked. MVS runtime validation uses registered
runtime symbol definitions, then retains executable compilation as the execution-support check.

#### MVS authoring API and checks

`@molstar/mvs-builder/molql` exports `MolScriptBuilder` and `compileScript`. This narrow authoring facade is a
deliberate exception to the no-convenience-re-export policy in §4.4. Direct defining-module imports from the language
package remain available; schema/serialization entry points do not import the facade or text parsers.

Unit tests cover syntax, nested callable names, array/map argument keys, custom vocabularies, MVS trees/snapshots, and
all four translators, including their existing example corpora. `scripts/workspace/query-language-check.mjs` verifies
source import boundaries and per-language parser independence. Packed smoke consumers exercise all languages and CLI
success/failure with only the standalone builder's dependency closure installed. Standalone-builder parity with
molviewspec-ts and coordinated JSR publication remain separate release tasks.

## 8. Builds and maintenance

### 8.1 Library and app builds

Use `tsc -b` with a root solution config and per-package references, `rootDir: src`, and `outDir: lib`. Cover all
published TypeScript packages, including extensions, CLI, servers, and Viewer. Copy assets per package and generate
plugin version source once for both library and app builds.

Apps use `scripts/esbuild/app.mjs` for shared plugins, source conditions, watch, and serve. Each app/example owns its
entry, output paths, dependencies, themes, and scripts. Resolve imports from the importing workspace package; do not
depend on a root `Apps` list or root dependency hoisting.

```text
pnpm --filter @molstar/viewer dev
pnpm -r --filter "./apps/**" --filter "./examples/**" build
```

Keep SCSS, copied HTML/images/icons, version injection, IIFE globals, and existing CDN output names. Shaders remain
`.glsl.ts` strings. The deploy script consumes the app outputs after its ESM/tooling conversion. Build and stage both
apps into `molstar/build/viewer/` and `molstar/build/mvs-stories/` when packing the root CDN package, including their
JS/CSS, HTML, and supporting assets. The root package's file allowlist must retain both directories; scoped-package
moves must not change these URLs.

### 8.2 Tests, skills, and docs

Rename `_spec/*.spec.ts` to `_test/*.test.ts`; keep tests colocated. Reserve “spec” for future agent/behavior
specifications. Choose a workspace runner with working ESM and export-condition resolution, remove Jest's
`moduleDirectories: ["lib"]`, and retain browser/WebGL coverage. Install headless dependencies on the packages that run
those tests.

Add `.agents/` skills for `add-extension`, `add-example-app`, `add-app`, `add-format`, `add-representation`,
`add-server`, and `update-dependencies`. Root `AGENTS.md` points to them. Each skill covers package ownership, imports,
registry entries, tests, and build wiring. The dependency skill covers catalog/version updates, lockfile refresh,
validation, and advisory checks. Update skills when architecture changes.

Rewrite mkdocs installation, plugin, examples, formats, extensions, MVS, and clone/build instructions. Add package
boundaries, composition, import conventions, migration, app setup, and adding-code guidance. Update source links and
branch references. Docs and skills ship with the implementation.

### 8.3 CI and release checks

- Use pnpm with a frozen lockfile, a store cache, supported checkout/setup actions, and the declared minimum Node
  version plus the release LTS used for validation.
- Run typechecking, lint, package-cycle/value-cycle checks, unit tests, and app/example builds. Reject aggregation
  barrels and enforce module boundaries between lean entry points/runtime leaves and defaults/full catalogs, even within
  one package. Include extensions, servers, and CLI in boundary checks.
- Verify source-based esbuild app builds and compiled JS consumption. Install tarballs in clean consumers to check
  exports, declarations, direct dependencies, CSS/assets, and CLI bins. Validate defining-module imports and their
  processed graphs/production output, with a Vite development smoke/profile case for downstream use.
- Verify the headless/CLI split: plugin and MVS browser imports do not reach headless/native/video modules; builder
  validation/schema bins run without native render dependencies. Test the packed rendering CLI's image, snapshot, and
  MP4 outputs on supported headless environments, plus basic headless capture without the MP4 integration.
- For each JSR package, run publication dry runs with `--allow-slow-types` and test its source/dependency graph with the
  pinned Deno version. Preserve public type precision through normal TypeScript checks and native npm
  declaration/consumer checks; fast-type compliance is not a v6 release gate.
- Run the slim-plugin acceptance example and the full Viewer; test snapshots with their required transformers imported
  and providers registered.
- Before advancing `molstar@latest`, test existing Viewer/MVS HTML fixtures against the packed root package, routing
  their unchanged CDN URLs to the candidate assets. Verify classic script loading, globals/API calls, custom-element
  registration, CSS/assets, MVSJ/MVSX loading, and independent named story contexts. Inspect both CDN directories in the
  tarball; file presence alone does not prove API compatibility.
- Check advisories with dependency review plus `pnpm audit --prod` or OSV; fail high/critical production findings. Track
  any justified exceptions explicitly.
- Build mkdocs for documentation changes. Before stable release, run the migrator and smoke-test `pdbe-molstar` and
  `rcsb-molstar`.

Use a root release script or configured Changesets workflow to version public packages together. On npm, publish
`6.0.0-dev.N` under `dev`, then stable `6.0.0` under `latest`; scoped packages use public access. Publish matching
versions of the validated JSR packages through the coordinated workflow in §5.6. Verify package-name availability and
access on both registries before the first prerelease.

## 9. Migration from 5.x

### 9.1 Import map

These are default prefix mappings; the dependency-relocation audit supplies explicit exceptions. Root/index imports need
export-map entries as well as prefix rewrites.

| 5.x prefix or API                                                        | 6.0 target                                                                                 |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ |
| `molstar/lib/mol-util`                                                   | `@molstar/core/util`                                                                       |
| `molstar/lib/mol-task`                                                   | `@molstar/core/task`                                                                       |
| `molstar/lib/mol-data`                                                   | `@molstar/core/data`                                                                       |
| `molstar/lib/mol-math`                                                   | `@molstar/core/math`                                                                       |
| `molstar/lib/mol-state`                                                  | `@molstar/core/state`                                                                      |
| `molstar/lib/mol-io`                                                     | `@molstar/io`                                                                              |
| `molstar/lib/mol-model`                                                  | `@molstar/model`                                                                           |
| `molstar/lib/mol-model-formats`                                          | `@molstar/model/formats`                                                                   |
| `molstar/lib/mol-model-props`                                            | `@molstar/model/props`                                                                     |
| `molstar/lib/mol-script`                                                 | `@molstar/query-language` for language/text; `@molstar/model/script` for molecular runtime |
| `molstar/lib/mol-gl`, `mol-geo`, `mol-theme`, `mol-repr`, `mol-canvas3d` | Corresponding `@molstar/graphics/gl`, `geo`, `theme`, `repr`, `canvas3d`                   |
| `molstar/lib/mol-plugin`                                                 | `@molstar/plugin`                                                                          |
| `molstar/lib/mol-plugin/headless-plugin-context`                         | `@molstar/plugin-headless/headless-plugin-context`                                         |
| `molstar/lib/mol-plugin/util/headless-screenshot`                        | `@molstar/plugin-headless/util/headless-screenshot`                                        |
| `molstar/lib/mol-plugin-state`                                           | `@molstar/plugin/state`                                                                    |
| `molstar/lib/mol-plugin-ui`                                              | `@molstar/plugin-ui`                                                                       |
| `DefaultPluginSpec`                                                      | `@molstar/plugin/default-spec`                                                             |
| `DefaultPluginUISpec`                                                    | `@molstar/plugin-ui/default-spec`                                                          |
| `StateTransforms`                                                        | Individual transformer imports from defining modules                                       |
| `molstar/lib/extensions/mvs`                                             | `@molstar/mvs` for runtime; `@molstar/mvs-builder` for builder/schema APIs                 |
| `molstar/lib/extensions/<name>`                                          | `@molstar/<name>-extension`                                                                |
| `molstar/lib/apps/viewer/app`                                            | `@molstar/viewer`                                                                          |

### 9.2 Migration tool

Ship `@molstar/migrate-6-cli`, exposing `molstar-migrate-6`, with dry-run output and a report of unresolved/manual
changes. It should:

1. Rewrite imports/re-exports using resolved paths and symbol-aware exceptions, including mixed imports of base spec
   helpers and default specs.
2. Replace in-repo cross-layer relatives with public package subpaths. Normalize `.js` suffixes on legacy package
   imports to the new extensionless exports.
3. In repository mode, resolve relative imports to their emitted `.js` paths, including directory indexes, and apply the
   new compiler settings. For downstream projects, respect their compiler/bundler convention; do not blindly apply the
   repository's relative-extension policy to every consumer.
4. Convert imports used solely as types to `import type` where analysis can establish that safely.
5. Update Mol* package dependencies and migrate removed barrel/facade imports to defining modules where symbol
   resolution is unambiguous. Report dynamic `StateTransforms` access and re-export side effects for manual migration.
   Flag unsupported CommonJS usage, implicit default specs, relocated APIs, and unintended full-catalog imports.
6. Migrate plugin specs: rewrite removed `actions`/`animations`/`customFormats` fields into `registry` entries, add
   `DefaultRegistry` to literal specs that relied on the implicit catalog, filter the matching default entry in
   spread-and-override specs, keep the full structure tools when a spec replaces `components`, and report specs and
   script languages it cannot resolve. The
   [plugin-composition plan](../plans/plugin-composition.md#5-migration-tool-additions) lists the rules.

No promise of a fully automatic upgrade. Validate idempotence and representative transformations, then run the tool
against `pdbe-molstar` and `rcsb-molstar` and smoke-test the resulting apps.

### 9.3 Compatibility contract

Keep transformer identifiers, snapshot JSON, provider name strings, `PluginSpec.Action`/`Behavior`, and the
`createPluginUI`/`Viewer.create` entry-point names. Snapshots still require their transformers to be imported; a
representation or theme name that is not registered falls back to the registry default with a warning. Explicit specs
and registration replace implicit catalog loading in base entry points.

Remove CommonJS and legacy deep import paths. Existing applications using library imports must migrate imports and
dependencies; installing `molstar@6` alone does not upgrade those consumers. There is no library compatibility-shim
package or later shim sunset.

The root `molstar` package preserves the existing CDN app contract:

| App         | Retained paths                                                                                  | Browser API                                      |
| ----------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------ |
| Viewer      | `build/viewer/`, including existing JS/CSS and supporting files                                 | Existing `molstar` global and Viewer API (below) |
| MVS Stories | `build/mvs-stories/`, including `mvs-stories.js`, `mvs-stories.css`, HTML, and supporting files | Existing `mvsStories` global and custom elements |

Accepted carve-out for the Viewer global: the library values under `molstar.lib.plugin` (`DefaultPluginSpec`,
`DefaultPluginUISpec`, `StateActions`, and the added `DefaultRegistry`) follow the 6.0 library API. `StateTransforms`
stays on the global with its 5.x keys and members, assembled from the split modules. Hand-built specs that use the
removed `actions`/`animations`/`customFormats` fields throw, and literal specs get empty registries. `Viewer.create` and
its options, including `customFormats`, and the Viewer's methods are unchanged
([plugin-composition §11](plugin-composition.md#11-extensions-mvs-and-apps)).

For MVS Stories, preserve the exports in the [app entry point](../src/apps/mvs-stories/index.tsx): `getContext`,
`loadFromURL`, `loadFromData`, `loadFromID`, `downloadCurrentStory`, and `MVSData`, including their existing arguments,
options, return behavior, and exposed context API. Retain `mvs-stories-viewer` and `mvs-stories-snapshot-markdown`,
their attributes (`context-name`, viewer `name`, and markdown `viewer-name`), and automatic registration when the
classic script loads.

Existing HTML using `https://cdn.jsdelivr.net/npm/molstar@latest/build/mvs-stories/mvs-stories.js` and the corresponding
CSS must work unchanged when `latest` advances to v6. Preserve equivalent package paths on other CDNs. Do not require
`type="module"`, scoped-package imports, new initialization calls, or edits to generated HTML. Internal plugin
composition and dependency packaging may change behind these app APIs. This compatibility promise is separate from
publishing a new stories library.

## 10. Implementation order

### 10.1 Release gate

The original rollout gate before default-branch integration was:

1. Land the planned 5.x work: particles, MVS changes, and bond-order perception.
2. Release the final 5.x feature minor.
3. Cut `v5` for any continuing fixes; freeze new features on that line.
4. Rename `master` to `main`, updating CI branch filters, docs, and clone instructions in the same change.

The workspace and composition changes have already landed on `master`; 5.13.0 and 5.13.1 are recorded in the release
history and `v5` exists for continuing fixes. The branch rename and matching CI/docs changes remain release-coordination
work; current workflows still target `master`. Keep later 5.x fixes on `v5` and forward-port them using the migration
records. Publish `dev` prereleases, then release stable after the acceptance checks.

### 10.2 Technical phases

Keep the major workstreams in separate PRs. Each phase ends with a working build; validate compiled packages and
source-based app bundles throughout.

Coordinate the [rendering-backend workstream](webgpu.md#6-minimal-implementation-sequence) with dependency cleanup and
packaging, and validate its contracts before freezing the v6 graphics API.

Keep the [fast-types workstream](fasttypes.md#5-effort-and-adoption) deferred for possible v7 adoption. For v6,
dual-registry release checks belong with packaging and CI, with slow types allowed on JSR and no associated annotation
or generator migration phase.

| Phase | Work and exit condition                                                                                                                                                                                                                                                                                             |
| ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0     | Verify package names; introduce pnpm, the workspace skeleton, and matching CI                                                                                                                                                                                                                                       |
| 1     | Audit and relocate reverse dependencies in §3.2; break plugin cycles and add explicit type imports. Record final API locations and verify the proposed package graph, including declaration edges. Retain existing module settings and override `verbatimModuleSyntax: false` for the temporary CJS build if needed |
| 2     | Split transform/provider modules and default-spec entry points; remove convenience barrels and the transform facade, split mixed implementation modules, migrate consumers to defining-module imports, and remove lazy getters; rename tests                                                                        |
| 3     | Introduce empty reference-counted registries, `spec.registry` entries, explicit base specs, preset lookup by id or alias, and presets that import what they run; prove the slim example and full default composition ([plugin-composition plan](../plans/plugin-composition.md) steps 2–4)                          |
| 4     | Drop CJS; set `type: module`, `NodeNext`, and `verbatimModuleSyntax`; use `.js` relative specifiers in source. Convert CommonJS globals/tooling and smoke-test emitted bins                                                                                                                                         |
| 5     | Move into grouped workspace packages; add exports, project references, direct dependencies, and per-app esbuild. Verify clean builds and packed consumers; stage the CDN-only `molstar` package                                                                                                                     |
| 6     | Finish standalone MVS builder/runtime, headless library, extension and server/CLI packaging; verify npm/JSR builder parity, headless dependency isolation, and command installation/migration instructions                                                                                                          |
| 7     | Ship the migrator, skills, mkdocs updates, advisory checks, and downstream smoke tests; publish `6-dev`, then stable when ready                                                                                                                                                                                     |

Maintain the migration map, docs, skills, and relevant checks as each phase lands; phase 7 closes remaining release
work. No legacy import shims are introduced at any phase.
