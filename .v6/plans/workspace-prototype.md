# Mol* v6 workspace prototype implementation plan

Status: implemented and locally verified structural prototype (2026-10-01).
Hosted Linux/Xvfb unit tests and native headless capture have passed.
Consumer smoke checks run locally for now; CI retains unit tests, builds, and package checks.

See [workspace-usage.md](workspace-usage.md) for commands and consumer examples,
and [migration-map.json](migration-map.json) for source ownership changes.

This plan covers a packaging-first prototype of the repository architecture in
[architecture.md](../designs/architecture.md) and [summary.md](../designs/summary.md).
It records the planning decisions: distribution packages live under
`distributions/`, plugin packages are grouped by family, all public packages use
lockstep versioning, and `smoke/` verifies consumer-facing ESM and packaging.

## 1. Outcome and scope

Produce a working pnpm workspace with physically separated packages, explicit
dependencies, compiled ESM library exports, source-based app builds, and a local
CDN distribution package. Preserve current default plugin behavior, rendering,
transformer identifiers, and snapshot data.

The prototype establishes repository and package layout. It does not redesign the
application UI or implement the full v6 plugin-composition proposal.

In scope:

- Package ownership, physical source moves, and an import/API migration map.
- Dependency cleanup required to build an acyclic package graph, including types
  exposed in declarations.
- pnpm workspace/catalog setup and direct dependency declarations.
- ESM package exports, TypeScript project references, and per-package assets.
- Shared app tooling with app-owned configuration, dependencies, and scripts.
- Separate headless, MVS, extension, CLI, and server ownership.
- Colocated unit tests named `_test/<name>.test.ts`, with matching runner and
  package exclusions.
- Local packaging and consumer checks; shared release-version synchronization.
- Direct ESM consumption without consumer compilation/bundling, including a
  browser-ready module distribution and root `smoke/` fixtures.

Deferred:

- Move the workspace to TypeScript 7 and replace ESLint with Biome as the next
  tooling step. See [section 10](#10-next-step-typescript-7-and-biome).
- Rendering-backend extraction, GL resource/pass/readback redesign, WebGPU, and
  Blender integration. See [webgpu.md](../designs/webgpu.md).
- `PluginFeature`, empty registries, explicit base specs, registry-aware presets,
  slim-plugin bundle guarantees, and comprehensive transformer/catalog splitting.
- Comprehensive convenience-barrel removal and `StateTransforms` facade removal.
  Remove or split existing modules when needed for package boundaries; track the
  remaining work. Do not introduce new convenience barrels or compatibility shims.
- A broad test-runner migration and the full maintainer skills/documentation
  rewrite. Adapt existing checks where necessary for ESM.
- JSR publication, migration CLI implementation, downstream migration validation,
  automated publishing, and stable-release readiness.
- Fast types, `isolatedDeclarations`, blanket API annotations, and syntax rewrites.
  See [fasttypes.md](../designs/fasttypes.md).

This deliberately stages packaging before composition, unlike the full release
sequence in the architecture design. Completion means a working structural
prototype, not completion of every v6 release requirement.

## 2. Target repository layout

```text
packages/
  core/
  io/
  model/
  graphics/
  plugin/
    core/                   # @molstar/plugin
    ui/                     # @molstar/plugin-ui
    headless/               # @molstar/plugin-headless
  mvs/
    builder/
    runtime/
extensions/<name>/
apps/
  viewer/
  docking-viewer/
  mesoscale-explorer/
  mvs-stories/
examples/<name>/
servers/<name>/
cli/<name>/
distributions/
  molstar/
    package.json
    build/                  # generated, ignored by Git
      viewer/
      mvs-stories/
      esm/                  # browser-ready ESM modules, shared chunks, CSS
smoke/                      # private consumer fixtures and validation tooling
scripts/
.agents/                    # expanded maintainer skills are follow-up work
version.json                # canonical public release version
pnpm-workspace.yaml
tsconfig.base.json
tsconfig.json               # solution references
package.json                # private repository tooling
```

Each TypeScript library owns `src/`, generated `lib/`, `package.json`, and a
composite `tsconfig.json`. Apps additionally own HTML, themes, assets, and build
configuration. Generated build output is disposable; manifests and build/staging
scripts are checked in. Reserve `dist/`, if used, for generated output.

`packages/plugin/` and `packages/mvs/` are grouping directories without their own
package manifests. Plugin core, UI, and headless support remain separate packages
at `packages/plugin/{core,ui,headless}/`, with the npm names shown above. Plugin
core is the plugin runtime; `packages/core/` owns general utilities and primitives.
Grouping does not change dependency direction: UI and headless depend on plugin
core, and plugin core depends on neither.

Workspace globs:

```yaml
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

`distributions/molstar/` publishes as `molstar`. Its repository location does not
change installed `molstar/build/viewer/` or `molstar/build/mvs-stories/` paths.
Additional distributions get sibling directories, manifests, and staging scripts.
They assemble app outputs rather than maintain duplicate application source.

### Direct ESM consumption

Mol* builds release artifacts once; consumers must be able to execute the published
JavaScript without their own compilation or bundling step. Raw TypeScript/TSX
execution is not the contract.

- Installed Node consumers import compiled library subpaths using normal package
  exports, without TypeScript loaders, source conditions, or path aliases.
- Browser consumers load ready-to-serve JavaScript modules and precompiled CSS
  from the distribution using native `script type="module"`. Provide documented
  URL entry points and an optional generated import map for supported library
  subpaths. npm export maps alone do not resolve browser imports.
- Generate browser modules from an explicit supported entry inventory. Resolve
  third-party dependencies into browser-compatible modules and reuse shared chunks
  so different entry points share Mol* state/transformer identities. Define a
  single React instance policy for UI entry points, including any external peers.
- Generate exact import-map entries from the export/output inventory; prefix maps
  alone cannot add `.js` to extensionless package subpaths. Do not rely on a remote
  CDN service to transform npm packages at request time for the acceptance checks.
- Preserve the classic-script distributions alongside this additive ESM output.
  Keep browser-only and Node/headless entry points explicit; ESM syntax does not
  make a module portable across environments.

The initial browser inventory must cover Viewer and a small application assembled
from supported library entry points. A Viewer-only bundle does not prove modular
library consumption. Inventory expansion can follow package migration; avoid
independently bundling each package into copies of its shared dependencies.

## 3. Ownership and dependency rules

| Package | Initial source ownership |
| --- | --- |
| `@molstar/core` | `mol-util`, `mol-task`, `mol-data`, CPU/general `mol-math`, `mol-state` |
| `@molstar/io` | `mol-io`, after separating model-aware writer integration |
| `@molstar/model` | `mol-model`, `mol-model-formats`, `mol-model-props`, `mol-script`, excluding rendering/plugin integration |
| `@molstar/graphics` | `mol-gl`, `mol-geo`, `mol-theme`, `mol-repr`, `mol-canvas3d`, and relocated rendering/GPU helpers |
| `@molstar/plugin` | `mol-plugin` and `mol-plugin-state`, excluding reusable Node headless support |
| `@molstar/plugin-ui` | `mol-plugin-ui` |
| `@molstar/plugin-headless` | Headless context, capture helpers, and Node output handling |
| `@molstar/mvs-builder` | Standalone schema, builder, serialization, validation, and validation/schema commands |
| `@molstar/mvs` | MVS loading, annotations, and plugin integration |
| `@molstar/<name>-extension` | Individual extensions and their direct dependencies |
| `@molstar/viewer` | Published Viewer API and application |
| `@molstar/<name>-cli` / `@molstar/<name>-server` | Supported commands and servers |
| `molstar` | Assembled Viewer and MVS Stories CDN assets |

The layer direction is `plugin -> graphics -> model -> io -> core`, where an arrow
means "depends on". Direct dependencies on any lower layer are allowed. UI,
headless support, extensions, and applications depend on the layers they use;
those layers must not depend back on their consumers.

Cross-layer-folder imports use public extensionless package subpaths, including
self-imports inside a package. Relative imports stay within a layer folder and
use emitted `.js` specifiers. Package exports must handle TSX, cohesive directory
entry points, and CSS/assets explicitly. Source mappings must describe the same
public API as compiled exports; they cannot hide dependencies on unowned source.

Every direct npm importer declares its dependencies. Shared external versions
live in the pnpm catalog. React/React DOM are UI peers, with an explicit policy for
consumer-facing React types. Declare `tslib` where emitted helpers require it.
Headless injected native modules remain optional peers. Native `gl` and `canvas`
are opt-in peers in the rendering CLI and Node examples too; the CLI declares its
pure JavaScript codecs directly. Native test setup lives in ignored `.cache/native`
and is invoked explicitly in CI. Browser plugin/MVS imports must
not acquire headless or MP4 dependencies.

### Initial boundary issues to resolve

The planning scan examined static imports/re-exports in the current source tree,
including type edges and excluding tests. It is a starting inventory, not a full
runtime, declaration, or cycle audit.

| Crossing | Examples | Resolution direction |
| --- | --- | --- |
| Core -> IO | `mol-data/db/column.ts`, `mol-util/data-source.ts` | Move general numeric/string/byte primitives into core; resolve tokenizer type ownership |
| Core -> graphics | GPU Gaussian density and related math APIs | Move GPU helpers into graphics; keep their WebGL implementation |
| Core -> model | `mol-util/param-definition.ts` uses query/script types | Extract a suitable lower-level contract or move the domain integration |
| IO -> model | CCP4 writer, ligand encoder, molecule writer integration | Separate low-level encoding contracts from model-aware adapters |
| Model -> graphics | Shape, particle transforms, locations, property visuals/themes | Separate domain data/computation from graphics integration |
| Model -> plugin | Custom element property label-provider integration | Move plugin integration above model or extract a lower-level contract |
| Graphics -> plugin | External-structure and related themes | Separate graphics provider logic from plugin state/query integration |

Converting an import to `import type` does not resolve a package-level declaration
cycle. Record exact destination APIs and imports in the migration map before
rewriting their consumers. Avoid moving entire cohesive APIs merely to silence a
graph check when a small contract split would preserve useful ownership.

## 4. Versioning

Use one version across every public Mol* library, extension, Viewer, CLI, server,
and Mol* distribution package. Start the prototype at `6.0.0-dev.0`; subsequent
development releases use `6.0.0-dev.N`. Unchanged public packages advance with the
release. Independent package versioning is deferred.

- `version.json` is the canonical release version, independent of the private root
  manifest. Each public manifest contains the synchronized concrete version.
- Internal dependencies use `workspace:*` in source manifests. Validate exact
  matching release versions in the packed manifests.
- A local synchronization command updates public manifests and generates the
  runtime version shared by library and app builds. Generate the build timestamp
  once per build invocation; version checks must not create timestamp-only churn.
- A read-only check rejects inconsistent public versions, internal ranges, or
  missing public-package inventory entries. Include peer/optional internal edges.
- Private apps/examples/tooling are excluded from public version synchronization;
  their dependencies still use workspace references.
- Maintain one release changelog identifying affected packages. Public publishing
  and cross-registry release coordination remain separate follow-up work.

Lockstep versions and exact internal dependencies reduce accidental mixing; they
do not prevent consumers from explicitly installing conflicting release versions.

## 5. Implementation checkpoints

Use separate reviewable changes with working builds at each completed checkpoint.
Do boundary cleanup in the existing tree before physical moves where practical.
Migration code may temporarily remain under `src/`, but each moved package must
have complete ownership and must not reach back into that legacy source tree.
Update remaining callers in the same change as each move; do not add old-path
re-export shims. Couple moves that cannot otherwise produce a buildable checkpoint.

### A. Baseline, inventory, and workspace foundation

- [ ] Capture a representative pre-migration Viewer render as a comparison baseline.
  Baseline typecheck, lint, and unit tests were captured; a pre-migration render
  was not captured, so the current browser fixtures provide functional evidence
  rather than a visual before/after comparison.
- [x] Create a reproducible source/package graph inventory covering imports,
  re-exports, dynamic loads, declaration dependencies, external imports, and assets.
- [x] Create the source/API migration map and explicit public/private package list.
- [x] Introduce pnpm, a pinned package-manager version, catalog, and lockfile; adapt
  current CI/install commands without requiring public publishing.
- [x] Add shared compiler configuration and local version sync/check tooling.
- [x] Establish the private `smoke/` harness and baseline fixture data; distinguish
  planned checks from checks that run against functional v6 packages.

Exit: a clean pnpm install reproduces current behavior, the audit is reproducible,
and the first independent package slice has an agreed ownership map. Do not count
empty package directories as a functioning prototype.

### B. First functioning slice: core and IO

- [x] Resolve core reverse dependencies, including general parsing primitives,
  query parameter types, and extraction of GPU-facing math.
- [x] Resolve model-aware writer dependencies so IO can compile below model.
- [x] Move core and IO into their packages, rewrite all affected callers, and add
  direct dependencies, exports, assets, and project references.
- [x] Establish compiled ESM and declaration consumption plus source-condition
  bundling for these packages. Convert relative imports and affected barrels as
  required to make the slice independently usable.
- [x] Keep not-yet-moved code buildable using the packages' supported imports.

Exit: core and IO build without remaining `src/` dependencies, the package graph
has no upward edges, and a clean packed consumer exercises utility/task and parsing
APIs through supported subpaths. Existing affected behavior checks pass.
Run the Node ESM, declaration, and packed-consumer smoke fixtures for this slice.

### C. Model, graphics, plugin, UI, and headless packages

- [x] Resolve remaining reverse dependencies and move packages in dependency order.
- [x] Preserve WebGL algorithms/shaders, Canvas3D behavior, representations, and
  current default registration. Limit graphics changes to ownership/import fixes.
- [x] Separate Node headless support and explicit MP4 integration from browser
  runtime ownership; preserve injected native-module support.
- [x] Keep default specs/catalogs owned by their packages. Track composition/barrel
  debt without promising slim bundles in this prototype.
- [x] Finish scoped imports, TSX exports, UI style exports/assets, and type dependencies.
- [x] Build moved packages with `NodeNext`, `verbatimModuleSyntax`, declarations,
  declaration maps, and the existing strictness settings; retain Node 22+.

During migration, retain existing compiler settings for unmoved code. If shared
`verbatimModuleSyntax` reaches the temporary CJS build, override it there to false.
Remove the CJS build when all its callers/bins have migrated. Preserve supported
namespaces, enums, parameter properties, and current inferred API precision.

Exit: the library graph and emitted declarations are acyclic across packages;
compiled plugin/UI consumers work; representative structure loading, rendering,
selection, and snapshot restoration retain baseline behavior. Browser imports
exclude headless/native/video modules.

### D. Per-app builds and distribution staging

- [x] Extract `scripts/esbuild/app.mjs` with source-condition resolution, SCSS,
  asset handling, version injection, production builds, watch, and serve support.
- [x] Move Viewer first; give it its own dependencies, entry points, themes, output
  configuration, and scripts. Preserve the published Viewer API.
- [x] Move the remaining apps/examples and remove dependence on the root `Apps`
  list and root dependency hoisting.
- [x] Create `distributions/molstar/` with an explicit staging/packing command that
  assembles freshly built Viewer and MVS Stories assets.
- [x] Preserve classic-script globals, output filenames, custom elements, HTML,
  CSS, themes, images, icons, and both installed CDN directory paths.
- [x] Generate browser-ready ESM modules/shared chunks, precompiled CSS, and the
  import map from the supported entry inventory; add no-build Viewer and modular
  library pages under `smoke/browser/`.

Exit: Viewer dev/build runs from its workspace package, app bundles build using
source exports after deleting library output, and a packed `molstar` distribution
passes existing Viewer/MVS Stories API and asset smoke checks. The repository root
is not the published CDN package.
Both browser smoke pages execute from packed artifacts on a static HTTP server
without a bundler, dev-server transformation, or consumer compilation.

### E. MVS, extensions, commands, servers, and complete workspace checks

- [x] Separate MVS builder/runtime ownership. Keep the builder independent of all
  Mol* packages and the runtime dependent on the builder rather than a copied API.
- [x] Move individual extensions and declare their direct imports/dependencies.
- [x] Package supported CLI/server entry points; assign shared server helpers and
  development-only generators explicit ownership rather than cross-package relatives.
- [x] Preserve executable names. Keep validation/schema commands with MVS builder;
  isolate rendering/native/video dependencies in the rendering CLI.
- [x] Convert affected `require`, path globals, scripts, and asset lookups to working
  compiled ESM execution. Never use typechecking alone as a bin-runtime check.
- [x] Finish root solution references and private tooling setup, remove obsolete
  monolithic build/CJS machinery, and verify all public package inventories.
- [x] Add workspace CI checks and short build/migration instructions. Record remaining
  v6 composition, publishing, migrator, and documentation work explicitly.

Exit: clean install/build/pack checks cover the complete workspace; source-built
apps and compiled consumers both work; relevant unit/browser/headless checks pass;
supported bins execute from installed tarballs with their required assets.

## 6. Prototype acceptance

- [x] Root standard `smoke/` checks run locally against freshly built/packed output;
  a successful import alone is not the complete runtime assertion. Hosted checks
  passed once; consumer smoke checks now remain local to avoid Chromium setup in CI.
- [x] No cross-package relative imports or dependencies on unowned legacy source.
- [x] No package cycles, including dependencies needed by emitted declarations.
- [x] Every public package declares its direct runtime and consumer-facing type
  dependencies; every export resolves in its packed artifact.
- [x] Clean `tsc -b` builds produce ESM JS/declarations without stale output.
- [x] Source-condition app builds work without prebuilt library `lib/` directories.
- [x] Plain Node `.mjs` consumers execute installed package exports with no loader.
- [x] Plain browser HTML/JavaScript consumers exercise Viewer and modular library
  ESM without consumer builds, with successful module/assets requests and no
  uncaught runtime errors or duplicated shared module instances.
- [x] Full default Viewer behavior survives representative loading/rendering,
  selection, and snapshot checks; no backend redesign is needed.
- [x] Native headless capture passed on hosted Linux/Xvfb. The capture fixture
  remains available locally; this macOS host returns no GL context. MVS validation passes without native rendering dependencies, and
  the base headless package has no automatic MP4 integration.
- [x] Packed distribution preserves Viewer/MVS Stories browser APIs and paths.
- [x] Public versions and packed internal dependency versions match `version.json`.
- [x] Required HTML, CSS, shader modules, images, and command data are present;
  tests/fixtures and unintended development files are excluded from tarballs.
- [x] Deferred architecture requirements remain documented as follow-up work.

Use existing meaningful tests and targeted consumer/build/browser fixtures. Do not
expand this into a separate test-framework project or add tests that merely mirror
directory names.

## 7. Smoke-test organization

Use root `smoke/` for small integration fixtures that act like external consumers,
not for moving the existing unit/browser test suites. See
[smoke/README.md](../../smoke/README.md) for the fixture layout and commands.
The harness is a private workspace package; temporary installed-consumer projects
are created outside the repository and must not inherit workspace dependencies.

| Fixture | Executes against | Required evidence |
| --- | --- | --- |
| Node ESM | Installed compiled exports | Task/utility execution, parsing/model results, supported subpath resolution, shared API identity |
| Types | Packed declarations in a clean TS consumer | Normal consumer typecheck without source aliases or ambient root dependencies |
| Browser ESM | Statically served packed browser distribution | Viewer and modular plugin creation, deterministic local structure loading/rendering, selection/state behavior, CSS/assets, shared instances |
| Source bundler | Package source conditions | A small bundle builds without library output and resolves TSX/assets correctly |
| CLI/headless | Installed command/headless packages | Validation/schema command output and capture with declared/injected dependencies; native checks in a suitable dedicated CI job |
| Classic browser compatibility | Statically served packed classic outputs | Existing Viewer/MVS Stories globals, APIs, custom elements, and assets |

Commands are `pnpm --dir smoke smoke` for the standard suite and `pnpm --dir smoke smoke:node`,
`pnpm --dir smoke smoke:types`, `pnpm --dir smoke smoke:browser`, `pnpm --dir smoke smoke:source`, and
`pnpm --dir smoke smoke:headless` for focused checks. `pnpm --dir smoke smoke` prepares fresh package/distribution artifacts,
installs the packed dependency closure, then runs the standard fixtures; focused
commands validate their prerequisites rather than silently accepting stale output.

Install tarballs into isolated temporary consumers, with every relevant internal
dependency resolved from the same packed release rather than the registry or
workspace. Run Node/types checks there so missing direct/type dependencies cannot
be satisfied by the root install. Serve browser fixtures and unpacked distribution
assets using a plain static server, with local pinned dependencies and fixture data.

Use meaningful result assertions, bounded timeouts, and nonzero failure exits. For
browser runs, capture runtime errors and failed module/asset requests, wait for a
real loaded scene, and verify rendering/state behavior. Save diagnostic logs/images
on failure, clean up servers/temp consumers, and exclude generated smoke artifacts
from Git and public tarballs. Keep package-level unit tests where they belong.

## 8. Implementation evidence and release separation

The workspace inventory contains 69 packages, including 43 public packages.
Production source has moved out of the original monolithic `src/` tree. The
structural work in checkpoints B–E is implemented; verification evidence and
remaining host-dependent checks are recorded below.

- Full lint and TypeScript solution builds pass.
- Jest: 131 suites and 1,477 tests pass; 11 suites and 14 tests skip.
- The source-condition app build succeeds with all 67 library output directories
  absent.
- Isolated packed Node ESM, strict declaration, source-bundler, and MVS CLI
  consumer checks pass.
- Workspace import/declaration graph and version checks pass for all 69 packages.
- Native unit tests and packed headless PNG capture passed on hosted Linux/Xvfb
  in run 36869833592. This macOS host cannot create a GL context. CI now retains
  native unit tests, while consumer smoke checks run locally.
- All 43 public tarballs pass version, dependency-range, bin, asset, and export
  target checks, including source/Sass conditions. Build metadata and source
  test fixtures are excluded.
- Compiled cif2bcif/cifschema and model/volume server command usage checks pass.
  Native rendering CLI startup is host-dependent and remains unverified locally.
- Packed browser fixtures pass for native ESM Viewer, modular plugin UI, classic
  Viewer, and classic MVS Stories. They exercise local structure loading and
  representations, selection/snapshot round trips, shared module identities,
  custom elements, and successful assets without runtime errors.
- No public publishing, deployment, or release branch changes were performed.

The full architecture's final 5.x release, `v5` maintenance branch, default-branch
rename, registry-name checks, and publishing gates still apply before v6 lands or
ships as a release. They are not actions performed by this prototype plan. Treat
default-branch integration and public publishing as separate release work.

The `distributions/molstar/` location supersedes the top-level `molstar/` location,
and `packages/plugin/{core,ui,headless}/` supersedes the three sibling plugin
directories in the design documents for this plan. When implementing, update
related design references and build examples to the chosen locations. Use the
actual filenames in `.v6/designs/` when updating links; several existing links
still use `v6-*` names.

## 9. Implementation orchestration

Requested profile: Sol with medium reasoning as the primary orchestrator, using
GPT-6 Luna subagents for bounded implementation tasks. The primary model/effort is
selected in the app; subagent model selection is explicit when delegating work.
This records the requested setup, not verification of the current primary model.

The orchestrator owns package-boundary decisions, API contract changes, migration
ordering, shared configuration, cross-package integration, and checkpoint acceptance.
Delegate substantial independent work to Luna once its inputs and acceptance
criteria are clear, including:

- Read-only dependency, external-import, export, and asset inventories.
- Package manifests/configuration within an agreed ownership/export scheme.
- Mechanical moves and import rewrites for explicitly assigned file sets.
- Smoke fixtures, version checks, asset-copy/staging helpers, and documentation.

Each task receives the relevant plan/contract, exact file ownership, dependencies,
required checks, and a request to report changed files and verification results.
Give workers focused context rather than the full conversation where possible.
Do not let concurrent workers rewrite overlapping files or independently invent
package boundaries. Stage broad source moves before parallel consumer rewrites.

Use up to three workers concurrently with the orchestrator under this session's
four-agent limit, only where tasks are independent. The orchestrator reviews all
changes, resolves shared-file edits, runs integration checks, and updates checkpoint
status. Reassign ambiguous work to the orchestrator rather than repeatedly issuing
mechanical fixes without understanding the dependency problem.

## 10. Next step: TypeScript 7 and Biome

After this structural prototype, migrate the workspace compiler to TypeScript 7
and replace ESLint with Biome. These changes are planned follow-up work, not
implemented in the current PR.

- [ ] Replace TypeScript 6.0.3 and JavaScript `tsc` with the TypeScript 7 native
  compiler. Update the shared catalog, package scripts, project-reference builds,
  declaration checks, and CI commands. Preserve strict type checking, ESM/source
  conditions, declaration output, package exports, and incremental rebuilds.
- [ ] Replace ESLint and its TypeScript parser/plugins with Biome. Map the current
  lint rules and exclusions, document unsupported rules and their replacements,
  and update `pnpm lint` and CI. Configure formatting to preserve repository
  conventions and keep any broad formatting changes separate from the tooling
  migration. Remove the superseded ESLint configuration and dependencies.
- [ ] Verify a clean install, lint, unit tests, library/app/distribution builds,
  workspace/version checks, all public tarballs, and local consumer smoke checks.
  Compare fresh and incremental build times against the current baseline.

Baseline: a clean build of 67 TypeScript projects took about 2m 24s locally; the
incremental rerun took 0.52s. Hosted CI spent about 5m 41s on the library build,
8s on apps/examples, and 1s on distribution (run 36869833592). Record the new
measurements after migration rather than assuming a particular speedup.
