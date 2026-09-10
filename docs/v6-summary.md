# Mol\* 6.0: proposal summary

Mol* 6.0 moves to ESM packages grouped by layer, explicit plugin composition, and per-app builds. This is a proposal; see the [architecture and implementation plan](v6-architecture.md) for details.

## Packages

Use a pnpm workspace with one release version across public packages.

| Package | Responsibility |
| --- | --- |
| `@molstar/core` | General util, task, data, math, and state code |
| `@molstar/io` | Readers and writers |
| `@molstar/model` | Domain models, formats, properties, and queries |
| `@molstar/graphics` | GL, geometry, themes, representations, and canvas |
| `@molstar/plugin` | Plugin and plugin-state runtime, with explicit default-spec/catalog entry points |
| `@molstar/plugin-ui` | React UI and an explicit default UI-spec entry point |
| `@molstar/mvs-builder` / `@molstar/mvs` | Standalone MVS builder / Mol* runtime |
| `@molstar/<name>-extension` | Extensions and their dependencies |
| `@molstar/viewer` | Published Viewer API and app |
| Server and CLI packages | Existing tools plus `@molstar/migrate-6` |
| `molstar` | CDN viewer assets only |

The target dependency direction is `plugin → graphics → model → io → core` (“depends on”). The current folders do not satisfy it: IO helpers used by core, GPU math, and graphics-dependent model APIs must be relocated before packaging. Preserve familiar leaf paths where possible and record exceptions in the migration map.

Internal Mol* dependencies use exact release versions through `workspace:*`. Every package declares its direct npm imports; shared versions live in the pnpm catalog. Root is private tooling. React/React DOM are UI peers; headless/native and cloud dependencies remain optional where appropriate. Published declarations must declare the type dependencies consumers need.

## Plugin composition

Registries start empty. A `PluginFeature` contributes formats, transforms/actions, representations, themes, and behaviors without registering them on import.

```ts
import { PluginSpec } from '@molstar/plugin/spec';
import { Core } from '@molstar/plugin/features/core';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { BallAndStick } from '@molstar/graphics/repr/structure/representation/ball-and-stick';

const spec = PluginSpec.fromFeatures(Core, Sdf, BallAndStick, {
    canvas3d: { /* custom overrides */ },
});
```

The optional last argument supplies settings and extra actions/behaviors. Base entry points, including `createPluginUI`, take an explicit spec. Full defaults come from `DefaultPluginSpec` in `@molstar/plugin/default-spec` or `DefaultPluginUISpec` in `@molstar/plugin-ui/default-spec`. Viewer extensions remain app choices.

Defaults and catalogs stay in their owning packages behind explicit entry points. Lean roots and runtime leaves must not import or re-export them. Split transform modules and keep the full `StateTransforms` facade at `@molstar/plugin/state/transforms` for explicit consumer use; internal code imports individual transformers. Delete the lazy getters and enforce these module boundaries in CI. Keep built-in name types for completion without introducing upward package dependencies. Presets choose applicable registered representations and themes.

The SDF/ball-and-stick example must render while excluding unrelated parsers, cartoon/volume representations, and MP4 export. Verify the import graph, bundle, and rendering; empty registries alone are insufficient.

## Imports, ESM, and builds

- **Package imports:** `molstar/lib/mol-util/color` becomes `@molstar/core/util/color`. Use the same extensionless package subpaths inside the repository, including hops between layer folders in one package.
- **Relative source imports:** stay within a layer folder and use emitted `.js` paths. TypeScript and esbuild resolve them to source during builds; JS output retains them.
- **Library:** `tsc -b` with `NodeNext`, declarations, and ESM-only `lib/`. This output is generated locally but ships in the npm package. Publish `src/` for bundlers/debugging too. Retain the Node 22+ baseline unless runtime/tooling requires more.
- **Execution:** Node runs compiled JavaScript, including CLI/server bins. Native Node TypeScript execution is out of scope; `molstar-src` is for bundlers. Validated JSR packages expose source for Deno/compatible tooling. See the [execution contract](v6-architecture.md#51-supported-execution-modes).
- **Fast types:** defer consideration to v7. No v6 `isolatedDeclarations` requirement or broad API annotation migration; retain the [analysis](v6-fasttypes.md) for future planning.
- **npm and JSR:** keep compiled npm artifacts and try JSR source publication with `--allow-slow-types`, starting with the MVS builder. Coordinate matching versions for validated packages. Accept slower JSR consumer checking and potentially incomplete generated docs/types; native npm declarations still come from `tsc`. See the [release design](v6-fasttypes.md#7-publishing-to-npm-and-jsr-together).
- **Apps/examples:** each owns dependencies, entry, output, and scripts. A shared esbuild helper bundles source via `molstar-src`, with SCSS/assets, watch, and serve. Viewer is published; other apps and examples are private workspace packages.
- **Rendering backends:** establish a boundary for scenes, passes, GPU resources/operations, and readback within `@molstar/graphics`, retaining WebGL. WebGPU implementation and parity come later; see the [blast-radius analysis and minimal design](v6-webgpu.md).
- **Future offline rendering:** a [Blender extension](v6-webgpu.md#71-offline-rendering-with-blender) could consume portable scene snapshots through a separate asynchronous render-job interface, while WebGL/WebGPU provides interactive preview. Reuse geometry export/readback; Blender integration and effect translation are later work.

Keep compiler-supported TypeScript syntax, including namespaces, enums, and parameter properties. No `erasableSyntaxOnly` requirement or blanket syntax rewrite. Retain hot `const enum`s under existing compiler constraints; benchmark affected paths when necessary refactors change their use or emit.

Keep explicit type imports and enable `verbatimModuleSyntax` for ESM. If enabled before CJS removal, override it to `false` in the temporary CJS config. The ESM phase also converts `require` and path globals in CLI/server code and root scripts, with compiled-JS smoke tests.

## MolViewSpec

`@molstar/mvs-builder` owns the schema, builder, MVSJ/MVSX serialization/validation, and validation/schema CLIs. It has no Mol* package dependency and replaces molviewspec-ts / JSR `@molstar/molviewspec` after parity checks, publishing to npm and JSR.

`@molstar/mvs` depends on the builder and owns loading, plugin integration, annotations, and `mvs-render`. Python remains in mol-view-spec. A published stories library and reconciliation with MolViewStories need a separate plan.

## Migration and maintenance

Ship **`@molstar/migrate-6`**, with dry-run output and a manual-work report. It rewrites imports and dependencies, handles relocated APIs, and flags CommonJS, implicit default specs, and full-catalog imports. Respect downstream compiler conventions when changing relative extensions. Validate the tool on `pdbe-molstar` and `rcsb-molstar` before stable release.

**No compatibility import shims.** `molstar@6` retains the CDN viewer paths, not `lib/mol-*` or CJS. Library consumers migrate to scoped packages. Keep transformer identifiers and snapshot JSON; restoring a snapshot requires its features to be loaded.

Rename tests to `_test/**/*.test.ts`. Add `.agents/` maintainer skills for extensions, formats, representations, apps/examples, servers, and dependency updates, referenced by root `AGENTS.md`. Rewrite mkdocs for packages, composition, builds, migration, and adding code.

CI covers pnpm, typechecking, JSR publication checks with slow types allowed, package/value cycles, tests, source-based app builds, packed consumers, compiled CLI/server smoke tests, docs, and dependency advisories. Fail high/critical production advisories with explicit exceptions where justified.

## Release order

First land the planned 5.x work (particles, MVS changes, bond-order perception), release the final feature minor, cut `v5` for continuing fixes, and rename `master` to `main` with matching CI/docs changes.

Then use separate, buildable phases:

1. pnpm/workspace preparation.
2. Dependency audit, relocations, explicit type imports, and cycle removal.
3. Split transform/catalog modules, introduce features and empty registries, and validate default/slim compositions.
4. ESM-only runtime and `.js` relative imports, including scripts and bins.
5. Physical package moves, exports, direct dependencies, per-app builds, and packed-consumer checks.
6. Finish MVS, extension/server/CLI packaging, migration tool, skills, mkdocs, CI, and downstream validation.

Maintain docs and the migration map throughout. Publish `6.0.0-dev.N` under the `dev` tag on the way to stable `6.0.0`. The [detailed phases](v6-architecture.md#102-technical-phases) define the implementation checkpoints.
