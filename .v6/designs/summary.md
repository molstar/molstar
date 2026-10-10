# Mol\* 6.0: proposal summary

Mol* 6.0 moves to ESM packages grouped by layer, explicit plugin composition, and per-app builds. The workspace/ESM
split, TypeScript 7/Biome tooling, formatting, and plugin composition are implemented. This design also includes work
still planned; the [checklist](../plans/checklist.md) records current status and the [architecture](architecture.md)
describes the contracts. Alex owns the rendering-backend workstream.

## Packages

Use a pnpm workspace with one release version across public packages.

| Package                                 | Responsibility                                                                   |
| --------------------------------------- | -------------------------------------------------------------------------------- |
| `@molstar/core`                         | General util, task, data, math, and state code                                   |
| `@molstar/io`                           | Readers and writers                                                              |
| `@molstar/model`                        | Domain models, formats, properties, and queries                                  |
| `@molstar/graphics`                     | GL, geometry, themes, representations, and canvas                                |
| `@molstar/plugin`                       | Plugin and plugin-state runtime, with explicit default-spec/catalog entry points |
| `@molstar/plugin-ui`                    | React UI and an explicit default UI-spec entry point                             |
| `@molstar/plugin-headless`              | Reusable Node headless context, screenshots, and output handling                 |
| `@molstar/mvs-builder` / `@molstar/mvs` | Standalone MVS builder / Mol* runtime                                            |
| `@molstar/<name>-extension`             | Extensions and their dependencies                                                |
| `@molstar/viewer`                       | Published Viewer API and app                                                     |
| `@molstar/<name>-cli`                   | Command-focused packages, including `mvs-render-cli` and `migrate-6-cli`         |
| `@molstar/<name>-server`                | Server packages and their related commands                                       |
| `molstar`                               | Viewer and MVS Stories CDN apps, retaining paths and browser APIs                |

The implemented dependency direction is `plugin → graphics → model → io → core` (“depends on”). Packaging relocated
reverse dependencies, including IO helpers used by core, GPU math, and graphics-dependent model APIs. Familiar leaf
paths are preserved where possible and exceptions are recorded in the migration map.

Internal Mol* dependencies use exact release versions through `workspace:*`. Every package declares its direct npm
imports; shared versions live in the pnpm catalog. Root is private tooling. React/React DOM are UI peers. The headless
library retains module injection and optional native peers; the rendering CLI directly declares its native
modules/codecs. Plugin and browser MVS packages do not depend on headless support. Published declarations must declare
the type dependencies consumers need.

Headless support depends on plugin and graphics, with MP4 integration in an explicit extension module. Command-focused
packages use `-cli`: `@molstar/mvs-render-cli`, `@molstar/cif2bcif-cli`, `@molstar/cifschema-cli`, and
`@molstar/migrate-6-cli`. Existing command names stay unchanged; the new migration command is `molstar-migrate-6`.
Libraries with ancillary bins and server packages retain their domain names. See the
[command/package map](architecture.md#71-cli-packages-and-executable-names).

## Plugin composition

See the [plugin-composition design](plugin-composition.md) and its
[implementation plan](../plans/plugin-composition.md).

Registries start empty. `PluginSpec.registry` lists declarative entries (formats, representations, themes, presets,
selection queries, and similar providers) that `PluginContext` registers in order, without dependency resolution.

```ts
const spec: PluginSpec = {
  registry: [Sdf, DefaultHierarchyPreset, BallAndStickPreset],
  behaviors: [/* ... */],
  config: [[PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-ball-and-stick']],
};
// ... create the plugin, download and parse an SDF file
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
```

What is included is what the spec imports and lists; single-file builds tree-shake the same way. Behaviors are unchanged
and can reuse the entry format through `plugin.register(entry)`. Presets import the providers they run and are applied
by id or alias; choices such as a format's default preset go through config. Unregistered representation and theme names
fall back to the registry default, with a warning from the builders, param helpers, and snapshot restore. PyMOL, VMD,
and Jmol scripts work when the app imports their transpiler modules. Base entry points, including `createPluginUI`, take
an explicit spec. Full defaults come from `DefaultPluginSpec` in `@molstar/plugin/default-spec` or `DefaultPluginUISpec`
in `@molstar/plugin-ui/default-spec`. Viewer extensions remain app choices.

The `StateTransforms` facade is removed and transform modules are split by functionality, with format-specific
transformers next to their providers. Transformers stay globally registered on import; snapshot loading checks their ids
before changing state. `BuiltIn*` name types remain available through type imports from catalog modules. Module
boundaries are enforced in workspace checks.

The SDF/ball-and-stick slim example renders while excluding unrelated parsers, representations, presets, and MP4 export.
Import-graph, bundle, and rendering acceptance checks are recorded in the plugin-composition plan.

## Imports, ESM, and builds

- **Package imports:** `molstar/lib/mol-util/color` becomes `@molstar/core/util/color`. Use the same extensionless
  package subpaths inside the repository, including hops between layer folders in one package.
- **No barrel files:** expose defining modules directly through package subpaths and use them internally and in
  examples. Remove convenience re-export modules and aggregate facades; keep default specs/catalogs restricted to
  deliberate composition. Enforce the [policy](architecture.md#44-no-barrel-files) in CI.
- **Relative source imports:** stay within a layer folder and use emitted `.js` paths. TypeScript and esbuild resolve
  them to source during builds; JS output retains them.
- **Library:** `tsc -b` with `NodeNext`, declarations, and ESM-only `lib/`. This output is generated locally but ships
  in the npm package. Publish `src/` for bundlers/debugging too. Retain the Node 22+ baseline unless runtime/tooling
  requires more.
- **Execution:** Node runs compiled JavaScript, including CLI/server bins. Native Node TypeScript execution is out of
  scope; `molstar-src` is for bundlers. Validated JSR packages expose source for Deno/compatible tooling. See the
  [execution contract](architecture.md#51-supported-execution-modes).
- **Fast types:** defer consideration to v7. No v6 `isolatedDeclarations` requirement or broad API annotation migration;
  retain the [analysis](fasttypes.md) for future planning.
- **npm and JSR:** keep compiled npm artifacts and try JSR source publication with `--allow-slow-types`, starting with
  the MVS builder. Coordinate matching versions for validated packages. Accept slower JSR consumer checking and
  potentially incomplete generated docs/types; native npm declarations still come from `tsc`. See the
  [release design](fasttypes.md#7-publishing-to-npm-and-jsr-together).
- **Apps/examples:** each owns dependencies, entry, output, and scripts. A shared esbuild helper bundles source via
  `molstar-src`, with SCSS/assets, watch, and serve. Viewer is published; other apps and examples are private workspace
  packages.
- **Rendering backends:** establish a boundary for scenes, passes, GPU resources/operations, and readback within
  `@molstar/graphics`, retaining WebGL. WebGPU implementation and parity come later; see the
  [blast-radius analysis and minimal design](webgpu.md).
- **Future offline rendering:** a [Blender extension](webgpu.md#71-offline-rendering-with-blender) could consume
  portable scene snapshots through a separate asynchronous render-job interface, while WebGL/WebGPU provides interactive
  preview. Reuse geometry export/readback; Blender integration and effect translation are later work.

Keep compiler-supported TypeScript syntax, including namespaces, enums, and parameter properties. No
`erasableSyntaxOnly` requirement or blanket syntax rewrite. Retain hot `const enum`s under existing compiler
constraints; benchmark affected paths when necessary refactors change their use or emit.

Keep explicit type imports and enable `verbatimModuleSyntax` for ESM. If enabled before CJS removal, override it to
`false` in the temporary CJS config. The ESM phase also converts `require` and path globals in CLI/server code and root
scripts, with compiled-JS smoke tests.

## MolViewSpec

`@molstar/mvs-builder` owns the schema, builder, MVSJ/MVSX serialization/validation, and validation/schema CLIs. It
depends on the standalone `@molstar/query-language` package. `@molstar/mvs-builder/molql` exposes `MolScriptBuilder` and
`compileScript` for MolScript, PyMOL, VMD, and Jmol. Validation checks syntax and callable/argument names against symbol
tables without the molecular query runtime; custom vocabularies use an explicit lookup. See the
[MolQL design](architecture.md#72-molql-builder-and-validation-design). Replacing molviewspec-ts / JSR
`@molstar/molviewspec` still requires parity checks and coordinated npm/JSR publication.

`@molstar/mvs` depends on the builder and owns loading, plugin integration, and annotations. `mvs-render` ships
separately in `@molstar/mvs-render-cli`, composing MVS, headless support, and MP4 export. Validation/schema commands
stay with the builder without native render dependencies. Python remains in mol-view-spec. The existing MVS Stories app
still ships in the root `molstar` package. A new stories library and reconciliation with MolViewStories need a separate
plan.

## Migration and maintenance

Ship **`@molstar/migrate-6-cli`** with the `molstar-migrate-6` command, dry-run output, and a manual-work report. It
rewrites imports and dependencies, handles relocated APIs, and flags CommonJS, implicit default specs, and full-catalog
imports. Respect downstream compiler conventions when changing relative extensions. Validate the tool on `pdbe-molstar`
and `rcsb-molstar` before stable release.

**No library compatibility import shims.** `molstar@6` retains `build/viewer/` and `build/mvs-stories/`, including their
classic-script globals, APIs, CSS/assets, and custom elements. Existing MVS HTML viewers importing `molstar@latest` from
a CDN must work without edits; verify against the packed candidate before advancing `latest`. See the
[browser compatibility contract](architecture.md#93-compatibility-contract). Library consumers migrate from
`lib/mol-*`/CJS to scoped packages. Keep transformer identifiers and snapshot JSON; restoring a snapshot requires its
transformers to be imported, and unregistered representation or theme names fall back to the registry default with a
warning.

Rename tests to `_test/**/*.test.ts`. Add `.agents/` maintainer skills for extensions, formats, representations,
apps/examples, servers, and dependency updates, referenced by root `AGENTS.md`. Rewrite mkdocs for packages,
composition, builds, migration, and adding code.

CI covers pnpm, typechecking, JSR publication checks with slow types allowed, package/value cycles, tests, source-based
app builds, packed consumers, compiled CLI/server smoke tests, docs, and dependency advisories. Fail high/critical
production advisories with explicit exceptions where justified.

## Release order

The original rollout called for finishing 5.x feature work, releasing a final feature minor, cutting `v5`, and renaming
`master` to `main` before v6 integration. The workspace and composition changes have already landed on `master`; the 5.x
release history includes 5.13.0 and 5.13.1, and `v5` exists for continuing fixes. The branch rename and corresponding
CI/docs changes remain release-coordination work; current workflows still target `master`.

The technical phases were planned as separate, buildable changes:

1. pnpm/workspace preparation.
2. Dependency audit, relocations, explicit type imports, and cycle removal.
3. Split transform/catalog modules, introduce empty registries and `spec.registry`, and validate default/slim
   compositions.
4. ESM-only runtime and `.js` relative imports, including scripts and bins.
5. Physical package moves, exports, direct dependencies, per-app builds, and packed-consumer checks.
6. Finish MVS, extension/server/CLI packaging, migration tool, skills, mkdocs, CI, and downstream validation.

Workspace preparation, dependency relocation, plugin composition, ESM conversion, and physical package moves are
implemented. Remaining MVS validation/parity, JSR publication, migrator, maintenance, and release work is tracked in the
checklist. Maintain docs and the migration map throughout. Publish `6.0.0-dev.N` under the `dev` tag on the way to
stable `6.0.0`. The [detailed phases](architecture.md#102-technical-phases) define the implementation checkpoints.
