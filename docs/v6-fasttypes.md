# Mol\* 6.0: npm/JSR distribution and deferred fast types

Analysis against `molstar@5.11.0`, commit `bbebb6082`, using TypeScript 6.0.3 and Deno 2.9.6 (also TypeScript 6.0.3) on September 10, 2026. This complements the [architecture](v6-architecture.md) and [summary](v6-summary.md). These are proposed changes; source edits and publishing configuration used for the pilot exist only in temporary audit copies.

## 1. Recommendation

**For v6, try JSR source publication with `--allow-slow-types` and defer fast types to consideration for v7.** Do not add an `isolatedDeclarations` requirement, broad annotation migration, or schema/factory redesign to v6. Preserve current API precision and use `tsc` to generate native npm declarations.

The audit indicates that fast types are feasible, often without runtime or declaration changes, but the schema/factory work could create pressure for public API changes. Its cost and review burden are outside the v6 scope. Sections 2–5 retain the measurements, effort estimates, and agent-assisted approach as future research; they are not v6 requirements or a commitment to implement fast types in v7.

Keep one source tree and two distribution outputs:

- **npm:** compiled ESM JavaScript and declarations in `lib/`, plus `src/` for opt-in bundlers and debugging.
- **JSR:** TypeScript source and source entry points, allowing slow types. Start with the standalone MVS builder and extend coverage package by package after validation.
- **Release:** one version and commit, separate publication steps and validation for each registry. Registry availability does not imply support for every runtime.

`lib/` is generated and disposable in the checkout, but its runtime files remain in the npm tarball. JSR source publication does not need that directory. Fast types do not require erasable syntax, native Node TypeScript execution, or replacing the existing build system.

Allowing slow types can make JSR consumer checking slower and leave JSR-generated documentation or npm-compatibility declarations incomplete or imprecise. Accept those limitations for this v6 experiment; native npm consumers keep the full `tsc` declaration output. The flag does not waive normal typechecking, import resolution, or other publication rules. See [JSR's slow-types limitations](https://jsr.io/docs/about-slow-types).

## 2. What fast types require

JSR derives documentation, declarations, and fast-check information from exported source APIs without full TypeScript inference. Public function returns, variables, and class members therefore need explicit or simply inferable types. Implementation-local inference can remain. JSR checks the API reachable from package exports; TypeScript's `isolatedDeclarations` checks module exports more broadly and is sometimes stricter. JSR also restricts augmentation, CommonJS export forms, and expando properties. See [JSR's rules and comparison with isolated declarations](https://jsr.io/docs/about-slow-types).

If fast-type migration is adopted for v7, use two complementary checks:

1. Gradually enable `isolatedDeclarations` for source modules across the codebase, beginning with published libraries. It provides an editor/compiler check before packaging. It requires declaration/composite mode and does not itself replace `tsc` with a faster emitter. See [TypeScript's design](https://www.typescriptlang.org/docs/handbook/release-notes/typescript-5-5.html#isolated-declarations).
2. Run `deno publish --dry-run` on each actual JSR package with its final export graph. A future fast-type guarantee requires this to pass without `--allow-slow-types`; the compiler check alone is not proof that publication succeeds. V6 explicitly permits the flag instead.

“Across the codebase” means explicit module API contracts, not annotating every callback or local variable. Private apps/tests need not become JSR packages. For a uniform compiler policy, check their exported helpers in a separate declaration-enabled configuration as well; keep implementation inference where it does not escape.

This is independent of `erasableSyntaxOnly`. Deno supports runtime-producing TypeScript constructs, including namespaces, enums, and parameter properties; the [TypeScript runtime guide](https://docs.deno.com/runtime/fundamentals/typescript/) distinguishes execution from checking. Keep these constructs unless a specific declaration or publication diagnostic requires a change. Likewise, fast-check compliance does not guarantee that a consumer's complex generic instantiations become cheap.

## 3. Measured blast radius

### 3.1 Method and limits

The existing `tsc --noEmit` check passed. For the audit, the TypeScript compiler API loaded the existing `tsconfig.json`, set `isolatedDeclarations: true`, and collected `getDeclarationDiagnostics()` without writing output. It covered 1,592 configured roots, including two declaration files, and 1,590 implementation source files. The current `_spec/` exclusion remains; apps, examples, servers, and browser/performance tests already included by this config remain in scope.

The follow-up Deno audit used temporary copies and a synthetic source export map; see §3.4. The proposed workspace export maps do not exist yet. Counts below describe the current compiler graph, not the final package API. Several diagnostics can belong to one object literal, and correcting its annotation can remove many at once; resolving a factory can also expose further work.

Pin Deno and reconcile its checking configuration with the package's supported TypeScript settings and browser/Node type libraries during the pilot. Deno's defaults and bundled compiler can produce additional errors; those are separate from the measured isolated-declaration diagnostics. Include namespace merging, enum emit, and TSX in source-consumer validation.

Reproduce the diagnostic collection from the repository root after installing its locked dependencies:

```js
const ts = require('typescript');
const config = ts.readConfigFile('tsconfig.json', ts.sys.readFile);
const parsed = ts.parseJsonConfigFileContent(config.config, ts.sys, process.cwd());
const program = ts.createProgram(parsed.fileNames, {
    ...parsed.options,
    isolatedDeclarations: true,
    incremental: false,
    emitDeclarationOnly: true,
});
const diagnostics = program.getDeclarationDiagnostics();
console.log(diagnostics.length);
console.log(new Set(diagnostics.map(d => d.file?.fileName)).size);
```

This is a current-tree CommonJS audit snippet; it is not a proposed published runtime entry point.

### 3.2 Distribution by proposed package

Grouping follows today's folders before the architecture's dependency relocations. File counts are affected files, not total files in each layer.

| Proposed owner / current scope | Diagnostics | Files |
| --- | ---: | ---: |
| Core: util, data, math, task, state | 1,125 | 137 |
| IO | 1,441 | 73 |
| Model: model, formats, properties, script | 951 | 175 |
| Graphics: GL, geometry, representations, themes, canvas | 1,690 | 270 |
| Plugin and plugin-state | 1,080 | 112 |
| Plugin UI | 704 | 50 |
| Extensions, including MVS | 717 | 167 |
| Apps, including published Viewer | 368 | 26 |
| CLI and servers | 165 | 50 |
| Examples and browser/performance tests | 166 | 33 |
| **Total** | **8,407** | **1,093** |

Broad leaf exports make a large fraction of the implementation a supported public API. Narrowing the package root alone will not remove the obligations of those leaf entry points. Keep the planned import surface; remove an export only when it is intentionally private and record any migration consequence.

### 3.3 Diagnostic categories

| Category | TS codes | Count |
| --- | --- | ---: |
| Function/method returns and accessors | 9007–9009 | 3,953 |
| Expression, spread, shorthand, and array inference | 9013, 9015–9017 | 3,495 |
| Variable and property annotations | 9010, 9012 | 749 |
| Parameters and implicit `undefined` | 9011, 9025 | 110 |
| Enum initializers, class factories/heritage, computed keys, private names | 9020–9022, 9038–9039 | 100 |

Most edits should be annotations, but the smaller final category and inference-heavy schemas carry disproportionate design risk.

### 3.4 Deno verification

Deno 2.9.6 ran `deno lint` with only `no-slow-types` enabled against a temporary copy of all 1,590 implementation modules from the TypeScript audit. The temporary manifest exported every module explicitly and its lint include list matched those files. No repository configuration was changed, no diagnostics were suppressed, and the lint result contained no processing errors.

| Deno diagnostic | Count |
| --- | ---: |
| Missing explicit return type | 3,141 |
| Missing explicit public type | 1,266 |
| Complex superclass expression | 60 |
| Unsupported public destructuring | 2 |
| Unresolved public type reference | 1 |
| **Total, across 957 files** | **4,470** |

This is a real Deno rule result for an intentionally broad module export surface, not proof about the future workspace packages. It does not replace a future isolated-declaration check or validate assets, package dependencies, and runtime behavior. The lower diagnostic count does not imply proportionally less work: the checkers report inference problems at different granularities. Neither diagnostic count is a v6 release gate.

The `Vec3` annotation pilot also reduced its Deno rule result from **62 to zero**. A `deno publish --dry-run --no-remote --no-lock` attempt on the temporary single-entry package first failed to resolve the current extensionless imports. Adding `--sloppy-imports` allowed typechecking and reached the publication fast-type check, which reported **116 missing-return-type diagnostics in referenced modules and none in the edited `Vec3` file**. No package was uploaded and no slow-types escape was used. A file passing lint is therefore a useful batch result, but the complete package graph still needs to pass publication.

After choosing the reduced v6 scope, the original, unannotated `Vec3` source was restored in the temporary package. `deno publish --dry-run --allow-slow-types --sloppy-imports --no-remote --no-lock` then passed with normal typechecking enabled. This verifies the allowance on one source graph, not the full workspace, registry upload, or generated documentation/types. Nothing was published.

## 4. Where the work is

### 4.1 Routine contracts

Add return types to exported functions and methods, types to public properties/accessors, and named contracts for factory results. Use the compiler's existing inferred declarations to suggest annotations, then review them. Preserve overloads, predicates, generics, literal unions, optionality, and mutability. Avoid copying giant anonymous declaration structures into source when an existing named interface expresses the contract.

Annotation-only changes should leave JavaScript behavior intact. Some annotations change contextual typing, however, so a clean declaration check is insufficient: run the normal typecheck and compare public type behavior too.

### 4.2 Parameter schemas, providers, and render schemas

[Ball-and-stick parameters](../src/mol-repr/structure/representation/ball-and-stick.ts) illustrate the main pattern: an exported object combines spreads and `PD.*` calls, and its inferred shape drives `PD.Values`, provider generics, and user completion. Similar patterns occur throughout themes, visuals, transforms, and [renderable schemas](../src/mol-gl/renderable/schema.ts), which alone produce 152 diagnostics.

Give these objects precise named definition types. Preserve every key and select-option union; compose types with explicit override semantics when later object spreads replace earlier fields. A blanket `PD.Params`, `Record<string, unknown>`, or `any` annotation would discard useful API information. `satisfies` checks compatibility but leaves inference in place, so it is not a general replacement for an annotation.

Avoid circular fixes such as annotating a value with a type alias defined as `typeof` that same value. Define the shape independently; other APIs can then safely refer to it. Type-level mapped/conditional utilities can remain when their inputs and exported contracts are explicit.

### 4.3 Generated schemas

Files under [CIF schemas](../src/mol-io/reader/cif/schema) account for **1,138 diagnostics**, including 805 in `mmcif.ts`. Update the [schema generator](../src/cli/cifschema/util/generate.ts) to emit precise schema types together with annotated values, then regenerate from pinned dictionaries. Preserve category/column names, aliases, value types, and documentation. Do not hand-annotate thousands of generated fields or broaden the result to a generic database schema.

The earlier decision against syntax-driven schema regeneration still stands; this regeneration would be specifically for public type contracts. Audit other schema generators for the same issue.

### 4.4 Class factories and namespace APIs

[StateObject.create](../src/mol-state/object.ts) returns a class expression, and [plugin state objects](../src/mol-plugin-state/objects.ts) extend factory results. The audit reports 61 expression-based `extends` diagnostics and two inferred class-expression diagnostics. Introduce precise constructor/static-side contracts and named, typed base bindings before subclassing. Preserve constructor arguments, instance data, static guards, and nominal/structural relationships; the existing broad `StateObject.Ctor` is not automatically a lossless replacement.

Namespace APIs can stay. Annotate their exported members and handle the 17 enum-initializer diagnostics individually without changing numeric values or hot-path behavior. An alias such as `import Schema = Column.Schema` is a namespace alias, not a CommonJS `require` import; do not conflate the two.

### 4.5 UI, ambient types, and dependencies

UI work includes component returns, state/property types, and inferred exported icon/component values. Retain React's public type precision and the architecture's consumer-facing React type dependency policy.

The local ambient scan found [background image typings](../src/extensions/backgrounds/typings.d.ts) declaring `*.jpg` with `export =`. Treat this as a bundler-asset boundary to resolve for JSR, not an indication that changing one declaration makes image imports portable. Audit WebXR/Node/React types and other dependency declarations in the actual publish graph; their ambient requirements are additional to this local diagnostic count. Keep browser assets, optional native code, and runtime-specific adapters explicitly owned.

## 5. Effort and adoption

**Deferred planning for possible v7 work.** These estimates and the sequence below are retained for comparison, not scheduled for v6.

These are engineering estimates for someone familiar with Mol*, using compiler-assisted edits and generators. They are not measured migration throughput and exclude the existing package moves, plugin redesign, and rendering-backend work.

| Work | Estimated engineer-days |
| --- | ---: |
| Pilot, real JSR checks, and type-regression fixtures | 3–5 |
| Schema generator changes and regeneration | 3–6 |
| Routine annotations and review across layers | 8–15 |
| Parameter/provider contracts, factories, and UI | 10–20 |
| Declaration/API validation and CI integration | 5–10 |
| **Total** | **29–56, approximately 6–12 weeks** |

For the pilot, include a math namespace, one generated CIF category, a representation parameter/provider pair, a state-object subclass, a React component, and an MVS builder entry point. The current MVS extension has 167 diagnostics, but includes runtime code and is not the future builder's isolated count. Do not assume the builder is already fast-type compliant.

Land the work in this order:

1. Establish the export inventory and record baseline declarations/type-consumer behavior. Validate the pilot with the pinned Deno publisher.
2. Fix shared contracts and generators, then annotate core/IO/model, graphics, plugin/UI, and extensions in dependency order. Coordinate relocations so annotations do not import types upward.
3. Enable `isolatedDeclarations` per migrated project and prevent regressions. Extend the compiler policy to private workspace code without publishing that code.
4. Claim fast-type support for a package only after its entire intended public surface passes without slow types. Until then, retain the v6 allowance when publishing that package to JSR.
5. Finish with clean npm and JSR consumers: validate schema keys, parameter defaults/options, transformer data types, constructor/static behavior, and declaration completeness. Measure consumer checking separately; no runtime speedup is promised.

Budget a further **1–2 engineer-weeks** for an initial dual-registry pipeline and portable-package fixtures once the types and package layout are ready. Broad graphics/UI/server coverage needs package-specific asset, dependency, and runtime validation; it cannot be estimated from declaration counts alone.

### 5.1 Agent-assisted migration

**Feasible, with compiler-assisted edits and explicit acceptance checks.** Agents can drive diagnostic collection, suggest named contracts, run codemods, verify batches, and repair failures. Give the compiler responsibility for determining existing inferred types wherever possible; use agents for selecting understandable type names and resolving the cases that need design.

The follow-up pilot used an in-memory copy of [Vec3](../src/mol-math/linear-algebra/3d/vec3.ts). It queried TypeScript's inferred return types for exported function declarations lacking annotations and inserted those types. No repository source was edited.

| Pilot check | Result |
| --- | --- |
| Return annotations added | 62 |
| Isolated-declaration diagnostics in the module | 63 → 0 |
| Deno `no-slow-types` diagnostics in the module | 62 → 0 |
| Full-project semantic diagnostics with the candidate | 0 |
| Emitted module JavaScript compared with baseline | Byte-identical |
| Emitted module declaration compared with baseline | Byte-identical |

This demonstrates a tractable mechanical case, not a migration success rate for the whole repository. The pilot did not run autonomous worker agents or fix a parameter schema/class factory. Deno verified the local annotations, while publication remained blocked by referenced modules (§3.4). These results support feasibility but do not measure agent speed, total cost, or how many files can be migrated without review.

| Work | Agent suitability | Required review |
| --- | --- | --- |
| Routine returns, accessors, and properties | High; compiler-derived candidates and repeatable checks | Public type/declaration changes and unexpected contextual typing |
| Generated CIF schemas | High once the output type pattern is agreed | Generator, precise keys/value types, and reproducible regeneration |
| Parameter/provider and render schemas | Medium; agents can expand an established pattern | Literal unions, overrides, defaults, and generic relationships |
| Class factories and merged namespace APIs | Medium; few shared decisions affect many callers | Constructor/static contracts, narrowing, and API compatibility |
| Publish graphs, assets, optional dependencies | Medium; suitable for fixture-driven investigation | Real runtime and registry behavior |

Use a lead maintainer/integrating agent and **2–4 workers with disjoint ownership**. More workers are not automatically faster: shared type definitions, generated output, and import moves create a serial integration path. Use separate worktrees or equivalent isolated patches, and avoid having multiple workers change the same foundational contract.

1. **Pilot and harness:** cover the representative cases in §5 with actual Deno checks, baseline declarations, and positive/negative type-consumer fixtures. Record accepted batches, maintainer review time, retries, and remaining diagnostic groups. Revise the schedule after this evidence.
2. **Shared patterns:** agree on parameter-schema composition, typed class-factory bases, and generator output before parallel expansion. Give one owner each shared contract and generator; other agents consume the accepted result.
3. **Bounded batches:** start around 10–25 related files, or one generator/factory family. Each task states its owned paths, allowed API changes, baseline, and required checks. Integrate completed batches frequently before assigning dependent work.
4. **Independent validation:** a worker's “checks pass” report is evidence to inspect, not approval by itself. Re-run the relevant checks on the combined branch and review every public declaration change. A separate reviewing agent can assist; the maintainer owns deliberate API decisions.

Do not let workers satisfy diagnostics by adding `any`, broadening schemas, suppressing errors, hiding intended public exports, or introducing circular `typeof` annotations. Existing `any` need not be redesigned as part of this migration, but new loss of precision requires an explicit API decision. Keep runtime refactors and package relocation out of annotation batches.

For ordinary annotation batches, seek unchanged emitted JavaScript and declarations, as in the pilot. Named-type extraction can legitimately alter declaration text: use type-consumer fixtures for accepted/rejected calls, schema keys, optional properties, literal unions, and constructor/static behavior. Mutual assignability alone does not establish equivalence. Run isolated-declaration checks, normal typechecking, and the actual JSR gate for the affected package; use behavior tests when factories or runtime code change.

For planning, allow **2–4 calendar weeks** with one available maintainer and 2–4 workers, and **8–15 maintainer-days** spread across that period. Routine batches should become parallel; shared design and final validation remain serial. Allow **4–6 weeks or more** if factory/schema contracts require redesign, source-publishing issues emerge, or this overlaps package moves. These ranges cover fast-type migration and its validation, not full dual-registry rollout or general Deno runtime compatibility.

The earlier 6–12 engineer-week estimate describes the broader implementation/review effort with compiler assistance; the agent schedule is an alternative execution scenario, not an additional workstream or a measured speedup. Do not estimate token spend or a fixed automation percentage from the diagnostic count. Use the pilot's measured accepted work and review burden to decide whether to expand concurrency.

## 6. TypeScript distribution through npm

Shipping `.ts` files and executing them by default are different promises. npm can carry source files; ordinary Node consumers should still resolve to compiled ESM. Node explicitly refuses native type stripping under `node_modules`, so source-only npm exports would break the proposed Node contract. See [Node's dependency restriction](https://nodejs.org/api/typescript.html#type-stripping-in-dependencies).

The linked [Deno TypeScript guide](https://docs.deno.com/runtime/fundamentals/typescript/) describes Deno's source execution and checking. It does not make a TypeScript-only npm package portable to Node. Deno has also historically restricted TypeScript imports from npm packages; its [tracking issue](https://github.com/denoland/deno/issues/24093) is closed as not planned. Do not base the release on such imports: use compiled npm entry points or JSR source entry points, and verify both with the supported Deno version.

Retain the [npm conditional export map](v6-architecture.md#55-exports-and-package-contents): `types` selects declarations, `import` selects JavaScript, and the explicit `molstar-src` condition selects source for configured bundlers. No default `.ts` export or consumer install-time compilation is needed. Source files and declaration maps improve navigation, but source consumers also require compatible TypeScript, dependencies, and TSX/asset handling.

| Distribution | Library payload | Consumer path |
| --- | --- | --- |
| Native npm package | `lib/**/*.js`, `.d.ts`, declaration maps, selected `src/**`, assets | Compiled JS by default; source only with an opted-in bundler |
| Native JSR package | Source modules and required assets | JSR source resolution and runtime/bundler transpilation |
| JSR npm compatibility package | JavaScript and potentially incomplete declarations generated by JSR | npm-compatible tools download from JSR's registry; native npm is preferred for complete types |
| CDN Viewer | Bundled JS/CSS/assets | Existing script-tag contract |

JSR's [npm compatibility layer](https://jsr.io/docs/npm-compatibility) lets npm-compatible tools install a JSR package and provides generated JavaScript. It does **not** publish the native `@molstar/*` package to npmjs.com. Keep native npm releases so existing users do not need JSR registry setup.

## 7. Publishing to npm and JSR together

Yes: use the same release version and source commit, with separate manifests and publication steps. JSR accepts ESM projects using `package.json`, including npm dependencies. Its publisher supports resolving `.js` relative imports to `.ts` files in that mode, so the current source import convention need not be reversed just for JSR. See [JSR publishing rules](https://jsr.io/docs/publishing-packages).

For v6, use `deno publish --dry-run --allow-slow-types` during validation and `deno publish --allow-slow-types` for actual JSR publication. Keep normal typechecking enabled. Validate source consumers and document generated documentation/type limitations; do not require JSR's generated declarations to match native npm's `tsc` output as a release gate. Other import, asset, dependency, and runtime failures still need to be resolved before publishing the affected package.

Generate registry metadata from one package/export inventory:

- npm keeps its conditional exports, `workspace:*`/catalog conversion, peers, bins, and file allowlist.
- JSR gets explicit source entry points, the same version, and a source/asset allowlist excluding `lib/`. Its [manifest](https://jsr.io/docs/package-configuration) has a different export shape; expand the npm wildcard/index/TSX cases into explicit source entries rather than copying the npm map.
- Resolve same-package imports to the package's source graph and cross-package imports to exact published versions. Prefer JSR-to-JSR dependencies for packages mirrored there and `npm:` dependencies for external/npm-only packages. Generate the appropriate import mappings; do not leak `workspace:*`, `catalog:`, or checkout paths into the upload.

Prove self-subpath resolution, `.js`-to-source normalization, TSX, and dependency mappings with a two-package dry run. Use publisher-supported normalization first; if metadata is insufficient, stage a deterministic specifier-only rewrite from the resolved export inventory. Do not maintain a second hand-edited source tree. Published JSR imports must resolve independently of the original checkout and its npm build output.

Publishing to JSR is not certification for Deno execution: graphics/UI may require a browser, and headless/server tools may need Node or native dependencies. Mirror the MVS builder first, then core/IO/model and the remaining public libraries as their checks pass. Record intentional runtime and distribution differences, especially CSS/assets, optional peers, and CLI bins. Keep the CDN-only package and private apps/examples out of JSR library publication.

Release sequence:

1. Verify names/access on both registries, version all selected packages once, and freeze the source commit and dependency map.
2. Build and inspect native npm tarballs; run JSR dry runs with `--allow-slow-types`. Test the packed Node/TypeScript/bundler consumers and Deno/JSR source graph.
3. Publish each registry's packages in dependency order, using exact internal versions. Keep one coherent Mol* dependency graph per consumer; mixing native npm and JSR copies can duplicate state and class identities.
4. Record per-package/per-registry completion. Publication is not atomic across registries: retry missing uploads from the same artifacts/commit and never overwrite an already released version with changed code.
5. Verify registry consumers before announcing the release or advancing npm's stable tag. Treat npm dist-tags separately from JSR's version selection.

### 7.1 Could `deno pack` replace the npm build?

Potentially for suitable libraries. The documented [`deno pack`](https://docs.deno.com/runtime/reference/cli/pack/) pipeline emits npm JavaScript/declarations and a generated manifest from Deno metadata. It does not reuse the existing `package.json`, synthesize CLI `bin` entries, or run npm lifecycle hooks. JSR dependencies become `@jsr/*` dependencies requiring the corresponding registry configuration.

That is a different packaging architecture from the proposed pnpm/`tsc -b` workspace. Keep this alternative deferred: v6 retains its existing compiler pipeline and full native npm declarations. Any later evaluation must preserve peers, bins, assets, declaration maps, `molstar-src`, and native npm-to-npm Mol* dependencies, as well as declaration completeness.

## 8. Is `lib/` only temporary?

**Generated locally, shipped on npm.** In the recommended layout, each package's `lib/` contains the JavaScript and declarations referenced by its published exports and bins. It is ignored by Git, rebuilt before packing, and may be deleted from the checkout afterward. It must still be present inside the installed npm package; packing does not replace compilation. Run the build explicitly or through a configured lifecycle step before collecting the files.

A separate packaging directory could place those same compiled files at the package root and rewrite manifests/maps accordingly. Then the checkout's `lib/` would be merely a staging input, but compiled output would still ship. This provides no required capability for v6, so retain `lib/` to avoid a second layout transformation. JSR's source artifact omits it, and source-based app builds continue to avoid depending on stale library output.
