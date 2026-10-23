# v6 work checklist

This is the consolidated list of remaining work. Detailed architecture and prototype evidence remain in
[workspace-prototype.md](workspace-prototype.md). Keep this checklist current as work is completed or new issues are
found. Accepted API changes belong in [breaking-v6-changes.md](breaking-v6-changes.md).

## Current PR review follow-up

Source: [PR #1951 review](https://github.com/molstar/molstar/pull/1951#pullrequestreview-5406406229).

- [x] Fix Node MP4 encoding: use the default `h264-mp4-encoder` import. Cover native Node ESM as well as browser
      bundles.
- [x] Fix the pre-existing headless MP4 behavior check: check the `VideoExport` transformer rather than its unqualified
      name. Verify that animation export reaches encoding; keep native GL installation optional.
- [x] Support caller-controlled MVSX ZIP options under the third `createMVSX` parameter's `zip` property. Explicit
      `mtime` enables identical bytes across different clock times; the default uses the current time. Verify both
      behaviors and identical archive contents; preserve archive loading compatibility without requiring v5-identical
      bytes.
- [x] Eliminate clean-install missing-bin warnings with checked-in executable launchers and correct packed-file
      inclusion. Audit all workspace bin entries, including `mvs-validate` and `mvs-print-schema`. Do not build during
      install.
- [x] Restore original formatting in `location-iterator.ts`. Preserve the required `PositionLocation` extraction and
      re-export.
- [x] Keep the extracted `Viewport` implementation and camera utility formatting close to the original source, including
      attribution. Preserve the package split. A follow-up formatting commit will not change the historical move commit.
- [x] Optimize clean TypeScript builds using measured compiler diagnostics. Evaluate declaration checking, build scope,
      and compiler migration; preserve package boundaries, source checks, declaration output, and consumer checks.
- [x] Audit and record small public API changes caused by the refactor, including classic globals, shape
      signatures/helpers, headless integration, and MVS validation. The ledger records known shape, bounds, theme
      registry, headless, declaration, and distribution changes. The 2026-10-04 checkpoint compares 1,547 mapped
      production modules and 17 extracted/new modules with v5 master. Keep updating the ledger as the remaining v6 work
      proceeds.
- [x] Restore cache-first `MVSData.toMVSX` asset fetching: a supplied cache must avoid requests and allow cached exports
      when the source is unavailable.
- [x] Define Node local-file asset loading for standalone MVS export through optional `toMVSX` `options.fetch`. Callers
      can provide a file-aware fetch adapter; the default remains platform `fetch`. Explicit `assets` also bypasses
      fetching. Keep the builder independent of plugin/rendering code.

## Deferred validation and tooling

- [x] Bring back `external-structure` and `external-volume` color themes during registry composition work. Done in
      plugin composition: the `ExternalColorThemes` entry (`@molstar/plugin/themes/external`) is part of
      `DefaultRegistry`; `PluginContext` has no special registration logic.

- [x] Add MolQL syntax and field-name validation to `mvs-validate` and the builder, against symbol and
      argument-definition tables without the molecular query runtime. Reject malformed expressions, unknown calls, and
      bad argument names; include file/expression paths and a nonzero CLI exit status. Custom symbols use an explicit
      definition lookup. Argument type, required-argument, and selection-result checks are outside this pass. See the
      [MolQL design](../designs/architecture.md#72-molql-builder-and-validation-design).
- [x] Expose `MolScriptBuilder` and `compileScript` through `@molstar/mvs-builder/molql`, backed by the standalone
      `@molstar/query-language` package. Keep executable molecular compilation/evaluation in model.
- [x] Provide text-to-MolQL translation for MolScript, PyMOL, VMD, and Jmol. Preserve explicit per-language and
      all-languages registration for plugin use; keep text parsers out of JSON builder/validation imports. Verify the
      standalone dependency boundary and packed authoring/CLI consumers.
- [ ] Verify standalone-builder API/schema/serialization parity before replacing `molviewspec-ts` and JSR
      `@molstar/molviewspec`. Existing packed smoke checks cover basic MVSJ/MVSX round trips; replacement parity and JSR
      source publication remain separate release tasks.
- [x] Migrate compilation to TypeScript 7, including project references, declaration checks, ESM/source conditions,
      incremental builds, and CI. Retain the separate TypeScript 6 compatibility API for AST/config parsing only.
- [x] Replace ESLint with Biome, documenting rule differences and configuring optional formatting.
- [x] Reformat the repository after the TypeScript 7 migration, then enable formatting checks in CI.
- [ ] Revisit build identification metadata beyond the release version, including whether to expose a Git commit
      fingerprint and how to support Windows and source archives. For now, generate and display only the version from
      `version.json`.
- [x] Compare clean and incremental builds after tooling changes and rerun the install, lint, test, build, workspace,
      version, tarball, and local smoke checks.
- [x] Move the repository's tests to Vitest, ESM throughout (this replaces "Move Jest to ESM"): one run for every
      package (`pnpm test:unit`, `vitest.config.ts`), including the AI packages (`pnpm jest:ai` and its config are gone)
      and the workspace tooling tests (`scripts/workspace/_test`, formerly `node --test`); `molstar-src` through
      `resolve.conditions` and `ssr.resolve.conditions`, so a new package needs no mapping; explicit `vitest` imports
      (no globals); `jest.*` to `vi.*`; `@jest-environment` docblocks to `@vitest-environment`. Four io tests that
      imported `@molstar/model` moved to the model package, which the Jest mapper had hidden. The sequence packages
      (`packages/sequences/_spec`) have joined this run (the `sequences` projects, the quick lane only; their other tiers are each
      package's own scripts).

Tooling verification (2026-10-05): native TypeScript 7.0.2 forced compilation of all 67 projects took 7.28 s, compared
with 22.01 s for TypeScript 6.0.3 using the same `skipLibCheck: true` settings. Incremental compiler runs took 0.45 s
and 0.48 s respectively. A clean `pnpm build:lib` (including assets/version staging) took 6.93 s. Fresh frozen-lockfile
installation, Biome lint, 1,504 unit tests (14 optional native tests skipped), six workspace tooling tests, the full
native declaration check, app/distribution builds, workspace/version checks, all 43 public tarballs, and local
Node/types/source/CLI/browser smokes passed. The native declaration gate also runs in CI. A CIF writer namespace alias
needed explicit qualification to preserve TS6 runtime behavior under TS7 emit; public exports are unchanged.
Repository-wide formatting is complete; CI enforces Biome/Prettier formatting.

Formatting checkpoint (2026-10-05): settings were committed separately from the 1,940-file formatting pass. Source and
configuration use two spaces; MkDocs Markdown retains four-space nested blocks for Python-Markdown compatibility.
Formatting and lint checks, 1,504 unit tests (14 optional native tests skipped), library/app/distribution builds, the
full native declaration check, workspace/version checks, and all 43 public tarballs passed. A comparison of all 1,777
changed code files found only equivalent syntax and JSX text changes. The existing toast JSX `@ts-ignore` was reattached
to its expression after wrapping. The formatting commit is recorded in `.git-blame-ignore-revs`.

Build evidence (2026-10-04): a forced rebuild of all 67 compiler projects took 149.10 s: aggregate checking 137.45 s,
emitting 6.64 s. The MP4 project loads 5 own source files, 433 workspace declaration files, and 364 external declaration
files. Separate no-emit checks took 3.41 s normally and 0.63 s with `skipLibCheck`. With `skipLibCheck` enabled, the
same forced workspace rebuild took 22.34 s, including 11.02 s checking and 6.82 s emitting (about 6.7× faster overall).
Normal builds skip declaration checking; the required `pnpm check:publish` gate forces it on for every project and
includes packed consumer validation. `pnpm check:types:full` also runs that declaration check independently.

Review-fix verification (2026-10-04): clean offline source installation linked all CLI bins before compiling and
produced no missing-bin warnings. Lint and 1,492 tests passed (14 optional native tests skipped); the full declaration
check, workspace/version checks, app/distribution builds, all 43 public tarballs, and local
Node/types/source/CLI/browser smokes passed. The Node smoke encodes real MP4 bytes without native GL and checks the
headless registration guard. MVSX bytes with explicit `mtime` match across two mocked clock times. Local date fields
allow the caller to choose timestamps that are reproducible across time zones.

## Remaining architecture and release work

- [ ] Seq\*, the sequence viewer, ported from Stories Untold into `packages/sequences` ([spec](../../packages/sequences/_spec/README.md)):
  workspace packages, the fast tier in the root Vitest run, Mol\*'s msgpack and zip, the Mol\* adapter on v6 (folded into the Prism viewer, `sequence/molstar/`), its Storybook
  and lint rules; then Prism over it ([spec](../../packages/prism/viewer/_spec/sequence.md)).
- [ ] Rendering-backend extraction and GL resource/pass/readback redesign (owner: Alex; separate workstream).
- [x] Plugin composition per the [design](../designs/plugin-composition.md) and [plan](plugin-composition.md): empty
      registries, `spec.registry` entries, explicit base specs, presets that import what they run, transformer/catalog
      splitting, snapshot pre-validation, and slim-plugin guarantees.
- [ ] Remaining convenience-barrel cleanup. The `StateTransforms` facade was removed in plugin composition step 1 (the
      classic Viewer global keeps an app-level object).
- [ ] Broader test-runner migration and maintainer skills/documentation rewrite.
- [ ] Add versioned documentation with a v5/v6 version switch. Preserve the existing v5 docs and update the v6 pages
      separately, including the video-export extension and API changes.
- [ ] In the documentation rewrite, cover plugin composition. `docs/` is deployed from `master` and still describes the
      published 5.x `molstar/lib/...` paths, so these pages keep their 5.x imports until then:
      `docs/docs/plugin/custom-library.md` and `instance.md` (`DefaultPluginSpec`/`DefaultPluginUISpec` now come from
      `@molstar/plugin/default-spec` and `@molstar/plugin-ui/default-spec`; literal specs add `DefaultRegistry`),
      `docs/docs/extensions/tunnels.md` (`StateTransforms` is gone; import `ShapeRepresentation3D` directly), and
      `docs/docs/plugin/selections.md` (selection queries are split under `@molstar/plugin/state/queries/structure/*`).
- [ ] JSR publication, migration CLI, and downstream migration validation.
- [ ] Publishing automation, release coordination, and stable-release readiness.
- [ ] Write the v6 `CHANGELOG.md` section from the [changelog draft](changelog.md).
- [ ] Before v6 is ready, forward-port every bug fix from the `v5` branch (`git log master..origin/v5`). Most files have
      moved, so apply fixes by hand using `migration-map.json` rather than by merging. Ported through 5.13.1
      (2026-10-08): molstar/molstar#1964, server path interpolation and dependency updates.
- [ ] Capture a representative pre-migration Viewer render comparison baseline.
- [ ] Bring the insertion codes of ranges to MolViewSpec selectors: `beg_auth_ins_code` and `end_auth_ins_code` are in
      `StructureElement.SchemaItem` (molstar/molstar#1972) but not in the MVS selector and annotation schemas
      (`packages/mvs/runtime/src/helpers/schemas.ts`, `MVSAnnotationRow`, the `auth_residue_range` annotation schema, the
      builder and the JSON schema). Settle the semantics in the MolViewSpec spec first: with a code on a bound, a range
      is a range of residues in the order of the input file, in each auth chain, not of `auth_seq_id` values (see
      [the ledger](breaking-v6-changes.md#mvs-selectors-through-model-indices-the-molql-generator-element-schema)).
- [ ] Give bond providers a MolViewSpec schema field once the interface has settled. Until then the structure node
      takes them through `custom.molstar_bond_provider` only ([bond-providers plan §4.2](bond-providers.md)). The
      field is the structure node's `bonds` (`{ name, params? }`) in the schemas, the builder and the JSON schema. The
      MVS registry test (`mvs/runtime/src/_test/registry.test.ts`) then checks the names a file emits.
- [ ] Add the periodic box keys to the Python `molviewspec` package (the builder's `custom` helpers and docs):
      `custom.molstar_periodic` on a structure node (`'whole' | 'box' | 'nojump' | 'off'`, or `{ mode, center }` with
      center `'box'`, `'polymers'` or `{ expression, language? }`) and `custom.molstar_boxes` on a coordinates node
      (`{ url }`, a NAMD `.xsc`/`.xst`, relative to the view). Their semantics: the MVS runtime's load.md §6.

## Deferred beyond v6

- Full WebGPU rendering/parity and Blender integration. The v6 scope is the rendering-backend boundary with preserved
  WebGL behavior; see the [rendering design](../designs/webgpu.md).
- Fast types, isolated declarations, and the corresponding API annotations. Consider these for v7; v6 permits slow types
  for validated JSR packages. See the [distribution and fast-types design](../designs/fasttypes.md).

## Other review findings to revisit

These were reported as existing on master, rather than refactor regressions.

- [x] Investigate the `state-docs` crash in `getOrientationParticlesParams`. Fixed in plugin composition step 4:
      data-less param getters (orientation particles, particle target, operator-hkl theme) accept missing data.
- [ ] Investigate the membrane-orientation server's `--bcifSource` handling.
- [ ] Recheck reported Node startup and ball-and-stick timing differences with repeated runs before making performance
      changes.
- [ ] Fix Viewer URL parameter parsing: `getParam` in `apps/viewer/src/index.html` matches substrings (the pattern is
      `name=` without a boundary). For example, `structure-url=` is also read as `url=` (and `structure-url-format=` as
      `url-format=`), so the file can be loaded twice. Anchor the pattern to `[?&]` or use `URLSearchParams`.
- [x] Fix two GLSL clipping quirks: the clip functions test the bit `1 << i` of `Clipping.Groups.Flag` (was the
      mask `i + 1`, so groups were exempted from the wrong objects from the third clip object on), in GLSL and WGSL;
      `assign-clipping-varying.glsl.ts` reads groupInstance clipping under `dClippingType_groupInstance` (was
      `dMarkerType_groupInstance`). Visual case `renderer/clipping/groups`; no other case changes.
