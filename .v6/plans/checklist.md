# v6 work checklist

This is the consolidated list of remaining work. Detailed architecture and
prototype evidence remain in [workspace-prototype.md](workspace-prototype.md).
Keep this checklist current as work is completed or new issues are found.
Accepted API changes belong in [breaking-v6-changes.md](breaking-v6-changes.md).

## Current PR review follow-up

Source: [PR #1951 review](https://github.com/molstar/molstar/pull/1951#pullrequestreview-5406406229).

- [x] Fix Node MP4 encoding: use the default `h264-mp4-encoder` import. Cover
  native Node ESM as well as browser bundles.
- [x] Fix the pre-existing headless MP4 behavior check: check the `Mp4Export`
  transformer rather than its unqualified name. Verify that animation export
  reaches encoding; keep native GL installation optional.
- [x] Support caller-controlled MVSX ZIP options under the third `createMVSX` parameter's `zip` property.
  Explicit `mtime` enables identical bytes across different clock times; the default
  uses the current time. Verify both behaviors and identical archive contents;
  preserve archive loading compatibility without requiring v5-identical bytes.
- [x] Eliminate clean-install missing-bin warnings with checked-in executable
  launchers and correct packed-file inclusion. Audit all workspace bin entries,
  including `mvs-validate` and `mvs-print-schema`. Do not build during install.
- [x] Restore original formatting in `location-iterator.ts`. Preserve the required
  `PositionLocation` extraction and re-export.
- [x] Keep the extracted `Viewport` implementation and camera utility formatting
  close to the original source, including attribution. Preserve the package split.
  A follow-up formatting commit will not change the historical move commit.
- [x] Optimize clean TypeScript builds using measured compiler diagnostics.
  Evaluate declaration checking, build scope, and compiler migration; preserve
  package boundaries, source checks, declaration output, and consumer checks.
- [x] Audit and record small public API changes caused by the refactor, including
  classic globals, shape signatures/helpers, headless integration, and MVS
  validation. The ledger records known shape, bounds, theme registry, headless, declaration,
  and distribution changes. The 2026-10-04 checkpoint compares 1,547 mapped
  production modules and 17 extracted/new modules with v5 master. Keep updating
  the ledger as the remaining v6 work proceeds.
- [ ] Restore cache-first `MVSData.toMVSX` asset fetching: a supplied cache must
  avoid requests and allow cached exports when the source is unavailable.
- [ ] Define/restore Node local-file asset loading for standalone MVS export.
  Platform `fetch` does not support the old `file://` path; explicit `assets` is
  the current workaround. Keep the builder independent of plugin/rendering code.

## Deferred validation and tooling

- [ ] Bring back `external-structure` and `external-volume` color themes during
  registry composition work. Their implementations depend on `PluginStateObject`
  and, for structures, the plugin backbone selection query. They remain in
  `@molstar/plugin/themes/*` but are intentionally unregistered in this prototype.
  Resolve those dependencies and define explicit composition; do not add special
  registration logic to `PluginContext`.

- [ ] Restore full MolQL validation in `mvs-validate` by importing mol-script.
  Decide how to expose the compiler dependency without making the standalone
  builder depend on plugin/rendering code. Then reject unknown symbols and invalid
  expressions with a nonzero CLI exit status. Until this is implemented, the
  builder validates expression shape and the runtime performs compiler validation.
  This is deferred to settle the integration design, not to remove full validation.
- [ ] Migrate TypeScript 6 to TypeScript 7, including project references,
  declaration checks, ESM/source conditions, incremental builds, and CI.
- [ ] Replace ESLint with Biome, documenting rule differences and preserving
  repository formatting. Keep broad formatting changes separate.
- [ ] Compare clean and incremental builds after tooling changes and rerun the
  install, lint, test, build, workspace, version, tarball, and local smoke checks.

Build evidence (2026-10-04): a forced rebuild of all 67 compiler projects took
149.10 s: aggregate checking 137.45 s, emitting 6.64 s. The MP4 project loads
5 own source files, 433 workspace declaration files, and 364 external declaration
files. Separate no-emit checks took 3.41 s normally and 0.63 s with `skipLibCheck`.
With `skipLibCheck` enabled, the same forced workspace rebuild took 22.34 s,
including 11.02 s checking and 6.82 s emitting (about 6.7× faster overall).
Normal builds skip declaration checking; the required `pnpm check:publish` gate
forces it on for every project and includes packed consumer validation.
`pnpm check:types:full` also runs that declaration check independently.

Review-fix verification (2026-10-04): clean offline source installation linked
all CLI bins before compiling and produced no missing-bin warnings. Lint and
1,492 tests passed (14 optional native tests skipped); the full declaration check,
workspace/version checks, app/distribution builds, all 43 public tarballs, and
local Node/types/source/CLI/browser smokes passed. The Node smoke encodes real
MP4 bytes without native GL and checks the headless registration guard. MVSX
bytes with explicit `mtime` match across two mocked clock times. Local date fields
allow the caller to choose timestamps that are reproducible across time zones.

## Remaining architecture and release work

- [ ] Rendering-backend extraction and GL resource/pass/readback redesign.
- [ ] WebGPU and Blender integration.
- [ ] Plugin features, empty registries, explicit base specs, registry-aware presets,
  slim-plugin guarantees, and transformer/catalog splitting.
- [ ] Remaining convenience-barrel cleanup and `StateTransforms` facade removal.
- [ ] Broader test-runner migration and maintainer skills/documentation rewrite.
- [ ] Fast types, isolated declarations, and the corresponding API annotations.
- [ ] JSR publication, migration CLI, and downstream migration validation.
- [ ] Publishing automation, release coordination, and stable-release readiness.
- [ ] Capture a representative pre-migration Viewer render comparison baseline.

## Other review findings to revisit

These were reported as existing on master, rather than refactor regressions.

- [ ] Investigate the `state-docs` crash in `getOrientationParticlesParams`.
- [ ] Investigate the membrane-orientation server's `--bcifSource` handling.
- [ ] Recheck reported Node startup and ball-and-stick timing differences with
  repeated runs before making performance changes.
