# v6 work checklist

This is the consolidated list of remaining work. Detailed architecture and prototype evidence remain in
[workspace-prototype.md](workspace-prototype.md). Keep this checklist current as work is completed or new issues are
found. Accepted API changes belong in [breaking-v6-changes.md](breaking-v6-changes.md).

## Current PR review follow-up

Source: [PR #1951 review](https://github.com/molstar/molstar/pull/1951#pullrequestreview-5406406229).

- [x] Fix Node MP4 encoding: use the default `h264-mp4-encoder` import. Cover native Node ESM as well as browser
      bundles.
- [x] Fix the pre-existing headless MP4 behavior check: check the `Mp4Export` transformer rather than its unqualified
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

- [ ] Add MolQL syntax and field-name validation to `mvs-validate` and the builder. Design a validation pass against
      symbol and argument-definition tables without loading the molecular query runtime; reject malformed expressions,
      unknown callable symbols, and invalid argument names with a nonzero CLI exit status. Argument type checking,
      required-argument checks, and selection-result checking are outside this pass. Until implemented, the builder
      checks expression shape and the runtime performs compiler validation. See the
      [MolQL design discussion](../designs/architecture.md#72-molql-builder-and-validation-design).
- [ ] Expose the MolQL expression builder to `mvs-builder` consumers. Decide the mol-script package boundary and public
      entry point together with validation; keep the builder free of plugin/rendering dependencies. This work is in
      design, with no validation or MolQL implementation changes yet.
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

- [ ] Rendering-backend extraction and GL resource/pass/readback redesign (owner: Alex; separate workstream).
- [x] Plugin composition per the [design](../designs/plugin-composition.md) and [plan](plugin-composition.md): empty
      registries, `spec.registry` entries, explicit base specs, presets that import what they run, transformer/catalog
      splitting, snapshot pre-validation, and slim-plugin guarantees.
- [ ] Remaining convenience-barrel cleanup. The `StateTransforms` facade was removed in plugin composition step 1 (the
      classic Viewer global keeps an app-level object).
- [ ] Broader test-runner migration and maintainer skills/documentation rewrite.
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
