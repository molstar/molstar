# Consumer smoke harness

From the repository root, run `pnpm --dir smoke smoke` to build, pack, and run the standard suite. Focused checks use
commands such as `pnpm --dir smoke smoke:node` and require current build output. Smoke scripts live in this package
rather than the root manifest.

The runner and workspace packer use `cross-spawn` to launch npm/pnpm on macOS, Linux, and Windows, including Windows
`.cmd` shims and paths containing spaces. The CLI smoke uses the installed `.cmd` wrappers on Windows. `tar` must be
available on `PATH` (included with current Windows versions and macOS).

`node smoke/run.mjs <node|types|browser|source|slim|cli|headless|all>` exercises package artifacts as an external
consumer. `all` runs Node, declarations, source bundling, the slim plugin, the packed MVS validation command, and
browser checks. `--prepare` invokes the repository's `build:workspace` and `pack:workspace` scripts before running
checks. Focused runs require a current `scripts/workspace/inventory.json` and report missing artifacts as failures.

The harness packs each requested public package and its internal dependency closure into temporary tarballs. Consumer
projects live in the OS temporary directory and install those tarballs through local file references, so internal
imports cannot resolve through workspace links or a registry. Temporary directories are removed on exit.

The node fixture checks plain Node ESM imports, task/color/CIF behavior, plugin initialization and PDB parsing into a
structure, and an MVSJ-to-MVSX-to-MVSJ round trip with an embedded local PDB. The types fixture runs strict TypeScript
`NodeNext` checks over core, IO, plugin, React UI/TSX, headless, and MVS builder APIs using packed declarations. The
source fixture builds a small esbuild bundle with the `molstar-src` condition and checks that its inputs resolve under
package `src/` paths. The slim fixture (`smoke/slim/`) installs the packed `@molstar/plugin`, `@molstar/core`,
`@molstar/plugin-ui`, and `@molstar/model` with their internal closure, copies the `examples/slim-plugin` sources
unchanged into the consumer project, and bundles them with esbuild and the `molstar-src` condition (one IIFE file;
styles, the page, and the ligand are not bundled) together with a small harness that exposes plugin state. It checks
that the bundle resolves installed package sources and contains no excluded module, then serves the bundle, the page,
and the SDF ligand and drives them with Playwright (same browser detection as the browser checks). It asserts no page
errors, console errors, or failed requests; a visible `ball-and-stick` representation with render objects after the
ligand loads; exactly the slim providers (structure representation `ball-and-stick`, format `sdf`, the default hierarchy
and ball-and-stick presets, only the `mol-script` script language); that restoring a snapshot whose structure
representation type is `cartoon` logs the unregistered-name warning and leaves a `ball-and-stick` representation with
render objects; and that evaluating a PyMOL script fails with the script-language error. On failure it saves the console
log and a screenshot to `build/smoke-diagnostics/slim/`. The CLI fixture installs packed `@molstar/mvs-builder` and
validates a local MVS document. Browser checks first pack and extract the `molstar` distribution, then serve only those
extracted files while exercising native ESM Viewer and modular plugin pages plus classic Viewer and MVS Stories globals
and custom elements. All browser pages load local PDB/MVS fixtures over a static server. Playwright is required to claim
a browser pass; set `MOLSTAR_SMOKE_BROWSER` to a browser executable when auto-detection is unsuitable. The separate
headless mode installs packed `@molstar/plugin-headless`, loads the local PDB, captures and decodes a PNG, and reports
missing native bindings as an explicit failure. Supply native module entry paths with `MOLSTAR_SMOKE_GL` and
`MOLSTAR_SMOKE_PNGJS` when they are not resolvable in the consumer or repository install.

Missing package inventory, outputs, or browser automation are surfaced as failures rather than treated as successful
checks.

## View the browser pages

Build first, then keep the packed-artifact server running:

```sh
pnpm build
pnpm --dir smoke smoke:browser --serve
```

Open <http://127.0.0.1:1339/viewer/> for the ESM Viewer or <http://127.0.0.1:1339/library/> for the modular ESM plugin
page. Classic pages are at `/classic/viewer.html` and `/classic/mvs-stories.html`. This mode serves the same extracted
distribution and import map as the automated check; it needs no Playwright or consumer build. Stop it with Ctrl+C. Use
`--port=1340` to change the port, or add `--prepare` to build/pack first.

## Optional native capture

Native dependencies are not installed by a normal workspace install. Set them up explicitly with `pnpm native:install`,
then run `pnpm native:run -- pnpm --dir smoke smoke:headless` from the repository root. The backend is stored in
`.cache/native`; the separate native test command is `pnpm test:native`. Use `pnpm native:install -- --canvas` for the
rendering CLI.

Consumer smoke checks run locally for now; CI does not install Chromium or run the smoke suite. CI retains lint,
unit/native tests, builds, and package checks.
