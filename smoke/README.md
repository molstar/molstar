# Consumer smoke harness

From the repository root, run `pnpm --dir smoke smoke` to build, pack, and run the
standard suite. Focused checks use commands such as
`pnpm --dir smoke smoke:node` and require current build output. Smoke scripts live
in this package rather than the root manifest.

`node smoke/run.mjs <node|types|browser|source|cli|headless|all>` exercises package artifacts as an external consumer. `all` runs Node, declarations, source bundling, the packed MVS validation command, and browser checks. `--prepare` invokes the repository's `build:workspace` and `pack:workspace` scripts before running checks. Focused runs require a current `scripts/workspace/inventory.json` and report missing artifacts as failures.

The harness packs each requested public package and its internal dependency closure into temporary tarballs. Consumer projects live in the OS temporary directory and install those tarballs through local file references, so internal imports cannot resolve through workspace links or a registry. Temporary directories are removed on exit.

The node fixture checks plain Node ESM imports, task/color/CIF behavior, plugin initialization and PDB parsing into a structure, and an MVSJ-to-MVSX-to-MVSJ round trip with an embedded local PDB. The types fixture runs strict TypeScript `NodeNext` checks over core, IO, plugin, React UI/TSX, headless, and MVS builder APIs using packed declarations. The source fixture builds a small esbuild bundle with the `molstar-src` condition and checks that its inputs resolve under package `src/` paths. The CLI fixture installs packed `@molstar/mvs-builder` and validates a local MVS document. Browser checks first pack and extract the `molstar` distribution, then serve only those extracted files while exercising native ESM Viewer and modular plugin pages plus classic Viewer and MVS Stories globals and custom elements. All browser pages load local PDB/MVS fixtures over a static server. Playwright is required to claim a browser pass; set `MOLSTAR_SMOKE_BROWSER` to a browser executable when auto-detection is unsuitable. The separate headless mode installs packed `@molstar/plugin-headless`, loads the local PDB, captures and decodes a PNG, and reports missing native bindings as an explicit failure. Supply native module entry paths with `MOLSTAR_SMOKE_GL` and `MOLSTAR_SMOKE_PNGJS` when they are not resolvable in the consumer or repository install.

Missing package inventory, outputs, or browser automation are surfaced as failures rather than treated as successful checks.

## View the browser pages

Build first, then keep the packed-artifact server running:

```sh
pnpm build
pnpm --dir smoke smoke:browser --serve
```

Open <http://127.0.0.1:1339/viewer/> for the ESM Viewer or
<http://127.0.0.1:1339/library/> for the modular ESM plugin page.
Classic pages are at `/classic/viewer.html` and `/classic/mvs-stories.html`.
This mode serves the same extracted distribution and import map as the automated
check; it needs no Playwright or consumer build. Stop it with Ctrl+C. Use
`--port=1340` to change the port, or add `--prepare` to build/pack first.

## Optional native capture

Native dependencies are not installed by a normal workspace install. Set them up
explicitly with `pnpm native:install`, then run
`pnpm native:run -- pnpm --dir smoke smoke:headless` from the repository root.
The backend is stored in `.cache/native`; the separate native test command is
`pnpm test:native`. Use `pnpm native:install -- --canvas` for the rendering CLI.

Consumer smoke checks run locally for now; CI does not install Chromium or run
the smoke suite. CI retains lint, unit/native tests, builds, and package checks.
