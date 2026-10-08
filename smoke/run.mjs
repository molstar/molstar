import assert from 'node:assert/strict';
import spawn from 'cross-spawn';
import { cp, mkdir, mkdtemp, readFile, readdir, rm, writeFile, access } from 'node:fs/promises';
import { createServer } from 'node:http';
import { tmpdir } from 'node:os';
import { dirname, extname, isAbsolute, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { createRequire } from 'node:module';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, '..');
const rootRequire = createRequire(import.meta.url);
const args = new Set(process.argv.slice(2));
const mode = [...args].find((arg) => !arg.startsWith('--')) ?? 'all';
const validModes = ['node', 'types', 'browser', 'source', 'slim', 'cli', 'headless', 'all'];
if (args.has('--serve') && mode !== 'browser') throw new Error('--serve requires browser mode.');
if (!validModes.includes(mode)) throw new Error(`Unknown smoke mode '${mode}'. Choose ${validModes.join(', ')}.`);
const run = (command, commandArgs, options = {}) =>
  new Promise((resolveRun, reject) => {
    const child = spawn(command, commandArgs, {
      cwd: options.cwd ?? root,
      stdio: 'inherit',
      env: process.env,
      ...options,
    });
    child.once('error', reject);
    child.once('exit', (code) =>
      code === 0 ? resolveRun() : reject(new Error(`${command} ${commandArgs.join(' ')} exited with ${code}`)),
    );
  });
const capture = (command, commandArgs, options = {}) =>
  new Promise((resolveRun, reject) => {
    const child = spawn(command, commandArgs, {
      cwd: options.cwd ?? root,
      env: process.env,
      ...options,
      stdio: ['ignore', 'pipe', 'pipe'],
    });
    let stdout = '',
      stderr = '';
    child.stdout.setEncoding('utf8').on('data', (chunk) => (stdout += chunk));
    child.stderr.setEncoding('utf8').on('data', (chunk) => (stderr += chunk));
    child.once('error', reject);
    child.once('exit', (code) =>
      code === 0
        ? resolveRun({ stdout, stderr })
        : reject(new Error(`${command} ${commandArgs.join(' ')} exited with ${code}\n${stderr}`)),
    );
  });
const exists = async (path) =>
  access(path).then(
    () => true,
    () => false,
  );
const temp = await mkdtemp(join(tmpdir(), 'molstar-smoke-'));
const cleanup = () => rm(temp, { recursive: true, force: true });
process.once('exit', () => {
  void cleanup();
});
const fail = (message) => {
  throw new Error(message);
};

async function prepare() {
  const rootPackage = JSON.parse(await readFile(join(root, 'package.json'), 'utf8'));
  for (const script of ['build:workspace', 'pack:workspace']) {
    if (!rootPackage.scripts?.[script]) fail(`--prepare requires the root '${script}' script.`);
    await run('npm', ['run', script]);
  }
}

async function inventory() {
  const path = join(root, 'scripts/workspace/inventory.json');
  if (!(await exists(path))) fail(`Missing ${relative(root, path)}. Generate the workspace package inventory first.`);
  const data = JSON.parse(await readFile(path, 'utf8'));
  assert(Array.isArray(data.packages), 'inventory.json must contain packages[]');
  return data.packages;
}

async function catalogVersion(packageName) {
  const workspace = await readFile(join(root, 'pnpm-workspace.yaml'), 'utf8');
  const escaped = packageName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const match = workspace.match(new RegExp(`^\\s+'?${escaped}'?:\\s*['\"]?([^'\\\"\\s]+)['\"]?\\s*$`, 'm'));
  if (!match) fail(`Could not read catalog version for ${packageName} from pnpm-workspace.yaml.`);
  return match[1];
}

async function unpackDistribution() {
  const packages = await inventory();
  const entry = packages.find((pkg) => pkg.name === 'molstar');
  if (!entry) fail('Distribution smoke requires the molstar package in inventory.json.');
  const packDir = join(temp, 'distribution-tarball');
  await mkdir(packDir, { recursive: true });
  await run('pnpm', ['pack', '--pack-destination', packDir], { cwd: resolve(root, entry.path) });
  const archives = (await readdir(packDir)).filter((name) => name.endsWith('.tgz'));
  if (archives.length !== 1) fail(`Expected one packed molstar distribution tarball; found ${archives.length}.`);
  const unpackDir = join(temp, 'packed-distribution');
  await mkdir(unpackDir, { recursive: true });
  await run('tar', ['-xzf', join(packDir, archives[0]), '-C', unpackDir]);
  const packageDir = join(unpackDir, 'package');
  if (!(await exists(join(packageDir, 'build')))) fail('Packed molstar distribution has no build directory.');
  return packageDir;
}

async function packClosure(packages, roots) {
  const byName = new Map(packages.map((p) => [p.name, { ...p, path: resolve(root, p.path) }]));
  const selected = new Set();
  const visit = async (name) => {
    const entry = byName.get(name);
    if (!entry) fail(`Required internal package '${name}' is missing from inventory.json.`);
    if (selected.has(name)) return;
    selected.add(name);
    const manifestPath = join(entry.path, 'package.json');
    const manifest = JSON.parse(await readFile(manifestPath, 'utf8'));
    for (const dep of Object.keys({
      ...manifest.dependencies,
      ...manifest.optionalDependencies,
      ...manifest.peerDependencies,
    })) {
      if (byName.has(dep) && !byName.get(dep).private) await visit(dep);
    }
  };
  for (const name of roots) await visit(name);
  const packDir = join(temp, 'tarballs');
  await (await import('node:fs/promises')).mkdir(packDir, { recursive: true });
  const tarballs = new Map();
  for (const name of selected) {
    const entry = byName.get(name);
    const output = await new Promise((resolvePack, reject) => {
      const child = spawn('pnpm', ['pack', '--json', '--pack-destination', packDir], {
        cwd: entry.path,
        stdio: ['ignore', 'pipe', 'inherit'],
      });
      let stdout = '';
      child.stdout.setEncoding('utf8').on('data', (chunk) => (stdout += chunk));
      child.once('error', reject);
      child.once('exit', (code) =>
        code === 0 ? resolvePack(stdout) : reject(new Error(`pnpm pack failed for ${name} (${code})`)),
      );
    });
    const packed = JSON.parse(output);
    const filename = typeof packed === 'string' ? packed : (packed.filename ?? packed[0]?.filename);
    if (!filename) fail(`pnpm pack did not report the tarball filename for ${name}: ${output}`);
    const tarball = isAbsolute(filename) ? filename : resolve(packDir, filename);
    const packedManifest = await new Promise((resolveTar, reject) => {
      const child = spawn('tar', ['-xOf', tarball, 'package/package.json'], { stdio: ['ignore', 'pipe', 'inherit'] });
      let stdout = '';
      child.stdout.setEncoding('utf8').on('data', (chunk) => (stdout += chunk));
      child.once('error', reject);
      child.once('exit', (code) =>
        code === 0 ? resolveTar(stdout) : reject(new Error(`Could not inspect packed manifest for ${name}`)),
      );
    });
    const consumerManifest = JSON.parse(packedManifest);
    for (const dep of Object.keys({
      ...consumerManifest.dependencies,
      ...consumerManifest.optionalDependencies,
      ...consumerManifest.peerDependencies,
    })) {
      if (selected.has(dep)) {
        const internal = JSON.parse(await readFile(join(byName.get(dep).path, 'package.json'), 'utf8'));
        if (consumerManifest.dependencies?.[dep]?.startsWith('workspace:'))
          fail(`${name} retains workspace protocol for internal dependency ${dep}.`);
        const spec = consumerManifest.dependencies?.[dep];
        if (spec && spec.replace(/^[~^]/, '') !== internal.version)
          fail(`${name} pins ${dep} to ${spec}, but packed ${dep} is ${internal.version}.`);
      }
    }
    tarballs.set(name, tarball);
  }
  return { byName, selected, tarballs };
}

async function consumer(roots, callback, options = {}) {
  const packages = await inventory();
  const packed = await packClosure(packages, roots);
  const dir = join(temp, options.name ?? 'consumer');
  await (await import('node:fs/promises')).mkdir(dir, { recursive: true });
  const internalSpecs = Object.fromEntries(
    [...packed.selected].map((name) => [name, `file:${relative(dir, packed.tarballs.get(name))}`]),
  );
  const dependencies = { ...internalSpecs, ...options.extraDependencies };
  const manifest = {
    name: 'molstar-smoke-consumer',
    version: '1.0.0',
    private: true,
    type: 'module',
    dependencies,
    overrides: internalSpecs,
  };
  await writeFile(join(dir, 'package.json'), JSON.stringify(manifest, null, 2));
  await run(
    'npm',
    [
      'install',
      '--no-audit',
      '--no-fund',
      '--fetch-retries=0',
      '--fetch-timeout=20000',
      '--cache',
      join(temp, 'npm-cache'),
    ],
    { cwd: dir },
  );
  return callback({ dir, packed });
}

async function nodeCheck() {
  await consumer(
    ['@molstar/core', '@molstar/io', '@molstar/plugin', '@molstar/mvs-builder', '@molstar/mp4-export-extension'],
    async ({ dir }) => {
      const require = createRequire(join(dir, 'package.json'));
      for (const native of ['gl', 'canvas']) {
        assert.throws(
          () => require.resolve(native),
          { code: 'MODULE_NOT_FOUND' },
          `Default consumers must not install ${native}.`,
        );
      }
      const fixture = join(dir, 'runtime.mjs');
      await cp(join(here, 'node/runtime.mjs'), fixture);
      await cp(join(here, 'fixtures/tiny.pdb'), join(dir, 'tiny.pdb'));
      await cp(join(here, 'fixtures/tiny.mvsj'), join(dir, 'tiny.mvsj'));
      await run(process.execPath, [fixture], { cwd: dir });
      const mp4Fixture = join(dir, 'mp4.mjs');
      await cp(join(here, 'node/mp4.mjs'), mp4Fixture);
      await run(process.execPath, [mp4Fixture], { cwd: dir });
    },
    { name: 'node-consumer', extraDependencies: { fflate: await catalogVersion('fflate') } },
  );
}

async function typesCheck() {
  await consumer(
    [
      '@molstar/core',
      '@molstar/io',
      '@molstar/plugin',
      '@molstar/plugin-ui',
      '@molstar/plugin-headless',
      '@molstar/mvs-builder',
    ],
    async ({ dir }) => {
      const ts = resolve(root, 'node_modules/typescript/bin/tsc');
      if (!(await exists(ts)))
        fail('TypeScript 7 is not installed in the repository; cannot run the packed declarations consumer.');
      const source = join(dir, 'consumer.tsx');
      await cp(join(here, 'types/consumer.tsx'), source);
      await run(
        process.execPath,
        [
          ts,
          '--noEmit',
          '--strict',
          '--module',
          'NodeNext',
          '--moduleResolution',
          'NodeNext',
          '--target',
          'ES2022',
          '--jsx',
          'react-jsx',
          '--types',
          'node,webxr',
          '--typeRoots',
          join(dir, 'node_modules/@types'),
          source,
        ],
        { cwd: dir },
      );
    },
    {
      name: 'types-consumer',
      extraDependencies: {
        '@types/node': await catalogVersion('@types/node'),
        '@types/webxr': await catalogVersion('@types/webxr'),
        react: await catalogVersion('react'),
        'react-dom': await catalogVersion('react-dom'),
      },
    },
  );
}

async function sourceCheck() {
  const esbuildPath = resolve(root, 'node_modules/esbuild/lib/main.js');
  if (!(await exists(esbuildPath))) fail('esbuild is not installed; source-condition smoke cannot run.');
  const esbuildModule = await import(pathToFileURL(esbuildPath));
  const esbuild = esbuildModule.default ?? esbuildModule;
  const output = join(temp, 'source-bundle.mjs');
  const result = await esbuild.build({
    entryPoints: [join(here, 'source/entry.ts')],
    bundle: true,
    format: 'esm',
    outfile: output,
    conditions: ['molstar-src'],
    metafile: true,
  });
  const inputs = Object.keys(result.metafile.inputs);
  assert(
    inputs.some((path) => path.includes('/src/')),
    `Source condition resolved no package source inputs: ${inputs.join(', ')}`,
  );
  assert(
    inputs.every((path) => !path.includes('/lib/')),
    `Source condition unexpectedly consumed built library output: ${inputs.filter((path) => path.includes('/lib/')).join(', ')}`,
  );
  assert(await exists(output), 'esbuild did not produce the source-condition bundle');
  console.log('Source-condition bundle passed');
}

async function slimCheck() {
  const esbuildPath = resolve(root, 'node_modules/esbuild/lib/main.js');
  if (!(await exists(esbuildPath))) fail('esbuild is not installed; slim smoke cannot run.');
  const esbuildModule = await import(pathToFileURL(esbuildPath));
  const esbuild = esbuildModule.default ?? esbuildModule;
  const started = Date.now();
  const timings = {};
  const diagnostics = resolve(root, 'build/smoke-diagnostics/slim');
  await rm(diagnostics, { recursive: true, force: true });
  await consumer(
    // the packages the slim example imports (@molstar/plugin, @molstar/core, @molstar/plugin-ui) and @molstar/model for
    // `Script`; their internal dependency closure is installed from the packed tarballs
    ['@molstar/plugin', '@molstar/core', '@molstar/plugin-ui', '@molstar/model'],
    async ({ dir }) => {
      timings.install = Date.now() - started;
      // The example sources are copied unchanged; the harness (smoke/slim/harness.js) imports the entry, so the bundle
      // resolves every package from the consumer's node_modules.
      const exampleSource = join(root, 'examples/slim-plugin/src');
      await cp(exampleSource, join(dir, 'slim-example'), { recursive: true });
      await cp(join(here, 'slim/harness.js'), join(dir, 'harness.js'));
      const site = join(dir, 'site');
      await mkdir(site, { recursive: true });
      await cp(join(here, 'slim/index.html'), join(site, 'index.html'));
      await cp(join(exampleSource, 'ligand.sdf'), join(site, 'ligand.sdf'));
      const buildStart = Date.now();
      const result = await esbuild.build({
        absWorkingDir: dir,
        entryPoints: [join(dir, 'harness.js')],
        bundle: true,
        platform: 'browser',
        format: 'iife',
        outfile: join(site, 'slim.js'),
        conditions: ['molstar-src', 'import', 'default'],
        external: ['crypto', 'fs', 'path', 'stream'],
        // The example imports its page, ligand, and skin; the page and the ligand are served from `site/` and the
        // skin is not needed to check the plugin state.
        loader: { '.html': 'empty', '.sdf': 'empty', '.scss': 'empty', '.css': 'empty' },
        tsconfigRaw: {
          compilerOptions: {
            target: 'ES2022',
            jsx: 'react-jsx',
            useDefineForClassFields: false,
            verbatimModuleSyntax: true,
          },
        },
        plugins: [
          {
            // sources import `./x.js` for `./x.ts`, as in scripts/esbuild/app.mjs
            name: 'molstar-source-js-extension',
            setup(build) {
              build.onResolve({ filter: /^\./ }, async (args) => {
                if (!args.path.endsWith('.js')) return null;
                const file = resolve(args.resolveDir, args.path);
                if (await exists(file)) return null;
                for (const ext of ['.ts', '.tsx']) {
                  const source = file.slice(0, -'.js'.length) + ext;
                  if (await exists(source)) return { path: source };
                }
                return null;
              });
            },
          },
        ],
        define: {
          'process.env.NODE_ENV': JSON.stringify('production'),
          'process.env.DEBUG': JSON.stringify(false),
          __MOLSTAR_PLUGIN_VERSION__: JSON.stringify('smoke'),
        },
        minify: true,
        minifyIdentifiers: false,
        metafile: true,
        logLevel: 'error',
      });
      timings.bundle = Date.now() - buildStart;
      const inputs = Object.keys(result.metafile.inputs).map((path) => path.split('\\').join('/'));
      assert(
        inputs.some((path) => /node_modules\/@molstar\/plugin\/src\//.test(path)),
        `The slim bundle resolved no installed package source: ${inputs.slice(0, 10).join(', ')}`,
      );
      const outside = inputs.filter((path) => path.startsWith('..'));
      assert.deepEqual(
        outside,
        [],
        `The slim bundle consumed files outside the consumer project: ${outside.join(', ')}`,
      );
      const excluded = inputs.filter((path) =>
        /\/representation\/cartoon\.ts$|\/formats\/structure\/mmcif\.ts$|\/reader\/ccp4\/|\/state\/formats\/volume\/ccp4\.ts$|\/transpilers\/(pymol|vmd|jmol)\//.test(
          path,
        ),
      );
      assert.deepEqual(excluded, [], `The slim bundle contains excluded modules: ${excluded.join(', ')}`);

      const types = {
        '.html': 'text/html',
        '.js': 'text/javascript',
        '.css': 'text/css',
        '.sdf': 'chemical/x-mdl-sdfile',
      };
      const server = createServer(async (req, res) => {
        try {
          let pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
          if (pathname === '/') pathname = '/index.html';
          if (pathname === '/favicon.ico') {
            res.writeHead(204);
            res.end();
            return;
          }
          const file = resolve(site, `.${pathname}`);
          if (!file.startsWith(site + sep)) {
            res.writeHead(403);
            res.end('forbidden');
            return;
          }
          const body = await readFile(file);
          res.writeHead(200, { 'content-type': types[extname(file)] ?? 'application/octet-stream' });
          res.end(body);
        } catch (error) {
          res.writeHead(error.code === 'ENOENT' ? 404 : 500);
          res.end(error.code === 'ENOENT' ? 'not found' : error.message);
        }
      });
      await new Promise((resolveListen, rejectListen) => {
        server.once('error', rejectListen);
        server.listen(0, '127.0.0.1', () => {
          server.removeListener('error', rejectListen);
          resolveListen();
        });
      });
      const { port } = server.address();
      const consoleLog = [];
      const browserErrors = [];
      let browser;
      let page;
      try {
        browser = await launchBrowser();
        page = await browser.newPage({ viewport: { width: 800, height: 600 } });
        page.on('pageerror', (error) => browserErrors.push(`pageerror: ${error.stack ?? error.message}`));
        page.on('requestfailed', (request) =>
          browserErrors.push(`requestfailed: ${request.url()} ${request.failure()?.errorText}`),
        );
        page.on('response', (response) => {
          if (response.status() >= 400) browserErrors.push(`HTTP ${response.status()} ${response.url()}`);
        });
        page.on('console', (message) => {
          consoleLog.push(`[${message.type()}] ${message.text()}`);
          if (message.type() === 'error') browserErrors.push(`console.error: ${message.text()}`);
        });
        const pageStart = Date.now();
        await page.goto(`http://127.0.0.1:${port}/`, { waitUntil: 'load', timeout: 30000 });
        await page.waitForFunction(() => Boolean(window.slimPluginReady), undefined, { timeout: 30000 });
        await page.evaluate(() => window.slimPluginReady.then(() => undefined));
        // The representation reaches the scene on the commit after the state update resolves.
        await page.waitForFunction(() => window.slimSmoke.describe().reprCount > 0, undefined, { timeout: 30000 });
        timings.load = Date.now() - pageStart;

        // 1. The SDF ligand renders as ball-and-stick.
        const loaded = await page.evaluate(() => window.slimSmoke.describe());
        assert.equal(loaded.structures, 1, `Expected one structure: ${JSON.stringify(loaded)}`);
        const ballAndStick = loaded.structureReprs.filter((r) => r.type === 'ball-and-stick');
        assert(ballAndStick.length > 0, `No ball-and-stick representation: ${JSON.stringify(loaded.structureReprs)}`);
        assert.equal(
          loaded.structureReprs.length,
          ballAndStick.length,
          `Unexpected representation types: ${JSON.stringify(loaded.structureReprs)}`,
        );
        for (const repr of ballAndStick) {
          assert.equal(repr.status, 'ok', `Representation ${repr.ref} is not ok: ${JSON.stringify(repr)}`);
          assert(repr.visible, `Representation ${repr.ref} is not visible`);
          assert(repr.renderObjectsWithGeometry > 0, `Representation ${repr.ref} has no render objects with geometry`);
        }
        assert(loaded.reprCount > 0, `canvas3d has no representations: ${JSON.stringify(loaded)}`);
        assert(loaded.sceneRenderObjects > 0, `canvas3d has no render objects: ${JSON.stringify(loaded)}`);
        assert((await page.locator('canvas').count()) > 0, 'The slim plugin did not create a render canvas');

        // 2. Only the slim set of providers is registered.
        assert.deepEqual(loaded.registered, {
          structureRepresentations: ['ball-and-stick'],
          formats: ['sdf'],
          hierarchyPresets: ['preset-trajectory-default'],
          representationPresets: ['preset-structure-representation-ball-and-stick'],
          scriptLanguages: ['mol-script'],
        });

        // 3. A snapshot naming `cartoon` reports it as unregistered and renders the registry default.
        const restored = await page.evaluate(() => window.slimSmoke.restoreCartoonSnapshot());
        assert(restored.changed > 0, 'The snapshot has no structure representation to rewrite');
        const afterRestore = await page.evaluate(() => window.slimSmoke.describe());
        const warnings = await page.evaluate(() => window.slimSmoke.warnings);
        assert(
          warnings.includes(
            "Structure representation 'cartoon' is not registered in this plugin; the registry default is used",
          ),
          `Missing the unregistered-cartoon warning; got ${JSON.stringify(warnings)}`,
        );
        assert.equal(afterRestore.structureReprs.length, restored.changed, JSON.stringify(afterRestore.structureReprs));
        for (const repr of afterRestore.structureReprs) {
          assert.equal(repr.type, 'ball-and-stick', `Restored representation is ${repr.type}: ${JSON.stringify(repr)}`);
          assert.equal(repr.status, 'ok', `Restored representation ${repr.ref} is not ok: ${JSON.stringify(repr)}`);
          assert(repr.renderObjectsWithGeometry > 0, `Restored representation ${repr.ref} has no render objects`);
        }
        assert(afterRestore.sceneRenderObjects > 0, 'canvas3d has no render objects after restoring the snapshot');

        // 4. A PyMOL script fails; the language is not enabled in the slim build.
        assert.equal(
          await page.evaluate(() => window.slimSmoke.evaluatePyMol()),
          "Script language 'pymol' is not available in this build",
        );

        assert.deepEqual(browserErrors, [], `Browser runtime or asset errors: ${browserErrors.join('; ')}`);
      } catch (error) {
        await mkdir(diagnostics, { recursive: true });
        await writeFile(join(diagnostics, 'console.log'), `${consoleLog.join('\n')}\n\n${browserErrors.join('\n')}\n`);
        await page?.screenshot({ path: join(diagnostics, 'screenshot.png') }).catch(() => undefined);
        console.error(`Slim smoke diagnostics saved to ${relative(root, diagnostics)}`);
        throw error;
      } finally {
        await browser?.close();
        await new Promise((resolveClose) => server.close(resolveClose));
      }
    },
    {
      name: 'slim-consumer',
      extraDependencies: { react: await catalogVersion('react'), 'react-dom': await catalogVersion('react-dom') },
    },
  );
  console.log(
    `Slim plugin renders the SDF ligand with only the slim providers (install ${timings.install} ms, bundle ${timings.bundle} ms, page ${timings.load} ms, total ${Date.now() - started} ms)`,
  );
}

/** Launches headless Chromium through the workspace Playwright (`MOLSTAR_SMOKE_BROWSER` overrides the executable). */
async function launchBrowser() {
  const playwrightPath = resolve(root, 'node_modules/playwright/index.mjs');
  if (!(await exists(playwrightPath)))
    fail(
      'Browser automation unavailable: install Playwright in the workspace to execute browser smoke pages. Pages and local fixture are ready, but browser behavior is unverified.',
    );
  const { chromium } = await import(pathToFileURL(playwrightPath));
  const browserPath =
    process.env.MOLSTAR_SMOKE_BROWSER ||
    ((await exists('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'))
      ? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
      : undefined);
  return chromium.launch({
    headless: true,
    ...(browserPath ? { executablePath: browserPath } : {}),
    args: ['--use-angle=swiftshader', '--enable-webgl', '--ignore-gpu-blocklist'],
  });
}

async function browserCheck() {
  const packedDistribution = await unpackDistribution();
  const dist = join(packedDistribution, 'build/esm');
  for (const artifact of ['viewer.js', 'viewer.css', 'import-map.json'])
    if (!(await exists(join(dist, artifact))))
      fail(`Browser smoke prerequisite missing: ${relative(root, join(dist, artifact))}`);
  for (const artifact of [
    'viewer/molstar.js',
    'viewer/molstar.css',
    'mvs-stories/mvs-stories.js',
    'mvs-stories/mvs-stories.css',
  ]) {
    if (!(await exists(join(packedDistribution, 'build', artifact))))
      fail(`Classic browser smoke prerequisite missing in the packed tarball: build/${artifact}`);
  }
  const importMap = await readFile(join(dist, 'import-map.json'), 'utf8');
  const server = createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      let file;
      if (pathname === '/favicon.ico') file = join(packedDistribution, 'build/viewer/favicon.ico');
      else if (pathname === '/fixtures/tiny.pdb') file = join(here, 'fixtures/tiny.pdb');
      else if (pathname === '/fixtures/tiny.mvsj') file = join(here, 'fixtures/tiny.mvsj');
      else if (pathname === '/viewer/') file = join(here, 'browser/viewer/index.html');
      else if (pathname === '/library/') file = join(here, 'browser/library/index.html');
      else if (pathname.startsWith('/build/esm/')) file = join(packedDistribution, pathname.slice(1));
      else if (pathname === '/classic/viewer.html') file = join(here, 'browser/classic/viewer.html');
      else if (pathname === '/classic/mvs-stories.html') file = join(here, 'browser/classic/mvs-stories.html');
      else if (pathname.startsWith('/classic/viewer/'))
        file = join(packedDistribution, 'build/viewer', pathname.slice('/classic/viewer/'.length));
      else if (pathname.startsWith('/classic/mvs-stories/'))
        file = join(packedDistribution, 'build/mvs-stories', pathname.slice('/classic/mvs-stories/'.length));
      else {
        res.writeHead(404);
        res.end('not found');
        return;
      }
      let body = await readFile(file);
      if (pathname === '/library/' || pathname === '/viewer/')
        body = body
          .toString('utf8')
          .replace(
            '<script type="importmap" id="smoke-importmap"></script>',
            `<script type="importmap">${importMap}</script>`,
          );
      const types = {
        '.html': 'text/html',
        '.js': 'text/javascript',
        '.mjs': 'text/javascript',
        '.css': 'text/css',
        '.json': 'application/json',
        '.pdb': 'chemical/x-pdb',
      };
      res.writeHead(200, { 'content-type': types[extname(file)] ?? 'application/octet-stream' });
      res.end(body);
    } catch (error) {
      res.writeHead(500);
      res.end(error.message);
    }
  });
  const serve = args.has('--serve');
  const portOption = [...args].find((arg) => arg.startsWith('--port='));
  const port = serve ? Number(portOption?.slice('--port='.length) ?? 1339) : 0;
  if (!Number.isInteger(port) || port < 0 || port > 65535) throw new Error('Use --port=<0–65535>.');
  await new Promise((resolveListen, rejectListen) => {
    server.once('error', rejectListen);
    server.listen(port, '127.0.0.1', () => {
      server.removeListener('error', rejectListen);
      resolveListen();
    });
  });
  const address = server.address();
  if (serve) {
    console.log('Packed browser smoke pages (Ctrl+C to stop):');
    for (const page of ['/viewer/', '/library/', '/classic/viewer.html', '/classic/mvs-stories.html']) {
      console.log(`  http://127.0.0.1:${address.port}${page}`);
    }
    try {
      await new Promise((resolveStop) => {
        const stop = () => {
          process.removeListener('SIGINT', stop);
          process.removeListener('SIGTERM', stop);
          resolveStop();
        };
        process.once('SIGINT', stop);
        process.once('SIGTERM', stop);
      });
    } finally {
      await new Promise((resolveClose) => server.close(resolveClose));
    }
    return;
  }
  const browserErrors = [];
  let browser;
  try {
    browser = await launchBrowser();
    for (const pagePath of ['/viewer/', '/library/', '/classic/viewer.html', '/classic/mvs-stories.html']) {
      const page = await browser.newPage();
      page.on('pageerror', (error) => browserErrors.push(`${pagePath}: ${error.stack ?? error.message}`));
      page.on('requestfailed', (request) =>
        browserErrors.push(`${pagePath}: ${request.url()} ${request.failure()?.errorText}`),
      );
      page.on('response', (response) => {
        if (response.status() >= 400) browserErrors.push(`${pagePath}: HTTP ${response.status()} ${response.url()}`);
      });
      page.on('console', (message) => {
        if (message.type() === 'error') browserErrors.push(`${pagePath}: ${message.text()}`);
      });
      await page.goto(`http://127.0.0.1:${address.port}${pagePath}`, { waitUntil: 'networkidle', timeout: 30000 });
      await page.waitForFunction(() => Boolean(window.smokeReady), undefined, { timeout: 45000 });
      try {
        await page.evaluate(() => window.smokeReady);
      } catch (error) {
        throw new Error(`${error.message}\nBrowser errors: ${browserErrors.join('; ')}`);
      }
      assert.equal(
        await page.locator('html').getAttribute('data-smoke'),
        'pass',
        `${pagePath} fixture assertions failed: ${JSON.stringify(await page.evaluate(() => window.smokeResult))}`,
      );
      if (pagePath === '/viewer/') {
        assert((await page.locator('canvas').count()) > 0, 'Viewer did not create a render canvas');
      }
      if (pagePath === '/viewer/' || pagePath === '/library/' || pagePath.startsWith('/classic/')) {
        assert(
          await page.evaluate(() => (window.smokeResult?.representations ?? 0) > 0),
          `${pagePath} loaded data but created no structure representations`,
        );
      }
      await page.close();
    }
    assert.deepEqual(browserErrors, [], `Browser runtime or asset errors: ${browserErrors.join('; ')}`);
  } finally {
    await browser?.close();
    await new Promise((resolveClose) => server.close(resolveClose));
  }
  console.log('Browser ESM and classic compatibility pages passed');
}

async function cliCheck() {
  await consumer(
    ['@molstar/mvs-builder', '@molstar/cifschema-cli'],
    async ({ dir }) => {
      const bin = join(dir, 'node_modules/.bin', process.platform === 'win32' ? 'mvs-validate.cmd' : 'mvs-validate');
      if (!(await exists(bin))) fail('Packed @molstar/mvs-builder does not expose mvs-validate.');
      const fixture = join(dir, 'tiny.mvsj');
      await cp(join(here, 'fixtures/tiny.mvsj'), fixture);
      const result = await capture(bin, [fixture], { cwd: dir });
      assert.match(result.stdout.trim(), /^OK\s+.*tiny\.mvsj$/, `Unexpected mvs-validate output: ${result.stdout}`);
      console.log('Packed mvs-validate accepted the local MVS fixture');
      const schemaBin = join(dir, 'node_modules/.bin', process.platform === 'win32' ? 'cifschema.cmd' : 'cifschema');
      const dictionary = join(dir, 'tiny.dic');
      const schema = join(dir, 'schema.ts');
      await cp(join(here, 'fixtures/tiny.dic'), dictionary);
      await capture(schemaBin, ['--preset', 'mmCIF', '--dicPath', dictionary, '--out', schema], { cwd: dir });
      assert.match(await readFile(schema, 'utf8'), /atom_site/);
      for (const name of await readdir(join(root, 'data/cif-field-names'))) {
        assert.deepEqual(
          await readFile(join(dir, 'node_modules/@molstar/cifschema-cli/lib/data/cif-field-names', name)),
          await readFile(join(root, 'data/cif-field-names', name)),
        );
      }
      assert(
        !(await exists(join(dir, 'node_modules/@molstar/cifschema-cli/data'))),
        'Packed CLI must use staged build assets.',
      );
      console.log('Packed cifschema generated a schema using filters staged from canonical root data');
    },
    { name: 'cli-consumer' },
  );
}

async function headlessCheck() {
  await consumer(
    ['@molstar/plugin-headless'],
    async ({ dir, packed }) => {
      const packageNames = [...packed.selected];
      assert(
        !packageNames.includes('@molstar/mp4-export-extension'),
        'Packed headless dependency closure unexpectedly includes MP4 export.',
      );
      for (const name of packageNames) {
        const manifest = JSON.parse(await readFile(join(packed.byName.get(name).path, 'package.json'), 'utf8'));
        const deps = { ...manifest.dependencies, ...manifest.optionalDependencies, ...manifest.peerDependencies };
        if (Object.hasOwn(deps, 'h264-mp4-encoder')) fail(`${name} unexpectedly depends on h264-mp4-encoder.`);
      }
      const resolver = createRequire(join(dir, 'headless-consumer.mjs'));
      const nativePath = (name, envName) => {
        if (process.env[envName]) return resolve(process.env[envName]);
        try {
          return resolver.resolve(name);
        } catch (consumerError) {
          try {
            return rootRequire.resolve(name);
          } catch (rootError) {
            if (rootError.code === 'MODULE_NOT_FOUND' && rootError.message.includes(`'${name}'`)) {
              fail(
                `Headless smoke unavailable: optional native peer '${name}' is not installed. Install it in the consumer or set ${envName} to its package entry path.`,
              );
            }
            throw consumerError;
          }
        }
      };
      const env = {
        ...process.env,
        MOLSTAR_SMOKE_GL: nativePath('gl', 'MOLSTAR_SMOKE_GL'),
        MOLSTAR_SMOKE_PNGJS: nativePath('pngjs', 'MOLSTAR_SMOKE_PNGJS'),
        MOLSTAR_SMOKE_INTERNAL_PACKAGES: JSON.stringify(packageNames),
      };
      const fixture = join(dir, 'capture.mjs');
      await cp(join(here, 'headless/capture.mjs'), fixture);
      await cp(join(here, 'fixtures/tiny.pdb'), join(dir, 'tiny.pdb'));
      await run(process.execPath, [fixture], { cwd: dir, env });
    },
    { name: 'headless-consumer' },
  );
}

if (args.has('--prepare')) await prepare();
const checks =
  mode === 'all'
    ? [
        ['node', nodeCheck],
        ['types', typesCheck],
        ['source', sourceCheck],
        ['slim', slimCheck],
        ['cli', cliCheck],
        ['browser', browserCheck],
      ]
    : [
        [
          mode,
          {
            node: nodeCheck,
            types: typesCheck,
            source: sourceCheck,
            slim: slimCheck,
            cli: cliCheck,
            browser: browserCheck,
            headless: headlessCheck,
          }[mode],
        ],
      ];
const errors = [];
for (const [name, check] of checks) {
  try {
    console.log(`\n== smoke:${name} ==`);
    await check();
  } catch (error) {
    errors.push(`${name}: ${error.message}`);
    console.error(`FAILED ${name}: ${error.message}`);
  }
}
if (errors.length) {
  console.error(`\n${errors.length} smoke check(s) failed:\n${errors.map((e) => `- ${e}`).join('\n')}`);
  process.exitCode = 1;
} else if (!args.has('--serve')) console.log('\nAll requested smoke checks passed.');
await cleanup();
