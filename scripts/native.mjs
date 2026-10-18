import fs from 'node:fs/promises';
import path from 'node:path';
import { delimiter } from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { spawn } from 'node:child_process';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const nativeDir = path.join(root, '.cache/native');
const manifestPath = path.join(nativeDir, 'package.json');
const require = createRequire(manifestPath);
const [action, ...rawArgs] = process.argv.slice(2);
const args = rawArgs[0] === '--' ? rawArgs.slice(1) : rawArgs;
const env = {
  ...process.env,
  NODE_PATH: [path.join(nativeDir, 'node_modules'), process.env.NODE_PATH].filter(Boolean).join(delimiter),
};
const command = (name) => (process.platform === 'win32' ? `${name}.cmd` : name);
// Node refuses to spawn .cmd/.bat files without a shell on Windows (EINVAL); the shell joins the arguments unquoted
const isBatch = (name) => process.platform === 'win32' && /\.(cmd|bat)$/i.test(name);
const quote = (arg) => (/[\s"&|<>^]/.test(arg) ? `"${arg.replaceAll('"', '""')}"` : arg);
const run = (name, commandArgs, cwd = root) =>
  new Promise((resolve, reject) => {
    const shell = isBatch(name);
    const child = shell
      ? spawn(quote(name), commandArgs.map(quote), { cwd, env, stdio: 'inherit', shell })
      : spawn(name, commandArgs, { cwd, env, stdio: 'inherit' });
    child.once('error', reject);
    child.once('exit', (code, signal) =>
      code === 0 ? resolve() : reject(new Error(`${name} failed (${signal ?? code}).`)),
    );
  });

const isInstalled = async (name) => {
  try {
    await fs.access(path.join(nativeDir, 'node_modules', name, 'package.json'));
    return true;
  } catch {
    return false;
  }
};

async function requireSetup() {
  if ((await isInstalled('gl')) || (await isInstalled('webgpu'))) return;
  throw new Error(
    'Native rendering dependencies are not installed. Run pnpm native:install first (add -- --canvas for the rendering CLI, --webgpu for the headless WebGPU backend).',
  );
}

async function main() {
  const help = args[0] === '--help';
  if (help || !['install', 'test', 'run'].includes(action)) {
    console.log(
      'pnpm native:install [-- --canvas] [--webgpu]  Install gl (and optionally canvas, webgpu) in .cache/native',
    );
    console.log('pnpm test:native                  Run the unit tests (Vitest) with the installed native GL backend');
    console.log('pnpm native:run -- <command>      Run a command with the native modules available');
    if (!help) process.exitCode = 2;
    return;
  }
  if (action === 'install') {
    if (args.some((arg) => arg !== '--canvas' && arg !== '--webgpu'))
      throw new Error('native:install accepts only --canvas, --webgpu or --help.');
    let previous;
    try {
      const manifest = JSON.parse(await fs.readFile(manifestPath, 'utf8'));
      previous = { dependencies: { ...manifest.dependencies, ...manifest.optionalDependencies } };
    } catch (error) {
      if (error.code !== 'ENOENT') throw error;
    }
    const dependencies = { gl: '8.1.6' };
    if (args.includes('--canvas') || previous?.dependencies?.canvas) dependencies.canvas = '2.11.2';
    // Dawn's Node.js bindings (prebuilt), for the headless WebGPU backend
    if (args.includes('--webgpu') || previous?.dependencies?.webgpu) dependencies.webgpu = '0.6.2';
    await fs.mkdir(nativeDir, { recursive: true });
    // optional, so that one module failing to build (e.g. gl without prebuilt binaries for the Node version) doesn't
    // keep the others from installing
    await fs.writeFile(
      manifestPath,
      JSON.stringify({ name: 'molstar-native-tools', private: true, optionalDependencies: dependencies }, null, 2) +
        '\n',
    );
    console.log(
      `Installing opt-in native modules in ${path.relative(root, nativeDir)}; workspace manifests and lockfile stay untouched.`,
    );
    await run(command('npm'), ['install', '--prefix', nativeDir, '--no-audit', '--no-fund'], nativeDir);
    for (const name of Object.keys(dependencies)) {
      if (!(await isInstalled(name))) console.warn(`${name} could not be installed (see the npm output above).`);
    }
    return;
  }
  await requireSetup();
  if (action === 'test') {
    if (args.length) throw new Error('Use pnpm native:run -- pnpm test:unit <args> for focused native tests.');
    if (await isInstalled('gl')) {
      const gl = require('gl');
      const context = gl(1, 1);
      if (!context)
        throw new Error('Native gl is installed but could not create a WebGL context. On Linux, run under Xvfb.');
      context.getExtension('STACKGL_destroy_context')?.destroy();
    }
    await run(command('pnpm'), ['test:unit']);
  } else {
    if (!args.length)
      throw new Error(
        'native:run requires a command, for example: pnpm native:run -- node cli/mvs-render/lib/mvs-render.js --help',
      );
    const [name, ...commandArgs] = args;
    await run(['npm', 'pnpm'].includes(name) ? command(name) : name, commandArgs);
  }
}

try {
  await main();
} catch (error) {
  console.error(error instanceof Error ? error.message : String(error));
  process.exitCode = 1;
}
