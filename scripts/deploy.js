/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { simpleGit } from 'simple-git';
import fsExtra from 'fs-extra';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const VERSION = JSON.parse(fs.readFileSync(path.join(root, 'version.json'), 'utf8')).version;
const storiesVersionSource = fs.readFileSync(path.join(root, 'apps/mvs-stories/src/version.ts'), 'utf8');
const storiesVersion = storiesVersionSource.match(/export\s+const\s+VERSION\s*=\s*([^;\n]+)/u)?.[1]?.trim();
if (!storiesVersion) throw new Error('Unable to read MVS Stories VERSION from apps/mvs-stories/src/version.ts');
const MVS_STORIES_VERSION = storiesVersion.replace(/^['"]|['"]$/gu, '');

const remoteUrl = 'https://github.com/molstar/molstar.github.io.git';
const dataDir = path.join(root, 'data');
const deployDir = path.join(root, 'deploy');
const localPath = path.join(deployDir, 'data');
const repositoryPath = path.join(deployDir, 'molstar.github.io');
const buildPaths = {
  viewer: path.join(root, 'distributions/molstar/build/viewer'),
  prismViewer: path.join(root, 'distributions/molstar/build/prism-viewer'),
  stories: path.join(root, 'distributions/molstar/build/mvs-stories'),
  mesoscale: path.join(root, 'apps/mesoscale-explorer/build'),
  storybook: path.join(root, 'storybook'),
  example: (name) => path.join(root, 'examples', name, 'build'),
};

/**
 * Where each app goes in the deploy (under deploy/data, and so on molstar.org): `--<app>-prefix <folder>` moves one
 * (`--prism-prefix prism-preview`). An app's pages and chunks address each other relatively, so any folder serves it.
 */
const prefixes = {
  viewer: 'viewer',
  prism: 'prism',
  storybook: 'ui-storybook',
  mesoscale: 'me/viewer',
  stories: `stories-viewer/v${MVS_STORIES_VERSION}`,
  demos: 'demos',
};

const analyticsTag = /<!-- __MOLSTAR_ANALYTICS__ -->/g;
const analyticsCode = `<!-- Cloudflare Web Analytics --><script defer src='https://static.cloudflareinsights.com/beacon.min.js' data-cf-beacon='{"token": "c414cbae2d284ea995171a81e4a3e721"}'></script><!-- End Cloudflare Web Analytics --><iframe src="https://web3dsurvey.com/collector-iframe.html" style="width: 1px; height: 1px;"></iframe>`;
const manifestTag = /<!-- __MOLSTAR_MANIFEST__ -->/g;
const manifestCode = `<link rel="manifest" href="./manifest.webmanifest">`;
const pwaTag = /<!-- __MOLSTAR_PWA__ -->/g;
const pwaCode = `<script src='./pwa.js'></script>`;

function replaceFile(filePath, pattern, replacement) {
  const data = fs.readFileSync(filePath, 'utf8');
  fs.writeFileSync(filePath, data.replace(pattern, replacement), 'utf8');
}
const addAnalytics = (filePath) => replaceFile(filePath, analyticsTag, analyticsCode);
const addManifest = (filePath) => replaceFile(filePath, manifestTag, manifestCode);
const addPwa = (filePath) => replaceFile(filePath, pwaTag, pwaCode);
const addVersion = (filePath) => replaceFile(filePath, '__MOLSTAR_VERSION__', VERSION);

function copyViewer() {
  console.log('\n### copy viewer files');
  const viewerDeployPath = path.join(localPath, prefixes.viewer);
  fsExtra.copySync(buildPaths.viewer, viewerDeployPath, { overwrite: true });
  addAnalytics(path.join(viewerDeployPath, 'index.html'));
  addManifest(path.join(viewerDeployPath, 'index.html'));
  addPwa(path.join(viewerDeployPath, 'index.html'));

  fsExtra.copySync(path.join(dataDir, 'pwa'), viewerDeployPath, { overwrite: true });
  addVersion(path.join(viewerDeployPath, 'sw.js'));
}
/** The Prism viewer (v6), beside the classic one: its page at `prism/` (`--prism-prefix`). */
function copyPrismViewer() {
  console.log(`\n### copy Prism viewer files (${prefixes.prism}/)`);
  const target = path.join(localPath, prefixes.prism);
  fsExtra.copySync(buildPaths.prismViewer, target, { overwrite: true });
  addAnalytics(path.join(target, 'index.html'));
}

/** The UI Storybook (`storybook/`), its static build at `ui-storybook/` (`--storybook-prefix`); its paths are relative. */
function buildStorybook() {
  console.log(`\n### build UI Storybook (${prefixes.storybook}/)`);
  const target = path.join(localPath, prefixes.storybook);
  execFileSync('pnpm', ['exec', 'storybook', 'build', '--output-dir', target, '--quiet'], {
    cwd: buildPaths.storybook,
    stdio: 'inherit',
  });
}

function copyMesoscale() {
  console.log('\n### copy mesoscale explorer files');
  const target = path.join(localPath, prefixes.mesoscale);
  fsExtra.copySync(buildPaths.mesoscale, target, { overwrite: true });
  addAnalytics(path.join(target, 'index.html'));
}
function copyMVSStories() {
  console.log('\n### copy MVS stories files');
  const target = path.join(localPath, prefixes.stories);
  fsExtra.copySync(buildPaths.stories, target, { overwrite: true });
  addAnalytics(path.join(target, 'index.html'));
}
function copyDemo(name) {
  console.log('\n### copy demo files for', name);
  const target = path.join(localPath, prefixes.demos, name);
  fsExtra.copySync(buildPaths.example(name), target, { overwrite: true });
  addAnalytics(path.join(target, 'index.html'));
}
function copyFiles() {
  console.log('\n### copy apps and demos');
  copyViewer();
  copyPrismViewer();
  copyMesoscale();
  copyMVSStories();
  copyDemo('lighting');
  copyDemo('alpha-orbitals');
  copyDemo('mvs-stories');
  buildStorybook();
}
function copyToRepository() {
  console.log('\n### copy repository files');
  fsExtra.copySync(localPath, repositoryPath, { overwrite: true });
}
function log(command, stdout, stderr) {
  if (!command) return;
  console.log('\n###', command);
  stdout.pipe(process.stdout);
  stderr.pipe(process.stderr);
}

async function syncRepository() {
  const exists = fs.existsSync(path.join(repositoryPath, '.git'));
  const git = exists ? simpleGit(repositoryPath) : simpleGit();
  git.outputHandler(log);
  if (!exists) {
    console.log('\n### clone repository');
    await git.clone(remoteUrl, repositoryPath);
  }
  const repoGit = simpleGit(repositoryPath).outputHandler(log);
  await repoGit.fetch(['--all']);
  if (exists) {
    console.log('\n### update repository');
    await repoGit.reset(['--hard', 'origin/master']);
  }
  copyToRepository();
  console.log('\n### commit changes');
  await repoGit.add(['-A']);
  await repoGit.commit(
    `Updated Apps and Demos\n- Mol* version: ${VERSION}\n- MVS Stories version: ${MVS_STORIES_VERSION}`,
  );
  await repoGit.push();
}

const usage = [
  'Usage: node scripts/deploy.js [--local | --dry-run] [--<app>-prefix <folder> ...]',
  '  --local    prepare deploy/data only (serve it to check: npx http-server deploy/data)',
  `  --dry-run  prepare deploy/data and print what would be pushed to ${remoteUrl}`,
  `  (none)     prepare deploy/data, commit it into deploy/molstar.github.io and push it`,
  '  --<app>-prefix <folder>  where an app goes (relative, no ..): ' +
    Object.entries(prefixes)
      .map(([app, folder]) => `${app} (${folder})`)
      .join(', '),
  '             e.g. pnpm deploy:local -- --prism-prefix prism-preview',
].join('\n');

const flags = new Set();
const argv = process.argv.slice(2);
for (let i = 0; i < argv.length; i++) {
  const arg = argv[i];
  const prefix = /^--([a-z]+)-prefix(?:=(.*))?$/.exec(arg);
  if (prefix) {
    const app = prefix[1];
    if (!(app in prefixes)) throw new Error(`Unknown app in ${arg}: one of ${Object.keys(prefixes).join(', ')}`);
    const value = prefix[2] ?? argv[++i];
    const folder = value?.replace(/^\/+|\/+$/g, '');
    if (!folder || folder.split('/').some((part) => !part || part === '.' || part === '..')) {
      throw new Error(`${arg} needs a relative folder without . or .. (got ${JSON.stringify(value)})`);
    }
    prefixes[app] = folder;
  } else flags.add(arg);
}
if (flags.has('--help') || flags.has('-h')) {
  console.log(usage);
  process.exit(0);
}
const unknown = [...flags].filter((arg) => arg !== '--local' && arg !== '--dry-run');
if (unknown.length) throw new Error(`Unknown deploy option(s): ${unknown.join(', ')}\n${usage}`);
const args = flags;
fs.mkdirSync(localPath, { recursive: true });
copyFiles();
if (args.has('--local')) process.exit(0);
if (args.has('--dry-run')) {
  console.log(`\n### dry run: would sync ${path.relative(root, localPath)} into ${remoteUrl} (master) and push`);
  for (const entry of fs.readdirSync(localPath).sort()) console.log(`  ${entry}/`);
  process.exit(0);
}
fs.mkdirSync(repositoryPath, { recursive: true });
await syncRepository();
