/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import simpleGit from 'simple-git';
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
    stories: path.join(root, 'distributions/molstar/build/mvs-stories'),
    mesoscale: path.join(root, 'apps/mesoscale-explorer/build'),
    example: name => path.join(root, 'examples', name, 'build'),
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
const addAnalytics = filePath => replaceFile(filePath, analyticsTag, analyticsCode);
const addManifest = filePath => replaceFile(filePath, manifestTag, manifestCode);
const addPwa = filePath => replaceFile(filePath, pwaTag, pwaCode);
const addVersion = filePath => replaceFile(filePath, '__MOLSTAR_VERSION__', VERSION);

function copyViewer() {
    console.log('\n### copy viewer files');
    const viewerDeployPath = path.join(localPath, 'viewer');
    fsExtra.copySync(buildPaths.viewer, viewerDeployPath, { overwrite: true });
    addAnalytics(path.join(viewerDeployPath, 'index.html'));
    addManifest(path.join(viewerDeployPath, 'index.html'));
    addPwa(path.join(viewerDeployPath, 'index.html'));

    fsExtra.copySync(path.join(dataDir, 'pwa'), viewerDeployPath, { overwrite: true });
    addVersion(path.join(viewerDeployPath, 'sw.js'));
}
function copyMesoscale() {
    console.log('\n### copy mesoscale explorer files');
    const target = path.join(localPath, 'me/viewer');
    fsExtra.copySync(buildPaths.mesoscale, target, { overwrite: true });
    addAnalytics(path.join(target, 'index.html'));
}
function copyMVSStories() {
    console.log('\n### copy MVS stories files');
    const target = path.join(localPath, `stories-viewer/v${MVS_STORIES_VERSION}`);
    fsExtra.copySync(buildPaths.stories, target, { overwrite: true });
    addAnalytics(path.join(target, 'index.html'));
}
function copyDemo(name) {
    console.log('\n### copy demo files for', name);
    const target = path.join(localPath, 'demos', name);
    fsExtra.copySync(buildPaths.example(name), target, { overwrite: true });
    addAnalytics(path.join(target, 'index.html'));
}
function copyFiles() {
    console.log('\n### copy apps and demos');
    copyViewer();
    copyMesoscale();
    copyMVSStories();
    copyDemo('lighting');
    copyDemo('alpha-orbitals');
    copyDemo('mvs-stories');
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
    await repoGit.commit(`Updated Apps and Demos\n- Mol* version: ${VERSION}\n- MVS Stories version: ${MVS_STORIES_VERSION}`);
    await repoGit.push();
}

const args = new Set(process.argv.slice(2));
if (args.has('--help') || args.has('-h')) {
    console.log('Usage: node scripts/deploy.js [--local]');
    process.exit(0);
}
const unknown = [...args].filter(arg => arg !== '--local');
if (unknown.length) throw new Error(`Unknown deploy option(s): ${unknown.join(', ')}`);
fs.mkdirSync(localPath, { recursive: true });
copyFiles();
if (args.has('--local')) process.exit(0);
fs.mkdirSync(repositoryPath, { recursive: true });
await syncRepository();
