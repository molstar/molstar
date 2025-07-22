/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const git = require('simple-git');
const path = require('path');
const fs = require("fs");
const fse = require("fs-extra");
const argparse = require('argparse');

const VERSION = require(path.resolve(__dirname, '../package.json')).version;
const MVS_STORIES_VERSION = require(path.resolve(__dirname, '../src/apps/mvs-stories/version.ts')).VERSION;

const remoteUrl = "https://github.com/molstar/molstar.github.io.git";
const dataDir = path.resolve(__dirname, '../data/');
const buildDir = path.resolve(__dirname, '../build/');
const deployDir = path.resolve(__dirname, '../deploy/');
const localPath = path.resolve(deployDir, 'data/');
const repositoryPath = path.resolve(deployDir, 'molstar.github.io/');

const analyticsTag = /<!-- __MOLSTAR_ANALYTICS__ -->/g;
const analyticsCode = `<!-- Cloudflare Web Analytics --><script defer src='https://static.cloudflareinsights.com/beacon.min.js' data-cf-beacon='{"token": "c414cbae2d284ea995171a81e4a3e721"}'></script><!-- End Cloudflare Web Analytics --><iframe src="https://web3dsurvey.com/collector-iframe.html" style="width: 1px; height: 1px;"></iframe>`;

const manifestTag = /<!-- __MOLSTAR_MANIFEST__ -->/g;
const manifestCode = `<link rel="manifest" href="./manifest.webmanifest">`;

const pwaTag = /<!-- __MOLSTAR_PWA__ -->/g;
const pwaCode = `<script src='./pwa.js'></script>`;

function log(command, stdout, stderr) {
    if (command) {
        console.log('\n###', command);
        stdout.pipe(process.stdout);
        stderr.pipe(process.stderr);
    }
}

function addAnalytics(path) {
    const data = fs.readFileSync(path, 'utf8');
    const result = data.replace(analyticsTag, analyticsCode);
    fs.writeFileSync(path, result, 'utf8');
}

function addManifest(path) {
    const data = fs.readFileSync(path, 'utf8');
    const result = data.replace(manifestTag, manifestCode);
    fs.writeFileSync(path, result, 'utf8');
}

function addPwa(path) {
    const data = fs.readFileSync(path, 'utf8');
    const result = data.replace(pwaTag, pwaCode);
    fs.writeFileSync(path, result, 'utf8');
}

function addVersion(path) {
    const data = fs.readFileSync(path, 'utf8');
    const result = data.replace('__MOLSTAR_VERSION__', VERSION);
    fs.writeFileSync(path, result, 'utf8');
}

function copyViewer() {
    console.log('\n###', 'copy viewer files');
    const viewerBuildPath = path.resolve(buildDir, 'viewer/');
    const viewerDeployPath = path.resolve(localPath, 'viewer/');
    fse.copySync(viewerBuildPath, viewerDeployPath, { overwrite: true });
    addAnalytics(path.resolve(viewerDeployPath, 'index.html'));
    addManifest(path.resolve(viewerDeployPath, 'index.html'));
    addPwa(path.resolve(viewerDeployPath, 'index.html'));

    const pwaDataPath = path.resolve(dataDir, 'pwa/');
    fse.copySync(pwaDataPath, viewerDeployPath, { overwrite: true });
    addVersion(path.resolve(viewerDeployPath, 'sw.js'));
}

function copyMe() {
    console.log('\n###', 'copy me files');
    const meBuildPath = path.resolve(buildDir, 'mesoscale-explorer/');
    const meDeployPath = path.resolve(localPath, 'me/viewer/');
    fse.copySync(meBuildPath, meDeployPath, { overwrite: true });
    addAnalytics(path.resolve(meDeployPath, 'index.html'));
}

function copyMVSStories() {
    console.log('\n###', 'copy MVS stories files');
    const mvsStoriesBuildPath = path.resolve(buildDir, 'mvs-stories/');
    const mvsStoriesDeployPath = path.resolve(localPath, `stories-viewer/v${MVS_STORIES_VERSION}/`);
    fse.copySync(mvsStoriesBuildPath, mvsStoriesDeployPath, { overwrite: true });
    addAnalytics(path.resolve(mvsStoriesDeployPath, 'index.html'));

    // TODO: add PWA
    // addManifest(path.resolve(mvsStoriesDeployPath, 'index.html'));
    // addPwa(path.resolve(mvsStoriesDeployPath, 'index.html'));
}

function copyDemo(name) {
    console.log('\n###', `copy demo files for ${name}`);
    const demoBuildPath = path.resolve(buildDir, `examples/${name}/`);
    const demoDeployPath = path.resolve(localPath, `demos/${name}/`);
    fse.copySync(demoBuildPath, demoDeployPath, { overwrite: true });
    addAnalytics(path.resolve(demoDeployPath, 'index.html'));
}

function copyDemos() {
    console.log('\n###', 'copy demos files');
    copyDemo('lighting');
    copyDemo('alpha-orbitals');
    copyDemo('mvs-stories');
}

function copyFiles() {
    try {
        copyViewer();
        copyMe();
        copyMVSStories();
        copyDemos();
    } catch (e) {
        console.error(e);
    }
}

function copyToRepository() {
    console.log('\n###', 'copy repository files');
    fse.copySync(localPath, repositoryPath, { overwrite: true });
}

function syncRepository() {
    console.log('\n###', 'sync repository');
    if (!fs.existsSync(path.resolve(repositoryPath, '.git/'))) {
        console.log('\n###', 'clone repository');
        git()
            .outputHandler(log)
            .clone(remoteUrl, repositoryPath)
            .fetch(['--all'])
            .exec(copyToRepository);
    } else {
        console.log('\n###', 'update repository');
        git()
            .outputHandler(log)
            .fetch(['--all'])
            .reset(['--hard', 'origin/master'])
            .exec(copyToRepository);
    }
}

function commit() {
    console.log('\n###', 'commit changes');
    git()
        .outputHandler(log)
        .add(['-A'])
        .commit(`Updated Apps and Demos
- Mol* version: ${VERSION}
- MVS Stories version: ${MVS_STORIES_VERSION}`)
        .push();
}

if (!fs.existsSync(localPath)) {
    console.log('\n###', 'create localPath');
    fs.mkdirSync(localPath, { recursive: true });
}

const argParser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Mol* Deploy'
});
argParser.add_argument('--local',{
    help: 'Do not commit to remote repository.',
    required: false,
    action: 'store_true',
});
const args = argParser.parse_args();

copyFiles();

if (args.local) {
    process.exit(0);
}

if (!fs.existsSync(repositoryPath)) {
    console.log('\n###', 'create repositoryPath');
    fs.mkdirSync(repositoryPath, { recursive: true });
}

process.chdir(repositoryPath);
syncRepository();
commit();