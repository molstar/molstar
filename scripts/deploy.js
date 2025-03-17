/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const git = require('simple-git');
const path = require('path');
const fs = require("fs");
const fse = require("fs-extra");

const VERSION = require(path.resolve(__dirname, '../package.json')).version;

const remoteUrl = "https://github.com/molstar/molstar.github.io.git";
const dataDir = path.resolve(__dirname, '../data/');
const buildDir = path.resolve(__dirname, '../build/');
const deployDir = path.resolve(buildDir, 'deploy/');
const localPath = path.resolve(deployDir, 'molstar.github.io/');

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

function copyDemos() {
    console.log('\n###', 'copy demos files');
    const lightingBuildPath = path.resolve(buildDir, 'examples/lighting/');
    const lightingDeployPath = path.resolve(localPath, 'demos/lighting/');
    fse.copySync(lightingBuildPath, lightingDeployPath, { overwrite: true });
    addAnalytics(path.resolve(lightingDeployPath, 'index.html'));

    const orbitalsBuildPath = path.resolve(buildDir, 'examples/alpha-orbitals/');
    const orbitalsDeployPath = path.resolve(localPath, 'demos/alpha-orbitals/');
    fse.copySync(orbitalsBuildPath, orbitalsDeployPath, { overwrite: true });
    addAnalytics(path.resolve(orbitalsDeployPath, 'index.html'));
}

function copyFiles() {
    try {
        copyViewer();
        copyMe();
        copyDemos();
    } catch (e) {
        console.error(e);
    }
}

if (!fs.existsSync(localPath)) {
    console.log('\n###', 'create localPath');
    fs.mkdirSync(localPath, { recursive: true });
}

process.chdir(localPath);

if (!fs.existsSync(path.resolve(localPath, '.git/'))) {
    console.log('\n###', 'clone repository');
    git()
        .outputHandler(log)
        .clone(remoteUrl, localPath)
        .fetch(['--all'])
        .exec(copyFiles)
        .add(['-A'])
        .commit('updated viewer & demos')
        .push();
} else {
    console.log('\n###', 'update repository');
    git()
        .outputHandler(log)
        .fetch(['--all'])
        .reset(['--hard', 'origin/master'])
        .exec(copyFiles)
        .add(['-A'])
        .commit('updated viewer & demos')
        .push();
}