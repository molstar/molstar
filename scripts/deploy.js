/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const git = require('simple-git');
const path = require('path');
const fs = require("fs");
const fse = require("fs-extra");

const remoteUrl = "https://github.com/molstar/molstar.github.io.git";
const buildDir = path.resolve(__dirname, '../build/');
const deployDir = path.resolve(buildDir, 'deploy/');
const localPath = path.resolve(deployDir, 'molstar.github.io/');

const analyticsTag = /<!-- __MOLSTAR_ANALYTICS__ -->/g;
const analyticsCode = `<!-- Cloudflare Web Analytics --><script defer src='https://static.cloudflareinsights.com/beacon.min.js' data-cf-beacon='{"token": "c414cbae2d284ea995171a81e4a3e721"}'></script><!-- End Cloudflare Web Analytics -->`;

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

function copyViewer() {
    console.log('\n###', 'copy viewer files');
    const viewerBuildPath = path.resolve(buildDir, '../build/viewer/');
    const viewerDeployPath = path.resolve(localPath, 'viewer/');
    fse.copySync(viewerBuildPath, viewerDeployPath, { overwrite: true });
    addAnalytics(path.resolve(viewerDeployPath, 'index.html'));
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
        .exec(copyViewer)
        .add(['-A'])
        .commit('updated viewer')
        .push();
} else {
    console.log('\n###', 'update repository');
    git()
        .outputHandler(log)
        .fetch(['--all'])
        .reset(['--hard', 'origin/master'])
        .exec(copyViewer)
        .add(['-A'])
        .commit('updated viewer')
        .push();
}