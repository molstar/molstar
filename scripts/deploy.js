/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const git = require('simple-git')
const path = require('path')
const fs = require("fs")
const fse = require("fs-extra")

const remoteUrl = "https://github.com/molstar/molstar.github.io.git"
const buildDir = path.resolve(__dirname, '../build/')
const deployDir = path.resolve(buildDir, 'deploy/')
const localPath = path.resolve(deployDir, 'molstar.github.io/')

if (!fs.existsSync(localPath) || !git(localPath).checkIsRepo()) {
    fs.mkdirSync(localPath, { recursive: true })
    git(deployDir)
        .silent(false)
        .clone(remoteUrl)
        .fetch(['--all'])
} else {
    git(localPath)
        .silent(false)
        .fetch(['--all'])
        .reset(['--hard', 'origin/master'])
}

const viewerBuildPath = path.resolve(buildDir, '../build/viewer/')
const viewerDeployPath = path.resolve(localPath, 'viewer/')
fse.copySync(viewerBuildPath, viewerDeployPath, { overwrite: true })

git(localPath)
    .silent(false)
    .add(['-A'])
    .commit('updated viewer')
    .push()

// #!/usr/bin/env bash

// LEVEL=$1
// DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";

// cd ${DIR};
// mkdir -p ../build/;
// cd ../build/;

// if [ -d "arose.github.io" ]; then
// 	cd ./arose.github.io/;
// 	git fetch --all;
// 	git reset --hard origin/master;
// 	cd ../
// else
// 	git clone "https://github.com/arose/arose.github.io.git";
// fi

// if [ "$LEVEL" = "prerelease" ]; then
// 	cd ./arose.github.io/ngldev/;
// else
// 	cd ./arose.github.io/ngl/;
// fi

// cp -r ${DIR}/../data/. ./data/;
// cp -r ${DIR}/../examples/css/. ./css/;
// cp -r ${DIR}/../examples/fonts/. ./fonts/;
// cp -r ${DIR}/../examples/js/. ./js/;
// cp -r ${DIR}/../examples/scripts/. ./scripts/;
// cp -r ${DIR}/../build/docs/. ./api/;
// cp -r ${DIR}/../build/gallery/. ./gallery/;
// cp ${DIR}/../build/scriptsList.json ./scripts/list.json;

// cp ${DIR}/../dist/ngl.js ./js/ngl.js;

// cd ../;
// git add -A;
// if [ "$LEVEL" = "prerelease" ]; then
// 	git commit -m "ngldev update";
// else
// 	git commit -m "ngl update";
// fi
// git push;
