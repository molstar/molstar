import * as esbuild from 'esbuild';
import * as fs from 'fs';
import * as path from 'path';
import { sassPlugin } from 'esbuild-sass-plugin';
import packageJSON from './package.json' with { type: 'json'};

function mkDir(dir) {
    if (!fs.existsSync(dir)) {
        fs.mkdirSync(dir);
    }
}

function fileLoaderPlugin(options) {
    mkDir(options.out, { recursive: true });
    
    return {
        name: 'file-loader',
        setup(build) {
            build.onLoad({ filter: /\.jpg$/ }, async (args) => {
                const name = path.basename(args.path);
                mkDir(path.resolve(options.out, 'images'));
                fs.copyFileSync(args.path, path.resolve(options.out, 'images', name));
                return {
                    contents: `images/${name}`,
                    loader: 'text',
                }
            });
            build.onLoad({ filter: /\.(html|ico)$/ }, async (args) => {
                const name = path.basename(args.path);
                fs.copyFileSync(args.path, path.resolve(options.out, name));
                return {
                    contents: '',
                    loader: 'empty',
                }
            });
        },
    }
}

const versionPlugin = {
    name: 'version',
    setup(build) {
        build.onStart(() => {
            fs.writeFileSync(
                path.resolve('./lib/mol-plugin/version.js'),
                `export var PLUGIN_VERSION = '${packageJSON.version}';\nexport var PLUGIN_VERSION_DATE = new Date(typeof __MOLSTAR_DEBUG_TIMESTAMP__ !== 'undefined' ? __MOLSTAR_DEBUG_TIMESTAMP__ : ${new Date().valueOf()});`);
        });
    }
}

const start = Date.now();

await esbuild.build({
    entryPoints: ['./src/apps/viewer/index.ts'],
    tsconfig: './tsconfig.json',
    bundle: true,
    // minify: true,
    globalName: 'molstar',
    outfile: './build/viewer/molstar.js',
    plugins: [
        versionPlugin,
        fileLoaderPlugin({ out: `./build/viewer` }),
        sassPlugin({ type: 'css' }),
    ],
    external: ['crypto', 'fs', 'path', 'stream'],
    loader: {
    },
});

console.log(`Build time: ${Date.now() - start}ms`);