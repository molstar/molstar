import * as esbuild from 'esbuild';
import * as fs from 'fs';
import * as path from 'path';
import { sassPlugin } from 'esbuild-sass-plugin';
import { Logger } from 'sass';

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
                await fs.promises.copyFile(args.path, path.resolve(options.out, 'images', name));
                return {
                    contents: `images/${name}`,
                    loader: 'text',
                }
            });
            build.onLoad({ filter: /\.(html|ico)$/ }, async (args) => {
                const name = path.basename(args.path);
                await fs.promises.copyFile(args.path, path.resolve(options.out, name));
                return {
                    contents: '',
                    loader: 'empty',
                }
            });
            // build.onResolve({ filter: /\.(html|ico)$/ }, async (args) => {
            //     console.log(args);
            //     return { watchFiles: [args.path] };
            // });
        },
    }
}

const apps = ['viewer', 'docking-viewer', 'mesoscale-explorer'];
const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting', 'alpha-orbitals', 'alphafolddb-pae', 'mvs-kinase-story', 'ihm-restraints'];

async function build(name, kind) {
    const prefix = kind === 'app'
        ? `./build/${name}`
        : `./build/examples/${name}`;

    let entry = `./src/${kind}s/${name}/index.ts`;
    if (!fs.existsSync(entry)) {
        entry = `./src/${kind}s/${name}/index.tsx`;
    }

    const ctx = await esbuild.context({
        entryPoints: [entry],
        tsconfig: './tsconfig.json',
        bundle: true,
        globalName: 'molstar',
        outfile: kind === 'app'
            ? `./build/${name}/molstar.js`
            : `./build/examples/${name}/index.js`,
        plugins: [
            fileLoaderPlugin({ out: prefix }),
            sassPlugin({
                type: 'css',
                logger: {
                    warn: (msg) => console.warn(msg),
                    debug: () => { },
                }
            }),
        ],
        external: ['crypto', 'fs', 'path', 'stream'],
        loader: {
        },
    });

    ctx.watch();
}

build('viewer', 'app');

for (const app of apps) {
    // build(app, 'app');
}

const ctx = await esbuild.context({});
ctx.serve({
    servedir: './build',
    port: 5888,
});

console.log('Serving on http://localhost:5888');


// await Promise.all(apps.map(name => build(name, 'app')));
// await Promise.all(examples.map(name => build(name, 'example')));

// console.log(`Build time: ${Date.now() - start}ms`);