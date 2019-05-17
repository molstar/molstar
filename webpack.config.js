const path = require('path');
const webpack = require('webpack');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
// const CircularDependencyPlugin = require('circular-dependency-plugin');

const sharedConfig = {
    module: {
        rules: [
            {
                loader: 'file-loader',
                test: /\.(woff2?|ttf|otf|eot|svg|html)$/,
                include: [path.resolve(__dirname, 'build/src/')],
                options: {
                    name: '[name].[ext]'
                }
            },
            {
                test: /\.(s*)css$/,
                use: [MiniCssExtractPlugin.loader, 'css-loader', 'resolve-url-loader', 'sass-loader']
            }
        ]
    },
    plugins: [
        // new CircularDependencyPlugin({
        //     include: [ path.resolve(__dirname, 'build/src/') ],
        //     failOnError: false,
        //     cwd: process.cwd(),
        // }),
        new ExtraWatchWebpackPlugin({
            files: [
                './build/src/**/*.vert',
                './build/src/**/*.frag',
                './build/src/**/*.glsl',
                './build/src/**/*.scss',
                './build/src/**/*.html'
            ],
        }),
        new webpack.DefinePlugin({
            __PLUGIN_VERSION_TIMESTAMP__: webpack.DefinePlugin.runtimeValue(() => `${new Date().valueOf()}`, true),
            'process.env.NODE_ENV': JSON.stringify(process.env.NODE_ENV),
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG)
        }),
        new MiniCssExtractPlugin({ filename: 'app.css' })
    ],
    resolve: {
        modules: [
            'node_modules',
            path.resolve(__dirname, 'build/src/')
        ],
    }
}


function createEntry(src, outFolder, outFilename, isNode) {
    return {
        node: isNode ? void 0 : { fs: 'empty' }, // TODO find better solution? Currently used in file-handle.ts
        target: isNode ? 'node' : void 0,
        entry: path.resolve(__dirname, `build/src/${src}.js`),
        output: { filename: `${outFilename}.js`, path: path.resolve(__dirname, `build/${outFolder}`) },
        ...sharedConfig
    }
}

function createEntryPoint(name, dir, out) {
    return {
        node: { fs: 'empty' }, // TODO find better solution? Currently used in file-handle.ts
        entry: path.resolve(__dirname, `build/src/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
        ...sharedConfig
    }
}

function createNodeEntryPoint(name, dir, out) {
    return {
        target: 'node',
        entry: path.resolve(__dirname, `build/src/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
        ...sharedConfig
    }
}

function createApp(name) { return createEntryPoint('index', `apps/${name}`, name) }
function createBrowserTest(name) { return createEntryPoint(name, 'tests/browser', 'tests') }
function createNodeApp(name) { return createNodeEntryPoint('index', `apps/${name}`, name) }

module.exports = [
    createApp('viewer'),
    createApp('basic-wrapper'),
    createEntry('examples/proteopedia-wrapper/index', 'examples/proteopedia-wrapper', 'index'),
    createNodeApp('state-docs'),
    createNodeEntryPoint('preprocess', 'servers/model', 'model-server'),
    createApp('model-server-query'),

    createBrowserTest('font-atlas'),
    createBrowserTest('render-asa'),
    createBrowserTest('marching-cubes'),
    createBrowserTest('render-lines'),
    createBrowserTest('render-mesh'),
    createBrowserTest('render-shape'),
    createBrowserTest('render-spheres'),
    createBrowserTest('render-structure'),
    createBrowserTest('render-text'),
]