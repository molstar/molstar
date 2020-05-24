const path = require('path');
const webpack = require('webpack');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const VersionFile = require('webpack-version-file-plugin');

const sharedConfig = {
    module: {
        rules: [
            {
                test: /\.(html|ico)$/,
                use: [{
                    loader: 'file-loader',
                    options: { name: '[name].[ext]' }
                }]
            },
            {
                test: /\.(s*)css$/,
                use: [
                    MiniCssExtractPlugin.loader,
                    { loader: 'css-loader', options: { sourceMap: false } },
                    { loader: 'sass-loader', options: { sourceMap: false } },
                ]
            }
        ]
    },
    plugins: [
        new ExtraWatchWebpackPlugin({
            files: [
                './lib/**/*.scss',
                './lib/**/*.html'
            ],
        }),
        new webpack.DefinePlugin({
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG),
            '__MOLSTAR_DEBUG_TIMESTAMP__': webpack.DefinePlugin.runtimeValue(() => `${new Date().valueOf()}`, true)
        }),
        new MiniCssExtractPlugin({ filename: 'molstar.css',  }),
        new VersionFile({
            extras: { timestamp: `${new Date().valueOf()}` },
            packageFile: path.resolve(__dirname, 'package.json'),
            templateString: `export const PLUGIN_VERSION = '<%= package.version %>';\nexport const PLUGIN_VERSION_DATE = new Date(typeof __MOLSTAR_DEBUG_TIMESTAMP__ !== 'undefined' ? __MOLSTAR_DEBUG_TIMESTAMP__ : <%= extras.timestamp %>);`,
            outputFile: path.resolve(__dirname, 'lib/mol-plugin/version.js')
        })
    ],
    resolve: {
        modules: [
            'node_modules',
            path.resolve(__dirname, 'lib/')
        ],
    },
    watchOptions: {
        aggregateTimeout: 750
    },
    devtool: ''
}

function createEntry(src, outFolder, outFilename, isNode) {
    return {
        node: isNode ? void 0 : { fs: 'empty' }, // TODO find better solution? Currently used in file-handle.ts
        target: isNode ? 'node' : void 0,
        entry: path.resolve(__dirname, `lib/${src}.js`),
        output: { filename: `${outFilename}.js`, path: path.resolve(__dirname, `build/${outFolder}`) },
        ...sharedConfig
    }
}

function createEntryPoint(name, dir, out, library) {
    return {
        node: { fs: 'empty' }, // TODO find better solution? Currently used in file-handle.ts
        entry: path.resolve(__dirname, `lib/${dir}/${name}.js`),
        output: { filename: `${library || name}.js`, path: path.resolve(__dirname, `build/${out}`), library: library || out, libraryTarget: 'umd' },
        ...sharedConfig
    }
}

function createNodeEntryPoint(name, dir, out) {
    return {
        target: 'node',
        entry: path.resolve(__dirname, `lib/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
        externals: {
            argparse: 'require("argparse")',
            'node-fetch': 'require("node-fetch")',
            'util.promisify': 'require("util.promisify")',
            xhr2: 'require("xhr2")',
        },
        ...sharedConfig
    }
}

function createApp(name, library) { return createEntryPoint('index', `apps/${name}`, name, library) }
function createExample(name) { return createEntry(`examples/${name}/index`, `examples/${name}`, 'index') }
function createBrowserTest(name) { return createEntryPoint(name, 'tests/browser', 'tests') }
function createNodeApp(name) { return createNodeEntryPoint('index', `apps/${name}`, name) }

module.exports = {
    createApp,
    createEntry,
    createExample,
    createBrowserTest,
    createNodeEntryPoint,
    createNodeApp
}