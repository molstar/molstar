const path = require('path');
const webpack = require('webpack');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
// const CircularDependencyPlugin = require('circular-dependency-plugin');

const sharedConfig = {
    module: {
        rules: [
            {
                test: /\.(woff2?|ttf|otf|eot|svg|html|ico)$/,
                use: [{
                    loader: 'file-loader',
                    options: { name: '[name].[ext]' }
                }]
            },
            {
                test: /\.(s*)css$/,
                use: [
                    MiniCssExtractPlugin.loader,
                    'css-loader', 'resolve-url-loader', 'sass-loader'
                ]
            }
        ]
    },
    plugins: [
        // new CircularDependencyPlugin({
        //     include: [ path.resolve(__dirname, 'lib/') ],
        //     failOnError: false,
        //     cwd: process.cwd(),
        // }),
        new ExtraWatchWebpackPlugin({
            files: [
                './lib/**/*.scss',
                './lib/**/*.html'
            ],
        }),
        new webpack.DefinePlugin({
            __VERSION__: webpack.DefinePlugin.runtimeValue(() => JSON.stringify(require('./package.json').version), true),
            __VERSION_TIMESTAMP__: webpack.DefinePlugin.runtimeValue(() => `${new Date().valueOf()}`, true),
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG)
        }),
        new MiniCssExtractPlugin({ filename: 'app.css' })
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

function createEntryPoint(name, dir, out) {
    return {
        node: { fs: 'empty' }, // TODO find better solution? Currently used in file-handle.ts
        entry: path.resolve(__dirname, `lib/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
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

function createApp(name) { return createEntryPoint('index', `apps/${name}`, name) }
function createBrowserTest(name) { return createEntryPoint(name, 'tests/browser', 'tests') }
function createNodeApp(name) { return createNodeEntryPoint('index', `apps/${name}`, name) }

module.exports = {
    createApp,
    createEntry,
    createBrowserTest,
    createNodeEntryPoint,
    createNodeApp
}