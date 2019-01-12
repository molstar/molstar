const path = require('path');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
// const CircularDependencyPlugin = require('circular-dependency-plugin');

const sharedConfig = {
    module: {
        rules: [
            {
                loader: 'raw-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ],
            },
            {
                loader: 'glslify-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ]
            },

            {
                loader: 'file-loader',
                test: /\.(woff2?|ttf|otf|eot|svg|html)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ],
                options: {
                    name: '[name].[ext]'
                }
            },
            {
                test:/\.(s*)css$/,
                use: [ MiniCssExtractPlugin.loader, 'css-loader', 'resolve-url-loader', 'sass-loader' ]
            }
        ]
    },
    plugins: [
        // new CircularDependencyPlugin({
        //     include: [ path.resolve(__dirname, 'build/node_modules/') ],
        //     failOnError: false,
        //     cwd: process.cwd(),
        // }),
        new ExtraWatchWebpackPlugin({
            files: [
                './build/node_modules/**/*.vert',
                './build/node_modules/**/*.frag',
                './build/node_modules/**/*.glsl',
                './build/node_modules/**/*.scss',
                './build/node_modules/**/*.html'
            ],
        }),
        new MiniCssExtractPlugin({ filename: "app.css" })
    ]
}

function createEntryPoint(name, dir, out) {
    return {
        entry: path.resolve(__dirname, `build/node_modules/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
        ...sharedConfig
    }
}

function createApp(name) { return createEntryPoint('index', `apps/${name}`, name) }
function createBrowserTest(name) { return createEntryPoint(name, 'tests/browser', 'tests') }

module.exports = [
    createApp('viewer'),
    createApp('model-server-query'),

    createBrowserTest('text-atlas'),
    createBrowserTest('render-text'),
    createBrowserTest('render-spheres'),
    createBrowserTest('render-mesh')
]