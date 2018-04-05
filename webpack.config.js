const path = require('path');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
module.exports = {
    module: {
        rules: [
            {
                loader: 'raw-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, "build/node_modules/") ],
            },
            {
                loader: 'glslify-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, "build/node_modules/") ]
            }
        ]
    },
    plugins: [
        new ExtraWatchWebpackPlugin({
            files: [
                './build/node_modules/**/*.vert',
                './build/node_modules/**/*.frag',
                './build/node_modules/**/*.glsl'
            ],
        }),
    ],
}