module.exports = {
    module: {
        rules: [
            {
                test: /\.(glsl|frag|vert)$/,
                loader: 'raw-loader',
            },
            {
                test: /\.(glsl|frag|vert)$/,
                loader: 'glslify-loader',
            }
        ]
    }
}