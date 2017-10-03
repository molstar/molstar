// import buble from 'rollup-plugin-buble';
import json from 'rollup-plugin-json';
import resolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';

var path = require('path');
var pkg = require('./package.json');
var external = Object.keys(pkg.dependencies);

export default {
  input: 'build/js/src/index.js',
  plugins: [
    resolve({
      jsnext: true,
      main: true
    }),
    commonjs(),
    json(),
    // buble()
  ],
  output: [
    {
      file: "build/js/molio.dev.js",
      format: 'umd',
      name: 'MOLIO',
      sourcemap: false
    },
    // {
    //   file: "build/js/molio.esm.js",
    //   format: 'es',
    //   sourcemap: false
    // }
  ],
  external: external,
  sourcemap: false,
  onwarn(warning, warn) {
    if (warning.code === 'THIS_IS_UNDEFINED') return;
    warn(warning); // this requires Rollup 0.46
  }
};