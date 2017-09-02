import buble from 'rollup-plugin-buble';
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
    buble()
  ],
  output: [
    {
      file: "build/js/molio.dev.js",
      format: 'umd',
      name: 'MOLIO',
      sourcemap: true
    },
    {
      file: "build/js/molio.esm.js",
      format: 'es',
      sourcemap: true
    }
  ],
  external: external,
  sourcemap: true
};