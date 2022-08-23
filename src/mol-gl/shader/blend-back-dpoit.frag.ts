export const blendBackDpoit_frag = `
  precision highp float;

  uniform sampler2D tDpoitBackColor;
  uniform vec2 uTexSize;

  void main() {
      vec2 coords = gl_FragCoord.xy / uTexSize;
      gl_FragColor = texture2D(tDpoitBackColor, coords);
      if (gl_FragColor.a == 0.0) {
          discard;
      }
  }
`;
