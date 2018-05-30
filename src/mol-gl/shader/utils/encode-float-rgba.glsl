vec4 encodeFloatRGBA(float v) {
  vec4 enc = vec4(1.0, 255.0, 65025.0, 16581375.0) * v;
  enc = frac(enc);
  enc -= enc.yzww * float4(1.0/255.0,1.0/255.0,1.0/255.0,0.0);
  return enc;
}

#pragma glslify: export(encodeFloatRGBA)