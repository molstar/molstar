// (c) 2014 Mikola Lysenko. MIT License
// https://github.com/glslify/glsl-specular-phong

float phongSpecular(vec3 lightDirection, vec3 viewDirection, vec3 surfaceNormal, float shininess) {
    //Calculate Phong power
    vec3 R = -reflect(lightDirection, surfaceNormal);
    return pow(max(0.0, dot(viewDirection, R)), shininess);
}

#pragma glslify: export(phongSpecular)