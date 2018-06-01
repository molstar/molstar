// (c) 2014 Mikola Lysenko. MIT License
// https://github.com/glslify/glsl-specular-phong

float phongSpecular(const in vec3 lightDirection, const in vec3 viewDirection, const in vec3 surfaceNormal, const in float shininess) {
    //Calculate Phong power
    vec3 R = -reflect(lightDirection, surfaceNormal);
    return pow(max(0.0, dot(viewDirection, R)), shininess);
}

#pragma glslify: export(phongSpecular)