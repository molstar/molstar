// by Tom Madams
// Simple:
// https://imdoingitwrong.wordpress.com/2011/01/31/light-attenuation/
//
// Improved
// https://imdoingitwrong.wordpress.com/2011/02/10/improved-light-attenuation/
float attenuation(const in float r, const in float f, const in float d) {
    float denom = d / r + 1.0;
    float attenuation = 1.0 / (denom*denom);
    float t = (attenuation - f) / (1.0 - f);
    return max(t, 0.0);
}

#pragma glslify: export(attenuation)