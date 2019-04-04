#ifdef dGeoTexture
    // aGroup is used as a triangle index here and the group id is retirieved from the tGroup texture
    float group = readFromTexture(tGroup, aGroup, uGeoTexDim).a;
#else
    float group = aGroup;
#endif