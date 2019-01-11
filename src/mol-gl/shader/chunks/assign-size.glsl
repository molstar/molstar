#if defined(dSizeType_uniform)
    float size = uSize;
#elif defined(dSizeType_attribute)
    float size = aSize;
#elif defined(dSizeType_instance)
    float size = readFromTexture(tSize, aInstance, uSizeTexDim).r;
#elif defined(dSizeType_group)
    float size = readFromTexture(tSize, aGroup, uSizeTexDim).r;
#elif defined(dSizeType_groupInstance)
    float size = readFromTexture(tSize, aInstance * float(uGroupCount) + aGroup, uSizeTexDim).r;
#endif