#if defined(dSizeType_uniform)
    float size = uSize;
#elif defined(dSizeType_attribute)
    float size = aSize;
#elif defined(dSizeType_instance)
    float size = readFromTexture(tSize, aInstance, uSizeTexDim).a;
#elif defined(dSizeType_group)
    float size = readFromTexture(tSize, aGroup, uSizeTexDim).a;
#elif defined(dSizeType_groupInstance)
    float size = readFromTexture(tSize, aInstance * float(uGroupCount) + aGroup, uSizeTexDim).a;
#endif

#if defined(dSizeType_instance) || defined(dSizeType_group) || defined(dSizeType_groupInstance)
    size = decodeFloatLog(size);
#endif