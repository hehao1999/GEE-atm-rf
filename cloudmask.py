"""影像去云"""
def maskS2clouds(image):
    qa = image.select('QA60')

    # Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11

    # Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudBitMask).eq(0) and qa.bitwiseAnd(cirrusBitMask).eq(0)
    # Return the masked and scaled data, without the QA bands.
    return image.updateMask(mask)\
        .divide(10000).select("B.*")\
        # .copyProperties(image, ["system:time_start"])
