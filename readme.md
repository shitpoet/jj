# jj

Compact lossless optimizing recompression module for JPEG.

The goal of the module is to be used in Web front-end as WASM module. It can be helpful in cases of uploading not optimized images over slow connectios.  In this cases it is faster to reencode image than to send it as is.

24-bit 3-channel YCbCr JPEG images are suppoted only. Non power-of-two sampling factors are not supported (almost no images have they as far as I know), usual cofigurations (4:4:4, 4:2:2 etc) are supported.

In the current version all meta-information is ignored and stripped. Including colorpace and EXIF orientation.

The gain is from 0% to 40% depending of the image. Already optimized images (f.ex. by `jpegoptim` or `mozjpeg` are not compressable). But the resulting image is almost never bigger than the original.

Recompression is completely lossless: it does not do any IDCT/DCT or requantization.

The module is heavely stipped down version of `libgpeg-turbo 2.1.4` (compiled in `libjpeg v6` compatibility mode). To reduce the size most of the original code was removed. So the module only includes the very minimum for parsing of image, decoding DCT coefficients, optimizing huffman tables and reencoding image back.

The x86-64 binary size is around 10k instead of 150-200k for `libjpeg` library.

When in doubt it fails: corrupted blocks, mismatch of any information etc. So it does not guarantee to provide output at all. Please, use the original image if the optimization fails. Alrough I believe it will work on 99% of JPEG images as of 2022.

Performance should be comparable to `libjpeg`/`libjpeg-turbo`. (Needs testing.)

