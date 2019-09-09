
Early JPEG XL support for Python image processing libraries
-----------------------------------------------------------

Standardization of JPEG XL is still in progress. Nevertheless, we can
already now make the most widely used Python image processing
libraries behave as if they had support for reading one particular
type of JPEG XL images for which the binary format likely has already
stabilized in its final form: losslessly recompressed JPEG images.


How it works
------------

In Python code, if one imports the module 'jxl_library_patches'
*after* image-processing modules, this will scan the list of imported
modules for known image libraries and patch them up in such a way that
their image-reading functions behave as if they understood JPEG XL
losslessly recompressed JPEG files. This works by checking headers and
using an external process to convert JPEG XL to JPEG on-the-fly,
deleting the intermediate file once it is no longer needed.

The 'jxl_library_patches' module scans sys.modules and (mostly) monkey-patches
up image loading in other modules, in an idempotent way. It is
controlled by two environment variables:

  "JXL_DEBUG" (default: "0"), which, when set to "1", will make the
  module print debug-information on which other modules it patched as it
  goes along, and for some libraries also what image-conversions it
  performed.

  "BIN_JPEG_FROM_JXL" (default: "dbrunsli"), which names an external
  converter command with signature {converter} {infile} {outfile} that
  converts JPEG XL-recompressed-JPEG to JPEG.


Supported Libraries
-------------------

This module patches JPEG XL support into (Python2 and Python3
versions of) these libraries:

  - OpenCV: cv2.imread() can read jpeg-xl files.

  - imageio: A new "JPEG XL" format will be registered that can read
      individual jpeg-xl files, from any source.

  - PythonMagic: PythonMagick.Image(filename) can read JPEG XL files.
      (Constructing Images from in-memory sources is not yet
      supported.)

  - PIL(/Pillow): A new jpeg-xl file format will be registered
      (extension: '.jxl') for JPEG XL-recompressed-JPEG.

  - scipy.misc.imread() / matplotlib.pyplot.imread() /
    matplotlib.image.imread():
      These internally use PIL(/Pillow) and work automatically as well.

