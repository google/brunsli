# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.
"""Adds basic reading of JPEG-XL recompressed JPEG to the 'PIL' library.

This registers a file format that will call an external converter which will
(securely) convert JPEG-XL to a temporary JPEG file behind the scenes.
"""

import io
import tempfile

from . import jxl_utils
from PIL import Image
from PIL import ImageFile


# Extending Pillow:
# https://pillow.readthedocs.io/en/5.1.x/
#   handbook/writing-your-own-file-decoder.html
class JpegXLImageFile(ImageFile.ImageFile):
  """JPEG-XL Image File."""

  format = 'JPEG-XL'
  format_description = 'JPEG-XL raster image'

  def _open(self):
    header = self.fp.read(6)
    if not header.startswith(jxl_utils.JPEGXL_RECOMPRESSED_JPEG_HEADER):
      # 'SyntaxError' is technically wrong, but expected by Pillow here.
      raise SyntaxError('Not a JXL file.')

    def fetch_jpeg(filename):
      with open(filename, 'rb') as h:
        return io.BytesIO(h.read())
    # We first have to make sure that the data which we process has a presence
    # on the filesystem so that the converter can see it. This might not
    # have been the case if file data was e.g. provided as a io.BytesIO
    # in-memory stream.
    with tempfile.NamedTemporaryFile(suffix='.jxl') as h:
      data = header + self.fp.read()
      h.write(data)
      h.flush()
      jpg_io = jxl_utils.call_with_jxl_filename_arg1_replaced_by_temp_jpeg_file(
          fetch_jpeg, h.name)
    self._jpeg = Image.open(jpg_io)
    # Forward information from underlying ._jpeg.
    self.mode = self._jpeg.mode
    self.size = self._jpeg.size
    # Forward .tile and .fp from JPEG.
    # Caution: This is relying on not-quite-properly-documented-behavior here.
    self.tile = self._jpeg.tile
    self.fp = self._jpeg.fp


def register_jxl_support(pil_image_module):
  """Registers JPEG-XL with PIL/Pillow."""
  # We are making the implicit assumption here that the ImageFile module
  # imported above for the ImageFile.ImageFile parent class is compatible
  # with pil_image_module.
  pil_image_module.register_open(JpegXLImageFile.format, JpegXLImageFile)
  pil_image_module.register_extension(JpegXLImageFile.format, '.jxl')
