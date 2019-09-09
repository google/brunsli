# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.
"""Extends the 'imageio' to support basic reading of JPEG-XL recompressed JPEG.

This registers a file format that will call an external converter which will
(securely) convert JPEG-XL to a temporary JPEG file behind the scenes.
"""

from . import jxl_utils
import imageio
import imageio.core


# Extending imageio with plugins:
# https://imageio.readthedocs.io/en/stable/plugins.html
class _JpegXLFormat(imageio.core.Format):
  """Jpeg-XL image format. Currently only supports recompressed JPEG."""

  def _can_read(self, request):
    return (request.mode[1] in (self.modes + '?') and
            request.firstbytes.startswith(
                jxl_utils.JPEGXL_RECOMPRESSED_JPEG_HEADER))

  def _can_write(self, request):
    return False

  class Reader(imageio.core.Format.Reader):
    """Image format Reader implementation."""

    def _open(self):
      """Opens the reader."""
      self._filename = self.request.get_local_filename()

    def _close(self):
      """Closes the reader."""

    def _get_length(self):
      """Returns the number of images."""
      return 1

    def _get_data(self, index):
      """Returns the data and metadata for the image with given index."""
      if not 0 <= index < 1:
        raise IndexError('Image index out of range: %s' % (index,))
      metadata = {}
      data = jxl_utils.call_with_jxl_filename_arg1_replaced_by_temp_jpeg_file(
          imageio.imread, self._filename)
      return data, metadata

    def _get_meta_data(self, index):
      return {}  # This format does not yet support returning metadata.

  class Writer(imageio.core.Format.Writer):
    """Image format Writer implementation."""

    def _open(self):
      raise NotImplementedError('Not implemented.')

    def _close(self):
      """Closes the writer."""

    def _append_data(self, image_data, metadata):
      del image_data, metadata  # Unused.
      raise NotImplementedError('Not implemented.')

    def set_meta_data(self, metadata):
      del metadata  # Unused.
      raise NotImplementedError('Not implemented.')


def register_jxl_support(imageio_module):
  """Registers JPEG-XL with imageio."""
  # Above, we are using "import imageio.core" to get a parent class.
  # We should make sure that the `imageio` module implicitly referred like this
  # matches the imageio module which we are extending with Jpeg-XL support.
  # In principle, we could avoid this by importing from the module passed in,
  # and defining a closure-class, but the situations when this fails are so
  # obscure that this should not be an issue for easing the transition into
  # Jpeg-XL.
  assert id(imageio_module) == id(imageio), (
      'User-passed `imageio` module differs from `import imageio` result.')
  imageio_module.formats.add_format(
      _JpegXLFormat(
          'JPEG-XL',
          'JPEG-XL (partial support; reads recompressed JPEG only.)',
          '.jxl',
          'iI',
      ))
