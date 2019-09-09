# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.
"""Extends 'PythonMagick' to support basic reading of JPEG-XL recompressed JPEG.

This registers a file format that will call an external converter which will
(securely) convert JPEG-XL to a temporary JPEG file behind the scenes.
"""

from . import jxl_utils
import PythonMagick


_patched_functions_by_id = {}


def _get_pymagic_wrapped_image_init(orig_init_func):
  """Returns a wrapped PythonMagick.Image.__init__ that understands .JXL."""
  def wrapped_init(self, *args, **kwargs):
    """{Docstring will get replaced by original function's}."""
    is_jxl = None
    try:
      with open(args[0], 'rb') as h:
        header = h.read(len(jxl_utils.JPEGXL_RECOMPRESSED_JPEG_HEADER))
        is_jxl = header == jxl_utils.JPEGXL_RECOMPRESSED_JPEG_HEADER
        if not is_jxl:
          raise RuntimeError()
    except Exception:  # pylint:disable=broad-except
      # If *anything* went wrong here, i.e. no args[0], or we could not open it
      # as a file, or we could not read a 6-bytes header, or it is not
      # JPEG-XL-recompressed, we just fall back to the base class's behavior.
      #
      # Note that PythonMagic.Image() supports constructing an image from many
      # different things. Here, we only tune creating an image from a file
      # given by name, not from some in-memory representation.
      return orig_init_func(self, *args, **kwargs)
    if not is_jxl:
      # Fall back to base class behavior.
      return orig_init_func(self, *args, **kwargs)
    # Otherwise, we have a Jpeg-XL image.
    def init_super(*jpeg_transformed_args):
      return orig_init_func(self, *jpeg_transformed_args, **kwargs)
    return jxl_utils.call_with_jxl_filename_arg1_replaced_by_temp_jpeg_file(
        init_super, *args)
  wrapped_init.__doc__ = orig_init_func.__doc__  # Forward docstring.
  return wrapped_init


def add_jxl_support(pythonmagick_module):
  pymagick_image_init = pythonmagick_module.Image.__init__
  init_id = id(pymagick_image_init)
  if init_id in _patched_functions_by_id:
    # This was already patched up, so, ignore.
    return
  wrapped_init = _get_pymagic_wrapped_image_init(pymagick_image_init)
  _patched_functions_by_id[init_id] = wrapped_init
  _patched_functions_by_id[id(wrapped_init)] = wrapped_init
  PythonMagick.Image.__init__ = wrapped_init
