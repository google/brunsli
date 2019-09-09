# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.
"""Patches major Python image libraries to support JPEG-XL.

This allows us to pretend that OpenCV2 and other relevant Python image
libraries already support JPEG-XL even where such support has not yet been
productionized.

This will be superseded soon by proper JPEG-XL support in the libraries
covered here.
"""

import os
import shutil
import subprocess
import sys
import tempfile

from . import jxl_utils


# Default value for JXL_DEBUG debug mode env-variable.
# _DEFAULT_DEBUG_ENV = '1'
_DEFAULT_DEBUG_ENV = '0'


# Global dict holding all the patches that were applied by this module.
# Key: {PyObject ID of function}. Value: {wrapped function}.
#
# Whenever a id(f): f_wrapped entry is added here, another
# id(f_wrapped): f_wrapped entry is also added which actually makes
# the patching idempotent.
#
# Since this module is only a hack, we are not making this a
# 'weak pointer hash' that would allow dropping entries that have
# become unused.
_patches_by_pyobj_id = {}


# The expected {module}.__name__ names for image library modules.
_OPENCV_MODULE_SELF_NAME = 'cv2'
_IMAGEIO_MODULE_SELF_NAME = 'imageio'
_PIL_MODULE_SELF_NAME = 'PIL.Image'
_PYTHONMAGICK_MODULE_SELF_NAME = 'PythonMagick'
_MATPLOTLIB_IMAGE_SELF_NAME = 'matplotlib.image'


### OpenCV ###


def _get_opencv_wrapped_imread(cv2_imread, debug=False):
  """Returns a JXL wrapper-function for the input cv2.imread() function."""
  # Note: The proper way to add a decoder to OpenCV will be to push a decoder
  # onto its static vector<ImageDecoder> decoders inside
  # ImageCodecInitializer().
  #
  # Make this function idempotent.
  id_imread = id(cv2_imread)
  wrapper = _patches_by_pyobj_id.get(id_imread)
  if wrapper:
    return wrapper
  # Otherwise, there was not yet a wrapper-function. Create one.
  def imread(filename, **kwargs):  # Has same name as replaced function.
    """[Docstring will get replaced]."""
    # For now, this wrapper does nothing - it just wraps the original imread().
    if debug:
      sys.stderr.write(
          'DDD wrapped OpenCV imread(filename=%r, **kwargs=%r)\n' % (
              filename, kwargs))
    if jxl_utils.is_jpegxl_recompressed_jpeg_file(filename):
      return jxl_utils.call_with_jxl_filename_arg1_replaced_by_temp_jpeg_file(
          cv2_imread, filename, **kwargs)
    else:
      return cv2_imread(filename, **kwargs)
  # We even fake the docstring.
  imread.__doc__ = cv2_imread.__doc__
  # Make both the old and the new function point to the wrapper in the global
  # 'patches' table.
  _patches_by_pyobj_id[id_imread] = imread
  _patches_by_pyobj_id[id(imread)] = imread
  return imread


def _opencv_try_patch_modules(debug=False):
  """Looks for OpenCV in loaded modules and patches all that are found."""
  # The 'cv2' module might have been imported with a different name, as in:
  # 'import cv2 as opencv'. Hence, we actually have to inspect modules.
  for module_name, module in sorted(sys.modules.items()):
    if module is None:
      # Happens e.g. for module_name == 'encodings.__builtin__' pseudo-modules.
      continue
    if module.__name__ != _OPENCV_MODULE_SELF_NAME:
      continue  # This is not OpenCV.
    # Let us make (reasonably) sure that this candidate module, which calls
    # itself 'cv2', actually is OpenCV.
    mdict = module.__dict__
    if mdict.get('CV_32FC3') != 21 or mdict.get('imread') is None:
      continue
    # We think we really have OpenCV here. Build a cv2.imload() wrapper
    # that replaces the original one. We need to be careful if different
    # OpenCV modules have been loaded - each one has to get a wrapper
    # for its own cv2.imload().
    mdict['imread'] = _get_opencv_wrapped_imread(module.imread, debug=debug)
    if debug:
      print('Patched JPEG-XL support into OpenCV module imported as %r.' % (
          module_name,))


### imageio ###


def _imageio_try_patch_modules(debug=False):
  """Looks for imageio in loaded modules and patches all matches."""
  for module_name, module in sorted(sys.modules.items()):
    if module is None:
      # Happens e.g. for module_name == 'encodings.__builtin__' pseudo-modules.
      continue
    if module.__name__ != _IMAGEIO_MODULE_SELF_NAME:
      continue  # This is not imageio.
    # Let us make (reasonably) sure that this candidate module actually is
    # imageio.
    mdict = module.__dict__
    if mdict.get('imread') is None or mdict.get('formats') is None:
      continue
    # We think we really have imageio here.
    # We only load the patching code if the user actually has 'imageio'.
    # We could not even import-at-top this otherwise.
    from . import jxl_imageio  # pylint:disable=g-import-not-at-top
    jxl_imageio.register_jxl_support(module)
    if debug:
      print('Patched JPEG-XL support into imageio module imported as %r.' % (
          module_name,))


### PIL / Pillow ###


def _pil_pillow_try_patch_modules(debug=False):
  """Looks for PIL/Pillow in loaded modules and patches all matches."""
  for module_name, module in sorted(sys.modules.items()):
    if module is None:
      # Happens e.g. for module_name == 'encodings.__builtin__' pseudo-modules.
      continue
    if module.__name__ != _PIL_MODULE_SELF_NAME:
      continue  # This is neither PIL nor Pillow.
    # Let us make (reasonably) sure that this candidate module actually is
    # PIL/Pillow.
    mdict = module.__dict__
    if not any('Image' in p for p in mdict.get('_plugins', ())):
      continue
    # We think we really have PIL/Pillow here.
    # We only load the patching code if the user actually has PIL/Pillow.
    # We could not even import-at-top this otherwise.
    from . import jxl_pillow  # pylint:disable=g-import-not-at-top
    jxl_pillow.register_jxl_support(module)
    if debug:
      print('Patched JPEG-XL support into PIL/Pillow module imported as %r.' % (
          module_name,))


### PythonMagick ###


def _pythonmagick_try_patch_modules(debug=False):
  """Looks for PythonMagick in loaded modules and patches all matches."""
  for module_name, module in sorted(sys.modules.items()):
    if module is None:
      # Happens e.g. for module_name == 'encodings.__builtin__' pseudo-modules.
      continue
    if module.__name__ != _PYTHONMAGICK_MODULE_SELF_NAME:
      continue  # This is not PythonMagick.
    # Let us make (reasonably) sure that this candidate module actually is
    # PythonMagick.
    if 'DrawableCircle' not in module.__dict__:
      continue
    # We think we really have PythonMagick here.
    # We only load the patching code if the user actually has that module.
    # We could not even import-at-top this otherwise.
    from . import jxl_pythonmagick  # pylint:disable=g-import-not-at-top
    jxl_pythonmagick.add_jxl_support(module)
    if debug:
      print(
          'Patched JPEG-XL support into PythonMagick module imported as %r.' % (
              module_name,))


######


def _try_patch_all(debug=False):
  """Scans for image-processing modules and patches in JPEG XL support."""
  if not jxl_utils.can_call_converter():
    sys.exit('Cannot call the JPEG-XL-to-JPEG converter. '
             'Please setenv the BIN_JPEG_FROM_JXL variable.')
  # `matplotlib` needs special handling: matplotlib.pyplot imports PIL and works
  # as-is, but matplotlib.image has native support for PNG only and on-demand
  # imports PIL as needed, e.g. for jpeg. We address this by loading PIL
  # straightaway if we find that matplotlib.image has been loaded. Afterwards,
  # we patch up all modules (including PIL).
  if any(m.__name__ == _MATPLOTLIB_IMAGE_SELF_NAME
         for m in sys.modules.values()):
    # Just make sure the PIL image module gets entered into sys.modules.
    # It is actually OK to do this even if we are not quite sure that
    # the module named as if it were matplotlib.image actually is that module.
    import PIL.Image  # pylint:disable=g-import-not-at-top,unused-variable
  # Patch up modules.
  _opencv_try_patch_modules(debug=debug)
  _imageio_try_patch_modules(debug=debug)
  _pil_pillow_try_patch_modules(debug=debug)
  _pythonmagick_try_patch_modules(debug=debug)


# We patch everything at import time.
_try_patch_all(debug=int(os.getenv('JXL_DEBUG', _DEFAULT_DEBUG_ENV)))
